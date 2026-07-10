#' Inspect and visualise a qvar_ms variance wrapper
#'
#' @description Given a variance estimation wrapper produced by
#'   \code{qvar_ms()}, this function extracts the internal per-stage
#'   structure (\code{technical_data$stages} and \code{technical_data$calib})
#'   and produces a diagnostic summary and an HTML report of the actual
#'   variance estimation pipeline, stage by stage.
#'
#'   Unlike \code{inspect_wrapper()}, the pipeline is derived from the
#'   technical data rather than from the variance function code, so that
#'   only the steps that will actually be executed are displayed (e.g. no
#'   calibration step at stages without calibration, one sampling variance
#'   step per stage).
#'
#' @param wrapper A variance estimation wrapper produced by \code{qvar_ms()}.
#' @param plot    Logical. If \code{TRUE}, generate and open an HTML report.
#' @param file    Character. Path for the HTML output file. If \code{NULL},
#'   a temporary file is used.
#' @param title   Character. Title of the HTML report.
#'
#' @return A list (invisibly) with elements \code{reference},
#'   \code{stages}, \code{calibration} and \code{pipeline}.
#'
#' @keywords internal

inspect_wrapper_ms <- function(wrapper, plot = TRUE, file = NULL,
                               title = "\U0001F50D Multistage Variance Wrapper Inspection") {
  
  # ══════════════════════════════════════════════════════════════
  # Part 1: Extract and validate wrapper internals
  # ══════════════════════════════════════════════════════════════
  
  if (!is.function(wrapper))
    stop("wrapper should be a variance estimation wrapper (a function).")
  env  <- environment(wrapper)
  form <- formals(wrapper)
  
  td <- tryCatch(env$technical_data, error = function(e) NULL)
  is_ms <- is.list(td) && all(c("stages", "calib") %in% names(td)) &&
    is.list(td$stages) && length(td$stages) >= 1
  if (!is_ms) stop(
    "This wrapper does not look like a qvar_ms() wrapper: technical_data ",
    "should contain a 'stages' list and a 'calib' element. ",
    "Use inspect_wrapper() for single-stage qvar() wrappers."
  )
  
  definition_log <- attr(wrapper, "definition_log")
  definition_call <- attr(wrapper, "definition_call")
  ref_info <- extract_reference_info(env, form)
  ms       <- extract_ms_structure(td)
  pipeline <- build_ms_pipeline(ms)
  
  # Technical parameters
  vf_args <- names(formals(env$variance_function))
  tp_names <- setdiff(vf_args, c("y","stages", "calib"))
  tp_info  <- extract_technical_param_info(form,tp_names)
  
  # ══════════════════════════════════════════════════════════════
  # Part 2: Print console summary
  # ══════════════════════════════════════════════════════════════

  print_ms_summary(ref_info, ms, tp_info, pipeline)
  
  # ══════════════════════════════════════════════════════════════
  # Part 3: Build and display the HTML report
  # ══════════════════════════════════════════════════════════════
  
  if (plot) {
    html <- build_ms_html(title, definition_log, definition_call, 
                          ref_info, tp_info, ms, pipeline)
    if (is.null(file)) file <- tempfile(fileext = ".html")
    writeLines(html, file)
    message("Report written to: ", file)
    if (interactive()) utils::browseURL(file)
  }
  
  invisible(list(
    reference   = ref_info,
    stages      = ms$stages,
    calibration = ms$calib,
    pipeline    = pipeline
  ))
}


# ══════════════════════════════════════════════════════════════════
# Internal: extract the per-stage structure from technical_data
# ══════════════════════════════════════════════════════════════════

extract_ms_structure <- function(td) {
  
  stages <- td$stages
  calib  <- td$calib
  L <- length(stages)
  
  fmt_range <- function(x) {
    x <- suppressWarnings(range(x, na.rm = TRUE))
    if (any(!is.finite(x))) return(NA_character_)
    paste0("[", format(x[1], digits = 4), ", ", format(x[2], digits = 4), "]")
  }
  
  stages_info <- lapply(seq_len(L), function(k) {
    
    lv <- stages[[k]]
    precalc <- lv$samp$precalc
    
    n_units    <- length(lv$samp$id)
    n_excluded <- sum(lv$samp$exclude, na.rm = TRUE)
    n_strata     <- tryCatch(NROW(precalc$A), error = function(e) NA_integer_)
    n_exhaustive <- tryCatch(sum(precalc$exh), error = function(e) NA_integer_)
    pik_range    <- tryCatch(fmt_range(precalc$pik), error = function(e) NA_character_)
    
    nrc <- if (is.null(lv$nrc)) NULL else list(
      n_resp  = length(lv$nrc$id),
      p_range = fmt_range(lv$nrc$response_prob),
      p_mean  = base::mean(lv$nrc$response_prob, na.rm = TRUE)
    )
    
    list(
      stage        = k,
      is_last      = k == L,
      n_units      = n_units,
      n_excluded   = n_excluded,
      n_kept       = n_units - n_excluded,
      n_strata     = n_strata,
      n_exhaustive = n_exhaustive,
      pik_range    = pik_range,
      agg_range    = fmt_range(lv$agg_weight),
      upper_range  = fmt_range(lv$upper_weight),
      is_top       = all(abs(lv$upper_weight - 1) < 1e-12, na.rm = TRUE),
      n_parents    = if (!is.null(lv$parent)) length(unique(lv$parent)) else NA_integer_,
      nrc          = nrc,
      has_calib    = !is.null(calib) && isTRUE(calib$stage == k)
    )
  })
  
  calib_info <- if (is.null(calib)) NULL else list(
    stage   = calib$stage,
    n_units = length(calib$id),
    n_var   = tryCatch({
      pc <- calib$precalc
      if (!is.null(pc$x)) NCOL(pc$x) else NA_integer_
    }, error = function(e) NA_integer_)
  )
  
  list(stages = stages_info, calib = calib_info, L = L, raw = td)
}


# ══════════════════════════════════════════════════════════════════
# Internal: build the pipeline from the structure (NOT from the AST)
# ══════════════════════════════════════════════════════════════════
# Reproduit l'ordre exact de qvar_ms_variance_function() : boucle
# remontante du degre L au degre 1, avec pour chaque degre
# add_zero -> [res_cal] -> [var_pois + reponderation] -> var_srs -> [sum_by]

build_ms_pipeline <- function(ms) {
  
  steps <- list()
  i <- 0
  add_step <- function(...) {
    i <<- i + 1
    steps[[i]] <<- c(list(step = i), list(...))
  }
  suffix <- function(k) if (ms$L > 1) paste0("_stage", k) else ""
  
  variance_names <- character(0)
  
  for (k in rev(seq_len(ms$L))) {
    
    li <- ms$stages[[k]]
    
    add_step(
      stage = k, function_name = "add_zero", category = "transition",
      label = paste0("Expand y to all stage-", k, " sampled units (zeros elsewhere)"),
      formula = "",
      detail = paste0(li$n_units, " units")
    )
    
    if (li$has_calib) add_step(
      stage = k, function_name = "res_cal", category = "calibration",
      label = "Calibration residuals",
      formula = "y_i - x_i' \\hat{B}",
      detail = paste0(
        ms$calib$n_units, " calibrated units",
        if (!is.na(ms$calib$n_var)) paste0(", ", ms$calib$n_var, " calibration variables (discretized)") else ""
      )
    )
    
    if (!is.null(li$nrc)) {
      vn <- paste0("nr", suffix(k))
      variance_names <- c(variance_names, vn)
      add_step(
        stage = k, function_name = "var_pois", category = "variance_nr",
        label = paste0("Non-response variance (Poisson) \u2014 stage ", k),
        formula = "\\sum_{i \\in r_k} w^{sup}_i \\pi_{k,i}^{-1} (1 - \\hat{p}_i) \\left( \\frac{y_i}{\\hat{p}_i} \\right)^{2}",
        detail = paste0(
          li$nrc$n_resp, " responding units, response prob. in ", li$nrc$p_range,
          " \u2014 weight w = upper_weight \u00d7 conditional sampling weight"
        ),
        variance_name = vn
      )
      add_step(
        stage = k, function_name = "Diagonal", category = "transition",
        label = paste0("Non-response reweighting \u2014 stage ", k),
        formula = "y_i \\leftarrow y_i / \\hat{p}_i",
        detail = paste0("expansion of the ", li$nrc$n_resp, " responding units")
      )
    }
    
    vn <- paste0("sampling", suffix(k))
    variance_names <- c(variance_names, vn)
    add_step(
      stage = k, function_name = "var_srs", category = "variance_sampling",
      label = paste0("Sampling variance (stratified SRS) \u2014 stage ", k),
      formula = "\\sum_{h} w^{sup}_h N_h^2 \\left(1 - \\frac{n_h}{N_h}\\right) \\frac{s_h^2}{n_h}",
      detail = paste0(
        li$n_kept, " units",
        if (li$n_excluded > 0) paste0(" (", li$n_excluded, " excluded: single-unit strata)") else "",
        if (!is.na(li$n_strata)) paste0(", ", li$n_strata, " strata",
                                        if (li$stage > 1) " (parent \u00d7 declared strata)" else "") else "",
        if (li$is_top) " \u2014 w = 1 (top stage)" else " \u2014 w = upper_weight"
      ),
      variance_name = vn
    )
    
    if (k > 1) add_step(
      stage = k, function_name = "sum_by", category = "aggregation",
      label = paste0("Aggregate to stage ", k - 1, " (estimated totals per parent unit)"),
      formula = "\\hat{Y}_j = \\sum_{i \\in j} \\frac{y_i}{\\pi_{k,i}}",
      detail = paste0(li$n_parents, " parent units")
    )
  }
  
  add_step(
    stage = NA_integer_, function_name = "Reduce", category = "summation",
    label = "Sum variance components",
    formula = paste0("\\hat{V} = ", paste0("\\hat{V}_{\\text{", gsub("_", "\\\\_", variance_names), "}}", collapse = " + ")),
    detail = paste0(length(variance_names), " components: ", paste(variance_names, collapse = ", "))
  )
  
  steps
}


# ══════════════════════════════════════════════════════════════════
# Internal: recursive tree description of technical data
# ══════════════════════════════════════════════════════════════════
# Remplace la description a un seul niveau d'imbrication de
# inspect_wrapper() : descend recursivement dans les listes (stages ->
# stage k -> samp/nrc -> precalc -> ...) jusqu'a max_depth.

describe_tree <- function(obj, name, depth = 0, max_depth = 4) {
  
  if (is.null(obj)) {
    return(list(list(name = name, type = "NULL", dimension = "-",
                     first_values = NULL, depth = depth)))
  }
  
  info <- describe_object(obj, name)
  info$depth <- depth
  rows <- list(info)
  
  if (is.list(obj) && !is.data.frame(obj) && depth < max_depth) {
    nms <- names(obj)
    if (is.null(nms)) nms <- paste0("[[", seq_along(obj), "]]")
    nms[nms == ""] <- paste0("[[", which(nms == ""), "]]")
    for (j in seq_along(obj)) {
      rows <- c(rows, describe_tree(obj[[j]], nms[j], depth + 1, max_depth))
    }
  }
  rows
}


# ══════════════════════════════════════════════════════════════════
# Internal: console summary
# ══════════════════════════════════════════════════════════════════

print_ms_summary <- function(ref_info, ms, tp_info, pipeline) {
  
  cat("\n")
  cat("============================================================\n")
  cat("  MULTISTAGE VARIANCE WRAPPER INSPECTION (qvar_ms)\n")
  cat("============================================================\n\n")
  
  cat("-- Reference population --\n")
  cat("  Analysis id:    ", ref_info$default_id, "\n")
  cat("  N units:        ", ref_info$n_units, "\n")
  cat("  First 5 ids:    ", paste(ref_info$id_first_5, collapse = ", "), "\n")
  cat("  Positive weight range: [",
      format(ref_info$weight_range[1], digits = 4), ", ",
      format(ref_info$weight_range[2], digits = 4), "]\n", sep = "")
  cat("  Zero-weight:    ", ref_info$n_zero_weight, " units\n\n")
  
  cat("-- Sampling design: ", ms$L, " stage(s) --\n\n", sep = "")
  
  for (li in ms$stages) {
    cat("  Stage ", li$stage,
        if (li$stage == 1) " (primary units)" else if (li$is_last) " (final units)" else "",
        "\n", sep = "")
    cat("    Units:                 ", li$n_units,
        if (li$n_excluded > 0) paste0("  (", li$n_excluded, " excluded: single-unit strata)") else "",
        "\n", sep = "")
    if (!is.na(li$n_strata))
      cat("    Strata:                ", li$n_strata,
          if (li$stage > 1) "  (parent x declared strata)" else "", "\n", sep = "")
    if (!is.na(li$n_exhaustive) && li$n_exhaustive > 0)
      cat("    Exhaustive units:      ", li$n_exhaustive, " (pik = 1, discarded by varDT)\n", sep = "")
    cat("    Cond. sampling weight: ", li$agg_range, "\n", sep = "")
    cat("    Upper weight:    ",
        if (li$is_top) "1 (top stage)" else li$upper_range, "\n", sep = "")
    if (!is.null(li$nrc)) {
      cat("    Non-response corr.:     YES - ", li$nrc$n_resp, " respondents, ",
          "response prob. in ", li$nrc$p_range,
          " (mean ", format(li$nrc$p_mean, digits = 3), ")\n", sep = "")
    } else {
      cat("    Non-response corr.:     no\n")
    }
    if (li$has_calib) {
      cat("    Calibration:            YES - ", ms$calib$n_units, " units",
          if (!is.na(ms$calib$n_var)) paste0(", ", ms$calib$n_var, " variables") else "",
          "\n", sep = "")
    } else {
      cat("    Calibration:            no\n")
    }
    if (!is.na(li$n_parents))
      cat("    Aggregates to:          ", li$n_parents, " parent units (stage ",
          li$stage - 1, ")\n", sep = "")
    cat("\n")
  }
  
  if (length(tp_info) > 0){
    cat("-- Technical parameters --\n")
    for (tp in tp_info) {
      if (length(tp$type) == 0 || tp$type == "NULL") {
        cat("  ", tp$name, ": NULL\n", sep = "")
        next
      }
      if (tp$type == "list") {
        cat("  ", tp$name, " [list, ", tp$n_elements, " elements]\n", sep = "")
        for (el in tp$elements) {
          cat("    - ", el$name, " <", el$type, "> ",
              el$dimension, sep = "")
          if (!is.null(el$first_values)) {
            cat("  (", el$first_values, ")", sep = "")
          }
          cat("\n")
        }
      } else {
        cat("  ", tp$name, " <", tp$type, "> ", sep = "")
        if (!is.null(tp$dimension)) cat(tp$dimension)
        cat("\n")
      }
    }
    cat("\n")
  }
  
  cat("-- Variance estimation pipeline (execution order) --\n")
  for (step in pipeline) {
    icon <- switch(step$category,
                   "calibration"       = "\u2500[CAL]\u2500\u2500",
                   "variance_nr"       = "\u2500[V_NR]\u2500",
                   "variance_sampling" = "\u2500[V_SMP]\u2500",
                   "transition"        = "\u2500[y]\u2500\u2500\u2500",
                   "aggregation"       = "\u2500[AGG]\u2500\u2500",
                   "summation"         = "\u2500[SUM]\u2500\u2500",
                   "\u2500[???]\u2500"
    )
    cat("  ", sprintf("%2d", step$step),
        if (!is.na(step$stage)) sprintf(" L%d ", step$stage) else "    ",
        icon, " ", step$function_name, "()  ", step$label, "\n", sep = "")
    if (nchar(step$detail) > 0 && nchar(step$detail) < 150)
      cat("           ", step$detail, "\n", sep = "")
  }
  cat("\n")
}


# ══════════════════════════════════════════════════════════════════
# Internal: Mermaid diagram with one subgraph per stage
# ══════════════════════════════════════════════════════════════════

build_ms_mermaid <- function(pipeline, ms) {
  
  node_shape <- function(step, sid) {
    fn_esc    <- mermaid_esc(step$function_name)
    label_esc <- mermaid_esc(step$label)
    if (step$category %in% c("variance_nr", "variance_sampling")) {
      lbl <- paste0('<b>', fn_esc, '()</b><br/>', label_esc,
                    if (!is.null(step$variance_name))
                      paste0('<br/><i>', mermaid_esc(step$variance_name), '</i>') else "")
      paste0('    ', sid, '(["', lbl, '"])')
    } else {
      paste0('    ', sid, '["<b>', fn_esc, '()</b><br/>', label_esc, '"]')
    }
  }
  
  # ── Node declarations, grouped by stage subgraph ──
  lines <- c("graph TB",
             '  START(["<b>y</b><br/>Variable of interest<br/>(last stage, disseminated units)"])')
  
  stage_of  <- vapply(pipeline, function(s) ifelse(is.na(s$stage), -1L, as.integer(s$stage)), integer(1))
  for (k in rev(seq_len(ms$L))) {
    li <- ms$stages[[k]]
    sub_title <- paste0("Stage ", k,
                        if (k == 1) " \u2014 primary units" else if (li$is_last) " \u2014 final units" else "",
                        " \u00b7 ", li$n_units, " units")
    lines <- c(lines, paste0('  subgraph LVL', k, '["', mermaid_esc(sub_title), '"]'),
               '    direction LR')
    for (idx in which(stage_of == k)) {
      lines <- c(lines, node_shape(pipeline[[idx]], paste0("S", pipeline[[idx]]$step)))
    }
    lines <- c(lines, '  end')
  }
  
  # Summation and result nodes (outside subgraphs)
  sum_idx <- which(vapply(pipeline, function(s) s$category == "summation", logical(1)))
  sum_sid <- paste0("S", pipeline[[sum_idx]]$step)
  lines <- c(lines,
             paste0('  ', sum_sid, '["<b>Reduce</b><br/>Sum variance components"]'),
             '  RESULT(["<b>V&#40;Y&#770;&#41;</b>"])')
  
  # ── Edges: main chain through transforms, dotted branches to variances ──
  prev <- "START"
  variance_nodes <- character(0)
  for (step in pipeline) {
    sid <- paste0("S", step$step)
    if (step$category %in% c("variance_nr", "variance_sampling")) {
      lines <- c(lines, paste0('  ', prev, ' -.-> ', sid))
      variance_nodes <- c(variance_nodes, sid)
    } else if (step$category == "summation") {
      for (vn in variance_nodes) lines <- c(lines, paste0('  ', vn, ' --> ', sid))
      lines <- c(lines, paste0('  ', sid, ' --> RESULT'))
    } else {
      lines <- c(lines, paste0('  ', prev, ' --> ', sid))
      prev <- sid
    }
  }
  
  # ── Styles ──
  for (step in pipeline) {
    sid <- paste0("S", step$step)
    style <- switch(step$category,
                    "calibration"       = "fill:#e8f5e9,stroke:#4caf50,color:#1b5e20",
                    "variance_nr"       = "fill:#fff3e0,stroke:#ff9800,color:#e65100",
                    "variance_sampling" = "fill:#e3f2fd,stroke:#2196f3,color:#0d47a1",
                    "transition"        = "fill:#f3e5f5,stroke:#9c27b0,color:#4a148c",
                    "aggregation"       = "fill:#fce4ec,stroke:#e91e63,color:#880e4f",
                    "summation"         = "fill:#eceff1,stroke:#607d8b,color:#263238",
                    NULL)
    if (!is.null(style)) lines <- c(lines, paste0('  style ', sid, ' ', style))
  }
  for (k in seq_len(ms$L)) {
    lines <- c(lines, paste0('  style LVL', k,
                             ' fill:#f8f9fc,stroke:#c5cae9,stroke-width:1px'))
  }
  lines <- c(lines, '  style START fill:#e8eaf6,stroke:#3f51b5,color:#1a237e')
  lines <- c(lines, '  style RESULT fill:#e8eaf6,stroke:#3f51b5,color:#1a237e')
  
  paste(lines, collapse = "\n")
}


# ══════════════════════════════════════════════════════════════════
# Internal: HTML building blocks
# ══════════════════════════════════════════════════════════════════

build_stage_card <- function(li, calib_info) {
  
  row <- function(label, value) paste0(
    '<tr><td>', label, '</td><td>', value, '</td></tr>'
  )
  
  nrc_cell <- if (!is.null(li$nrc)) paste0(
    '<span class="badge badge-vnr">YES</span> ', li$nrc$n_resp,
    ' respondents \u00b7 response prob. in <code>', htmlesc(li$nrc$p_range), '</code>',
    ' (mean ', format(li$nrc$p_mean, digits = 3), ')'
  ) else '<span class="badge badge-off">no</span>'
  
  cal_cell <- if (li$has_calib) paste0(
    '<span class="badge badge-cal">YES</span> ', calib_info$n_units, ' units',
    if (!is.na(calib_info$n_var)) paste0(' \u00b7 ', calib_info$n_var, ' calibration variables') else ''
  ) else '<span class="badge badge-off">no</span>'
  
  paste0(
    '<div class="card">',
    '<h2><span class="stage-badge">Stage ', li$stage, '</span> ',
    if (li$stage == 1) 'Primary units' else if (li$is_last) 'Final units' else 'Intermediate units',
    '</h2>',
    '<table class="info-table">',
    row('Units', paste0(li$n_units,
                        if (li$n_excluded > 0) paste0(' <span class="muted">(', li$n_excluded,
                                                      ' excluded: single-unit strata)</span>') else '')),
    if (!is.na(li$n_strata)) row('Strata', paste0(li$n_strata,
                                                  if (li$stage > 1) ' <span class="muted">(parent \u00d7 declared strata)</span>' else '')) else '',
    if (!is.na(li$n_exhaustive) && li$n_exhaustive > 0)
      row('Exhaustive units (pik = 1)', li$n_exhaustive) else '',
    row('Conditional sampling weight', paste0('<code>', htmlesc(li$agg_range), '</code>')),
    row('Upper weight',
        if (li$is_top) '1 <span class="muted">(top stage)</span>'
        else paste0('<code>', htmlesc(li$upper_range), '</code>')),
    row('Non-response correction', nrc_cell),
    row('Calibration', cal_cell),
    if (!is.na(li$n_parents))
      row('Aggregates to', paste0(li$n_parents, ' parent units <span class="muted">(stage ',
                                  li$stage - 1, ')</span>')) else '',
    '</table>',
    '</div>'
  )
}


build_tree_table <- function(raw_td) {
  
  rows_data <- describe_tree(raw_td, "technical_data", depth = 0, max_depth = 4)
  
  rows <- vapply(rows_data, function(r) {
    indent <- paste0('padding-left: ', 0.5 + r$depth * 1.3, 'rem !important;')
    prefix <- if (r$depth > 0) '\u2514 ' else ''
    preview <- if (!is.null(r$first_values))
      paste0('<code>', htmlesc(r$first_values), '</code>') else ''
    paste0(
      '<tr><td style="', indent, '">', prefix,
      if (r$depth <= 1) paste0('<b>', htmlesc(r$name), '</b>') else htmlesc(r$name),
      '</td>',
      '<td><code>', htmlesc(r$type), '</code></td>',
      '<td>', htmlesc(r$dimension), '</td>',
      '<td>', preview, '</td></tr>'
    )
  }, character(1))
  
  paste0(
    '<div class="card">',
    '<details><summary><h2 style="display:inline">Raw technical data structure</h2>',
    ' <span class="muted">(click to expand)</span></summary>',
    '<table class="td-table" style="margin-top: 0.8rem;">',
    '<thead><tr><th>Name</th><th>Type</th><th>Dimension</th><th>Preview</th></tr></thead>',
    '<tbody>', paste(rows, collapse = "\n"), '</tbody>',
    '</table>',
    '</details>',
    '</div>'
  )
}


build_ms_html <- function(title, definition_log, definition_call, 
                          ref_info, tp_info, ms, pipeline) {
  
  mermaid <- build_ms_mermaid(pipeline, ms)
  
  # ── Reference card ──
  ref_html <- paste0(
    '<div class="card">',
    '<h2>Reference Population</h2>',
    '<table class="info-table">',
    '<tr><td>Analysis id</td><td><code>', htmlesc(ref_info$default_id), '</code></td></tr>',
    '<tr><td>Sampling stages</td><td>', ms$L, '</td></tr>',
    '<tr><td>Number of units</td><td>', ref_info$n_units, '</td></tr>',
    '<tr><td>First 5 ids</td><td><code>',
    htmlesc(paste(ref_info$id_first_5, collapse = ", ")), '</code></td></tr>',
    '<tr><td>Positive weight range</td><td>[',
    format(ref_info$weight_range[1], digits = 4), ', ',
    format(ref_info$weight_range[2], digits = 4), ']</td></tr>',
    '<tr><td>Zero-weight units</td><td>', ref_info$n_zero_weight, '</td></tr>',
    '</table>',
    '</div>'
  )
  
  notes_html <- build_notes_card(definition_log, definition_call)
  
  if (length(tp_info) > 0){
    tp_html <- build_td_table(tp_info, "Technical Parameters")
  }else{
    tp_html <- ""
  }
  
  # ── One card per stage (stage 1 first) ──
  stage_cards <- paste(
    vapply(ms$stages, build_stage_card, character(1), calib_info = ms$calib),
    collapse = "\n"
  )
  
  # ── Pipeline table (with a Stage column) ──
  pipe_rows <- vapply(pipeline, function(step) {
    badge_class <- switch(step$category,
                          "calibration"       = "badge-cal",
                          "variance_nr"       = "badge-vnr",
                          "variance_sampling" = "badge-vsmp",
                          "transition"        = "badge-trans",
                          "aggregation"       = "badge-agg",
                          "summation"         = "badge-sum",
                          "badge-other")
    formula_cell <- if (nchar(step$formula) > 0) paste0('\\(', step$formula, '\\)') else ''
    vname_cell <- if (!is.null(step$variance_name))
      paste0('<code>', htmlesc(step$variance_name), '</code>') else ''
    paste0(
      '<tr>',
      '<td>', step$step, '</td>',
      '<td>', if (!is.na(step$stage))
        paste0('<span class="stage-badge">L', step$stage, '</span>') else '', '</td>',
      '<td><span class="badge ', badge_class, '">', htmlesc(step$category), '</span></td>',
      '<td><code>', htmlesc(step$function_name), '()</code></td>',
      '<td>', htmlesc(step$label), '</td>',
      '<td>', vname_cell, '</td>',
      '<td>', formula_cell, '</td>',
      '<td class="args-cell">', htmlesc(step$detail), '</td>',
      '</tr>'
    )
  }, character(1))
  
  pipe_html <- paste0(
    '<div class="card">',
    '<h2>Pipeline Steps (execution order: last stage \u2192 first stage)</h2>',
    '<table class="pipe-table">',
    '<thead><tr><th>#</th><th>Stage</th><th>Type</th><th>Function</th>',
    '<th>Description</th><th>Component</th><th>Formula</th><th>Details</th></tr></thead>',
    '<tbody>', paste(pipe_rows, collapse = "\n"), '</tbody>',
    '</table>',
    '</div>'
  )
  
  tree_html <- build_tree_table(ms$raw)
  
  # ── Assemble full HTML ──
  paste0(
    '<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Multistage Variance Wrapper Inspection</title>
<script src="https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
<style>
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    background: #f5f7fa; color: #333; padding: 2rem; line-height: 1.5;
  }
  h1 { font-size: 1.6rem; margin-bottom: 1.5rem; color: #1a1a2e; }
  h2 { font-size: 1.15rem; margin-bottom: 0.8rem; color: #16213e; }
  .layout {
    display: grid;
    grid-template-columns: 1fr 2fr;
    gap: 1.2rem;
    align-items: start;
  }
  .layout-full {
    display: grid;
    grid-template-columns: 1fr;
    gap: 1.2rem;
    align-items: start;
  }
  .left-panel  { display: flex; flex-direction: column; gap: 1.2rem; }
  .right-panel {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
    gap: 1.2rem;
    align-items: start;
  }
  .card {
    background: #fff; border-radius: 8px; padding: 1.2rem;
    box-shadow: 0 1px 3px rgba(0,0,0,0.08); border: 1px solid #e8ecf1;
  }
  .info-table { width: 100%; border-collapse: collapse; }
  .info-table td { padding: 0.35rem 0.6rem; border-bottom: 1px solid #f0f0f0; font-size: 0.9rem; }
  .info-table td:first-child { font-weight: 600; width: 200px; color: #555; }
  code { background: #f0f3f8; padding: 0.15rem 0.4rem; border-radius: 3px; font-size: 0.85rem; }
  .muted { color: #999; font-size: 0.82rem; }
  .td-table, .pipe-table { width: 100%; border-collapse: collapse; font-size: 0.85rem; }
  .td-table th, .pipe-table th { text-align: left; padding: 0.4rem 0.5rem; background: #f8f9fb;
    border-bottom: 2px solid #e0e4ea; color: #555; font-weight: 600; }
  .td-table td, .pipe-table td { padding: 0.4rem 0.5rem; border-bottom: 1px solid #f0f0f0;
    vertical-align: top; }
  .td-table tr:hover, .pipe-table tr:hover { background: #fafbfd; }
  .args-cell { max-width: 350px; }
  .badge {
    display: inline-block; padding: 0.15rem 0.5rem; border-radius: 10px;
    font-size: 0.75rem; font-weight: 600; letter-spacing: 0.02em;
  }
  .badge-cal  { background: #e8f5e9; color: #2e7d32; }
  .badge-vnr  { background: #fff3e0; color: #e65100; }
  .badge-vsmp { background: #e3f2fd; color: #1565c0; }
  .badge-trans { background: #f3e5f5; color: #7b1fa2; }
  .badge-agg  { background: #fce4ec; color: #c62828; }
  .badge-sum  { background: #eceff1; color: #37474f; }
  .badge-off  { background: #f5f5f5; color: #9e9e9e; }
  .badge-other { background: #f5f5f5; color: #616161; }
  .log-item {
  display: flex; align-items: flex-start; gap: 0.6rem;
  padding: 0.55rem 0.7rem; margin-bottom: 0.5rem;
  border-radius: 6px; border-left: 3px solid;
  font-size: 0.85rem;
  }
  .log-item .badge { flex-shrink: 0; margin-top: 0.1rem; }
  .log-text { white-space: pre-wrap; min-width: 0; }
  .log-mono { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
    font-size: 0.8rem; }
  .log-info    { background: #f8f9fb; border-left-color: #607d8b; }
  .log-note    { background: #f4f9fe; border-left-color: #2196f3; }
  .log-warning { background: #fff8f0; border-left-color: #ff9800; }
  .log-call    { background: #f8f9fb; border-left-color: #9e9e9e; }
  .stage-badge {
    display: inline-block; padding: 0.1rem 0.55rem; border-radius: 4px;
    font-size: 0.78rem; font-weight: 700;
    background: #e8eaf6; color: #283593; border: 1px solid #c5cae9;
  }
  .mermaid-container { background: #fff; border-radius: 8px; padding: 1.5rem;
    box-shadow: 0 1px 3px rgba(0,0,0,0.08); border: 1px solid #e8ecf1;
    margin-bottom: 1.2rem; overflow-x: auto; }
  .mermaid { text-align: center; }
  details summary { cursor: pointer; }
</style>
</head>
<body>

<div class="layout">
  <div class="left-panel">
    <h1>', title, '</h1>
    ', ref_html, '
  </div>
  <div class="right-panel">
  ', build_notes_card(definition_log, definition_call), '
  </div>
</div>
<br/>
<div class="right-panel">
    ', stage_cards, '
</div>
<br/>
<div class="layout-full">
  ', pipe_html, '
  <div class="mermaid-container">
    <h2>Variance Estimation Pipeline</h2>
    <div class="mermaid">
      ', mermaid, '
    </div>
  </div>
  ', tp_html, '
  ', tree_html, '
</div>

<script>mermaid.initialize({ startOnLoad: true, theme: "neutral", flowchart: { useMaxWidth: true, htmlLabels: true } });</script>
</body>
</html>'
  )
}

# ══════════════════════════════════════════════════════════════════
# Internal: definition log card (notes/warnings emitted by qvar_ms)
# ══════════════════════════════════════════════════════════════════

build_notes_card <- function(definition_log, definition_call = NULL) {
  
  if (is.null(definition_log) || length(definition_log) == 0) return("")
  
  items <- vapply(definition_log, function(entry) {
    type <- if (startsWith(entry, "Warning: ")) "warning"
    else if (startsWith(entry, "Note: ")) "note"
    else "info"
    text  <- sub("^(Warning|Note): ", "", entry)
    badge <- switch(type,
                    "warning" = '<span class="badge badge-vnr">WARNING</span>',
                    "note"    = '<span class="badge badge-vsmp">NOTE</span>',
                    "info"    = '<span class="badge badge-sum">INFO</span>'
    )
    paste0(
      '<div class="log-item log-', type, '">', badge,
      '<div class="log-text">', htmlesc(text), '</div>',
      '</div>'
    )
  }, character(1))
  
  call_html <- if (!is.null(definition_call)) paste0(
    '<div class="log-item log-call">',
    '<span class="badge badge-other">CALL</span>',
    '<div class="log-text log-mono">',
    htmlesc(paste(deparse(definition_call, width.cutoff = 80), collapse = "\n")),
    '</div></div>'
  ) else ""
  
  paste0(
    '<div class="card">',
    '<h2>Methodological notes</h2>',
    call_html,
    paste(items, collapse = "\n"),
    '</div>'
  )
}