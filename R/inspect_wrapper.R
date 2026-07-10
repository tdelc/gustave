#' Inspect and visualise a gustave variance wrapper
#'
#' @description Given a variance estimation wrapper produced by
#'   \code{qvar()}, or \code{define_variance_wrapper()}, this function 
#'   extracts the internal structure and produces a diagnostic summary 
#'   and an HTML diagram of the variance estimation pipeline.
#'
#' @param wrapper A variance estimation wrapper function.
#' @param plot    Logical. If \code{TRUE}, generate and open an HTML
#'   diagram of the pipeline.
#' @param file    Character. Path for the HTML output file. If \code{NULL},
#'   a temporary file is used.
#'
#' @return A list (invisibly) with the following elements:
#'   \describe{
#'     \item{reference}{List with reference_id and reference_weight summary.}
#'     \item{technical_data}{List describing each technical data element.}
#'     \item{pipeline}{List describing the variance function steps.}
#'   }
#'
#' @export

inspect_wrapper <- function(wrapper, plot = TRUE, file = NULL,
                            title = "🔍 Variance Wrapper Inspection") {
  
  # ══════════════════════════════════════════════════════════════
  # Part 1: Extract wrapper internals
  # ══════════════════════════════════════════════════════════════
  
  if (!is.function(wrapper))
    stop("wrapper should be a variance estimation wrapper (a function).")
  env  <- environment(wrapper)
  form <- formals(wrapper)
  
  td <- tryCatch(env$technical_data, error = function(e) NULL)
  is_ms <- is.list(td) && all(c("stages", "calib") %in% names(td)) &&
    is.list(td$stages) && length(td$stages) >= 1
  if (is_ms){
    warning(
      "This wrapper does look like a qvar_ms() wrapper: technical_data 
      contain a 'stages' list and a 'calib' element."
    )
    return(inspect_wrapper_ms(wrapper, plot = plot, file = file, title = title))
  } 
  
  # Extract technical data and technical paramaters
  vf <- env$variance_function
  vf_args <- names(formals(vf))
  form_args <- names(form)
  td_names <- setdiff(vf_args, c("y",form_args)) # Remove y and technical parameters
  tp_names <- intersect(vf_args, form_args) # Only technical parameters
  
  ref_info <- extract_reference_info(env,form)
  td_info  <- extract_technical_data_info(env,td_names)
  tp_info  <- extract_technical_param_info(form,tp_names)
  vf_info  <- extract_variance_function_info(env)
  
  # ══════════════════════════════════════════════════════════════
  # Part 2: Print console summary
  # ══════════════════════════════════════════════════════════════
  
  print_summary(ref_info, td_info, tp_info, vf_info)
  
  # ══════════════════════════════════════════════════════════════
  # Part 3: Build and display the HTML diagram
  # ══════════════════════════════════════════════════════════════
  
  if (plot) {
    html <- build_html_report(title, ref_info, td_info, tp_info, vf_info)
    if (is.null(file)) {
      file <- tempfile(fileext = ".html")
    }
    writeLines(html, file)
    message("Report written to: ", file)
    if (interactive()) utils::browseURL(file)
  }
  
  invisible(list(
    reference      = ref_info,
    technical_data = td_info,
    pipeline       = vf_info
  ))
}


# ══════════════════════════════════════════════════════════════════
# Internal: extract reference info
# ══════════════════════════════════════════════════════════════════

extract_reference_info <- function(env, form) {
  
  ref_id <- env$reference_id
  ref_w  <- env$reference_weight
  
  list(
    n_units        = length(ref_id),
    id_first_5     = utils::head(ref_id, 5),
    weight_range   = range(ref_w[ref_w != 0], na.rm = TRUE),
    weight_first_5 = utils::head(ref_w, 5),
    n_zero_weight  = sum(ref_w == 0, na.rm = TRUE),
    default_id     = tryCatch(form$id, error = function(e) NA_character_)
  )
}


# ══════════════════════════════════════════════════════════════════
# Internal: extract technical data inventory
# ══════════════════════════════════════════════════════════════════

extract_technical_data_info <- function(env, td_names) {
  
  lapply(stats::setNames(td_names, td_names), function(nm) {
    obj <- tryCatch(env$technical_data[[nm]], error = function(e) NULL)
    if (is.null(obj)) {
      return(list(name = nm, type = "NULL", detail = "Not provided"))
    }
    
    if (is.list(obj) && !is.data.frame(obj)) {
      # Nested list (typical for psu, hh, mb, calib)
      elements <- lapply(names(obj), function(el_name) {
        el <- obj[[el_name]]
        describe_object(el, el_name)
      })
      list(
        name          = nm,
        type          = "list",
        n_elements    = length(obj),
        element_names = names(obj),
        elements      = elements
      )
    } else if (is.data.frame(obj)) {
      list(
        name      = nm,
        type      = "data.frame",
        dimension = paste(nrow(obj), "rows x", ncol(obj), "cols"),
        col_names = names(obj)
      )
    } else {
      describe_object(obj, nm)
    }
  })
}


extract_technical_param_info <- function(form, tp_names) {
  
  lapply(stats::setNames(tp_names, tp_names), function(nm) {
    obj <- tryCatch(form[[nm]], error = function(e) NULL)
    if (is.null(obj)) {
      return(list(name = nm, type = "NULL", detail = "Not provided"))
    }
    
    if (is.list(obj) && !is.data.frame(obj)) {
      # Nested list (typical for psu, hh, mb, calib)
      elements <- lapply(names(obj), function(el_name) {
        el <- obj[[el_name]]
        describe_object(el, el_name)
      })
      list(
        name          = nm,
        type          = "list",
        n_elements    = length(obj),
        element_names = names(obj),
        elements      = elements
      )
    } else if (is.data.frame(obj)) {
      list(
        name      = nm,
        type      = "data.frame",
        dimension = paste(nrow(obj), "rows x", ncol(obj), "cols"),
        col_names = names(obj)
      )
    } else {
      describe_object(obj, nm)
    }
  })
}

#' Describe a single R object for the summary
#' @keywords internal

describe_object <- function(obj, name) {
  
  obj_class <- paste(class(obj), collapse = "/")
  
  # Dimension string
  if (is.data.frame(obj)) {
    dim_str <- paste(nrow(obj), "rows x", ncol(obj), "cols")
  } else if (is.matrix(obj) || inherits(obj, "Matrix")) {
    d <- dim(obj)
    dim_str <- paste(d[1], "x", d[2])
  } else if (is.atomic(obj) || is.factor(obj)) {
    dim_str <- paste("length", length(obj))
  } else if (is.list(obj)) {
    dim_str <- paste(length(obj), "elements")
  } else if (is.function(obj)) {
    dim_str <- "function"
  } else {
    dim_str <- "?"
  }
  
  # First values preview
  first_vals <- NULL
  if (is.atomic(obj) && length(obj) > 0) {
    fv <- utils::head(obj, 3)
    first_vals <- paste(
      if (!is.null(names(fv))) paste0(names(fv), "=", format(fv, digits = 4)) else format(fv, digits = 4),
      collapse = ", "
    )
  } else if (is.factor(obj) && length(obj) > 0) {
    fv <- utils::head(as.character(obj), 3)
    first_vals <- paste(fv, collapse = ", ")
  }
  
  list(
    name        = name,
    type        = obj_class,
    dimension   = dim_str,
    first_values = first_vals
  )
}


# ══════════════════════════════════════════════════════════════════
# Internal: parse variance function AST
# ══════════════════════════════════════════════════════════════════

extract_variance_function_info <- function(env) {
  
  vf <- env$variance_function
  vf_body <- body(vf)
  
  target_fns <- c("var_srs", "var_pois", "res_cal", "add_zero",
                  "sum_by", "varDT", "Reduce")
  
  steps <- list()
  step_counter <- 0
  seen_calls <- character(0)
  
  # Recursive AST walker
  walk_ast <- function(expr) {
    if (!is.call(expr)) return()
    
    fn_name <- deparse(expr[[1]])
    fn_name <- gsub("var_srs_safe","var_srs",fn_name)
    
    # Handle simple names and namespaced calls (pkg::fn)
    if (fn_name %in% target_fns) {
      # Deduplicate: if the exact same call appears inside a sapply/map,
      # we still only record it once
      call_str <- paste(deparse(expr, width.cutoff = 120), collapse = " ")
      signature <- paste0(fn_name, ":", call_str)
      
      if (!signature %in% seen_calls) {
        seen_calls <<- c(seen_calls, signature)
        step_counter <<- step_counter + 1
        steps[[step_counter]] <<- parse_step(expr, fn_name, step_counter)
      }
    }
    
    # Recurse into sub-expressions
    for (i in seq_along(expr)[-1]) {
      if (!is.null(expr[[i]])) walk_ast(expr[[i]])
    }
  }
  
  parse_step <- function(expr, fn_name, step_num) {
    args <- as.list(expr)[-1]
    arg_names <- names(args)
    if (is.null(arg_names)) arg_names <- rep("", length(args))
    
    list(
      step          = step_num,
      function_name = fn_name,
      category      = categorise_function(fn_name),
      label         = label_function(fn_name),
      formula       = formula_function(fn_name),
      args_summary  = summarise_args(args, arg_names),
      raw_call      = paste(deparse(expr, width.cutoff = 80), collapse = " ")
    )
  }
  
  walk_ast(vf_body)
  steps
}


categorise_function <- function(fn_name) {
  switch(fn_name,
         "res_cal"  = "calibration",
         "var_pois" = "variance_nr",
         "var_srs"  = "variance_sampling",
         "varDT"    = "variance_sampling",
         "add_zero" = "transition",
         "sum_by"   = "aggregation",
         "Reduce"   = "summation",
         "other"
  )
}

label_function <- function(fn_name) {
  switch(fn_name,
         "res_cal"  = "Calibration residuals",
         "var_pois" = "Non-response variance (Poisson)",
         "var_srs"  = "Sampling variance (SRS)",
         "varDT"    = "Sampling variance (Deville-Till\u00e9)",
         "add_zero" = "Reintroduce zeros",
         "sum_by"   = "Aggregate",
         "Reduce"   = "Sum variance components",
         fn_name
  )
}

formula_function <- function(fn_name) {
  switch(fn_name,
         "res_cal"  = "y_k - x_k' \\hat{\\beta}",
         "var_pois" = "\\sum_{i \\in \\text{resp}} \\frac{1-p_i}{p_i} \\left(\\frac{y_i}{\\pi_i}\\right)^2",
         "var_srs"  = "\\sum_{h=1}^{H} N_h^2 \\left(1 - \\frac{n_h}{N_h}\\right) \\frac{s_h^2}{n_h}",
         "varDT"    = "\\text{Deville-Till\\'{e} approximation}",
         "add_zero" = "",
         "sum_by"   = "\\sum_{h \\in j} \\frac{y_h}{\\pi_{h|j} \\cdot p_h}",
         "Reduce"   = "",
         ""
  )
}

summarise_args <- function(args, arg_names) {
  if (length(args) == 0) return("")
  summaries <- mapply(function(arg, nm) {
    arg_str <- paste(deparse(arg, width.cutoff = 50), collapse = " ")
    arg_str <- gsub(", , drop = FALSE","",arg_str)
    if (nchar(arg_str) > 50) arg_str <- paste0(substr(arg_str, 1, 47), "...")
    if (nm != "" && !is.null(nm)) paste0(nm, " = ", arg_str) else arg_str
  }, args, arg_names, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  paste(summaries, collapse = ", ")
}


# ══════════════════════════════════════════════════════════════════
# Internal: console summary
# ══════════════════════════════════════════════════════════════════

print_summary <- function(ref_info, td_info, tp_info, vf_info) {
  
  cat("\n")
  cat("============================================================\n")
  cat("  VARIANCE WRAPPER INSPECTION\n")
  cat("============================================================\n\n")
  
  # Reference info
  cat("-- Reference population --\n")
  cat("  Analysis id:    ", ref_info$default_id, "\n")
  cat("  N units:        ", ref_info$n_units, "\n")
  cat("  First 5 ids:    ", paste(ref_info$id_first_5, collapse = ", "), "\n")
  cat("  Positive Weight range:    [",
      format(ref_info$weight_range[1], digits = 4), ", ",
      format(ref_info$weight_range[2], digits = 4), "]\n", sep = "")
  cat("  Zero-weight:    ", ref_info$n_zero_weight, " units\n\n")
  
  # Technical data
  cat("-- Technical data --\n")
  for (td in td_info) {
    if (length(td$type) == 0 || td$type == "NULL") {
      cat("  ", td$name, ": NULL\n", sep = "")
      next
    }
    if (td$type == "list") {
      cat("  ", td$name, " [list, ", td$n_elements, " elements]\n", sep = "")
      for (el in td$elements) {
        cat("    - ", el$name, " <", el$type, "> ",
            el$dimension, sep = "")
        if (!is.null(el$first_values)) {
          cat("  (", el$first_values, ")", sep = "")
        }
        cat("\n")
      }
    } else {
      cat("  ", td$name, " <", td$type, "> ", sep = "")
      if (!is.null(td$dimension)) cat(td$dimension)
      cat("\n")
    }
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
  }
  
  
  # Pipeline
  cat("\n-- Variance estimation pipeline --\n")
  for (step in vf_info) {
    icon <- switch(step$category,
                   "calibration"       = "\u2500[CAL]\u2500",
                   "variance_nr"       = "\u2500[V_NR]\u2500",
                   "variance_sampling" = "\u2500[V_SMP]\u2500",
                   "transition"        = "\u2500[+0]\u2500\u2500",
                   "aggregation"       = "\u2500[AGG]\u2500",
                   "summation"         = "\u2500[SUM]\u2500",
                   "\u2500[???]\u2500"
    )
    cat("  ", sprintf("%2d", step$step), " ", icon, " ",
        step$function_name, "()\n", sep = "")
    cat("       ", step$label, "\n", sep = "")
    if (nchar(step$args_summary) > 0 && nchar(step$args_summary) < 120) {
      cat("        args: ", step$args_summary, "\n", sep = "")
    }
  }
  cat("\n")
}


# ══════════════════════════════════════════════════════════════════
# Internal: build HTML report with Mermaid diagram
# ══════════════════════════════════════════════════════════════════

build_html_report <- function(title, ref_info, td_info, tp_info, vf_info) {
  
  # ── Build mermaid diagram ──
  mermaid <- build_mermaid(vf_info)
  
  # ── Build technical data HTML table ──
  td_html <- build_td_table(td_info)
  
  if (length(tp_info) > 0){
    tp_html <- build_td_table(tp_info, "Technical Parameters")
  }else{
    tp_html <- ""
  }
  
  # ── Build reference info HTML ──
  ref_html <- paste0(
    '<div class="card">',
    '<h2>Reference Population</h2>',
    '<table class="info-table">',
    '<tr><td>Analysis id</td><td><code>', htmlesc(ref_info$default_id), '</code></td></tr>',
    '<tr><td>Number of units</td><td>', ref_info$n_units, '</td></tr>',
    '<tr><td>First 5 ids</td><td><code>',
    htmlesc(paste(ref_info$id_first_5, collapse = ", ")), '</code></td></tr>',
    '<tr><td>Positive Weight range</td><td>[',
    format(ref_info$weight_range[1], digits = 4), ', ',
    format(ref_info$weight_range[2], digits = 4), ']</td></tr>',
    '<tr><td>First 5 weights</td><td><code>',
    htmlesc(paste(format(ref_info$weight_first_5, digits = 4), collapse = ", ")),
    '</code></td></tr>',
    '<tr><td>Zero-weight units</td><td>', ref_info$n_zero_weight, '</td></tr>',
    '</table>',
    '</div>'
  )
  
  # ── Build pipeline detail table ──
  pipe_rows <- vapply(vf_info, function(step) {
    badge_class <- switch(step$category,
                          "calibration"       = "badge-cal",
                          "variance_nr"       = "badge-vnr",
                          "variance_sampling" = "badge-vsmp",
                          "transition"        = "badge-trans",
                          "aggregation"       = "badge-agg",
                          "summation"         = "badge-sum",
                          "badge-other"
    )
    
    # Formula cell: wrap in \( \) for MathJax inline rendering
    formula_cell <- if (nchar(step$formula) > 0) {
      paste0('\\(', step$formula, '\\)')
    } else {
      ''
    }
    
    paste0(
      '<tr>',
      '<td>', step$step, '</td>',
      '<td><span class="badge ', badge_class, '">', htmlesc(step$category), '</span></td>',
      '<td><code>', htmlesc(step$function_name), '()</code></td>',
      '<td>', htmlesc(step$label), '</td>',
      '<td>', formula_cell, '</td>',
      '<td class="args-cell"><code>', htmlesc(step$args_summary), '</code></td>',
      '</tr>'
    )
  }, character(1))
  
  pipe_html <- paste0(
    '<div class="card">',
    '<h2>Pipeline Steps</h2>',
    '<table class="pipe-table">',
    '<thead><tr><th>#</th><th>Type</th><th>Function</th><th>Description</th><th>Formula</th><th>Arguments</th></tr></thead>',
    '<tbody>', paste(pipe_rows, collapse = "\n"), '</tbody>',
    '</table>',
    '</div>'
  )
  
  # ── Assemble full HTML ──
  html <- paste0(
    '<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Variance Wrapper Inspection</title>
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
  .left-panel {
  display: flex;
  flex-direction: column;          /* empile les cards verticalement */
  gap: 1.2rem;                     /* espacement entre les cards */
  }
  .right-panel {
    display: flex;
    flex-direction: column;
    gap: 1.2rem;
  }
  .card {
    background: #fff; border-radius: 8px; padding: 1.2rem;
    box-shadow: 0 1px 3px rgba(0,0,0,0.08); border: 1px solid #e8ecf1;
  }
  .card-full { grid-column: 1 / -1; }
  .info-table { width: 100%; border-collapse: collapse; }
  .info-table td { padding: 0.35rem 0.6rem; border-bottom: 1px solid #f0f0f0; font-size: 0.9rem; }
  .info-table td:first-child { font-weight: 600; width: 180px; color: #555; }
  code { background: #f0f3f8; padding: 0.15rem 0.4rem; border-radius: 3px; font-size: 0.85rem; }
  .td-table { width: 100%; border-collapse: collapse; font-size: 0.85rem; }
  .td-table th { text-align: left; padding: 0.4rem 0.5rem; background: #f8f9fb;
    border-bottom: 2px solid #e0e4ea; color: #555; font-weight: 600; }
  .td-table td { padding: 0.35rem 0.5rem; border-bottom: 1px solid #f0f0f0; }
  .td-table tr:hover { background: #fafbfd; }
  .td-indent { padding-left: 1.8rem !important; color: #666; }
  .pipe-table { width: 100%; border-collapse: collapse; font-size: 0.85rem; }
  .pipe-table th { text-align: left; padding: 0.4rem 0.5rem; background: #f8f9fb;
    border-bottom: 2px solid #e0e4ea; color: #555; font-weight: 600; }
  .pipe-table td { padding: 0.4rem 0.5rem; border-bottom: 1px solid #f0f0f0;
    vertical-align: top; }
  .pipe-table tr:hover { background: #fafbfd; }
  .args-cell { max-width: 350px; overflow-x: auto; white-space: nowrap; }
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
  .badge-other { background: #f5f5f5; color: #616161; }
  .mermaid-container { background: #fff; border-radius: 8px; padding: 1.5rem;
    box-shadow: 0 1px 3px rgba(0,0,0,0.08); border: 1px solid #e8ecf1;
    margin-bottom: 1.2rem; overflow-x: auto; }
  .mermaid { text-align: center; }
</style>
</head>
<body>

<div class="layout">
  <div class="left-panel">
  <h1>',title,'</h1>
    ', ref_html, '
  </div>
  <div class="right-panel">
    ', td_html, '
    ', tp_html, '
  </div>
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
  </div>
</div>

<script>mermaid.initialize({ startOnLoad: true, theme: "neutral", flowchart: { useMaxWidth: true, htmlLabels: true } });</script>
</body>
</html>'
  )
  
  html
}


# ══════════════════════════════════════════════════════════════════
# Internal: build Mermaid diagram from pipeline steps
# ══════════════════════════════════════════════════════════════════

build_mermaid <- function(vf_info) {
  
  lines <- "graph TB"
  
  # Start node
  lines <- c(lines, '  START(["<b>y</b><br/>Variable of interest"])')
  
  # Collect variance nodes to connect them to a final summation
  variance_nodes <- character(0)
  prev_transform_node <- "START"
  
  for (step in vf_info) {
    sid <- paste0("S", step$step)
    fn_esc <- mermaid_esc(step$function_name)
    label_esc <- mermaid_esc(step$label)
    
    if (step$category %in% c("variance_nr", "variance_sampling")) {
      # Variance nodes: rounded shape, branch off from the current transform
      lines <- c(lines, paste0(
        '  ', sid, '(["<b>', fn_esc, '()</b><br/>', label_esc, '"])'
      ))
      lines <- c(lines, paste0('  ', prev_transform_node, ' -.-> ', sid))
      variance_nodes <- c(variance_nodes, sid)
      
    } else if (step$category == "summation") {
      # Reduce node: all variance nodes feed in
      lines <- c(lines, paste0(
        '  ', sid, '["<b>', fn_esc, '</b><br/>', label_esc, '"]'
      ))
      for (vn in variance_nodes) {
        lines <- c(lines, paste0('  ', vn, ' --> ', sid))
      }
      lines <- c(lines, paste0('  ', sid, ' --> RESULT(["<b>V(Y&#770;)</b>"])'))
      
    } else {
      # Transform nodes: rectangular, part of the main chain
      lines <- c(lines, paste0(
        '  ', sid, '["<b>', fn_esc, '()</b><br/>', label_esc, '"]'
      ))
      lines <- c(lines, paste0('  ', prev_transform_node, ' --> ', sid))
      prev_transform_node <- sid
    }
  }
  
  # Style nodes by category
  for (step in vf_info) {
    sid <- paste0("S", step$step)
    style <- switch(step$category,
                    "calibration"       = "fill:#e8f5e9,stroke:#4caf50,color:#1b5e20",
                    "variance_nr"       = "fill:#fff3e0,stroke:#ff9800,color:#e65100",
                    "variance_sampling" = "fill:#e3f2fd,stroke:#2196f3,color:#0d47a1",
                    "transition"        = "fill:#f3e5f5,stroke:#9c27b0,color:#4a148c",
                    "aggregation"       = "fill:#fce4ec,stroke:#e91e63,color:#880e4f",
                    "summation"         = "fill:#eceff1,stroke:#607d8b,color:#263238",
                    NULL
    )
    if (!is.null(style)) {
      lines <- c(lines, paste0('  style ', sid, ' ', style))
    }
  }
  
  lines <- c(lines, '  style START fill:#e8eaf6,stroke:#3f51b5,color:#1a237e')
  lines <- c(lines, '  style RESULT fill:#e8eaf6,stroke:#3f51b5,color:#1a237e')
  
  paste(lines, collapse = "\n")
}


# ══════════════════════════════════════════════════════════════════
# Internal: build technical data HTML table
# ══════════════════════════════════════════════════════════════════

build_td_table <- function(td_info, td_name = "Technical Data") {
  
  rows <- character(0)
  
  for (td in td_info) {
    if (length(td$type) == 0 || td$type == "NULL") {
      rows <- c(rows, paste0(
        '<tr><td><b>', htmlesc(td$name), '</b></td>',
        '<td>NULL</td><td>-</td><td>Not provided</td></tr>'
      ))
      next
    }
    
    if (td$type == "list") {
      # Parent row
      rows <- c(rows, paste0(
        '<tr><td><b>', htmlesc(td$name), '</b></td>',
        '<td>list</td><td>', td$n_elements, ' elements</td>',
        '<td><code>', htmlesc(paste(td$element_names, collapse = ", ")), '</code></td></tr>'
      ))
      # Child rows
      for (el in td$elements) {
        preview <- if (!is.null(el$first_values)) {
          paste0('<code>', htmlesc(el$first_values), '</code>')
        } else { "" }
        rows <- c(rows, paste0(
          '<tr><td class="td-indent">\u2514 ', htmlesc(el$name), '</td>',
          '<td><code>', htmlesc(el$type), '</code></td>',
          '<td>', htmlesc(el$dimension), '</td>',
          '<td>', preview, '</td></tr>'
        ))
      }
    } else {
      dim_str <- if (!is.null(td$dimension)) td$dimension else ""
      preview <- if (!is.null(td$first_values)) {
        paste0('<code>', htmlesc(td$first_values), '</code>')
      } else { "" }
      rows <- c(rows, paste0(
        '<tr><td><b>', htmlesc(td$name), '</b></td>',
        '<td><code>', htmlesc(td$type), '</code></td>',
        '<td>', htmlesc(dim_str), '</td>',
        '<td>', preview, '</td></tr>'
      ))
    }
  }
  
  paste0(
    '<div class="card">',
    '<h2>',td_name,'</h2>',
    '<table class="td-table">',
    '<thead><tr><th>Name</th><th>Type</th><th>Dimension</th><th>Preview</th></tr></thead>',
    '<tbody>', paste(rows, collapse = "\n"), '</tbody>',
    '</table>',
    '</div>'
  )
}


# ══════════════════════════════════════════════════════════════════
# Internal: HTML/Mermaid escaping utilities
# ══════════════════════════════════════════════════════════════════

htmlesc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

mermaid_esc <- function(x) {
  # Mermaid is sensitive to special chars in labels
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x <- gsub("(", "&#40;", x, fixed = TRUE)
  x <- gsub(")", "&#41;", x, fixed = TRUE)
  x
}
