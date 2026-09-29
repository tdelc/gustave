#' Call a gustave variance wrapper for a single indicator and optional grouping
#'
#' This helper builds and evaluates the wrapper call programmatically,
#' handling single and crossed breakdowns.
#'
#' @param df        Data frame passed to the wrapper.
#' @param indic     Character, name of the variable of interest.
#' @param var_by    Character vector of grouping variable(s), or NULL.
#' @param wrapper   Character, name of the variance wrapper function.
#' @param fun_name  Character, name of the statistic wrapper ("total", "mean", etc.).
#'
#' @return A data frame with variance estimation results and an `indic` column.
#' @export
call_gustave <- function(df, indic, var_by = NULL, wrapper, fun_name) {
  
  # Build the command string depending on the number of grouping variables
  if (is.null(var_by)) {
    command <- sprintf("%s(df, %s(%s))", wrapper, fun_name, indic)
  } else if (length(var_by) == 1) {
    command <- sprintf("%s(df, %s(%s), by = %s)", wrapper, fun_name, indic, var_by)
  } else {
    # For crossed breakdowns: create a temporary interaction variable
    df <- df %>% mutate(.grp = interaction(!!!syms(var_by), sep = "|||"))
    command <- sprintf("%s(df, %s(%s), by = .grp)", wrapper, fun_name, indic)
  }
  
  result <- eval(parse(text = command)) %>% mutate(indic = indic)
  
  # Split the interaction column back into individual breakdown variables
  if (length(var_by) >= 1 && "by" %in% colnames(result)) {
    result <- result %>% tidyr::separate(by, into = var_by, sep = "\\|\\|\\|")
  }
  
  result
}


#' Batch variance estimation across indicators, groups, and domains
#'
#' Runs a gustave variance wrapper over all combinations of indicators,
#' grouping variables, and domain variables. Returns a single tidy data frame
#' with a human-readable `pop` column describing the subpopulation.
#'
#' @param df           Data frame passed to the wrapper.
#' @param wrapper      Character, name of the variance wrapper function.
#' @param fun_name     Character, name of the statistic wrapper ("total", "mean", etc.).
#' @param vars_indics  Character vector of variable names to estimate.
#' @param vars_groups  Character vector of grouping variables (crossed), or NULL.
#' @param vars_domain  Character vector of domain variables (not crossed with each other), or NULL.
#'
#' @return A tidy data frame with columns: step, indic, pop, n, est, variance,
#'         std, cv, lower, upper, plus the original group/domain columns.
#' @export
batch_gustave <- function(df,
                          wrapper,
                          fun_name,
                          vars_indics,
                          vars_groups = NULL,
                          vars_domain = NULL) {
  
  # --- Input validation ---
  missing_indics <- vars_indics[!vars_indics %in% colnames(df)]
  missing_groups <- vars_groups[!vars_groups %in% colnames(df)]
  missing_domain <- vars_domain[!vars_domain %in% colnames(df)]
  
  if (length(missing_indics) > 0)
    stop("Indicator variables not found: ", paste(missing_indics, collapse = ", "))
  if (length(missing_groups) > 0)
    stop("Group variables not found: ", paste(missing_groups, collapse = ", "))
  if (length(missing_domain) > 0)
    stop("Domain variables not found: ", paste(missing_domain, collapse = ", "))
  
  # --- Step 1: Overall estimates (no breakdown) ---
  message("Step 1: Overall estimates (no breakdown)")
  step_1 <- purrr::map(vars_indics, function(v) {
    df %>% call_gustave(indic = v, wrapper = wrapper, fun_name = fun_name)
  }) %>% bind_rows() %>% mutate(step = "overall")
  
  step_1[, vars_groups] <- NA
  step_1[, vars_domain] <- NA
  step_1 <- step_1 %>%
    select(step, indic, !!!syms(vars_domain), !!!syms(vars_groups),
           n, est, variance, std, cv, lower, upper, call)
  
  results <- list(step1 = step_1)
  
  # --- Coerce group/domain variables to character, recode NA ---
  for (var in c(vars_groups, vars_domain)) {
    df[[var]] <- as.character(df[[var]])
    df[[var]][is.na(df[[var]])] <- "NA"
  }
  
  # --- Step 2: Marginal estimates (each group/domain variable separately) ---
  message("Step 2: Marginal estimates (each group and domain variable)")
  all_by_vars <- c(vars_groups, vars_domain)
  results[["step2"]] <- purrr::map(all_by_vars, function(b) {
    purrr::map(vars_indics, function(v) {
      df %>% call_gustave(v, b, wrapper = wrapper, fun_name = fun_name)
    }) %>% bind_rows()
  }) %>% bind_rows() %>% mutate(step = "marginal")
  
  # --- Step 3: Crossed group estimates (all group variables combined) ---
  if (!is.null(vars_groups)) {
    message("Step 3: Crossed group estimates")
    results[["step3"]] <- purrr::map(vars_indics, function(v) {
      df %>% call_gustave(v, vars_groups, wrapper = wrapper, fun_name = fun_name)
    }) %>% bind_rows() %>% mutate(step = "crossed_groups")
    
    # --- Step 4: Crossed groups within each domain ---
    if (!is.null(vars_domain)) {
      message("Step 4: Crossed groups by domain")
      results[["step4"]] <- purrr::map(vars_domain, function(d) {
        purrr::map(vars_indics, function(v) {
          df %>% call_gustave(v, c(vars_groups, d),
                              wrapper = wrapper, fun_name = fun_name)
        }) %>% bind_rows()
      }) %>% bind_rows() %>% mutate(step = "crossed_groups_by_domain")
    }
  }
  
  # --- Assemble and build the human-readable population label ---
  results <- purrr::map_df(results, rbind) %>% mutate(pop = "")
  
  for (var in c(vars_domain, vars_groups)) {
    results <- results %>%
      mutate(pop = if_else(
        is.na(!!sym(var)), pop,
        paste0(pop, " ; ", var, " = ", !!sym(var))
      ))
  }
  
  results <- results %>%
    mutate(
      pop = stringr::str_remove(pop, "^ *; *"),
      pop = stringr::str_remove(pop, " *; *$"),
      pop = stringr::str_replace(pop, "^$", "ALL")
    )
  
  results %>%
    remove_rownames() %>%
    select(step, indic, pop, n, est, variance, std, cv, lower, upper,
           all_of(c(vars_groups, vars_domain)), call)
}
