


#' Quickly perform a variance estimation for multistage designs
#'
#' @description \code{qvar_ms} extends \code{\link{qvar}} to multistage
#'   sampling designs, that is: \itemize{\item stratified simple random
#'   sampling at each stage (units of stage k being drawn within the units
#'   of stage k - 1) \item non-response correction (if any) through
#'   reweighting, possibly at each stage \item calibration (if any) at one
#'   given stage}
#'
#' Used with \code{define = TRUE}, it defines a variance wrapper to be
#' applied to the survey dataset of the last stage (see
#' \code{\link{define_variance_wrapper}}).
#'
#' @param data A \code{data.frame} (single-stage case, as in
#'   \code{\link{qvar}}) or a list of \code{data.frame} of length
#'   \code{sampling_stages}, ordered from the highest stage (e.g. primary
#'   sampling units) to the lowest one (final units). Each file should
#'   contain all the units sampled at the corresponding stage, including
#'   the out-of-scope and non-responding units.
#' @param ...,where,by,alpha,display,define,envir See \code{\link{qvar}}.
#' @param sampling_stages A numeric vector of length 1, the number of
#'   sampling stages (\code{1} by default). It should match the length of
#'   \code{data} when \code{data} is a list.
#' @param id A character vector (or list) of length \code{sampling_stages},
#'   the name of the identification variable at each stage.
#' @param parent_id A character vector of length \code{sampling_stages - 1}
#'   (or a list of length \code{sampling_stages} whose first element is
#'   \code{NULL}), the name, in each file from stage 2 onwards, of the
#'   variable identifying the unit of the stage above. If \code{NULL}
#'   (default), the id variable of the stage above is assumed to be
#'   present in the file.
#' @param dissemination_dummy,dissemination_weight,scope_dummy Character
#'   vectors of length 1 referring to variables of the LAST stage file
#'   (see \code{\link{qvar}}).
#' @param sampling_weight A character vector (or list) of length
#'   \code{sampling_stages}, the name of the CONDITIONAL sampling weight
#'   at each stage (inverse of the inclusion probability given the stages
#'   above).
#' @param strata A list of length \code{sampling_stages} (\code{NULL}
#'   elements meaning no stratification at the corresponding stage). From
#'   stage 2 onwards, the effective strata are the crossing of the parent
#'   unit and of the declared strata variable.
#' @param nrc_weight,response_dummy,nrc_dummy Lists of length
#'   \code{sampling_stages} (\code{NULL} elements meaning no non-response
#'   correction at the corresponding stage). A single variable name is
#'   assumed to refer to the LAST stage. \code{nrc_weight} is the
#'   CONDITIONAL weight after non-response correction at the stage.
#' @param calibration_stage A numeric vector of length 1, the stage at
#'   which calibration was performed (\code{sampling_stages} by default).
#' @param calibration_weight A character vector of length 1, the name, in
#'   the \code{calibration_stage} file, of the CUMULATED weight after
#'   calibration (the weight actually used in the calibration process).
#' @param calibration_dummy,calibration_var Character vectors referring to
#'   variables of the \code{calibration_stage} file (see
#'   \code{\link{qvar}}).
#'
#' @seealso \code{\link{qvar}}, \code{\link{define_variance_wrapper}}
#' 
#' @details \code{qvar_ms} extends the variance estimation methodology of
#'   \code{\link{qvar}} to multistage sampling designs. Like \code{qvar}, it
#'   performs not only technical but also methodological checks in order to
#'   ensure that the standard variance estimation methodology does apply at
#'   each stage (e.g. equal probability of inclusion within strata, number of
#'   units per stratum, consistency of the weights across stages).
#'
#'   \strong{Conventions}
#'
#'   \itemize{
#'   \item Stage 1 is the highest sampling stage (e.g. primary sampling
#'   units); stage \code{L = sampling_stages} is the final stage (surveyed
#'   units).
#'   \item \code{data} is a list of \code{L} \code{data.frame} (or a single
#'   \code{data.frame} when \code{L = 1}), ordered from stage 1 to stage
#'   \code{L}. Each file should contain all the units sampled at the
#'   corresponding stage, including the out-of-scope and non-responding
#'   units.
#'   \item Each \code{data.frame} from stage 2 onwards contains a linking
#'   variable towards the stage above (\code{parent_id} argument; by
#'   default, the \code{id} variable of stage \code{k - 1}).
#'   \item \code{sampling_weight} and \code{nrc_weight} are CONDITIONAL
#'   weights at the considered stage, that is the inverse of the inclusion
#'   (resp. response) probability given the stages above (Poulpe
#'   convention).
#'   \item \code{calibration_weight} is the CUMULATED weight at the
#'   calibration stage, that is the weight actually used in the calibration
#'   equations.
#'   \item \code{dissemination_dummy}, \code{dissemination_weight} and
#'   \code{scope_dummy} refer to the LAST stage only.
#'   \item Non-response correction may be declared at each stage
#'   (list-valued arguments); a single variable name is assigned to the
#'   last stage.
#'   \item Calibration is declared at one single stage
#'   (\code{calibration_stage} argument, last stage by default).
#'   }
#'
#'   \strong{Methodology}
#'
#'   The variance is estimated as the sum of the conditional variances of
#'   each stage, following the Rao (1975) formula as implemented in the SAS
#'   Poulpe macros and natively supported by the gustave package: at stage
#'   \code{k}, the variance contributions (\code{\link{var_srs}},
#'   \code{\link{var_pois}}) are LINEARLY weighted by the cumulated
#'   corrected weight of the stages above, through the \code{w} argument of
#'   the variance functions ("row weight used at the final summation step
#'   [...] applying the Rao (1975) formula", see \code{\link{varDT}} and
#'   \code{\link{var_pois}}). The first stage variance is computed on the
#'   estimated totals per primary sampling unit. This decomposition relies
#'   on the assumption that the sampling design is invariant and
#'   independent across stages.
#'
#' @examples \dontrun{
#' # Two-stage design: municipalities (PSU) then dwellings
#' precision_survey <- qvar_ms(
#'
#'   # One file per stage, from the highest to the lowest
#'   data = list(psu_sample, dwelling_sample),
#'   sampling_stages = 2,
#'   id = c("psu_id", "dwelling_id"),
#'   parent_id = "psu_id",
#'
#'   # Dissemination and scope (last stage only)
#'   dissemination_dummy = "dissemination",
#'   dissemination_weight = "w_final",
#'   scope_dummy = "scope",
#'
#'   # Sampling design (conditional weights)
#'   sampling_weight = c("w_psu", "w_dwelling"),
#'   strata = list("region", NULL),
#'
#'   # Non-response correction (last stage by default)
#'   nrc_weight = "w_nrc",
#'   response_dummy = "resp",
#'
#'   # Calibration (last stage by default, cumulated weight)
#'   calibration_weight = "w_final",
#'   calibration_var = c("age_group", "sex"),
#'
#'   define = TRUE
#' )
#' precision_survey(dwelling_survey, mean(income))
#' }
#'
#' @export

qvar_ms <- function(data, ..., by = NULL, where = NULL,
                    alpha = 0.05, display = TRUE,
                    sampling_stages = 1,
                    id, parent_id = NULL,
                    dissemination_dummy, dissemination_weight,
                    sampling_weight, strata = NULL,
                    scope_dummy = NULL,
                    nrc_weight = NULL, response_dummy = NULL, nrc_dummy = NULL,
                    calibration_stage = NULL,
                    calibration_weight = NULL, calibration_dummy = NULL, calibration_var = NULL,
                    define = FALSE, envir = parent.frame()
                    # TODO: Add objects_to_include
){
  
  # Step 1: Define the variance wrapper
  call <- as.list(match.call())[-1]
  call$envir <- envir
  qvar_variance_wrapper <- do.call(
    define_qvar_ms_variance_wrapper,
    call[names(call) %in% names(formals(define_qvar_ms_variance_wrapper))]
  )
  
  # Step 2: Export the variance wrapper
  if(define){
    note("As define = TRUE, a ready-to-use variance wrapper is (invisibly) returned.")
    return(invisible(qvar_variance_wrapper))
  }
  
  # Step 3: Estimate variance
  data_last <- if(is.data.frame(data)) data else data[[length(data)]]
  id_last <- if(is.list(id)) id[[length(id)]] else id[length(id)]
  qvar_data <- data_last[data_last[, id_last] %in% environment(qvar_variance_wrapper)$reference_id, ]
  call$data <- substitute(qvar_data)
  call$envir <- environment()
  do.call(
    qvar_variance_wrapper,
    call[names(call) == "" | names(call) %in% names(formals(qvar_variance_wrapper))]
  )
  
}



# Unexported (and undocumented) functions

define_qvar_ms_variance_wrapper <- function(data, sampling_stages = 1,
                                            id, parent_id = NULL,
                                            dissemination_dummy, dissemination_weight,
                                            sampling_weight, strata = NULL,
                                            scope_dummy = NULL,
                                            nrc_weight = NULL, response_dummy = NULL, nrc_dummy = NULL,
                                            calibration_stage = NULL,
                                            calibration_weight = NULL, calibration_var = NULL, calibration_dummy = NULL,
                                            envir = parent.frame()
){
  
  # Step 1: Control arguments consistency and display the welcome message ----
  
  # Step 1.1: Arguments consistency (missing arguments)
  is_missing <- c(
    data = missing(data),
    id = missing(id),
    dissemination_dummy = missing(dissemination_dummy),
    dissemination_weight = missing(dissemination_weight),
    sampling_weight = missing(sampling_weight)
  )
  if(any(is_missing)) stop(
    "The following arguments are missing: ",
    paste(names(which(is_missing)), collapse = ", "), "."
  )
  
  deparse_data <- deparse(substitute(data))
  data <- eval(substitute(data), envir = envir)
  if(is.data.frame(data)) data <- list(data)
  if(!is.list(data) || !all(sapply(data, is.data.frame))) stop(
    "data argument must refer to a data.frame or to a list of data.frame ",
    "(one per sampling stage, from the highest to the lowest one)."
  )
  if(missing(sampling_stages)){
    sampling_stages <- length(data)
    if(sampling_stages > 1) note(
      "As the sampling_stages argument is not provided, it is inferred from the length of data: ",
      sampling_stages, " sampling stages."
    )
  }
  if(length(data) != sampling_stages) stop(
    "The data argument is a list of ", length(data), " data.frame(s) but sampling_stages = ",
    sampling_stages, ". Both should match."
  )
  
  as_stage_list <- function(x, arg_name, default_last = FALSE){
    if(is.null(x)) return(vector("list", sampling_stages))
    if(is.list(x)){
      if(length(x) != sampling_stages) stop(
        "The ", arg_name, " argument is a list of length ", length(x),
        " but sampling_stages = ", sampling_stages,
        ". It should contain one element per sampling stage (possibly NULL)."
      )
      return(x)
    }
    if(is.character(x)){
      if(length(x) == sampling_stages) return(as.list(x))
      if(length(x) == 1 && default_last){
        out <- vector("list", sampling_stages)
        out[[sampling_stages]] <- x
        if(sampling_stages > 1) note(
          "The ", arg_name, " argument is of length 1: it is assumed to refer to the last sampling stage (stage ",
          sampling_stages, ")."
        )
        return(out)
      }
    }
    stop(
      "The ", arg_name, " argument should be a character vector or a list of length ",
      sampling_stages, " (one element per sampling stage)."
    )
  }
  
  id <- as_stage_list(id, "id")
  sampling_weight <- as_stage_list(sampling_weight, "sampling_weight")
  strata <- as_stage_list(strata, "strata")
  nrc_weight <- as_stage_list(nrc_weight, "nrc_weight", default_last = TRUE)
  response_dummy <- as_stage_list(response_dummy, "response_dummy", default_last = TRUE)
  nrc_dummy <- as_stage_list(nrc_dummy, "nrc_dummy", default_last = TRUE)
  
  if(is.null(parent_id)){
    parent_id <- c(list(NULL), id[-sampling_stages])
    if(sampling_stages > 1) note(
      "As the parent_id argument is NULL, each file from stage 2 onwards is assumed to contain ",
      "the identification variable of the stage above (",
      paste(unlist(id[-sampling_stages]), collapse = ", "), ")."
    )
  }else{
    if(!is.list(parent_id)) parent_id <- as.list(parent_id)
    if(length(parent_id) == sampling_stages - 1) parent_id <- c(list(NULL), parent_id)
    if(length(parent_id) != sampling_stages || !is.null(parent_id[[1]])) stop(
      "The parent_id argument should refer to one variable name per sampling stage but the first one ",
      "(either a character vector of length sampling_stages - 1 or a list of length sampling_stages ",
      "whose first element is NULL)."
    )
  }
  
  ## <<< ML: niveau du calage (dernier niveau par defaut)
  if(is.null(calibration_stage)) calibration_stage <- sampling_stages
  if(!is.numeric(calibration_stage) || length(calibration_stage) != 1 ||
     !(calibration_stage %in% seq_len(sampling_stages))) stop(
       "The calibration_stage argument should be a single integer between 1 and sampling_stages (",
       sampling_stages, ")."
     )
  
  # Step 1.3: Arguments consistency (inconsistent arguments)
  inconsistency <- list(
    nrc_weight_but_no_response_dummy = sapply(seq_len(sampling_stages), function(k)
      !is.null(nrc_weight[[k]]) && is.null(response_dummy[[k]])),
    resp_or_nrc_dummy_but_no_nrc_weight = sapply(seq_len(sampling_stages), function(k)
      is.null(nrc_weight[[k]]) && (!is.null(response_dummy[[k]]) || !is.null(nrc_dummy[[k]]))),
    calibration_weight_but_no_calibration_var = !is.null(calibration_weight) && is.null(calibration_var),
    calibration_or_calibration_var_but_no_calibration_weight = is.null(calibration_weight) && (!is.null(calibration_dummy) || !is.null(calibration_var))
  )
  if(any(unlist(inconsistency))) stop(
    "Some arguments are inconsistent:",
    if(any(inconsistency$nrc_weight_but_no_response_dummy)) paste0(
      "\n  - at stage(s) ", paste(which(inconsistency$nrc_weight_but_no_response_dummy), collapse = ", "),
      ": weights after non-response correction are provided (nrc_weight argument) but no variable indicating responding units (response_dummy argument)"
    ) else "",
    if(any(inconsistency$resp_or_nrc_dummy_but_no_nrc_weight)) paste0(
      "\n  - at stage(s) ", paste(which(inconsistency$resp_or_nrc_dummy_but_no_nrc_weight), collapse = ", "),
      ": a variable indicating responding units and/or a variable indicating the units taking part in the non-response correction process are provided (response_dummy and nrc_dummy argument) but no weights after non-response correction (nrc_weight argument)."
    ) else "",
    if(inconsistency$calibration_weight_but_no_calibration_var)
      "\n  - calibrated weights are provided (calibration_weight argument) but no calibration variables (calibration_var argument)" else "",
    if(inconsistency$calibration_or_calibration_var_but_no_calibration_weight)
      "\n  - a variable indicating the units taking part in a calibration process and/or calibration variables are provided (calibration_dummy and calibration_var arguments) but no calibrated weights (calibration_weight argument)" else ""
  )
  
  # Step 1.4: Welcome message
  message(
    "Survey variance estimation with the gustave package",
    "\n\nThe following features are taken into account:",
    if(sampling_stages > 1) paste0("\n  - ", sampling_stages, "-stage sampling design") else "",
    paste0(sapply(seq_len(sampling_stages), function(k) paste0(
      if(sampling_stages > 1) paste0("\n  - stage ", k, ": ") else "\n  - ",
      if(!is.null(strata[[k]])) "stratified simple random sampling" else
        "simple random sampling WITHOUT stratification",
      if(k > 1) " (within the units of the stage above)" else "",
      if(!is.null(nrc_weight[[k]])) ", with non-response correction through reweighting" else ""
    )), collapse = ""),
    if(!is.null(scope_dummy)) "\n  - out-of-scope units (last stage)" else "",
    if(!is.null(calibration_weight)) paste0(
      "\n  - calibration on margins",
      if(sampling_stages > 1) paste0(" at stage ", calibration_stage) else ""
    ) else "",
    "\n"
  )
  
  # Step 2: Control that arguments do exist and retrieve their value ----
  
  # Step 2.1: Evaluation of all arguments
  ## <<< ML: data a deja ete evalue au Step 1.2 ; arg est reconstruit a partir
  ## des arguments normalises (et non plus via match.call) pour que les listes
  ## par niveau soient directement indexables (arg$strata[[k]], etc.)
  arg <- list(
    id = id, parent_id = parent_id,
    dissemination_dummy = dissemination_dummy, dissemination_weight = dissemination_weight,
    sampling_weight = sampling_weight, strata = strata, scope_dummy = scope_dummy,
    nrc_weight = nrc_weight, response_dummy = response_dummy, nrc_dummy = nrc_dummy,
    calibration_weight = calibration_weight, calibration_dummy = calibration_dummy,
    calibration_var = calibration_var
  )
  
  # Step 2.2: Expected length
  should_be_single_variable_name_by_stage <- c(
    "id", "parent_id", "sampling_weight", "strata",
    "nrc_weight", "response_dummy", "nrc_dummy"
  )
  should_be_single_variable_name <- intersect(c(
    "dissemination_dummy", "dissemination_weight", "scope_dummy",
    "calibration_weight", "calibration_dummy"
  ), names(arg)[!sapply(arg, is.null)])
  should_be_variable_name_vector <- intersect(c("calibration_var"), names(arg)[!sapply(arg, is.null)])
  
  # Step 2.3: Check whether arguments are character vectors and
  # have the expected length
  is_single_variable_name_by_stage <- sapply(
    should_be_single_variable_name_by_stage,
    function(param) all(sapply(
      arg[[param]],
      function(x) is.null(x) || is_variable_name(x, max_length = 1)
    ))
  )
  if(any(!is_single_variable_name_by_stage)) stop(
    "The following arguments do not refer to one variable name (character vector of length 1) per sampling stage: ",
    paste(names(is_single_variable_name_by_stage)[!is_single_variable_name_by_stage], collapse = ", ")
  )
  is_single_variable_name <- sapply(
    arg[should_be_single_variable_name],
    function(x) is.null(x) || is_variable_name(x, max_length = 1)
  )
  if(any(!is_single_variable_name)) stop(
    "The following arguments do not refer to a variable name (character vector of length 1): ",
    names(is_single_variable_name)[!is_single_variable_name]
  )
  is_variable_name_vector <- sapply(
    arg[should_be_variable_name_vector],
    function(x) is.null(x) || is_variable_name(x, max_length = Inf)
  )
  if(any(!is_variable_name_vector)) stop(
    "The following arguments do not refer to a vector of variable names: ",
    names(is_variable_name_vector)[!is_variable_name_vector]
  )
  
  # Step 2.4: Check the presence of the variables in data
  is_not_in_data <- unlist(lapply(seq_len(sampling_stages), function(k){
    param_at_stage <- c(
      stats::setNames(
        lapply(should_be_single_variable_name_by_stage, function(param) arg[[param]][[k]]),
        should_be_single_variable_name_by_stage
      ),
      if(k == sampling_stages) list(
        dissemination_dummy = arg$dissemination_dummy,
        dissemination_weight = arg$dissemination_weight,
        scope_dummy = arg$scope_dummy
      ) else NULL,
      if(k == calibration_stage) list(
        calibration_weight = arg$calibration_weight,
        calibration_dummy = arg$calibration_dummy,
        calibration_var = arg$calibration_var
      ) else NULL
    )
    lapply(names(param_at_stage), function(param){
      if(is.null(param_at_stage[[param]])) return(NULL)
      tmp <- variable_not_in_data(var = param_at_stage[[param]], data = data[[k]])
      if(is.null(tmp)) return(NULL)
      paste0("\n  - stage ", k, ", ", param, " argument: ", paste0(tmp, collapse = " "))
    })
  }))
  if(length(is_not_in_data) > 0) stop(
    "Some variables do not exist in ", deparse_data, ": ",
    paste0(is_not_in_data, collapse = "")
  )
  
  # Step 2.5: Retrieve the value of the arguments
  lev <- lapply(seq_len(sampling_stages), function(k){
    d <- data[[k]]
    d <- d[order(d[[arg$id[[k]]]]), ]
    id_k <- as.character(d[[arg$id[[k]]]])
    get_var <- function(param){
      if(is.null(arg[[param]][[k]])) return(NULL)
      stats::setNames(d[[arg[[param]][[k]]]], id_k)
    }
    out <- list(
      id = stats::setNames(id_k, id_k),
      parent = if(k > 1) stats::setNames(as.character(d[[arg$parent_id[[k]]]]), id_k) else NULL,
      sampling_weight = get_var("sampling_weight"),
      strata = get_var("strata"),
      nrc_weight = get_var("nrc_weight"),
      response_dummy = get_var("response_dummy"),
      nrc_dummy = get_var("nrc_dummy")
    )
    if(k == sampling_stages){
      out$dissemination_dummy <- stats::setNames(d[[arg$dissemination_dummy]], id_k)
      out$dissemination_weight <- stats::setNames(d[[arg$dissemination_weight]], id_k)
      out$scope_dummy <- if(is.null(arg$scope_dummy)) NULL else stats::setNames(d[[arg$scope_dummy]], id_k)
    }
    if(k == calibration_stage && !is.null(arg$calibration_weight)){
      out$calibration_weight <- stats::setNames(d[[arg$calibration_weight]], id_k)
      out$calibration_dummy <- if(is.null(arg$calibration_dummy)) NULL else stats::setNames(d[[arg$calibration_dummy]], id_k)
      tmp <- d[arg$calibration_var]
      row.names(tmp) <- id_k
      out$calibration_var <- tmp
    }
    out
  })
  
  
  # Step 3: Control arguments value ----
  
  calibration_var_quanti <- NULL
  calibration_var_quali <- NULL
  
  for(k in seq_len(sampling_stages)){
    
    lv <- lev[[k]]
    lab <- if(sampling_stages > 1) paste0(" at stage ", k) else ""
    
    # id
    if(anyNA(lv$id))
      stop("The id variable (", arg$id[[k]], ")", lab, " should not contain any missing (NA) values.")
    if(any(duplicated(lv$id)))
      stop("The id variable (", arg$id[[k]], ")", lab, " should not contain any duplicated values.")
    
    if(k > 1){
      if(anyNA(lv$parent))
        stop("The parent id variable (", arg$parent_id[[k]], ")", lab, " should not contain any missing (NA) values.")
      parent_not_found <- unique(lv$parent[!(lv$parent %in% lev[[k - 1]]$id)])
      if(length(parent_not_found) > 0) stop(
        "The following parent units (", arg$parent_id[[k]], ")", lab,
        " do not exist in the stage ", k - 1, " file: ",
        display_only_n_first(parent_not_found), "."
      )
      if(!is.null(arg$nrc_weight[[k - 1]])){
        parent_nonresponding <- unique(lv$parent[!(lev[[k - 1]]$response_dummy[lv$parent] %in% TRUE)])
        if(length(parent_nonresponding) > 0) stop(
          "The following units of stage ", k - 1, " are non-responding (", arg$response_dummy[[k - 1]],
          ") but have units of stage ", k, " attached to them: ",
          display_only_n_first(parent_nonresponding), "."
        )
      }
    }
    
    # sampling_weight
    if(!is.numeric(lv$sampling_weight))
      stop("The sampling weights (", arg$sampling_weight[[k]], ")", lab, " should be numeric.")
    if(anyNA(lv$sampling_weight))
      stop("The sampling weights (", arg$sampling_weight[[k]], ")", lab, " should not contain any missing (NA) values.")
    
    # strata
    if(is.null(lv$strata)) lv$strata <- stats::setNames(factor(rep("1", length(lv$id))), lv$id)
    if(!is.null(lv$strata)){
      if(is.character(lv$strata)){
        note("The strata variable (", arg$strata[[k]], ")", lab, " is of type character. It is automatically coerced to factor.")
        lv$strata <- stats::setNames(factor(lv$strata), lv$id)
      }
      if(!is.factor(lv$strata))
        stop("The strata variable (", arg$strata[[k]], ")", lab, " should be of type factor or character.")
      if(anyNA(lv$strata))
        stop("The strata variable (", arg$strata[[k]], ")", lab, " should not contain any missing (NA) values.")
    }
    
    # scope_dummy
    if(k < sampling_stages || is.null(lv$scope_dummy)){
      lv$scope_dummy <- stats::setNames(rep(TRUE, length(lv$id)), lv$id)
    }else{
      if(is.numeric(lv$scope_dummy)){
        note("The scope dummy variable (", arg$scope_dummy, ") is of type numeric. It is automatically coerced to logical.")
        lv$scope_dummy <- stats::setNames(as.logical(lv$scope_dummy), lv$id)
      }
      if(!is.logical(lv$scope_dummy))
        stop("The scope dummy variable (", arg$scope_dummy, ") should be of type logical or numeric.")
      if(anyNA(lv$scope_dummy))
        stop("The scope dummy variable (", arg$scope_dummy, ") should not contain any missing (NA) values.")
    }
    
    # response_dummy
    if(is.null(lv$response_dummy)) lv$response_dummy <- lv$scope_dummy else{
      if(is.numeric(lv$response_dummy)){
        note("The response dummy variable (", arg$response_dummy[[k]], ")", lab, " is of type numeric. It is automatically coerced to logical.")
        lv$response_dummy <- stats::setNames(as.logical(lv$response_dummy), lv$id)
      }
      if(!is.logical(lv$response_dummy))
        stop("The response dummy variable (", arg$response_dummy[[k]], ")", lab, " should be of type logical or numeric.")
      if(anyNA(lv$response_dummy))
        stop("The response dummy variable (", arg$response_dummy[[k]], ")", lab, " should not contain any missing (NA) values.")
    }
    
    # nrc_dummy
    if(is.null(lv$nrc_dummy)) lv$nrc_dummy <- lv$scope_dummy else{
      if(is.numeric(lv$nrc_dummy)){
        note("The non-reponse correction dummy variable (", arg$nrc_dummy[[k]], ")", lab, " is of type numeric. It is automatically coerced to logical.")
        lv$nrc_dummy <- stats::setNames(as.logical(lv$nrc_dummy), lv$id)
      }
      if(!is.logical(lv$nrc_dummy))
        stop("The non-reponse correction dummy variable (", arg$nrc_dummy[[k]], ")", lab, " should be of type logical or numeric.")
      if(anyNA(lv$nrc_dummy))
        stop("The non-reponse correction dummy variable (", arg$nrc_dummy[[k]], ")", lab, " should not contain any missing (NA) values.")
    }
    
    # nrc_weight
    if(!is.null(lv$nrc_weight)){
      if(!is.numeric(lv$nrc_weight))
        stop("The weights after non-response correction (", arg$nrc_weight[[k]], ")", lab, " should be numeric.")
      if(anyNA(lv$nrc_weight[lv$response_dummy %in% TRUE & lv$nrc_dummy %in% TRUE])) stop(
        "The weights after non-response correction (", arg$nrc_weight[[k]], ")", lab, " should not contain any missing (NA) values ",
        "for responding units (", arg$response_dummy[[k]], ") having taken part in the non-reponse correction process (", arg$nrc_dummy[[k]], ")."
      )
    }
    
    if(k == sampling_stages){
      
      # dissemination_dummy
      if(is.numeric(lv$dissemination_dummy)){
        note("The dissemination dummy variable (", arg$dissemination_dummy, ") is of type numeric. It is automatically coerced to logical.")
        lv$dissemination_dummy <- stats::setNames(as.logical(lv$dissemination_dummy), lv$id)
      }
      if(!is.logical(lv$dissemination_dummy))
        stop("The dissemination dummy variable (", arg$dissemination_dummy, ") should be of type logical or numeric.")
      if(anyNA(lv$dissemination_dummy))
        stop("The dissemination dummy variable (", arg$dissemination_dummy, ") should not contain any missing (NA) values.")
      
      # dissemination_weight
      if(!is.numeric(lv$dissemination_weight))
        stop("The dissemination weights (", arg$dissemination_weight, ") should be numeric.")
      if(anyNA(lv$dissemination_weight[lv$dissemination_dummy])) stop(
        "The dissemination weights (", arg$dissemination_weight, ") should not contain ",
        "any missing (NA) values for disseminated units (", arg$dissemination_dummy, ")."
      )
      
      # disseminated out-of-scope units
      if(!is.null(arg$scope_dummy)){
        disseminated_out_of_scope <- lv$id[lv$dissemination_dummy & !lv$scope_dummy]
        if(length(disseminated_out_of_scope) > 0) stop(
          "The following units are out-of-scope (", arg$scope_dummy, ") but nonetheless disseminated (",
          arg$dissemination_dummy, "): ", display_only_n_first(disseminated_out_of_scope), "."
        )
      }
    }
    
    if(k == calibration_stage && !is.null(arg$calibration_weight)){
      
      # calibration_dummy
      if(is.null(lv$calibration_dummy)) lv$calibration_dummy <- lv$response_dummy
      if(is.numeric(lv$calibration_dummy)){
        note("The dummy variable indicating the units used in the calibation process (", arg$calibration_dummy, ") is of type numeric. It is automatically coerced to logical.")
        lv$calibration_dummy <- stats::setNames(as.logical(lv$calibration_dummy), lv$id)
      }
      if(!is.logical(lv$calibration_dummy))
        stop("The dummy variable indicating the units used in the calibation process (", arg$calibration_dummy, ") should be of type logical or numeric.")
      if(anyNA(lv$calibration_dummy))
        stop("The dummy variable indicating the units used in the calibation process (", arg$calibration_dummy, ") should not contain any missing (NA) values.")
      
      # calibration_weight
      if(!is.numeric(lv$calibration_weight))
        stop("The weights after calibration (", arg$calibration_weight, ") should be numeric.")
      if(anyNA(lv$calibration_weight[lv$calibration_dummy %in% TRUE]))
        stop("The weights after calibration (", arg$calibration_weight, ") should not contain any missing (NA) values for units used in the calibration process.")
      
      # calibration_var
      calibration_var_quanti <- names(which(sapply(lv$calibration_var, function(var) is.numeric(var) || is.logical(var))))
      calibration_var_quali <- names(which(sapply(lv$calibration_var, function(var) is.factor(var) || is.character(var))))
      calibration_var_pb_type <- setdiff(arg$calibration_var, c(calibration_var_quanti, calibration_var_quali))
      if(length(calibration_var_pb_type) > 0) stop(
        "The following calibration variables are neither quantitative (numeric, logical) nor qualitative (factor, character): ",
        display_only_n_first(calibration_var_pb_type), "."
      )
      if(length(calibration_var_quali) > 0) note(
        "Note: The following calibration variables are qualitative (factor, character): ",
        display_only_n_first(calibration_var_quali), ". They will be automatically discretized."
      )
      calibration_var_pb_NA <- names(which(sapply(lv$calibration_var, function(var) anyNA(var[lv$calibration_dummy %in% TRUE]))))
      if(length(calibration_var_pb_NA) > 0) stop(
        "The following calibration variables contain missing (NA) values for units used in the calibration process: ",
        display_only_n_first(calibration_var_pb_NA, collapse = " "), "."
      )
    }
    
    lev[[k]] <- lv
  }
  
  
  # Step 4: Define methodological quantities ----
  
  technical_stages <- vector("list", sampling_stages)
  cumulated_weight <- NULL
  calib <- NULL
  
  for(k in seq_len(sampling_stages)){
    
    lv <- lev[[k]]
    lab <- if(sampling_stages > 1) paste0(" at stage ", k) else ""
    samp_exclude <- stats::setNames(rep(FALSE, length(lv$id)), lv$id)
    
    # Logical controls
    inconsistency <- list(
      out_of_scope_and_responding = lv$id[!lv$scope_dummy & lv$response_dummy]
    )
    if(k == sampling_stages && !is.null(arg$scope_dummy) && !is.null(arg$response_dummy[[k]]) &&
       any(sapply(inconsistency, length) > 0)) stop(
         "Some arguments are inconsistent:",
         if(length(inconsistency$out_of_scope_and_responding) > 0) paste0(
           "\n  - the following units are classified both as out-of-scope units (",
           arg$scope_dummy, " variable) and as responding units (", arg$response_dummy[[k]],
           " variable): ", display_only_n_first(inconsistency$out_of_scope_and_responding), "."
         )
       )
    
    if(k == 1){
      strata_eff <- lv$strata
    }else{
      strata_eff <- stats::setNames(
        interaction(lv$parent, lv$strata, drop = TRUE, sep = ":"),
        lv$id
      )
      note(
        "At stage ", k, ", units are drawn within the units of the stage above: ",
        "the effective strata used for variance estimation are the crossing of the parent unit (",
        arg$parent_id[[k]], ")",
        if(!is.null(arg$strata[[k]])) paste0(" and of the strata variable (", arg$strata[[k]], ")") else "",
        "."
      )
    }
    
    # Exclude strata with only one sampled unit
    strata_with_one_sampled_unit <-
      names(which(tapply(lv$id[!samp_exclude], strata_eff[!samp_exclude], length) == 1))
    if(length(strata_with_one_sampled_unit) > 0){
      warn(
        "The following strata", lab, " contain less than two sampled units: ",
        display_only_n_first(strata_with_one_sampled_unit), ". ",
        "They are excluded from the variance estimation process (but kept for point estimates)."
      )
      samp_exclude <- samp_exclude | as.character(strata_eff) %in% strata_with_one_sampled_unit
    }
    
    # Enforce equal probabilities in each stratum
    sampling_weight_equal <- lv$sampling_weight
    strata_with_unequal_sampling_weight <-
      names(which(tapply(sampling_weight_equal[!samp_exclude], strata_eff[!samp_exclude], stats::sd) > 1e-6))
    if(length(strata_with_unequal_sampling_weight) > 0){
      # TODO: Enhance warning message when strata = NULL
      warn(
        "The following strata", lab, " contain units whose sampling weights are not exactly equal: ",
        display_only_n_first(strata_with_unequal_sampling_weight), ". ",
        "The mean weight per stratum is used instead."
      )
      sampling_weight_equal[!samp_exclude] <-
        tapply(sampling_weight_equal, strata_eff, base::mean)[as.character(strata_eff[!samp_exclude])]
    }
    
    corrected_weight <- lv$sampling_weight
    if(!is.null(lv$nrc_weight))
      corrected_weight[lv$response_dummy & lv$nrc_dummy] <-
      lv$nrc_weight[lv$response_dummy & lv$nrc_dummy]
    
    upper_weight <- if(k == 1) stats::setNames(rep(1, length(lv$id)), lv$id) else
      stats::setNames(cumulated_weight[lv$parent], lv$id)
    cumulated_weight <- upper_weight * corrected_weight
    
    if(k == calibration_stage && !is.null(arg$calibration_weight)){
      g_ratio <- lv$calibration_weight[lv$calibration_dummy] / cumulated_weight[lv$calibration_dummy]
      if(stats::median(abs(log(g_ratio)), na.rm = TRUE) > log(2)) warn(
        "The calibrated weights (", arg$calibration_weight, ") strongly differ from the cumulated weights ",
        "guessed from the survey description (median calibration ratio: ", round(stats::median(g_ratio, na.rm = TRUE), 2),
        " instead of about 1). This most often indicates that some sampling or non-response correction weights ",
        "are CUMULATED across stages whereas conditional weights are expected (see conventions at the top of the file)."
      )
      cumulated_weight[lv$calibration_dummy] <- lv$calibration_weight[lv$calibration_dummy]
    }
    
    if(k == sampling_stages){
      guessed_weight_not_matching_dissemination_weight <- lv$id[
        lv$dissemination_dummy &
          abs(cumulated_weight - lv$dissemination_weight) > 1e-6 * abs(lv$dissemination_weight)
      ]
      if(length(guessed_weight_not_matching_dissemination_weight)) stop(
        "The following units have a disseminated weight (", arg$dissemination_weight,
        ") that does not match the one guessed from the survey description: ",
        display_only_n_first(guessed_weight_not_matching_dissemination_weight), "."
      )
    }
    
    # Sampling
    samp <- list()
    samp$id <- lv$id
    samp$exclude <- samp_exclude[samp$id]
    samp$weight <- sampling_weight_equal[samp$id]
    samp$strata <- strata_eff[samp$id]
    samp$precalc <- suppressMessages(with(samp, var_srs(
      y = NULL, pik = 1 / weight[!exclude], strata = strata[!exclude]
    )))
    samp <- samp[c("id", "exclude", "precalc")]
    
    # Non-reponse
    if(!is.null(lv$nrc_weight)){
      nrc <- list()
      nrc$id <- lv$id[lv$response_dummy]
      nrc$sampling_weight <- lv$sampling_weight[nrc$id]
      nrc$response_prob <- (lv$sampling_weight / lv$nrc_weight)[nrc$id]
    }else nrc <- NULL
    
    # Calibration
    if(k == calibration_stage && !is.null(arg$calibration_weight)){
      calib <- list()
      calib$stage <- k
      calib$id <- lv$id[lv$response_dummy & lv$calibration_dummy]
      calib$weight <- lv$calibration_weight[calib$id]
      calib$var <- lv$calibration_var[calib$id, , drop = FALSE]
      calib$var[calibration_var_quanti] <-
        lapply(calib$var[calibration_var_quanti], Matrix)
      calib$var[calibration_var_quali] <-
        lapply(calib$var[calibration_var_quali], discretize_qualitative_var)
      calib$var <- do.call(cbind, calib$var)
      # TODO: Handle the node stack overflow problem
      calib$precalc <- res_cal(y = NULL, x = calib$var, w = calib$weight)
      calib <- calib[c("stage", "id", "precalc")]
    }
    
    technical_stages[[k]] <- list(
      samp = samp,
      nrc = nrc,
      upper_weight = upper_weight,
      agg_weight = lv$sampling_weight,
      parent = if(k > 1) lv$parent else NULL
    )
    
  }
  
  # Reference id and reference weight
  reference_id <- lev[[sampling_stages]]$id[lev[[sampling_stages]]$dissemination_dummy]
  reference_weight <- lev[[sampling_stages]]$dissemination_weight[reference_id]
  
  
  # Step 5: Define the variance wrapper ----
  qvar_variance_wrapper <- define_variance_wrapper(
    variance_function = qvar_ms_variance_function,
    reference_id = reference_id,
    reference_weight = reference_weight,
    default_id = arg$id[[sampling_stages]],
    technical_data = list(stages = technical_stages, calib = calib)
  )
  
  qvar_variance_wrapper
  
}


#' Variance function for multistage designs (internal, used by qvar_ms)
#'
#' @description \code{qvar_ms_variance_function} is the variance function
#'   embedded in the variance wrappers defined by \code{\link{qvar_ms}}. It
#'   implements a stage-by-stage variance decomposition for stratified
#'   multistage sampling designs with non-response correction through
#'   reweighting and calibration on margins, following the Rao (1975)
#'   formula as implemented in the SAS Poulpe macros (Caron, 1998).
#'
#' @param y A (sparse) numerical matrix of the variable(s) whose variance of
#'   their total is to be estimated. Its rows match the disseminated units
#'   of the last sampling stage (see \code{\link{define_variance_wrapper}}).
#' @param stages A list with one element per sampling stage, ordered from the
#'   highest stage (e.g. primary sampling units) to the lowest one (final
#'   units), as built by \code{define_qvar_ms_variance_wrapper}. Each element
#'   contains: \itemize{
#'   \item \code{samp}: data for the sampling variance estimation at this
#'   stage (\code{id}, \code{exclude}, \code{precalc});
#'   \item \code{nrc}: data for the non-response variance estimation at this
#'   stage (\code{id}, \code{sampling_weight}, \code{response_prob}), or
#'   \code{NULL} when no non-response correction applies;
#'   \item \code{upper_weight}: the cumulated corrected weight of the stages
#'   above (1 at the first stage), used as the Rao (1975) row weight;
#'   \item \code{agg_weight}: the conditional sampling weight of the stage,
#'   used to aggregate the estimated totals towards the stage above;
#'   \item \code{parent}: the identifiers of the parent units at the stage
#'   above (\code{NULL} at the first stage).}
#' @param calib Data for the calibration step (\code{stage}, \code{id},
#'   \code{precalc}), or \code{NULL} when no calibration applies.
#'
#' @details The loop over stages runs upwards, from the last stage to the
#'   first one. At each stage: \enumerate{
#'   \item the y matrix is expanded to all the units sampled at the stage
#'   (\code{add_zero}): units not already present (non-responding,
#'   out-of-scope or non-disseminated units) are given a zero value. This
#'   expansion does not include non-responding units in the calibration or
#'   non-response correction steps, as these steps only operate on the rows
#'   they explicitly select; it only guarantees their presence, with a zero
#'   value, at the sampling variance step, exactly like the \code{add_zero}
#'   call in the single-stage \code{\link{qvar}} case;
#'   \item y is replaced by the calibration residuals (\code{\link{res_cal}})
#'   if calibration took place at the stage;
#'   \item the non-response variance is estimated (\code{\link{var_pois}},
#'   Poisson approximation) and y is then expanded by the inverse of the
#'   estimated response probabilities;
#'   \item the sampling variance conditional on the stages above is
#'   estimated (\code{\link{var_srs}}, stratified simple random sampling,
#'   the strata being crossed with the parent units from stage 2 onwards);
#'   \item y is aggregated towards the stage above as estimated totals per
#'   parent unit, weighted by the conditional sampling weights
#'   (\code{sum_by}, which handles sparse and multi-column y matrices).}
#'
#'   Following the Rao (1975) formula, the variance contributions of steps 3
#'   and 4 are weighted by \code{upper_weight}, the cumulated corrected
#'   weight of the stages above. This weight enters \emph{linearly}, through
#'   the \code{w} argument of \code{\link{var_pois}} and
#'   \code{\link{var_srs}} ("row weight used at the final summation step
#'   [...] applying the Rao (1975) formula", see \code{\link{varDT}}). This
#'   is consistent both with the single-stage \code{\link{qvar}} case (where
#'   \code{var_pois} is called with \code{w = sampling_weight}) and with the
#'   SAS Poulpe macros. Scaling y itself by \code{upper_weight} would amount
#'   to a quadratic weighting (\code{upper_weight^2}), as the variance
#'   formulas are quadratic in y, and would therefore overestimate the
#'   variance components of stages 2 and beyond by a factor of about
#'   \code{upper_weight}. Note also that, as \code{\link{varDT}} internally
#'   discards exhaustive units (\code{pik = 1}), the Rao row weight is
#'   restricted accordingly before being passed to \code{\link{var_srs}}.
#'
#' @return The estimated variances as a numerical vector of size the number
#'   of columns of y (sum of the per-stage sampling and non-response
#'   components, plus the calibration effect through the residuals).
#'
#' @seealso \code{\link{qvar_ms}}, \code{\link{define_variance_wrapper}},
#'   \code{\link{varDT}}, \code{\link{var_pois}}, \code{\link{res_cal}}
#'
#' @references
#'   Caron N. (1998), "Le logiciel Poulpe : aspects méthodologiques",
#'   \emph{Actes des Journées de méthodologie statistique}
#'   \url{http://jms-insee.fr/jms1998s03_1/}
#'
#'   Rao, J.N.K. (1975), "Unbiased variance estimation for multistage
#'   designs", \emph{Sankhya}, C n°37
#'
#' @keywords internal

qvar_ms_variance_function <- function(y, stages, calib){
  
  var <- list()
  
  for(k in rev(seq_along(stages))){
    
    lev_k <- stages[[k]]
    
    y <- add_zero(y, rownames = lev_k$samp$id)
    upper_weight <- lev_k$upper_weight[rownames(y)]
    
    # Calibration
    if(!is.null(calib) && calib$stage == k){
      y[calib$id, ] <- res_cal(y = y[calib$id, , drop = FALSE], precalc = calib$precalc)
    }
    
    # Non-response
    if(!is.null(lev_k$nrc)){
      var[[paste0("nr", if(length(stages) > 1) paste0("_stage", k) else "")]] <- var_pois(
        y[lev_k$nrc$id, , drop = FALSE],
        pik = lev_k$nrc$response_prob,
        w = upper_weight[lev_k$nrc$id] * lev_k$nrc$sampling_weight
      )
      y[lev_k$nrc$id, ] <- as.matrix(Diagonal(x = 1 / lev_k$nrc$response_prob) %*% y[lev_k$nrc$id, , drop = FALSE])
    }
    
    # Sampling
    w_rao <- upper_weight[!lev_k$samp$exclude]
    if(!is.null(lev_k$samp$precalc$exh)) w_rao <- w_rao[!lev_k$samp$precalc$exh]
    var[[paste0("sampling", if(length(stages) > 1) paste0("_stage", k) else "")]] <- var_srs(
      y = y[!lev_k$samp$exclude, , drop = FALSE],
      precalc = lev_k$samp$precalc,
      w = w_rao
    )
    
    if(k > 1){
      y <- sum_by(
        y, by = lev_k$parent[rownames(y)],
        w = lev_k$agg_weight[rownames(y)]
      )
    }
    
  }
  
  print(var)
  
  # Final summation
  Reduce(`+`, var)
  
}