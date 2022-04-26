#' @title Calibrate weights from a primary survey to estimated totals from a control survey,
#' with replicate-weight adjustments that account for variance of the control totals
#' @description Calibrate the weights of a primary survey to match estimated totals from a control survey,
#' using adjustments to the replicate weights to account for the variance of the estimated control totals.
#' The adjustments to replicate weights are conducted using the method proposed by Fuller (1998).
#' This method can be used to implement general calibration as well as post-stratification or raking specifically
#' (see the details for the \code{calfun} parameter).
#' @details With the Fuller method, each of \code{k} randomly-selected replicate columns from the primary survey
#' are calibrated to control totals formed by perturbing the \code{k}-dimensional vector of
#' estimated control totals using a spectral decomposition of the variance-covariance matrix
#' of the estimated control totals. Other replicate columns are simply calibrated to the unperturbed control totals.
#' \cr
#'
#' Because the set of replicate columns whose control totals are perturbed should be random,
#' there are multiple ways to ensure that this matching is reproducible.
#' The user can either call \link[base]{set.seed} before using the function,
#' or supply a vector of randomly-selected column indices to the argument \code{perturbed_control_cols}.
#'
#' @param rep_design A replicate design object for the primary survey, created with either the \code{survey} or \code{srvyr} packages.
#' @param estimate A vector of estimated control totals.
#' The names of entries must match the names from calling \code{svytotal(x = cal_formula, design = rep_design)}.
#' @param vcov_estimate A variance-covariance matrix for the estimated control totals.
#' The column names and row names must match the names of \code{estimate}.
#' @param cal_formula A formula listing the variables to use for calibration.
#' All of these variables must be included in \code{rep_design}.
#' @param calfun A calibration function from the \code{survey} package,
#' such as \link[survey]{cal.linear}, \link[survey]{cal.raking}, or \link[survey]{cal.logit}.
#' Use \code{cal.linear} for ordinary post-stratification, and \code{cal.raking} for raking.
#' See \link[survey]{calibrate} for additional details.
#' @param bounds Parameter passed to \link[survey]{grake} for calibration. See \link[survey]{calibrate} for details.
#' @param verbose Parameter passed to \link[survey]{grake} for calibration. See \link[survey]{calibrate} for details.
#' @param maxit Parameter passed to \link[survey]{grake} for calibration. See \link[survey]{calibrate} for details.
#' @param epsilon Parameter passed to \link[survey]{grake} for calibration. See \link[survey]{calibrate} for details.
#' @param variance Parameter passed to \link[survey]{grake} for calibration. See \link[survey]{calibrate} for details.
#' @param perturbed_control_cols Optional parameter to determine which replicate columns
#' will have their control totals perturbed. If supplied, \code{col_selection} must be an integer vector
#' with length equal to the length of \code{estimate}.
#' @return A replicate design object, with full-sample weights calibrated to totals from \code{estimate},
#' and replicate weights adjusted to account for variance of the control totals.
#' The element \code{perturbed_control_cols} indicates, for each replicate column of the calibrated primary survey,
#' which column of replicate weights it was matched to from the control survey.
#' @references
#' Fuller, W.A. (1998).
#' "Replication variance estimation for two-phase samples."
#' \strong{Statistica Sinica}, \emph{8}: 1153-1164.
#' @export
#'
#' @examples
#'
#' # Load example data for primary survey ----
#'
#'   suppressPackageStartupMessages(library(survey))
#'   data(api)
#'
#'   primary_survey <- svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc) |>
#'     as.svrepdesign(type = "JK1")
#'
#' # Load example data for control survey ----
#'
#'   control_survey <- svydesign(id = ~ 1, fpc = ~fpc, data = apisrs) |>
#'     as.svrepdesign(type = "JK1")
#'
#' # Estimate control totals ----
#'
#'   estimated_controls <- svytotal(x = ~ stype + enroll,
#'                                  design = control_survey)
#'   control_point_estimates <- coef(estimated_controls)
#'   control_vcov_estimate <- vcov(estimated_controls)
#'
#' # Calibrate totals for one categorical variable and one numeric ----
#'
#'   calibrated_rep_design <- calibrate_to_estimate(
#'     rep_design = primary_survey,
#'     estimate = control_point_estimates,
#'     vcov_estimate = control_vcov_estimate,
#'     cal_formula = ~ stype + enroll
#'   )
#'
#' # Inspect estimates before and after calibration ----
#'
#'   ##_ For the calibration variables, estimates and standard errors
#'   ##_ from calibrated design will match those of the control survey
#'
#'     svytotal(x = ~ stype + enroll, design = primary_survey)
#'     svytotal(x = ~ stype + enroll, design = control_survey)
#'     svytotal(x = ~ stype + enroll, design = calibrated_rep_design)
#'
#'   ##_ Estimates from other variables will be changed as well
#'
#'     svymean(x = ~ api00 + api99, design = primary_survey)
#'     svymean(x = ~ api00 + api99, design = control_survey)
#'     svymean(x = ~ api00 + api99, design = calibrated_rep_design)
#'
#' # Inspect weights before and after calibration ----
#'
#'   summarize_rep_weights(primary_survey, type = 'overall')
#'   summarize_rep_weights(calibrated_rep_design, type = 'overall')
#'
#' # For reproducibility, specify which columns are randomly selected for Fuller method ----
#'
#'   column_selection <- calibrated_rep_design$perturbed_control_cols
#'   print(column_selection)
#'
#'   calibrated_rep_design <- calibrate_to_estimate(
#'     rep_design = primary_survey,
#'     estimate = control_point_estimates,
#'     vcov_estimate = control_vcov_estimate,
#'     cal_formula = ~ stype + enroll,
#'     col_selection = column_selection
#'   )

calibrate_to_estimate <- function(rep_design,
                                  estimate, vcov_estimate,
                                  cal_formula,
                                  calfun = survey::cal.linear,
                                  bounds = list('lower' = -Inf, 'upper' = Inf),
                                  verbose = FALSE, maxit = 50,
                                  epsilon = 1e-7, variance = NULL,
                                  col_selection = NULL) {

  if (!inherits(rep_design, "svyrep.design")) {
    stop("`rep_design` must be a replicate survey design object, with class `svyrep.design`")
  }

  # Get description of estimate ----

  k <- length(estimate)
  estimate_names <- names(estimate)

  if (is.null(estimate_names)) {
    stop("`estimate` must be a named vector, with names matching result from `svytotal(x = cal_formula, design = rep_design)`.")
  }
  if (length(estimate) > 1) {
    if (!is.matrix(vcov_estimate)) {
      stop("`vcov_estimate` must be a matrix.")
    }
    if (any(estimate_names != colnames(vcov_estimate))) {
      stop("`vcov_estimate` must have row names and column names exactly matching the names of `estimate`.")
    }
    if (!isSymmetric.matrix(vcov_estimate)) {
      stop("`vcov_estimate` must be a symmetric matrix.")
    }
  }

  # Determine parameters describing replicate designs ----
  R_primary <- ncol(rep_design$repweights)

  A_primary <- rep_design$scale

  rscales_primary <- rep_design$rscales

  # Use user-supplied matching selection of columns, otherwise select at random ----

  if (is.null(col_selection)) {
    col_selection <- sample(x = R_primary, size = k, replace = FALSE)
    matching_msg <- paste("Selection of replicate columns whose control totals will be perturbed will be done at random.",
                          "For tips on reproducible selection, see `help('calibrate_to_estimate')`",
                          sep = "\n")
    message(matching_msg)
  } else {
    if (length(col_selection) != k) {
      stop("`col_selection` must have the same length as `estimate`, with no duplicate entries.")
    }
    if (any(col_selection) != as.integer(col_selection)) {
      stop("`col_selection` must only contain integer entries.")
    }
    if (any(col_selection < 1) || any(col_selection > k)) {
      stop("`col_selection` must be an integer vector with entries whose value is between 1 and k, where k is the length of `estimate`.")
    }
    if (length(col_selection) != length(unique(col_selection))) {
      stop("`col_selection` must have k distinct entries, where k is the length of `estimate`.")
    }
  }

  # Calculate spectral decomposition ----
  eigen_decomposition <- eigen(x = vcov_estimate,
                               symmetric = TRUE)

  # Calculate matrix of replicate control totals ----
  v <- sapply(X = seq_along(eigen_decomposition$values),
              FUN = function(k) {
                truncated_eigenvalue <- ifelse(eigen_decomposition$values[k] < 0,
                                               0, eigen_decomposition$values[k])
                sqrt(truncated_eigenvalue) * eigen_decomposition$vectors[,k]
              })

  # Create needed data matrices for primary design ----

    ##_ Get dataframe with all of the variables mentioned in the formula
    mf <- model.frame(cal_formula, rep_design$variables, na.action=na.pass)

    ##__For each factor/character variable, obtain a matrix of dummy variables, one dummy per category.
    ##_ For numeric variables, simply return a matrix with the variable itself
    xx <- lapply(
      attr(terms(cal_formula),"variables")[-1],
      function(tt) model.matrix( eval(bquote( ~0 + .(tt))), mf)
    )

    ##_ Determine dimensions of matrix of all variables
    cols <- sapply(xx, NCOL)

    x <- matrix(nrow = NROW(xx[[1]]), ncol = sum(cols))

    scols <- c(0, cumsum(cols))

    ##_ Combine all of the separate matrices into a single matrix
    for (i in 1:length(xx)){
      x[,scols[i]+1:cols[i]] <- xx[[i]]
    }
    colnames(x) <- do.call("c",lapply(xx,colnames))

  # Check that all of the control totals have corresponding variables from primary survey ----

    if (any(!colnames(x) %in% estimate_names) || any(!estimate_names %in% colnames(x))) {
      stop("Using `svytotal(x = cal_formula, design = rep_design)` should yield estimates with the same names as `estimate`.")
    } else {
      x <- x[,estimate_names, drop = FALSE]
    }

  # Extract replicate weights matrix ----

  primary_replicate_weights <- weights(rep_design, type = 'analysis')

  # Generate replicate-specific control totals ----

    replicate_control_totals <- matrix(data = estimate,
                                       nrow = R_primary,
                                       ncol = k,
                                       byrow = TRUE)

    for (i in seq_len(k)) {
      i_star <- col_selection[i]
      replicate_control_totals[i_star,] <-  estimate + (sqrt(1/(A_primary * rscales_primary[i_star])) * v[,i])
    }

  # Calibrate the replicate weights ----

  adjusted_replicate_weights <- matrix(nrow = nrow(primary_replicate_weights),
                                       ncol = ncol(primary_replicate_weights))
  for (i in seq_len(R_primary)) {
    g_weights <- survey::grake(mm = x, ww = primary_replicate_weights[,i, drop = TRUE],
                               population = replicate_control_totals[i, ,drop = TRUE],
                               calfun = calfun,
                               bounds = bounds,
                               verbose = verbose, maxit = maxit,
                               epsilon = epsilon, variance = variance)
    adjusted_replicate_weights[,i] <- as.vector(primary_replicate_weights[,i]) * g_weights
  }

  # Calibrate the full-sample weights ----

  g_weights <- grake(mm = x, ww = as.vector(rep_design$pweights),
                     population = estimate,
                     calfun = calfun,
                     bounds = bounds,
                     verbose = verbose, maxit = maxit,
                     epsilon = epsilon, variance = variance)

  adjusted_fullsample_weights <- as.vector(rep_design$pweights) * g_weights
  attr(adjusted_fullsample_weights, 'eta') <- NULL

  # Assemble the updated replicate design object ----

  calibrated_rep_design <- rep_design

  calibrated_rep_design$pweights <- adjusted_fullsample_weights
  names(calibrated_rep_design$pweights) <- names(rep_design)

  calibrated_rep_design$repweights <- adjusted_replicate_weights
  class(calibrated_rep_design$repweights) <- 'repweights'
  calibrated_rep_design$combined.weights <- TRUE

  if (rep_design$type %in% c("JK1", "JKn", "JK2", "ACS", "successive-difference")) {
    rep_design_rho <- NULL
  } else {
    rep_design_rho <- rep_design$rho
  }

  rep_design_type <- ifelse(
    rep_design$type %in% c("bootstrap", "subbootstrap", "mrbbootstrap"),
    "bootstrap", rep_design$type
  )

  calibrated_rep_design <- survey::svrepdesign(
    data = rep_design$variables,
    repweights = adjusted_replicate_weights,
    weights = adjusted_fullsample_weights,
    type = rep_design_type,
    combined.weights = TRUE,
    rho = rep_design_rho,
    scale = rep_design$scale,
    rscales = rep_design$rscales,
    fpc = rep_design$fpc,
    fpctype = rep_design$fpctype,
    mse = TRUE
  )

  if (!rep_design$mse) {
    warning("Setting `mse` to TRUE; variance estimates will be centered around full-sample estimate, not mean of replicates.")
  }

  # Indicate which replicate columns correspond had their control totals perturbed using Fuller's method ----
  calibrated_rep_design$perturbed_control_cols <- col_selection

  # Return the result ----
  return(calibrated_rep_design)
}
