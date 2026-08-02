
#' @title Replicate Design Object for a Two-phase Sample
#' @description Derives a replicate design object for a two-phase sample,
#' based on a replicate design object created for the first phase sample.
#' The full-sample and replicate weights are adjusted to implement
#' the reweighted expansion estimator (REE).
#' @param design A replicate survey design object for the first phase sample.
#' @param phase_two_indicators A string giving the name of a variable in the data
#' that indicates which cases are selected for the phase two sample.
#' @param phase_two_strata A string giving the name of the stratification
#' variable for the second phase sample, or \code{NULL} (the default) if
#' the second-phase sampling was not stratified.
#' @param phase_two_probs A string giving the name of a variable in the data
#' that represents the inclusion probabilities for the second-phase sampling.
#'
#' @returns A replicate survey design object, containing only the observations
#' from the second phase sample. The full-sample and replicate weights
#' are adjusted using the approach developed by Kim and Yu (2011).
#'
#' @details
#' This function adjusts the full-sample and replicate weights so that they correspond
#' to the reweighted expansion estimator (REE), as described in Kim and Yu (2011).
#' Note that the REE gives different estimates than the double expansion estimator (DEE),
#' which is what is implemented by the function [survey::twophase()].
#' 
#' The REE weights for the second-phase sample are derived in two steps. 
#' First, the phase-one weights are divided by the second-phase conditional inclusion probabilities. 
#' Next, the resulting weights for the second-phase sample undergo a post-stratification 
#' adjustment, where the phase-two sampling strata are the post-strata. This
#' ensures that the weight sums for the second-phase sample in each second-phase stratum 
#' match the corresponding weight sums for the entire first-phase sample. 
#' 
#' See Kim and Yu (2011) for the underlying theory behind this method. Section 4 of 
#' Opsomer et al. (2016) provides a clear summary of the method, in the context of discussing
#' its application to successive difference replication weights. 
#' 
#' Opsomer et al. (2016) point out that special care is needed when the second-phase sampling method uses
#' systematic sampling. In this case, the user should identify any implicit strata created
#' by the the second-phase systematic sampling. The variable \code{phase_two_strata} should
#' be defined so as to represent the combination of explicit and implicit strata.
#' 
#' If the cube method was used for the second-phase sampling, then this function's approach
#' to variance estimation will overestimate variance due to second-phase sampling,
#' unless additional calibration is applied. After calling this function, use the
#' function [calibrate_to_sample()] to calibrate the totals for balancing variables.
#' 
#' @references
#' Kim, J. K., and Yu, C. L. (2011), "Replication Variance Estimation under Two-Phase Sampling."
#' Survey Methodology, 37, 67–74. https://www150.statcan.gc.ca/n1/en/pub/12-001-x/2011001/article/11448-eng.pdf
#'
#' Opsomer, J. D., Jay Breidt, F., White, M., and Li, Y. (2016). "Successive Difference Replication Variance Estimation in Two-Phase Sampling." Journal of Survey Statistics and Methodology, 4(1), 43–70. https://doi.org/10.1093/jssam/smv033
#' @export
#' @seealso Use [calibrate_to_sample()] to calibrate the second-phase sample
#' to match the first-phase sample for specific calibration variables.
#' 
#' The generalized replication functions [as_fays_gen_rep_design()] and [as_gen_boot_design()] 
#' can be used to create replicate weights for a two-phase design object created with [survey::twophase()].
#' @examples
#' # Create example data
#' phase_one_data <- data.frame(
#'   id = c(1, 2, 3, 4, 5, 6, 7, 8),
#'   x  = c(109, 95, 86, 102, 106, 105, 106, 115),
#'   y  = c(103, 60, 92, 76, 104, 132, 127, 88),
#'   PHASE_ONE_PSU     = c(1, 2, 3, 4, 5, 6, 7, 8),
#'   PHASE_ONE_WGT     = rep(100, times = 8),
#'   PHASE_TWO_STRATA  = rep(c(1,2), each = 4),
#'   PHASE_TWO_PROB    = rep(c(0.5, 0.75), each = 4),
#'   PHASE_TWO_SAMPLED = c(TRUE, TRUE, FALSE, FALSE, 
#'                         TRUE, TRUE, TRUE, FALSE)
#' )
#' 
#' # Create replicate weights for the first phase
#' phase_one_design <- svydesign(
#'   data    = phase_one_data,
#'   ids     = ~ PHASE_ONE_PSU,
#'   weights = ~ PHASE_ONE_WGT
#' )
#' 
#' phase_one_rep_design <- phase_one_design |>
#'   as.svrepdesign(type = "JK1", mse = TRUE)
#' 
#' # Derive a replicate design for the two-phase sample
#' phase_two_rep_design <- derive_twophase_rep_design(
#'   design               = phase_one_rep_design,
#'   phase_two_indicators = "PHASE_TWO_SAMPLED",
#'   phase_two_probs      = "PHASE_TWO_PROB",
#'   phase_two_strata     = "PHASE_TWO_STRATA"
#' )
#' 
#' # Check estimates (and standard errors)
#' svytotal(x = ~ x + y, design = phase_two_rep_design)
#' svytotal(x = ~ x + y, design = phase_one_rep_design)
#' 
#' # Additional calibration to improve precision
#' # and align estimates with first-phase sample
#' # (NOTE: important to correctly specify `control_col_matches`)
#' calibrated_phase_two_rep_design <- calibrate_to_sample(
#'   primary_rep_design  = phase_two_rep_design,
#'   control_rep_design  = phase_one_rep_design,
#'   cal_formula         = ~ x,
#'   control_col_matches = seq_len(ncol(phase_two_rep_design$repweights))
#' )
#' 
#' svytotal(x = ~ x + y, design = calibrated_phase_two_rep_design)
#' @md
derive_twophase_rep_design <- function(
  design,
  phase_two_indicators,
  phase_two_probs,
  phase_two_strata = NULL
) {
  UseMethod("derive_twophase_rep_design", design)
}
#' @export
derive_twophase_rep_design.svyrep.design <- function(
  design,
  phase_two_indicators,
  phase_two_probs,
  phase_two_strata = NULL
) {

  # Validate the arguments
  design_varnames <- colnames(design)
  if (!is.character(phase_two_indicators) || !all(phase_two_indicators %in% design_varnames)) {
    stop("`phase_two_indicators` must be a string giving the name of a variable in the data.")
  }
  if (!is.character(phase_two_probs) || !all(phase_two_probs %in% design_varnames)) {
    stop("`phase_two_probs` must be a string giving the name of a variable in the data.")
  }
  if (length(phase_two_indicators) != 1) {
    stop("`phase_two_indicators` can only be used to specify one variable in the data.")
  }
  if (length(phase_two_probs) != 1) {
    stop("`phase_two_probs` can only be used to specify one variable in the data.")
  }
  phase_two_indicators <- design[['variables']][[phase_two_indicators]]
  phase_two_probs <- design[['variables']][[phase_two_probs]]

  if (!is.null(phase_two_strata)) {
    if (!is.character(phase_two_strata) || !all(phase_two_strata %in% design_varnames)) {
      stop("`phase_two_strata` must be a string giving the name of a variable in the data.")
    }
  }

  stopifnot(
    'The `phase_two_indicators` variable must have only TRUE or FALSE values' = is.logical(phase_two_indicators) && !any(is.na(phase_two_indicators)),
    'The `phase_two_strata` variable cannot have any missing values' = !any(is.na(design[['variables']][[phase_two_strata]]))
  )
  phase_two_probs <- phase_two_probs[phase_two_indicators]
  if (any(is.na(phase_two_probs))) {
    stop("The `phase_two_probs` variable cannot have missing values for cases in the second phase sample.")
  }
  if (any((phase_two_probs > 1) | (phase_two_probs <= 0))) {
    stop("The `phase_two_probs` variable must have values between 0 and 1")
  }

  # Create the two-phase replicate design by subsetting the input replicate design
  phase_two_rep_design <- design[phase_two_indicators,]
  
  # Adjust the full-sample weights using the conditional phase-two probabilities
  phase_two_rep_design[['pweights']] <- phase_two_rep_design[['pweights']] / phase_two_probs
  if (phase_two_rep_design$combined.weights) {
    phase_two_rep_design$repweights <- weights(phase_two_rep_design, "analysis") / phase_two_probs
  }

  # Adjust full-sample and replicate weights within strata,
  # so that the second-phase weight sums match the first-phase weight sums
  n_replicates <- length(phase_two_rep_design[['rscales']])
  if (is.null(phase_two_strata)) {
    calibration_formula <- ~ 1
  } else {
    calibration_formula <- reformulate(sprintf("factor(%s)", phase_two_strata))
  }

  phase_two_rep_design <- svrep::calibrate_to_sample(
    primary_rep_design  = phase_two_rep_design,
    control_rep_design  = design,
    cal_formula         = calibration_formula,
    calfun              = survey::cal.linear, 
    control_col_matches = seq_len(n_replicates)
  )
  phase_two_rep_design$call <- sys.call(which = -1)

  return(phase_two_rep_design)
}
