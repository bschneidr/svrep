#' @title Retain only a random subset of the replicates in a design
#' @description Randomly subsamples the replicates of survey design object,
#' to keep only a subset. The scale factor used in estimation is increased
#' to account for the subsampling.
#' @param design A survey design object, created with either the \code{survey} or \code{srvyr} packages.
#' @param n The number of replicates to keep after subsampling
#'
#' @return An updated survey design object, where only a random selection
#' of the replicates has been retained. The overall 'scale' factor for the design
#' (accessed with \code{design$scale}) is increased to account for the sampling of replicates.
#' @export
#'
#' @examples
#' library(survey)
#' set.seed(2023)
#'
#' # Create an example survey design object
#'
#'   sample_data <- data.frame(
#'     STRATUM = c(1,1,1,1,2,2,2,2),
#'     PSU     = c(1,2,3,4,5,6,7,8)
#'   )
#'
#'   survey_design <- svydesign(
#'     data = sample_data,
#'     strata = ~ STRATUM,
#'     ids = ~ PSU,
#'     weights = ~ 1
#'   )
#'
#'   rep_design <- survey_design |>
#'     as_fays_gen_rep_design(variance_estimator = "Ultimate Cluster")
#'
#' # Inspect replicates before subsampling
#'
#'   rep_design |> getElement("repweights")
#'
#' # Inspect replicates after subsampling
#'
#'   rep_design |>
#'     subsample_replicates(n = 4) |>
#'     getElement("repweights")
subsample_replicates <- function(design, n) {

  if (!inherits(design, "svyrep.design")) {
    stop("`design` must be a replicate design object.")
  }
  is_compressed <- inherits(design$repweights, "repweights_compressed")

  if ((length(n) > 1) || (!is.numeric(n)) || (is.na(n))) {
    stop("`n` must be a single number")
  }

  # Count the number of replicates
  if (!is_compressed) {
    n_reps <- ncol(design$repweights)
  }
  if (is_compressed) {
    n_reps <- ncol(design$repweights[['weights']])
  }

  # Determine the new order of the replicates
  selected_subsample <- sample(
    x = seq_len(n_reps), size = n, replace = FALSE
  )
  subsample_rate <- n / n_reps
  scale_adjustment <- subsample_rate^(-1)

  # Retrieve the original 'scale' and 'rscales' attributes
  if (!is_compressed) {
    orig_scale_attribute <- scale_adjustment * attr(design$repweights, 'scale')
    orig_rscales_attribute <- attr(design$repweights, 'rscales')
  }

  # Update the matrix of replicate weights
  if (!is_compressed) {
    design$repweights <- design$repweights[,selected_subsample,drop=FALSE]
    attr(design$repweights, 'scale') <- scale_adjustment * orig_scale_attribute
    attr(design$repweights, 'rscales') <- orig_rscales_attribute[selected_subsample]
  }
  if (is_compressed) {
    design$repweights[['weights']] <- design$repweights[['weights']][,selected_subsample,drop=FALSE]
    attr(design$repweights[['weights']], 'scale') <- scale_adjustment * orig_scale_attribute
    attr(design$repweights[['weights']], 'rscales') <- orig_rscales_attribute[selected_subsample]
  }

  # Update the order of the replicate-specific scale factors
  design$rscales <- design$rscales[selected_subsample]

  return(design)
}
