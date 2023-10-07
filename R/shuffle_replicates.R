#' @title Shuffle the order of replicates in a survey design object
#' @description Shuffle the order of replicates in a survey design object.
#' In other words, the order of the columns of replicate weights is randomly permuted.
#' @param design A survey design object, created with either the \code{survey} or \code{srvyr} packages.
#'
#' @return An updated survey design object, where the order of the replicates
#' has been shuffled (i.e., the order has been randomly permuted).
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
#' # Inspect replicates before shuffling
#'
#'   rep_design |> getElement("repweights")
#'
#' # Inspect replicates after shuffling
#'
#'   rep_design |>
#'     shuffle_replicates() |>
#'     getElement("repweights")
shuffle_replicates <- function(design) {

  if (!inherits(design, "svyrep.design")) {
    stop("`design` must be a replicate design object.")
  }

  is_compressed <- inherits(design$repweights, "repweights_compressed")

  # Count the number of replicates
  if (!is_compressed) {
    n_reps <- ncol(design$repweights)
  }
  if (is_compressed) {
    n_reps <- ncol(design$repweights[['weights']])
  }

  # Determine the new order of the replicates
  new_order <- sample(x = seq_len(n_reps), size = n_reps, replace = FALSE)

  # Retrieve the original 'scale' and 'rscales' attributes
  if (!is_compressed) {
    orig_scale_attribute <- attr(design$repweights, 'scale')
    orig_rscales_attribute <- attr(design$repweights, 'rscales')
  }
  if (is_compressed) {
    orig_scale_attribute <- attr(design$repweights[['weights']], 'scale')
    orig_rscales_attribute <- attr(design$repweights[['weights']], 'rscales')
  }

  # Update the matrix of replicate weights
  if (!is_compressed) {
    design$repweights <- design$repweights[,new_order,drop=FALSE]
    attr(design$repweights, 'scale') <- orig_scale_attribute
    attr(design$repweights, 'rscales') <- orig_rscales_attribute[new_order]
  }
  if (is_compressed) {
    design$repweights[['weights']] <- design$repweights[['weights']][,new_order,drop=FALSE]
    attr(design$repweights[['weights']], 'scale') <- orig_scale_attribute
    attr(design$repweights[['weights']], 'rscales') <- orig_rscales_attribute[new_order]
  }

  # Update the order of the replicate-specific scale factors
  design$rscales <- design$rscales[new_order]

  return(design)
}
