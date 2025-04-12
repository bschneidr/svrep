#' @title Retain only a random subset of the replicates in a design
#' @description Randomly subsamples the replicates of a survey design object,
#' to keep only a subset. The scale factor used in estimation is increased
#' to account for the subsampling.
#' @param design A survey design object, created with either the \code{survey} or \code{srvyr} packages.
#' @param n_reps The number of replicates to keep after subsampling
#'
#' @return An updated survey design object, where only a random selection
#' of the replicates has been retained. The overall 'scale' factor for the design
#' (accessed with \code{design$scale}) is increased to account for the sampling of replicates.
#'
#' @section Statistical Details:
#'
#' Suppose the initial replicate design has \eqn{L} replicates, with
#' respective constants \eqn{c_k} for \eqn{k=1,\dots,L} used to estimate variance
#' with the formula
#' \deqn{v_{R} = \sum_{k=1}^L c_k\left(\hat{T}_y^{(k)}-\hat{T}_y\right)^2}
#'
#' With subsampling of replicates, \eqn{L_0} of the original \eqn{L} replicates
#' are randomly selected, and then variances are estimated using the formula:
#' \deqn{v_{R} = \frac{L}{L_0} \sum_{k=1}^{L_0} c_k\left(\hat{T}_y^{(k)}-\hat{T}_y\right)^2}
#'
#' This subsampling is suggested for certain replicate designs in Fay (1989).
#' Kim and Wu (2013) provide a detailed theoretical justification and
#' also propose alternative methods of subsampling replicates.
#'
#' @references
#'
#' Fay, Robert. 1989.
#' "Theory And Application Of Replicate Weighting For Variance Calculations."
#' In, 495-500. Alexandria, VA: American Statistical Association.
#' http://www.asasrms.org/Proceedings/papers/1989_033.pdf
#'
#' Kim, J.K. and Wu, C. 2013.
#' "Sparse and Efficient Replication Variance Estimation for Complex Surveys."
#' \strong{Survey Methodology}, Statistics Canada, 39(1), 91-120.
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
#'     subsample_replicates(n_reps = 4) |>
#'     getElement("repweights")
subsample_replicates <- function(design, n_reps) {

  if (!inherits(design, "svyrep.design")) {
    stop("`design` must be a replicate design object.")
  }
  is_compressed <- inherits(design$repweights, "repweights_compressed")

  if (is.null(n_reps) || (length(n_reps) > 1) || (!is.numeric(n_reps)) || (is.na(n_reps))) {
    stop("`n` must be a single number")
  }
  if ((n_reps < 2)) {
    stop("`n_reps` must be at least 2.")
  }

  # Count the number of replicates
  if (!is_compressed) {
    orig_n_reps <- ncol(design$repweights)
  }
  if (is_compressed) {
    orig_n_reps <- ncol(design$repweights[['weights']])
  }

  # Determine the new order of the replicates
  selected_subsample <- sample(
    x = seq_len(orig_n_reps), size = n_reps, replace = FALSE
  )

  # Update the overall scale factor
  subsample_rate <- n_reps / orig_n_reps
  scale_adjustment <- subsample_rate^(-1)

  design$scale <- scale_adjustment * design$scale

  # Retrieve the original 'scale' and 'rscales' attributes of the weights matrix
  if (!is_compressed) {
    orig_scale_attribute <- scale_adjustment * attr(design$repweights, 'scale')
    orig_rscales_attribute <- attr(design$repweights, 'rscales')
  }
  if (is_compressed) {
    orig_scale_attribute <- scale_adjustment * attr(design$repweights[['weights']], 'scale')
    orig_rscales_attribute <- attr(design$repweights[['weights']], 'rscales')
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
