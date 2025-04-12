#' @title Add inactive replicates to a survey design object
#' @description Adds inactive replicates to a survey design object. An inactive
#' replicate is a replicate that does not contribute to variance estimates but
#' adds to the matrix of replicate weights so that the matrix has the desired
#' number of columns. The new replicates' values are simply equal to the full-sample weights.
#' @param design A survey design object, created with either the \code{survey} or \code{srvyr} packages.
#' @param n_total The total number of replicates
#' that the result should contain. If the design already contains \code{n_total}
#' replicates (or more), then no update is made.
#' @param n_to_add The number of additional replicates to add.
#' Can only use the \code{n_total} argument OR the \code{n_to_add} argument,
#' not both.
#' @param location Either \code{"first"}, \code{"last"} (the default), or \code{"random"}.
#' Specifies where the columns of new replicates should be located in the
#' matrix of replicate weights. Use \code{"first"}
#' to place new replicates first (i.e., in the leftmost part of the matrix),
#' \code{"last"} to place the new replicates last (i.e., in the rightmost part
#' of the matrix). Use \code{"random"} to intersperse the new replicates
#' in random column locations of the matrix; the original replicates will
#' still be in their original order.
#' @return An updated survey design object, where the number of columns
#' of replicate weights has potentially increased. The increase only happens
#' if the user specifies the \code{n_to_add} argument instead of \code{n_total},
#' of if the user specifies \code{n_total} and \code{n_total} is less than the number
#' of columns of replicate weights that the design already had.
#' @section Statistical Details:
#' Inactive replicates are also sometimes referred to as "dead replicates",
#' for example in Ash (2014). The purpose of adding inactive replicates
#' is to increase the number of columns of replicate weights without impacting
#' variance estimates. This can be useful, for example, when combining data
#' from a survey across multiple years, where different years use different number
#' of replicates, but a consistent number of replicates is desired in the combined
#' data file.
#'
#' Suppose the initial replicate design has \eqn{L} replicates, with
#' respective constants \eqn{c_k} for \eqn{k=1,\dots,L} used to estimate variance
#' with the formula
#' \deqn{v_{R} = \sum_{k=1}^L c_k\left(\hat{T}_y^{(k)}-\hat{T}_y\right)^2}
#' where \eqn{\hat{T}_y} is the estimate produced using the full-sample weights
#' and \eqn{\hat{T}_y^{(k)}} is the estimate from replicate \eqn{k}.
#'
#' Inactive replicates are simply replicates that are exactly equal to the full sample:
#' that is, the replicate \eqn{k} is called "inactive" if its vector of replicate
#' weights exactly equals the full-sample weights. In this case, when using the formula
#' above to estimate variances, these replicates contribute nothing to the variance estimate.
#'
#' If the analyst uses the variant of the formula above where the full-sample estimate
#' \eqn{\hat{T}_y} is replaced by the average replicate estimate (i.e., \eqn{L^{-1}\sum_{k=1}^{L}\hat{T}_y^{(k)}}),
#' then variance estimates will differ before vs. after adding the inactive replicates.
#' For this reason, it is strongly recommend to explicitly specify \code{mse=TRUE}
#' when creating a replicate design object in R with functions such as \code{svrepdesign()},
#' \code{as_bootstrap_design()}, etc. If working with an already existing replicate design,
#' you can update the \code{mse} option to \code{TRUE} simply by using code such as
#' \code{my_design$mse <- TRUE}.
#' @references
#' Ash, S. (2014). "\emph{Using successive difference replication for estimating variances}."
#' \strong{Survey Methodology}, Statistics Canada, 40(1), 47-59.
#'
#' @export
#'
#' @examples
#' library(survey)
#' set.seed(2023)
#'
#' # Create an example survey design object
#'
#'   sample_data <- data.frame(
#'     PSU     = c(1,2,3)
#'   )
#'
#'   survey_design <- svydesign(
#'     data = sample_data,
#'     ids = ~ PSU,
#'     weights = ~ 1
#'   )
#'
#'   rep_design <- survey_design |>
#'     as.svrepdesign(type = "JK1", mse = TRUE)
#'
#' # Inspect replicates before subsampling
#'
#'   rep_design |> weights(type = "analysis")
#'
#' # Inspect replicates after adding inactive replicates
#'
#'   rep_design |>
#'     add_inactive_replicates(n_total = 5, location = "first") |>
#'     weights(type = "analysis")
#'
#'   rep_design |>
#'     add_inactive_replicates(n_to_add = 2, location = "last") |>
#'     weights(type = "analysis")
#'
#'   rep_design |>
#'     add_inactive_replicates(n_to_add = 5, location = "random") |>
#'     weights(type = "analysis")
#'
add_inactive_replicates <- function(design, n_total, n_to_add, location = "last") {

  # Check for invalid inputs
  if (!inherits(design, "svyrep.design")) {
    stop("`design` must be a replicate design object.")
  }
  if (!xor(missing(n_total), missing(n_to_add))) {
    stop("Must use either `n_total` or `n_to_add`. Cannot use both.")
  }
  if (!missing(n_total)) {
    if (is.null(n_total) || (length(n_total) > 1) || (!is.numeric(n_total)) || (is.na(n_total))) {
      stop("`n_total` must be a single number")
    }
  }
  if (!missing(n_to_add)) {
    if (is.null(n_to_add) || (length(n_to_add) > 1) || (!is.numeric(n_to_add)) || (is.na(n_to_add))) {
      stop("`n_to_add` must be a single number")
    }
  }
  if ((!is.character(location)) || (is.na(location)) || (length(location) > 1)) {
    stop("`location` must be either 'first', 'last', or 'random'")
  }
  if (!location %in% c("first", "last", "random")) {
    stop("`location` must be either 'first', 'last', or 'random'")
  }

  # Check whether the repweights are stored using compression
  is_compressed <- inherits(design$repweights, "repweights_compressed")

  # Count the number of replicates
  if (!is_compressed) {
    orig_n_reps <- ncol(design$repweights)
  }
  if (is_compressed) {
    orig_n_reps <- ncol(design$repweights[['weights']])
  }

  # Determine number of replicates to add
  if (missing(n_to_add)) {
    n_to_add <- n_total - orig_n_reps
  }
  if (missing(n_total)) {
    n_total <- orig_n_reps + n_to_add
  }

  # If there are already enough replicates, return the original design object
  if (orig_n_reps >= n_total) {
    return(design)
  }



  # Retrieve the (uncompressed) matrix of replicate weights
  orig_rep_wts <- stats::weights(design, type = "replication")

  # Create the inactive replicates
  if (design$combined.weights) {
    inactive_reps <- matrix(
      data = stats::weights(design, type = "sampling"),
      nrow = nrow(design), ncol = n_to_add,
      byrow = FALSE
    )
  } else if (!design$combined.weights) {
    inactive_reps <- matrix(
      data = 1,
      nrow = nrow(design), ncol = n_to_add,
      byrow = FALSE
    )
  }

  # Combine the original and the new replicates
  rep_positions <- seq_len(n_total)
  if (location == "first") {
    new_rep_positions <- seq_len(n_to_add)
    old_rep_positions <- seq(from = n_to_add + 1, to = n_total, by = 1)
  }
  if (location == "last") {
    old_rep_positions <- seq_len(orig_n_reps)
    new_rep_positions <- seq(from = orig_n_reps + 1, to = n_total, by = 1)
  }
  if (location == "random") {
    new_rep_positions <- sample(x = n_total, size = n_to_add, replace = FALSE)
    old_rep_positions <- setdiff(x = seq_len(n_total), y = new_rep_positions)
  }

  updated_rep_wts <- matrix(data = NA_real_, nrow = nrow(design), ncol = n_total)
  updated_rep_wts[,old_rep_positions] <- orig_rep_wts
  updated_rep_wts[,new_rep_positions] <- inactive_reps

  if (is_compressed) {
    design$repweights <- survey::compressWeights(updated_rep_wts)
  }
  if (!is_compressed) {
    design$repweights <- updated_rep_wts
  }

  # Add new 'rscales' elements, as necessary
  updated_rscales <- rep(1, times = n_total)
  updated_rscales[old_rep_positions] <- design$rscales
  design$rscales <- updated_rscales

  if (!design$mse) {
    paste0(
      "The design object has `mse = FALSE`: ",
      "variance estimates may differ before and after adding inactive replicates. ",
      "For details, see `help('add_inactive_replicates', package = 'svrep')`."
    ) |> warning()
  }

  return(design)
}
