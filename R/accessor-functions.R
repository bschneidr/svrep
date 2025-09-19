#' @title Access Type of Replication Method
#' @description Identify the type of replication
#' method used for the replicates in a 
#' replicate design object.
#' @param rep_design A replicate design object
#' @returns A character string giving the
#' type of replication method.
#' @export
#' @examples
#' data('scd', package = 'survey')
#' 
#' scd_design <- svydesign(
#'   data   = scd, 
#'   id     = ~ ambulance, 
#'   prob   = ~ 1,
#'   strata = ~ ESA,
#'   nest   = TRUE
#' )
#'
#' scd_design |>
#'   as_bootstrap_design(replicates = 5) |>
#'   get_rep_type()
#' 
#' scd_design |>
#'   as_fays_gen_rep_design(
#'     variance_estimator = "Ultimate Cluster"
#'   ) |>
#'   get_rep_type()
#' 
get_rep_type <- function(rep_design) {
  if (!inherits(rep_design, 'svyrep.design')) {
    stop("Must be a replicate design object with class 'svyrep.design'")
  }
  rep_design[['type']]
}

#' @title Access Replication Scale Coefficients
#' @description Get the scale coefficients
#' used for variance estimation in a replicate design object.
#' @inheritParams get_rep_type
#' @param type Either \code{'overall'}, \code{'specific'}, 
#' or \code{'combined'}. See the details section below.
#' The result for \code{'overall'} is the single overall
#' scale coefficient. The result for \code{'specific'} is the vector of
#' replicate-specific coefficients. The result for \code{'combined'}
#' is the product of the overall and replicate-specific coefficients.
#' @returns If \code{type = 'overall'}, the result is a single number.
#' Otherwise, the result is a vector with length matching
#' the number of replicates.
#' @details
#' For a statistic \eqn{\hat{\theta}},
#' replication methods estimate the sampling variance
#' using \eqn{R} replicate estimates, with the estimate
#' for the \eqn{r}-th replicate denoted \eqn{\hat{\theta}_r}.
#' 
#' The formula for the variance estimate is the following:
#' \deqn{
#'   v(\hat{\theta}) = C \sum_{r=1}^{R} c_r (\hat{\theta_r} - \hat{\theta})^2
#' }
#' 
#' The terms \eqn{C} and \eqn{c_r, r=1,\dots,R} are scale coefficients.
#' \eqn{C} is the overall coefficient,
#' and \eqn{c_r, r=1,\dots,R} are replicate-specific coefficients.
#' 
#' Specifying \code{get_rep_scale_coefs(type='overall')} returns
#' the overall coefficient \eqn{C}. Specifying \code{type='specific'}
#' returns the replicate-specific coefficients \eqn{c_r, r=1,\dots,R}.
#' 
#' Specifying \code{type='combined'} returns a vector with 
#' \eqn{R} elements,
#' where the \eqn{r}-th element is \eqn{C \times c_r}.
#' @export
#' @examples
#' data('api', package = 'survey')
#' 
#' api_design <- svydesign(
#'   data    = apistrat, 
#'   id      = ~ 1, 
#'   strata  = ~ stype,
#'   weights = ~ pw,
#'   nest    = TRUE
#' )
#' 
#' jk_design <- api_design |>
#'   as_random_group_jackknife_design(
#'     replicates = 12
#'   )
#'
#' jk_design |>
#'   get_rep_scale_coefs('overall')
#' 
#' jk_design |>
#'   get_rep_scale_coefs('specific')
#' 
#' jk_design |>
#'   get_rep_scale_coefs('combined')
#' 
get_rep_scale_coefs <- function(rep_design, type = "combined") {
  if (!inherits(rep_design, 'svyrep.design')) {
    stop("Must be a replicate design object with class 'svyrep.design'")
  }
  if (length(type) != 1 || (!type %in% c("overall", "specific", "combined"))) {
    stop("`type` must be 'overall', 'specific', or 'combined'")
  }
  result <- switch(type,
    'overall' = rep_design[['scale']],
    'specific' = rep_design[['rscales']],
    'combined' = rep_design[['scale']] * rep_design[['rscales']]
  )
  return(result)
}
