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
#' @description Get the scale coefficents
#' used for variance estimation in a replicate design object.
#' @param rep_design A replicate design object
#' @param output Either \code{'separate'} (which is the default)
#' or \code{'combined'}. If \code{'separate'}, the result is a
#' list with elements 'scale' and 'rscales', where 'scale' is the 
#' overall scale coefficient and 'rscales' is a vector of coefficients
#' for each replicate. If \code{'combined'}, the result is a vector
#' of coefficients for each replicate, which incorporates the overall
#' scale coefficient.
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
#' When calling \code{get_rep_scale_coefs(output = 'separate')},
#' the list object includes an element \code{'scale'} corresponding
#' to \eqn{C} and a vector named \code{'rscales'} corresponding
#' to \eqn{c_r, r=1,\dots,R}.
#' 
#' When calling \code{get_rep_scale_coefs(output = 'combined')},
#' the output is a vector with \eqn{R} elements,
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
#'   fpc     = ~ fpc,
#'   nest    = TRUE
#' )
#' 
#' jk_design <- api_design |>
#'   as_random_group_jackknife_design(replicates = 12)
#'
#' jk_design |>
#'   get_rep_scale_coefs('separate')
#' 
#' jk_design |>
#'   get_rep_scale_coefs('combined')
#' 
get_rep_scale_coefs <- function(rep_design, output = "separate") {
  if (!inherits(rep_design, 'svyrep.design')) {
    stop("Must be a replicate design object with class 'svyrep.design'")
  }
  if (output == "separate") {
    return(list('scale'  = rep_design[['scale']],
                'rscales' = rep_design[['rscales']]))
  }
  if (output == "combined") {
    return(rep_design[['scale']] * rep_design[['rscales']])
  }
}