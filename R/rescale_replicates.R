#' @title Rescale replicate factors
#' @description Rescale replicate factors.
#' The main application of this rescaling is to ensure
#' that all replicate weights are strictly positive.
#'
#' Note that this rescaling has no impact on variance estimates for totals (or other linear statistics),
#' but variance estimates for nonlinear statistics will be affected by the rescaling.
#'
#' @param x Either a replicate survey design object,
#' or a numeric matrix of replicate weights.
#' If \code{x} is a matrix, then it must have an attribute
#' \code{scale} specifying the scale factor for variance estimation.
#' @param new_scale Either a single positive number, or \code{NULL}.
#' If supplied, \code{new_scale} will be the new scale factor
#' used for estimating variances. For example, with \code{B=10} bootstrap replicates,
#' specifying \code{new_scale=0.2} will change the scale factor for variance estimation
#' from \eqn{B^{-1}=0.1} to \eqn{0.2}. \cr
#' If \code{new_scale=NULL} or is left unspecified, 
#' then the argument \code{min_wgt} should be used instead.
#' @param min_wgt Should only be used if \code{new_scale=NULL} or \code{new_scale} is left unspecified. 
#' Specifies the minimum acceptable value for the rescaled weights,
#' which will be used to automatically determine the new scale.
#' Must be at least zero and must be less than one.
#' @param digits Only used if the argument \code{min_wgt} is used. 
#' Specifies the number of decimal places
#' to use for the scale adjustment factor,
#' so that the ratio \code{new_scale/orig_scale}
#' will be a number with at most \code{digits} decimal places.
#' For example, let \code{orig_scale} denote the original
#' scale and \code{new_scale} denote the new scale.
#' Then if the user specifies \code{digits=2},
#' the ratio \code{new_scale/orig_scale} will
#' be a number such as 0.25 or 0.10.
#' This results in documentation that is easier to read.
#'
#' @return If the input is a numeric matrix, returns the rescaled matrix.
#' If the input is a replicate survey design object, returns an updated replicate survey design object.
#'
#' For a replicate survey design object, results depend on
#' whether the object has a matrix of replicate factors rather than
#' a matrix of replicate weights (which are the product of replicate factors and sampling weights).
#' If the design object has \code{combined.weights=FALSE},
#' then the replication factors are adjusted.
#' If the design object has \code{combined.weights=TRUE},
#' then the replicate weights are adjusted. It is strongly
#' recommended to only use the rescaling method for replication factors
#' rather than the weights.
#'
#' For a replicate survey design object, the \code{scale} element
#' of the design object will be updated appropriately.
#' @details
#' Let \eqn{\mathbf{A} = \left[ \mathbf{a}^{(1)} \cdots \mathbf{a}^{(b)} \cdots \mathbf{a}^{(B)} \right]} denote the \eqn{(n \times B)} matrix of replicate adjustment factors.
#' The overall scale factor \eqn{C} is used for estimating variances
#' with the formula \deqn{v(\hat{y}) = C \sum_{b=1}^{B} (\hat{y}_b - \hat{y})^2}
#' 
#' The scale factor is changed from \eqn{C} to \eqn{C^{\prime}}
#' by rescaling the replicate factor matrix \eqn{\mathbf{A}} using the transformation 
#' \deqn{
#'   1 + \sqrt{\frac{C}{C^{\prime}}}(\mathbf{A}-1)
#' }
#' @references
#' Rescaling was suggested by Fay (1989) for the specific application
#' of creating replicate factors using his generalized replication method.
#' This kind of rescaling is commony used in balanced repeated replication
#' to implement Fay's method of balanced repeated replication.
#' Beaumont and Patak (2012) provided an extended discussion of rescaling
#' methods in the context of rescaling generalized bootstrap replication factors
#' to avoid negative replicate weights.
#'
#' - Beaumont, Jean-François, and Zdenek Patak. 2012.
#' "On the Generalized Bootstrap for Sample Surveys with Special Attention to Poisson Sampling: Generalized Bootstrap for Sample Surveys."
#' International Statistical Review 80 (1): 127-48.
#' https://doi.org/10.1111/j.1751-5823.2011.00166.x.
#' \cr \cr
#' - Fay, Robert. 1989. "Theory And Application Of Replicate Weighting For Variance Calculations."
#' In, 495-500. Alexandria, VA: American Statistical Association.
#' http://www.asasrms.org/Proceedings/papers/1989_033.pdf
#'
#' @export
#'
#' @examples
#' # Example 1: Rescaling a matrix of replicate weights to avoid negative weights
#'
#'  rep_wgts <- matrix(
#'    c(1.69742746694909, -0.230761178913411, 1.53333377634192,
#'      0.0495043413294782, 1.81820367441039, 1.13229198793703,
#'      1.62482013925955, 1.0866133494029, 0.28856654131668,
#'      0.581930729719006, 0.91827012312825, 1.49979905894482,
#'      1.26281337410693, 1.99327362761477, -0.25608700039304),
#'    nrow = 3, ncol = 5
#'  )
#'  attr(rep_wgts, 'scale') <- 1/5
#'
#'  rescaled_wgts <- rescale_replicates(rep_wgts, min_wgt = 0.01)
#'
#'  print(rep_wgts)
#'  print(rescaled_wgts)
#'  
#'  # Example 2: Rescaling replicate weights with a specified value of 'tau'
#'  
#'  rescaled_wgts <- rescale_replicates(rep_wgts, new_scale = 1/10)
#'  print(rescaled_wgts)
#'
#'  # Example 3: Rescaling replicate weights of a survey design object
#'  set.seed(2023)
#'  library(survey)
#'  data('mu284', package = 'survey')
#'
#'  ## First create a bootstrap design object
#'  svy_design_object <- svydesign(
#'    data = mu284,
#'    ids = ~ id1 + id2,
#'    fpc = ~  n1 +  n2
#'  )
#'
#'  boot_design <- as_gen_boot_design(
#'    design = svy_design_object,
#'    variance_estimator = "Stratified Multistage SRS",
#'    replicates = 5
#'  )
#'
#'  ## Rescale the weights
#'  rescaled_boot_design <- boot_design |>
#'    rescale_replicates(min_wgt = 0.01)
#'
#'  boot_wgts <- weights(boot_design, "analysis")
#'  rescaled_boot_wgts <- weights(rescaled_boot_design, 'analysis')
#'
#'  print(boot_wgts)
#'  print(rescaled_boot_wgts)
rescale_replicates <- function(x, new_scale = NULL, min_wgt = 0.01, digits = 2) {
  
  if (!is.null(new_scale) && missing(min_wgt)) {
    min_wgt <- NULL
  }
  if (!xor(!is.null(new_scale), !is.null(min_wgt))) {
    stop("Can only specify `new_scale` or `min_wgt`: cannot use both arguments.")
  }
  if (!is.null(new_scale)) {
    if (!is.numeric(new_scale) || any(is.na(new_scale)) || (length(new_scale) > 1) || (new_scale < 0)) {
      stop("If `new_scale` is specified, it must be a single positive number.")
    }
    if (!missing(digits) && !is.null(digits)) {
      message('Arguments `digits` is ignored when using `tau` instead of `min_wgt`.')
    }
  }
  if ((!is.null(min_wgt)) && (!is.numeric(min_wgt) || length(min_wgt) != 1 || is.na(min_wgt) || min_wgt < 0 || min_wgt > 1)) {
    stop("When `tau=NULL`, the argument `min_wgt` must be at least 0 and less than 1.")
  }
  if (!missing(digits)) {
    if (!is.numeric(digits) || (length(digits) != 1) || is.na(digits) || (digits < 1)) {
      stop("`digits` must be an integer greater than or equal to 1.")
    }
  }
  if (!is.null(min_wgt) && (min_wgt != 0) && (round(min_wgt, digits) == 0)) {
    stop("round(min_wgt, digits) equals 0; increase either `min_wgt` or `digits`.")
  }
  
  UseMethod("rescale_replicates", x)
}

#' @export
rescale_replicates.matrix <- function(x, new_scale = NULL, min_wgt = 0.01, digits = 2) {
  
  rep_weights <- x
  orig_scale <- attr(rep_weights, 'scale')

  if (is.null(orig_scale)) {
    stop("The matrix `x` must have an attribute named 'scale', giving the overall scale factor for variance estimation.")
  }
  
  if ((is.null(new_scale)) && any(rep_weights < min_wgt)) {
    rescaling_constant <- min((1-rep_weights)/(min_wgt-1))
    rescaling_constant <- abs(rescaling_constant)
    rescaling_constant <- ceiling(rescaling_constant^2 * 10^(digits))/(10^(digits))
    new_scale <- orig_scale * rescaling_constant
  }
  if ((is.null(new_scale)) && !any(rep_weights < min_wgt)) {
    rescaling_constant   <- 1
    rescaled_rep_weights <- rep_weights
  }
  if (!is.null(new_scale)) {
    if (!is.numeric(orig_scale) || any(is.na(orig_scale)) || (length(orig_scale) > 1) || (orig_scale < 0)) {
      stop("The scale attribute of `x` must be a single positive number.")
    }
    rescaling_constant <- sqrt(orig_scale / new_scale)
  }

  rescaled_rep_weights <- 1 + (rep_weights - 1)*rescaling_constant
  
  attr(rescaled_rep_weights, 'scale') <- new_scale
  attr(rescaled_rep_weights, 'tau') <- NULL

  return(rescaled_rep_weights)
}

#' @export
rescale_replicates.svyrep.design <- function(x, new_scale = NULL, min_wgt = 0.01, digits = 2) {
  
  rep_weights <- weights(x, type = "replication")
  attr(rep_weights, 'scale') <- x$scale
  rescaled_rep_weights <- rescale_replicates.matrix(
    x = rep_weights,
    new_scale = new_scale, 
    min_wgt = min_wgt, digits = digits
  )
  x$repweights <- rescaled_rep_weights
  x$scale <- attr(rescaled_rep_weights, 'scale')
  attr(rescaled_rep_weights, 'tau') <- NULL
  return(x)
}
#' @title Rescale replicate factors
#' @description 
#' \emph{Deprecated}: Please use the function
#' \code{rescale_replicates} instead.
#' 
#' Rescale replicate factors.
#' The main application of this rescaling is to ensure
#' that all replicate weights are strictly positive.
#'
#' Note that this rescaling has no impact on variance estimates for totals (or other linear statistics),
#' but variance estimates for nonlinear statistics will be affected by the rescaling.
#'
#' @param x Either a replicate survey design object,
#' or a numeric matrix of replicate weights.
#' @param tau Either a single positive number, or \code{NULL}. 
#' This is the rescaling constant \eqn{\tau} used in the transformation 
#' \eqn{\frac{w + \tau - 1}{\tau}},
#' where \eqn{w} is the original weight. \cr
#' If \code{tau=NULL} or is left unspecified, 
#' then the argument \code{min_wgt} should be used instead,
#' in which case, \eqn{\tau} is automatically set to the smallest value needed to rescale
#' the replicate weights such that they are all at least \code{min_wgt}.
#' @param min_wgt Should only be used if \code{tau=NULL} or \code{tau} is left unspecified. 
#' Specifies the minimum acceptable value for the rescaled weights,
#' which will be used to automatically determine the value \eqn{\tau} used
#' in the transformation \eqn{\frac{w + \tau - 1}{\tau}},
#' where \eqn{w} is the original weight.
#' Must be at least zero and must be less than one.
#' @param digits Only used if the argument \code{min_wgt} is used. 
#' Specifies the number of decimal places
#' to use for choosing \code{tau}. Using a smaller number of \code{digits}
#' is useful simply for producing easier-to-read documentation.
#'
#' @return If the input is a numeric matrix, returns the rescaled matrix.
#' If the input is a replicate survey design object, returns an updated replicate survey design object.
#'
#' For a replicate survey design object, results depend on
#' whether the object has a matrix of replicate factors rather than
#' a matrix of replicate weights (which are the product of replicate factors and sampling weights).
#' If the design object has \code{combined.weights=FALSE},
#' then the replication factors are adjusted.
#' If the design object has \code{combined.weights=TRUE},
#' then the replicate weights are adjusted. It is strongly
#' recommended to only use the rescaling method for replication factors
#' rather than the weights.
#'
#' For a replicate survey design object, the \code{scale} element
#' of the design object will be updated appropriately,
#' and an element \code{tau} will also be added.
#' If the input is a matrix instead of a survey design object,
#' the result matrix will have an attribute named \code{tau}
#' which can be retrieved using \code{attr(x, 'tau')}.
#' @details
#' Let \eqn{\mathbf{A} = \left[ \mathbf{a}^{(1)} \cdots \mathbf{a}^{(b)} \cdots \mathbf{a}^{(B)} \right]} denote the \eqn{(n \times B)} matrix of replicate adjustment factors.
#' To eliminate negative adjustment factors, Beaumont and Patak (2012) propose forming a rescaled matrix of nonnegative replicate factors \eqn{\mathbf{A}^S} by rescaling each adjustment factor \eqn{a_k^{(b)}} as follows:
#' \deqn{
#'    a_k^{S,(b)} = \frac{a_k^{(b)} + \tau - 1}{\tau}
#'  }
#' where \eqn{\tau \geq 1 - a_k^{(b)} \geq 1} for all \eqn{k} in \eqn{\left\{ 1,\ldots,n \right\}} and all \eqn{b} in \eqn{\left\{1, \ldots, B\right\}}.
#'
#' The value of \eqn{\tau} can be set based on the realized adjustment factor matrix \eqn{\mathbf{A}} or by choosing \eqn{\tau} prior to generating the adjustment factor matrix \eqn{\mathbf{A}} so that \eqn{\tau} is likely to be large enough to prevent negative adjustment factors.
#'
#' If the adjustment factors are rescaled in this manner, it is important to adjust the scale factor used in estimating the variance with the bootstrap replicates.
#' For example, for bootstrap replicates, the adjustment factor becomes \eqn{\frac{\tau^2}{B}} instead of \eqn{\frac{1}{B}}.
#' \deqn{
#'  \textbf{Prior to rescaling: } v_B\left(\hat{T}_y\right) = \frac{1}{B}\sum_{b=1}^B\left(\hat{T}_y^{*(b)}-\hat{T}_y\right)^2
#'  }
#' \deqn{
#'  \textbf{After rescaling: } v_B\left(\hat{T}_y\right) = \frac{\tau^2}{B}\sum_{b=1}^B\left(\hat{T}_y^{S*(b)}-\hat{T}_y\right)^2
#' }
#' @references
#' This method was suggested by Fay (1989) for the specific application
#' of creating replicate factors using his generalized replication method.
#' Beaumont and Patak (2012) provided an extended discussion on this rescaling
#' method in the context of rescaling generalized bootstrap replication factors
#' to avoid negative replicate weights.
#' 
#' The notation used in this documentation is taken from Beaumont and Patak (2012).
#'
#' - Beaumont, Jean-François, and Zdenek Patak. 2012.
#' "On the Generalized Bootstrap for Sample Surveys with Special Attention to Poisson Sampling: Generalized Bootstrap for Sample Surveys."
#' International Statistical Review 80 (1): 127-48.
#' https://doi.org/10.1111/j.1751-5823.2011.00166.x.
#' \cr \cr
#' - Fay, Robert. 1989. "Theory And Application Of Replicate Weighting For Variance Calculations."
#' In, 495-500. Alexandria, VA: American Statistical Association.
#' http://www.asasrms.org/Proceedings/papers/1989_033.pdf
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Rescaling a matrix of replicate weights to avoid negative weights
#'
#'  rep_wgts <- matrix(
#'    c(1.69742746694909, -0.230761178913411, 1.53333377634192,
#'      0.0495043413294782, 1.81820367441039, 1.13229198793703,
#'      1.62482013925955, 1.0866133494029, 0.28856654131668,
#'      0.581930729719006, 0.91827012312825, 1.49979905894482,
#'      1.26281337410693, 1.99327362761477, -0.25608700039304),
#'    nrow = 3, ncol = 5
#'  )
#'
#'  rescaled_wgts <- rescale_reps(rep_wgts, min_wgt = 0.01)
#'
#'  print(rep_wgts)
#'  print(rescaled_wgts)
#'  
#'  # Example 2: Rescaling replicate weights with a specified value of 'tau'
#'  
#'  rescaled_wgts <- rescale_reps(rep_wgts, tau = 2)
#'  print(rescaled_wgts)
#'
#'  # Example 3: Rescaling replicate weights of a survey design object
#'  set.seed(2023)
#'  library(survey)
#'  data('mu284', package = 'survey')
#'
#'  ## First create a bootstrap design object
#'  svy_design_object <- svydesign(
#'    data = mu284,
#'    ids = ~ id1 + id2,
#'    fpc = ~ n1 + n2
#'  )
#'
#'  boot_design <- as_gen_boot_design(
#'    design = svy_design_object,
#'    variance_estimator = "Stratified Multistage SRS",
#'    replicates = 5, tau = 1
#'  )
#'
#'  ## Rescale the weights
#'  rescaled_boot_design <- boot_design |>
#'    rescale_reps(min_wgt = 0.01)
#'
#'  boot_wgts <- weights(boot_design, "analysis")
#'  rescaled_boot_wgts <- weights(rescaled_boot_design, 'analysis')
#'
#'  print(boot_wgts)
#'  print(rescaled_boot_wgts)
#' }
rescale_reps <- function(x, tau = NULL, min_wgt = 0.01, digits = 2) {
  
  warning(paste(
    "`rescale_reps()` is deprecated and will be removed in a future version of 'svrep'.",
    "Please use `rescale_replicates()` instead."
  ))

  if (!is.null(tau) && missing(min_wgt)) {
    min_wgt <- NULL
  }
  if (!xor(!is.null(tau), !is.null(min_wgt))) {
    stop("Can only specify `tau` or `min_wgt`: cannot use both arguments.")
  }
  if (!is.null(tau)) {
    if (!is.numeric(tau) || any(is.na(tau)) || (length(tau) > 1) || (tau < 0)) {
      stop("If `tau` is specified, it must be a single positive number.")
    }
    if (!missing(digits) && !is.null(digits)) {
      message('Arguments `digits` is ignored when using `tau` instead of `min_wgt`.')
    }
  }
  if ((!is.null(min_wgt)) && (!is.numeric(min_wgt) || length(min_wgt) != 1 || is.na(min_wgt) || min_wgt < 0 || min_wgt > 1)) {
    stop("When `tau=NULL`, the argument `min_wgt` must be at least 0 and less than 1.")
  }
  if (!missing(digits)) {
    if (!is.numeric(digits) || (length(digits) != 1) || is.na(digits) || (digits < 1)) {
      stop("`digits` must be an integer greater than or equal to 1.")
    }
  }
  if (!is.null(min_wgt) && (min_wgt != 0) && (round(min_wgt, digits) == 0)) {
    stop("round(min_wgt, digits) equals 0; increase either `min_wgt` or `digits`.")
  }
  
  UseMethod("rescale_reps", x)
}

#' @export
rescale_reps.matrix <- function(x, tau = NULL, min_wgt = 0.01, digits = 2) {
  rep_weights <- x
  
  if ((is.null(tau)) && any(rep_weights < min_wgt)) {
    rescaling_constant <- min((1-rep_weights)/(min_wgt-1))
    rescaling_constant <- abs(rescaling_constant)
    rescaling_constant <- ceiling(rescaling_constant * 10^(digits))/(10^(digits))
  }
  if ((is.null(tau)) && !any(rep_weights < min_wgt)) {
    rescaling_constant <- 1
    rescaled_rep_weights <- rep_weights
  }
  if (!is.null(tau)) {
    rescaling_constant <- tau
  }
  rescaled_rep_weights <- (rep_weights + (rescaling_constant-1))/rescaling_constant
  
  attr(rescaled_rep_weights, 'tau') <- rescaling_constant
  orig_scale <- attr(x, "scale")
  if (!is.null(orig_scale)) {
    attr(rescaled_rep_weights, 'scale') <- (rescaling_constant^2) * orig_scale
  }
  return(rescaled_rep_weights)
}

#' @export
rescale_reps.svyrep.design <- function(x, tau = NULL, min_wgt = 0.01, digits = 2) {
  
  rep_weights <- weights(x, type = "replication")
  rescaled_rep_weights <- rescale_reps.matrix(x = rep_weights,
                                              tau = tau, min_wgt = min_wgt,
                                              digits = digits)
  attr(rescaled_rep_weights, 'scale') <- (attr(rescaled_rep_weights, 'tau')^2) * x$scale
  x$scale <- attr(rescaled_rep_weights, 'scale')
  x$tau <- attr(rescaled_rep_weights, 'tau')
  x$repweights <- rescaled_rep_weights
  return(x)
}
