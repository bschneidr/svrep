#' @title Rescale replicate factors to have specified lower bound
#' @description Rescale replicate factors to ensure that they all exceed
#' a specified lower bound. The main use of this rescaling is to ensure
#' that all replicate weights are strictly positive.
#'
#' Note that this rescaling has no impact on variance estimates for totals,
#' but variance estimates for nonlinear statistics will be affected by the rescaling.
#'
#' @param x Either a replicate survey design object,
#' or a numeric matrix of replicate weights.
#' @param tau Either \code{"auto"}, or a single number. This is the rescaling constant
#' used to avoid negative weights through the transformation \eqn{\frac{w + \tau - 1}{\tau}},
#' where \eqn{w} is the original weight and \eqn{\tau} is the rescaling constant \code{tau}. \cr
#' If \code{tau="auto"}, the rescaling factor is determined automatically as follows:
#' if all of the adjustment factors exceed the minimum value \code{min_wgt},
#'  then \code{tau} is set equal to 1;
#' otherwise, \code{tau} is set to the smallest value needed to rescale
#' the adjustment factors such that they are all at least \code{min_wgt}.
#' @param min_wgt Only used if \code{tau='auto'}. Specifies the minimum acceptable value for the rescaled weights.
#' Must be at least zero and must be less than one.
#' @param digits Only used if \code{tau='auto'}. Specifies the number of decimal places
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
#' then the replicate weights are adjusted.
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
#'
#' @export
#'
#' @examples
#' # Example 1: Rescaling a matrix of replicate weights
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
#'  rescaled_wgts <- rescale_reps(rep_wgts, tau = 'auto', min_wgt = 0.01)
#'
#'  print(rep_wgts)
#'  print(rescaled_wgts)
#'
#'  # Example 2: Rescaling replicate weights of a survey design object
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
#'    rescale_reps(tau = 'auto', min_wgt = 0.01)
#'
#'  boot_wgts <- weights(boot_design, "analysis")
#'  rescaled_boot_wgts <- weights(rescaled_boot_design, 'analysis')
#'
#'  print(boot_wgts)
#'  print(rescaled_boot_wgts)
rescale_reps <- function(x, tau = "auto", min_wgt = 0.01, digits = 2) {

  if (length(tau) != 1 || is.na(tau) || (tau != "auto" & !is.numeric(tau)) || (is.numeric(tau) & tau < 0)) {
    stop("`tau` must be either 'auto' or a single positive number.")
  }
  if ((tau == "auto") && (min_wgt < 0 || min_wgt > 1)) {
    stop("When `tau='auto'`, the argument `min_wgt` must be at least 0 and less than 1.")
  }

  UseMethod("rescale_reps", x)
}

#' @export
rescale_reps.matrix <- function(x, tau = "auto", min_wgt = 0.01, digits = 2) {
  rep_weights <- x
  if (any(rep_weights < min_wgt)) {
    if (tau == "auto") {
      rescaling_constant <- min((1-rep_weights)/(min_wgt-1))
      rescaling_constant <- abs(rescaling_constant)
      rescaling_constant <- ceiling(rescaling_constant * 100)/100
    } else {
      rescaling_constant <- tau
    }
    rescaled_rep_weights <- (rep_weights + (rescaling_constant-1))/rescaling_constant
  } else {
    rescaling_constant <- 1
    rescaled_rep_weights <- rep_weights
  }
  attr(rescaled_rep_weights, 'tau') <- rescaling_constant
  orig_scale <- attr(x, "scale")
  if (!is.null(orig_scale)) {
    attr(rescaled_rep_weights, 'scale') <- (rescaling_constant^2) * orig_scale
  }
  return(rescaled_rep_weights)
}

#' @export
rescale_reps.svyrep.design <- function(x, tau = "auto", min_wgt = 0.01, digits = 2) {

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

#' @title Creates replicate factors for the generalized survey bootstrap
#' @description Creates replicate factors for the generalized survey bootstrap method.
#' The generalized survey bootstrap is a method for forming bootstrap replicate weights
#' from a textbook variance estimator, provided that the variance estimator
#' can be represented as a quadratic form whose matrix is positive semidefinite
#' (this covers a large class of variance estimators).
#' @param Sigma The matrix of the quadratic form used to represent the variance estimator.
#' Must be positive semidefinite.
#' @param num_replicates The number of bootstrap replicates to create.
#' @param tau Either \code{"auto"}, or a single number. This is the rescaling constant
#' used to avoid negative weights through the transformation \eqn{\frac{w + \tau - 1}{\tau}},
#' where \eqn{w} is the original weight and \eqn{\tau} is the rescaling constant \code{tau}. \cr
#' If \code{tau="auto"}, the rescaling factor is determined automatically as follows:
#' if all of the adjustment factors are nonnegative, then \code{tau} is set equal to 1;
#' otherwise, \code{tau} is set to the smallest value needed to rescale
#' the adjustment factors such that they are all at least \code{0.01}.
#' @param exact_vcov If \code{exact_vcov=TRUE}, the replicate factors will be generated
#' such that their variance-covariance matrix exactly matches the target variance estimator's
#' quadratic form (within numeric precision).
#' This is desirable as it causes variance estimates for totals to closely match
#' the values from the target variance estimator.
#' This requires that \code{num_replicates} exceeds the rank of \code{Sigma}.
#' The replicate factors are generated by applying PCA-whitening to a collection of draws
#' from a multivariate Normal distribution, then applying a coloring transformation
#' to the whitened collection of draws.
#'
#' @section Statistical Details:
#' Let \eqn{v( \hat{T_y})} be the textbook variance estimator for an estimated population total \eqn{\hat{T}_y} of some variable \eqn{y}.
#' The base weight for case \eqn{i} in our sample is \eqn{w_i}, and we let \eqn{\breve{y}_i} denote the weighted value \eqn{w_iy_i}.
#' Suppose we can represent our textbook variance estimator as a quadratic form: \eqn{v(\hat{T}_y) = \breve{y}\Sigma\breve{y}^T},
#' for some \eqn{n \times n} matrix \eqn{\Sigma}.
#' The only constraint on \eqn{\Sigma} is that, for our sample, it must be symmetric and positive semidefinite.
#'
#' The bootstrapping process creates \eqn{B} sets of replicate weights, where the \eqn{b}-th set of replicate weights is a vector of length \eqn{n} denoted \eqn{\mathbf{a}^{(b)}}, whose \eqn{k}-th value is denoted \eqn{a_k^{(b)}}.
#' This yields \eqn{B} replicate estimates of the population total, \eqn{\hat{T}_y^{*(b)}=\sum_{k \in s} a_k^{(b)} \breve{y}_k}, for \eqn{b=1, \ldots B}, which can be used to estimate sampling variance.
#'
#' \deqn{
#'   v_B\left(\hat{T}_y\right)=\frac{\sum_{b=1}^B\left(\hat{T}_y^{*(b)}-\hat{T}_y\right)^2}{B}
#' }
#'
#' This bootstrap variance estimator can be written as a quadratic form:
#'
#'   \deqn{
#'     v_B\left(\hat{T}_y\right) =\mathbf{\breve{y}}^{\prime}\Sigma_B \mathbf{\breve{y}}
#'   }
#'   where
#'   \deqn{
#'     \boldsymbol{\Sigma}_B = \frac{\sum_{b=1}^B\left(\mathbf{a}^{(b)}-\mathbf{1}_n\right)\left(\mathbf{a}^{(b)}-\mathbf{1}_n\right)^{\prime}}{B}
#'   }
#'
#' Note that if the vector of adjustment factors \eqn{\mathbf{a}^{(b)}} has expectation \eqn{\mathbf{1}_n} and variance-covariance matrix \eqn{\boldsymbol{\Sigma}},
#' then we have the bootstrap expectation \eqn{E_{*}\left( \boldsymbol{\Sigma}_B \right) = \boldsymbol{\Sigma}}. Since the bootstrap process takes the sample values \eqn{\breve{y}} as fixed, the bootstrap expectation of the variance estimator is \eqn{E_{*} \left( \mathbf{\breve{y}}^{\prime}\Sigma_B \mathbf{\breve{y}}\right)= \mathbf{\breve{y}}^{\prime}\Sigma \mathbf{\breve{y}}}.
#' Thus, we can produce a bootstrap variance estimator with the same expectation as the textbook variance estimator simply by randomly generating \eqn{\mathbf{a}^{(b)}} from a distribution with the following two conditions:
#' \cr \cr
#'     \strong{Condition 1}: \eqn{\quad \mathbf{E}_*(\mathbf{a})=\mathbf{1}_n}
#' \cr
#'     \strong{Condition 2}: \eqn{\quad \mathbf{E}_*\left(\mathbf{a}-\mathbf{1}_n\right)\left(\mathbf{a}-\mathbf{1}_n\right)^{\prime}=\mathbf{\Sigma}}
#' \cr \cr
#' While there are multiple ways to generate adjustment factors satisfying these conditions,
#' the simplest general method is to simulate from a multivariate normal distribution: \eqn{\mathbf{a} \sim MVN(\mathbf{1}_n, \boldsymbol{\Sigma})}.
#' This is the method used by this function.
#'
#' @section Details on Rescaling to Avoid Negative Adjustment Factors:
#' Let \eqn{\mathbf{A} = \left[ \mathbf{a}^{(1)} \cdots \mathbf{a}^{(b)} \cdots \mathbf{a}^{(B)} \right]} denote the \eqn{(n \times B)} matrix of bootstrap adjustment factors.
#' To eliminate negative adjustment factors, Beaumont and Patak (2012) propose forming a rescaled matrix of nonnegative replicate factors \eqn{\mathbf{A}^S} by rescaling each adjustment factor \eqn{a_k^{(b)}} as follows:
#' \deqn{
#'    a_k^{S,(b)} = \frac{a_k^{(b)} + \tau - 1}{\tau}
#'  }
#' where \eqn{\tau \geq 1 - a_k^{(b)} \geq 1} for all \eqn{k} in \eqn{\left\{ 1,\ldots,n \right\}} and all \eqn{b} in \eqn{\left\{1, \ldots, B\right\}}.
#'
#' The value of \eqn{\tau} can be set based on the realized adjustment factor matrix \eqn{\mathbf{A}} or by choosing \eqn{\tau} prior to generating the adjustment factor matrix \eqn{\mathbf{A}} so that \eqn{\tau} is likely to be large enough to prevent negative bootstrap weights.
#'
#' If the adjustment factors are rescaled in this manner, it is important to adjust the scale factor used in estimating the variance with the bootstrap replicates, which becomes \eqn{\frac{\tau^2}{B}} instead of \eqn{\frac{1}{B}}.
#' \deqn{
#'  \textbf{Prior to rescaling: } v_B\left(\hat{T}_y\right) = \frac{1}{B}\sum_{b=1}^B\left(\hat{T}_y^{*(b)}-\hat{T}_y\right)^2
#'  }
#' \deqn{
#'  \textbf{After rescaling: } v_B\left(\hat{T}_y\right) = \frac{\tau^2}{B}\sum_{b=1}^B\left(\hat{T}_y^{S*(b)}-\hat{T}_y\right)^2
#' }
#' When sharing a dataset that uses rescaled weights from a generalized survey bootstrap, the documentation for the dataset should instruct the user to use replication scale factor \eqn{\frac{\tau^2}{B}} rather than \eqn{\frac{1}{B}} when estimating sampling variances.
#'
#' @references
#' The generalized survey bootstrap was first proposed by Bertail and Combris (1997).
#' See Beaumont and Patak (2012) for a clear overview of the generalized survey bootstrap.
#' The generalized survey bootstrap represents one strategy for forming replication variance estimators
#' in the general framework proposed by Fay (1984) and Dippo, Fay, and Morganstein (1984).
#' \cr \cr
#' - Beaumont, Jean-François, and Zdenek Patak. 2012. “On the Generalized Bootstrap for Sample Surveys with Special Attention to Poisson Sampling: Generalized Bootstrap for Sample Surveys.” International Statistical Review 80 (1): 127–48. https://doi.org/10.1111/j.1751-5823.2011.00166.x.
#' \cr \cr
#' - Bertail, and Combris. 1997. “Bootstrap Généralisé d’un Sondage.” Annales d’Économie Et de Statistique, no. 46: 49. https://doi.org/10.2307/20076068.
#' \cr \cr
#' - Dippo, Cathryn, Robert Fay, and David Morganstein. 1984. “Computing Variances from Complex Samples with Replicate Weights.” In, 489–94. Alexandria, VA: American Statistical Association. http://www.asasrms.org/Proceedings/papers/1984_094.pdf.
#' \cr \cr
#' - Fay, Robert. 1984. “Some Properties of Estimates of Variance Based on Replication Methods.” In, 495–500. Alexandria, VA: American Statistical Association. http://www.asasrms.org/Proceedings/papers/1984_095.pdf.
#' \cr
#'
#' @seealso The function \code{\link[svrep]{make_quad_form_matrix}} can be used to
#' represent several common variance estimators as a quadratic form's matrix,
#' which can then be used as an input to \code{make_gen_boot_factors()}.
#'
#' @return A matrix with the same number of rows as \code{Sigma}, and the number of columns
#' equal to \code{num_replicates}. The object has an attribute named \code{tau} which can be retrieved
#' by calling \code{attr(which = 'tau')} on the object. The value \code{tau} is a rescaling factor
#' which was used to avoid negative weights.
#'
#' In addition, the object has attributes named \code{scale} and \code{rscales} which can be
#' passed directly to \link[survey]{svrepdesign}. Note that the value of \code{scale} is \eqn{\tau^2/B},
#' while the value of \code{rscales} is vector of length \eqn{B}, with every entry equal to \eqn{1}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   library(survey)
#'
#' # Load an example dataset that uses unequal probability sampling ----
#'   data('election', package = 'survey')
#'
#' # Create matrix to represent the Horvitz-Thompson estimator as a quadratic form ----
#'   n <- nrow(election_pps)
#'   pi <- election_jointprob
#'   horvitz_thompson_matrix <- matrix(nrow = n, ncol = n)
#'   for (i in seq_len(n)) {
#'     for (j in seq_len(n)) {
#'       horvitz_thompson_matrix[i,j] <- 1 - (pi[i,i] * pi[j,j])/pi[i,j]
#'     }
#'   }
#'
#'   ## Equivalently:
#'
#'   horvitz_thompson_matrix <- make_quad_form_matrix(
#'     variance_estimator = "Horvitz-Thompson",
#'     joint_probs = election_jointprob
#'   )
#'
#' # Make generalized bootstrap adjustment factors ----
#'
#'   bootstrap_adjustment_factors <- make_gen_boot_factors(
#'     Sigma = horvitz_thompson_matrix,
#'     num_replicates = 80,
#'     tau = 'auto'
#'   )
#'
#' # Determine replication scale factor for variance estimation ----
#'   tau <- attr(bootstrap_adjustment_factors, 'tau')
#'   B <- ncol(bootstrap_adjustment_factors)
#'   replication_scaling_constant <- tau^2 / B
#'
#' # Create a replicate design object ----
#'   election_pps_bootstrap_design <- svrepdesign(
#'     data = election_pps,
#'     weights = 1 / diag(election_jointprob),
#'     repweights = bootstrap_adjustment_factors,
#'     combined.weights = FALSE,
#'     type = "other",
#'     scale = attr(bootstrap_adjustment_factors, 'scale'),
#'     rscales = attr(bootstrap_adjustment_factors, 'rscales')
#'   )
#'
#' # Compare estimates to Horvitz-Thompson estimator ----
#'
#'   election_pps_ht_design <- svydesign(
#'     id = ~1,
#'     fpc = ~p,
#'     data = election_pps,
#'     pps = ppsmat(election_jointprob),
#'     variance = "HT"
#'   )
#'
#' svytotal(x = ~ Bush + Kerry,
#'          design = election_pps_bootstrap_design)
#' svytotal(x = ~ Bush + Kerry,
#'          design = election_pps_ht_design)
#' }
make_gen_boot_factors <- function(Sigma, num_replicates, tau = "auto", exact_vcov = FALSE) {

  n <- nrow(Sigma)

  # Generate replicate factors by simulating from a multivariate normal distribution
  if (!exact_vcov) {
    replicate_factors <- t(
      mvtnorm::rmvnorm(n = num_replicates,
                       mean = rep(1, times = n),
                       sigma = Sigma,
                       checkSymmetry = TRUE,
                       method = "eigen")
    )
  }
  if (exact_vcov) {
    qr_rank <- qr(Sigma)$rank
    if (num_replicates <= qr_rank) {
      err_msg <- paste(
        "`exact_vcov=TRUE` only works if `num_replicates` exceeds the rank of `Sigma`,",
        sprintf("which is %s", qr_rank)
      )
      stop(err_msg)
    }
    # Draw from standard multivariate normal
    draws <- mvtnorm::rmvnorm(n = num_replicates,
                              mean = rep(0, times = n),
                              sigma = diag(n),
                              checkSymmetry = FALSE,
                              method = "eigen")
    # Whiten the draws
    draws <- scale(draws, center = TRUE, scale = FALSE)
    draws <- draws %*% svd(draws, nu = 0, nv = n)$v
    draws <- scale(draws, FALSE, TRUE)
    # Color the draws
    replicate_factors <- with(svd(Sigma, nu = 0, nv = n), {
      (v %*% diag(sqrt(pmax(d, 0)), n)) %*% t(draws)
    })
    # Adjust the scale of the replicates,
    # since the generalized bootstrap uses the scale (1/B),
    # but the coloring process gives the exact result only for scale (1/(B-1))
    scale_adjustment <- sqrt((num_replicates)/(num_replicates-1))
    replicate_factors <- replicate_factors * scale_adjustment
    # Make sure that draws have the exact specified mean
    replicate_factors <- rep(1, times = n) + replicate_factors
  }

  # (Potentially) rescale to avoid negative weights
  rescaled_replicate_factors <- rescale_reps(
    x = replicate_factors,
    tau = tau,
    min_wgt = 0.01
  )

  selected_tau <- attr(rescaled_replicate_factors, 'tau')

  attr(rescaled_replicate_factors, 'scale') <- (selected_tau^2)/num_replicates
  attr(rescaled_replicate_factors, 'rscales') <- rep(1, times = num_replicates)

  # Set column names
  colnames(rescaled_replicate_factors) <- sprintf("REP_%s", seq_len(num_replicates))

  # Return result
  return(rescaled_replicate_factors)
}

#' @title Convert a survey design object to a generalized bootstrap replicate design
#' @description Converts a survey design object to a replicate design object
#' with replicate weights formed using the generalized bootstrap method.
#' The generalized survey bootstrap is a method for forming bootstrap replicate weights
#' from a textbook variance estimator, provided that the variance estimator
#' can be represented as a quadratic form whose matrix is positive semidefinite
#' (this covers a large class of variance estimators).
#' @param design A survey design object created using the 'survey' (or 'srvyr') package,
#' with class \code{'survey.design'} or \code{'svyimputationList'}.
#' @param variance_estimator The name of the variance estimator
#' whose quadratic form matrix should be created. See the section "Variance Estimators" below.
#' Options include:
#' \itemize{
#'   \item{\strong{"Yates-Grundy"}: }{The Yates-Grundy variance estimator based on
#'   first-order and second-order inclusion probabilities.}
#'   \item{\strong{"Horvitz-Thompson"}: }{The Horvitz-Thompson variance estimator based on
#'   first-order and second-order inclusion probabilities.}
#'   \item{\strong{"Poisson Horvitz-Thompson"}: }{The Horvitz-Thompson variance estimator
#'   based on assuming Poisson sampling, with first-order inclusion probabilities
#'   inferred from the sampling probabilities of the survey design object.}
#'   \item{\strong{"Stratified Multistage SRS"}: }{The usual stratified multistage variance estimator
#'   based on estimating the variance of cluster totals within strata at each stage.}
#'   \item{\strong{"Ultimate Cluster"}: }{The usual variance estimator based on estimating
#'   the variance of first-stage cluster totals within first-stage strata.}
#'   \item{\strong{"Deville-1"}: }{A variance estimator for unequal-probability
#'   sampling without replacement, described in Matei and Tillé (2005)
#'   as "Deville 1".}
#'   \item{\strong{"Deville-2"}: }{A variance estimator for unequal-probability
#'   sampling without replacement, described in Matei and Tillé (2005)
#'   as "Deville 2".}
#'   \item{\strong{"SD1"}: }{The non-circular successive-differences variance estimator described by Ash (2014),
#'   sometimes used for variance estimation for systematic sampling.}
#'   \item{\strong{"SD2"}: }{The circular successive-differences variance estimator described by Ash (2014).
#'   This estimator is the basis of the "successive-differences replication" estimator commonly used
#'   for variance estimation for systematic sampling.}
#' }
#' @param replicates Number of bootstrap replicates (should be as large as possible, given computer memory/storage limitations).
#' A commonly-recommended default is 500.
#' @param tau Either \code{"auto"}, or a single number. This is the rescaling constant
#' used to avoid negative weights through the transformation \eqn{\frac{w + \tau - 1}{\tau}},
#' where \eqn{w} is the original weight and \eqn{\tau} is the rescaling constant \code{tau}. \cr
#' If \code{tau="auto"}, the rescaling factor is determined automatically as follows:
#' if all of the adjustment factors are nonnegative, then \code{tau} is set equal to 1;
#' otherwise, \code{tau} is set to the smallest value needed to rescale
#' the adjustment factors such that they are all at least \code{0.01}.
#' @param exact_vcov If \code{exact_vcov=TRUE}, the replicate factors will be generated
#' such that variance estimates for totals exactly match the results from the target variance estimator.
#' This requires that \code{num_replicates} exceeds the rank of \code{Sigma}.
#' The replicate factors are generated by applying PCA-whitening to a collection of draws
#' from a multivariate Normal distribution, then applying a coloring transformation
#' to the whitened collection of draws.
#' @param psd_option Either \code{"warn"} (the default) or \code{"error"}.
#' This option specifies what will happen if the target variance estimator
#' has a quadratic form matrix which is not positive semidefinite. This
#' can occasionally happen, particularly for two-phase designs. \cr
#' If \code{psd_option="error"}, then an error message will be displayed. \cr
#' If \code{psd_option="warn"}, then a warning message will be displayed,
#' and the quadratic form matrix will be approximated by the most similar
#' positive semidefinite matrix.
#' This approximation was suggested by Beaumont and Patak (2012),
#' who note that this is conservative in the sense of producing
#' overestimates of variance.
#' Beaumont and Patak (2012) argue that this overestimation is expected to be
#' small in magnitude. See \code{\link[svrep]{get_nearest_psd_matrix}}
#' for details of the approximation.
#' @param compress This reduces the computer memory required to represent the replicate weights and has no
#' impact on estimates.
#' @param mse If \code{TRUE}, compute variances from sums of squares around the point estimate from the full-sample weights,
#' If \code{FALSE}, compute variances from sums of squares around the mean estimate from the replicate weights.
#' @return
#' A replicate design object, with class \code{svyrep.design}, which can be used with the usual functions,
#' such as \code{svymean()} or \code{svyglm()}. \cr
#' Use \code{weights(..., type = 'analysis')} to extract the matrix of replicate weights. \cr
#' Use \code{as_data_frame_with_weights()} to convert the design object to a data frame with columns
#' for the full-sample and replicate weights.
#' @export
#' @seealso Use \code{\link[svrep]{estimate_boot_reps_for_target_cv}} to help choose the number of bootstrap replicates. \cr
#' For greater customization of the method, \code{\link[svrep]{make_quad_form_matrix}} can be used to
#' represent several common variance estimators as a quadratic form's matrix,
#' which can then be used as an input to \code{\link[svrep]{make_gen_boot_factors}}.
#' The function \code{\link[svrep]{rescale_reps}} is used to implement
#' the rescaling of the bootstrap adjustment factors.
#' @section Statistical Details:
#' Let \eqn{v( \hat{T_y})} be the textbook variance estimator for an estimated population total \eqn{\hat{T}_y} of some variable \eqn{y}.
#' The base weight for case \eqn{i} in our sample is \eqn{w_i}, and we let \eqn{\breve{y}_i} denote the weighted value \eqn{w_iy_i}.
#' Suppose we can represent our textbook variance estimator as a quadratic form: \eqn{v(\hat{T}_y) = \breve{y}\Sigma\breve{y}^T},
#' for some \eqn{n \times n} matrix \eqn{\Sigma}.
#' The only constraint on \eqn{\Sigma} is that, for our sample, it must be symmetric and positive semidefinite.
#'
#' The bootstrapping process creates \eqn{B} sets of replicate weights, where the \eqn{b}-th set of replicate weights is a vector of length \eqn{n} denoted \eqn{\mathbf{a}^{(b)}}, whose \eqn{k}-th value is denoted \eqn{a_k^{(b)}}.
#' This yields \eqn{B} replicate estimates of the population total, \eqn{\hat{T}_y^{*(b)}=\sum_{k \in s} a_k^{(b)} \breve{y}_k}, for \eqn{b=1, \ldots B}, which can be used to estimate sampling variance.
#'
#' \deqn{
#'   v_B\left(\hat{T}_y\right)=\frac{\sum_{b=1}^B\left(\hat{T}_y^{*(b)}-\hat{T}_y\right)^2}{B}
#' }
#'
#' This bootstrap variance estimator can be written as a quadratic form:
#'
#'   \deqn{
#'     v_B\left(\hat{T}_y\right) =\mathbf{\breve{y}}^{\prime}\Sigma_B \mathbf{\breve{y}}
#'   }
#'   where
#'   \deqn{
#'     \boldsymbol{\Sigma}_B = \frac{\sum_{b=1}^B\left(\mathbf{a}^{(b)}-\mathbf{1}_n\right)\left(\mathbf{a}^{(b)}-\mathbf{1}_n\right)^{\prime}}{B}
#'   }
#'
#' Note that if the vector of adjustment factors \eqn{\mathbf{a}^{(b)}} has expectation \eqn{\mathbf{1}_n} and variance-covariance matrix \eqn{\boldsymbol{\Sigma}},
#' then we have the bootstrap expectation \eqn{E_{*}\left( \boldsymbol{\Sigma}_B \right) = \boldsymbol{\Sigma}}. Since the bootstrap process takes the sample values \eqn{\breve{y}} as fixed, the bootstrap expectation of the variance estimator is \eqn{E_{*} \left( \mathbf{\breve{y}}^{\prime}\Sigma_B \mathbf{\breve{y}}\right)= \mathbf{\breve{y}}^{\prime}\Sigma \mathbf{\breve{y}}}.
#' Thus, we can produce a bootstrap variance estimator with the same expectation as the textbook variance estimator simply by randomly generating \eqn{\mathbf{a}^{(b)}} from a distribution with the following two conditions:
#' \cr
#'     \strong{Condition 1}: \eqn{\quad \mathbf{E}_*(\mathbf{a})=\mathbf{1}_n}
#' \cr
#'     \strong{Condition 2}: \eqn{\quad \mathbf{E}_*\left(\mathbf{a}-\mathbf{1}_n\right)\left(\mathbf{a}-\mathbf{1}_n\right)^{\prime}=\mathbf{\Sigma}}
#' \cr \cr
#' While there are multiple ways to generate adjustment factors satisfying these conditions,
#' the simplest general method is to simulate from a multivariate normal distribution: \eqn{\mathbf{a} \sim MVN(\mathbf{1}_n, \boldsymbol{\Sigma})}.
#' This is the method used by this function.
#'
#' @section Details on Rescaling to Avoid Negative Adjustment Factors:
#' Let \eqn{\mathbf{A} = \left[ \mathbf{a}^{(1)} \cdots \mathbf{a}^{(b)} \cdots \mathbf{a}^{(B)} \right]} denote the \eqn{(n \times B)} matrix of bootstrap adjustment factors.
#' To eliminate negative adjustment factors, Beaumont and Patak (2012) propose forming a rescaled matrix of nonnegative replicate factors \eqn{\mathbf{A}^S} by rescaling each adjustment factor \eqn{a_k^{(b)}} as follows:
#' \deqn{
#'    a_k^{S,(b)} = \frac{a_k^{(b)} + \tau - 1}{\tau}
#'  }
#' where \eqn{\tau \geq 1 - a_k^{(b)} \geq 1} for all \eqn{k} in \eqn{\left\{ 1,\ldots,n \right\}} and all \eqn{b} in \eqn{\left\{1, \ldots, B\right\}}.
#'
#' The value of \eqn{\tau} can be set based on the realized adjustment factor matrix \eqn{\mathbf{A}} or by choosing \eqn{\tau} prior to generating the adjustment factor matrix \eqn{\mathbf{A}} so that \eqn{\tau} is likely to be large enough to prevent negative bootstrap weights.
#'
#' If the adjustment factors are rescaled in this manner, it is important to adjust the scale factor used in estimating the variance with the bootstrap replicates, which becomes \eqn{\frac{\tau^2}{B}} instead of \eqn{\frac{1}{B}}.
#' \deqn{
#'  \textbf{Prior to rescaling: } v_B\left(\hat{T}_y\right) = \frac{1}{B}\sum_{b=1}^B\left(\hat{T}_y^{*(b)}-\hat{T}_y\right)^2
#'  }
#' \deqn{
#'  \textbf{After rescaling: } v_B\left(\hat{T}_y\right) = \frac{\tau^2}{B}\sum_{b=1}^B\left(\hat{T}_y^{S*(b)}-\hat{T}_y\right)^2
#' }
#' When sharing a dataset that uses rescaled weights from a generalized survey bootstrap, the documentation for the dataset should instruct the user to use replication scale factor \eqn{\frac{\tau^2}{B}} rather than \eqn{\frac{1}{B}} when estimating sampling variances.
#'
#' @inheritSection make_quad_form_matrix Variance Estimators
#' @section Two-Phase Designs:
#' For a two-phase design, \code{variance_estimator} should be a list of variance estimators' names,
#' with two elements, such as \code{list('Ultimate Cluster', 'Poisson Horvitz-Thompson')}.
#' In two-phase designs, only the following estimators may be used for the second phase:
#' \itemize{
#'   \item "Ultimate Cluster"
#'   \item "Stratified Multistage SRS"
#'   \item "Poisson Horvitz-Thompson"
#' }
#' For statistical details on the handling of two-phase designs,
#' see the documentation for \link[svrep]{make_twophase_quad_form}.
#' @references
#' The generalized survey bootstrap was first proposed by Bertail and Combris (1997).
#' See Beaumont and Patak (2012) for a clear overview of the generalized survey bootstrap.
#' The generalized survey bootstrap represents one strategy for forming replication variance estimators
#' in the general framework proposed by Fay (1984) and Dippo, Fay, and Morganstein (1984).
#' \cr \cr
#' - Ash, S. (2014). "\emph{Using successive difference replication for estimating variances}."
#' \strong{Survey Methodology}, Statistics Canada, 40(1), 47–59.
#' \cr \cr
#' - Bellhouse, D.R. (1985). "\emph{Computing Methods for Variance Estimation in Complex Surveys}."
#' \strong{Journal of Official Statistics}, Vol.1, No.3.
#' \cr \cr
#' - Beaumont, Jean-François, and Zdenek Patak. 2012. “On the Generalized Bootstrap for Sample Surveys with Special Attention to Poisson Sampling: Generalized Bootstrap for Sample Surveys.” International Statistical Review 80 (1): 127–48. https://doi.org/10.1111/j.1751-5823.2011.00166.x.
#' \cr \cr
#' - Bertail, and Combris. 1997. “Bootstrap Généralisé d’un Sondage.” Annales d’Économie Et de Statistique, no. 46: 49. https://doi.org/10.2307/20076068.
#' \cr \cr
#' - Dippo, Cathryn, Robert Fay, and David Morganstein. 1984. “Computing Variances from Complex Samples with Replicate Weights.” In, 489–94. Alexandria, VA: American Statistical Association. http://www.asasrms.org/Proceedings/papers/1984_094.pdf.
#' \cr \cr
#' - Fay, Robert. 1984. “Some Properties of Estimates of Variance Based on Replication Methods.” In, 495–500. Alexandria, VA: American Statistical Association. http://www.asasrms.org/Proceedings/papers/1984_095.pdf.
#' \cr \cr
#' - Matei, Alina, and Yves Tillé. (2005).
#' “\emph{Evaluation of Variance Approximations and Estimators
#' in Maximum Entropy Sampling with Unequal Probability and Fixed Sample Size.}”
#' \strong{Journal of Official Statistics}, 21(4):543–70.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(survey)
#'
#'# Example 1: Bootstrap based on the Yates-Grundy estimator ----
#'    set.seed(2014)
#'
#'    data('election', package = 'survey')
#'
#'    ## Create survey design object
#'    pps_design_yg <- svydesign(
#'      data = election_pps,
#'      id = ~1, fpc = ~p,
#'      pps = ppsmat(election_jointprob),
#'      variance = "YG"
#'    )
#'
#'    ## Convert to generalized bootstrap replicate design
#'    gen_boot_design_yg <- pps_design_yg |>
#'      as_gen_boot_design(variance_estimator = "Yates-Grundy",
#'                         replicates = 1000, tau = "auto")
#'
#'    svytotal(x = ~ Bush + Kerry, design = pps_design_yg)
#'    svytotal(x = ~ Bush + Kerry, design = gen_boot_design_yg)
#'
#'# Example 2: Bootstrap based on the successive-difference estimator ----
#'
#'    data('library_stsys_sample', package = 'svrep')
#'
#'    ## First, ensure data are sorted in same order as was used in sampling
#'    library_stsys_sample <- library_stsys_sample[
#'      order(library_stsys_sample$SAMPLING_SORT_ORDER),
#'    ]
#'
#'    ## Create a survey design object
#'    design_obj <- svydesign(
#'      data = library_stsys_sample,
#'      strata = ~ SAMPLING_STRATUM,
#'      ids = ~ 1,
#'      fpc = ~ STRATUM_POP_SIZE
#'    )
#'
#'    ## Convert to generalized bootstrap replicate design
#'    gen_boot_design_sd2 <- as_gen_boot_design(
#'      design = design_obj,
#'      variance_estimator = "SD2",
#'      replicates = 2000
#'    )
#'
#'    ## Estimate sampling variances
#'    svytotal(x = ~ TOTSTAFF, na.rm = TRUE, design = gen_boot_design_sd2)
#'    svytotal(x = ~ TOTSTAFF, na.rm = TRUE, design = design_obj)
#'
#' # Example 3: Two-phase sample ----
#' # -- First stage is stratified systematic sampling,
#' # -- second stage is response/nonresponse modeled as Poisson sampling
#'
#'   nonresponse_model <- glm(
#'     data = library_stsys_sample,
#'     family = quasibinomial('logit'),
#'     formula = I(RESPONSE_STATUS == "Survey Respondent") ~ 1,
#'     weights = 1/library_stsys_sample$SAMPLING_PROB
#'   )
#'
#'   library_stsys_sample[['RESPONSE_PROPENSITY']] <- predict(
#'     nonresponse_model,
#'     newdata = library_stsys_sample,
#'     type = "response"
#'   )
#'
#'   twophase_design <- twophase(
#'     data = library_stsys_sample,
#'     # Identify cases included in second phase sample
#'     subset = ~ I(RESPONSE_STATUS == "Survey Respondent"),
#'     strata = list(~ SAMPLING_STRATUM, NULL),
#'     id = list(~ 1, ~ 1),
#'     probs = list(NULL, ~ RESPONSE_PROPENSITY),
#'     fpc = list(~ STRATUM_POP_SIZE, NULL),
#'   )
#'
#'   twophase_boot_design <- as_gen_boot_design(
#'     design = twophase_design,
#'     variance_estimator = list(
#'       "SD2", "Poisson Horvitz-Thompson"
#'     )
#'   )
#'
#'   svytotal(x = ~ LIBRARIA, design = twophase_boot_design)
#'
#' }
as_gen_boot_design <- function(design, variance_estimator = NULL,
                                replicates = 500, tau = "auto", exact_vcov = FALSE,
                                psd_option = "warn",
                                mse = getOption("survey.replicates.mse"),
                                compress = TRUE) {
  UseMethod("as_gen_boot_design", design)
}

#' @export
as_gen_boot_design.twophase2 <- function(design, variance_estimator = NULL,
                                         replicates = 500, tau = "auto",
                                         exact_vcov = FALSE, psd_option = "warn",
                                         mse = getOption("survey.replicates.mse"),
                                         compress = TRUE) {

  Sigma <- get_design_quad_form(design, variance_estimator)

  if (!is_psd_matrix(Sigma)) {
    problem_msg <- paste0(
      "The sample quadratic form matrix",
      " for this design and variance estimator",
      " is not positive semidefinite."
    )
    if (psd_option == "warn") {

      warning_msg <- paste0(
        problem_msg,
        " It will be approximated by the nearest",
        " positive semidefinite matrix."
      )
      warning(warning_msg)
      Sigma <- get_nearest_psd_matrix(Sigma)

    } else {
      error_msg <- paste0(
        problem_msg,
        " This can be handled using the `psd_option` argument."
      )
      stop(error_msg)
    }
  }

  adjustment_factors <- make_gen_boot_factors(
    Sigma = Sigma,
    num_replicates = replicates,
    tau = tau,
    exact_vcov = exact_vcov
  )

  rep_design <- survey::svrepdesign(
    variables = design$phase1$full$variables[design$subset,,drop=FALSE],
    weights = stats::weights(design, type = "sampling"),
    repweights = adjustment_factors,
    combined.weights = FALSE,
    compress = compress, mse = mse,
    scale = attr(adjustment_factors, 'scale'),
    rscales = attr(adjustment_factors, 'rscales'),
    type = "other"
  )

  rep_design$tau <- attr(adjustment_factors, 'tau')

  rep_design$call <- sys.call(which = -1)

  return(rep_design)
}

#' @export
as_gen_boot_design.survey.design <- function(design, variance_estimator = NULL,
                                             replicates = 500, tau = "auto", exact_vcov = FALSE,
                                             psd_option = 'warn',
                                             mse = getOption("survey.replicates.mse"),
                                             compress = TRUE) {

  Sigma <- get_design_quad_form(design, variance_estimator)

  if (!is_psd_matrix(Sigma)) {
    problem_msg <- "The sample quadratic form matrix for this design and variance estimator is not positive semidefinite."
    if (psd_option == "warn") {

      warning_msg <- paste0(
        problem_msg,
        " It will be approximated by the nearest positive semidefinite matrix."
      )
      warning(warning_msg)
      Sigma <- get_nearest_psd_matrix(Sigma)

    } else {
      error_msg <- paste0(
        problem_msg, " This can be handled using the `psd_option` argument."
      )
      stop(error_msg)
    }
  }

  adjustment_factors <- make_gen_boot_factors(
    Sigma = Sigma,
    num_replicates = replicates,
    tau = tau,
    exact_vcov = exact_vcov
  )

  rep_design <- survey::svrepdesign(
    variables = design$variables,
    weights = stats::weights(design, type = "sampling"),
    repweights = adjustment_factors,
    combined.weights = FALSE,
    compress = compress, mse = mse,
    scale = attr(adjustment_factors, 'scale'),
    rscales = attr(adjustment_factors, 'rscales'),
    type = "other"
  )

  rep_design$tau <- attr(adjustment_factors, 'tau')

  rep_design$call <- sys.call(which = -1)

  return(rep_design)
}
