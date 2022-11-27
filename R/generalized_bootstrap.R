#' @title Creates replicate factors for the generalized survey bootstrap
#' @description Creates replicate factors for the generalized survey bootstrap method.
#' The generalized survey bootstrap is a method for forming bootstrap replicate weights
#' from a textbook variance estimator, provided that the variance estimator
#' can be represented as a quadratic form whose matrix is positive semi-definite
#' (this covers a large class of variance estimators).
#' @param Sigma The matrix of the quadratic form used to represent the variance estimator.
#' Must be positive semi-definite.
#' @param num_replicates The number of bootstrap replicates to create.
#' @param tau Either \code{"auto"}, or a single number. This is the rescaling constant
#' used to avoid negative weights through the transformation \eqn{\frac{w + \tau - 1}{\tau}},
#' where \eqn{w} is the original weight and \eqn{\tau} is the rescaling constant \code{tau}. \cr
#' If \code{tau="auto"}, the rescaling factor is determined automatically as follows:
#' if all of the adjustment factors are nonnegative, then \code{tau} is set equal to 1;
#' otherwise, \code{tau} is set to the smallest value needed to rescale
#' the adjustment factors such that they are all at least \code{0.01}.
#'
#' @section Statistical Details:
#' Let \eqn{v( \hat{T_y})} be the textbook variance estimator for an estimated population total \eqn{\hat{T}_y} of some variable \eqn{y}.
#' The base weight for case \eqn{i} in our sample is \eqn{w_i}, and we let \eqn{\breve{y}_i} denote the weighted value \eqn{w_iy_i}.
#' Suppose we can represent our textbook variance estimator as a quadratic form: \eqn{v(\hat{T}_y) = \breve{y}\Sigma\breve{y}^T},
#' for some \eqn{n \times n} matrix \eqn{\Sigma}.
#' The only constraint on \eqn{\Sigma} is that, for our sample, it must be symmetric and positive semi-definite.
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
#'     \begin{aligned}
#'     v_B\left(\hat{T}_y\right) &=\mathbf{\breve{y}}^{\prime}\Sigma_B \mathbf{\breve{y}} \\
#'     \textit{where}& \\
#'     \boldsymbol{\Sigma}_B &= \frac{\sum_{b=1}^B\left(\mathbf{a}^{(b)}-\mathbf{1}_n\right)\left(\mathbf{a}^{(b)}-\mathbf{1}_n\right)^{\prime}}{B}
#'     \end{aligned}
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
#'  \begin{aligned}
#'  a_k^{S,(b)} &= \frac{a_k^{(b)} + \tau - 1}{\tau} \\
#'  \textit{where } \tau &\geq 1 - a_k^{(b)} \geq 1 \\
#'  &\textit{for all }k \textit{ in } \left\{ 1,\ldots,n \right\} \\
#'  &\textit{and all }b \textit{ in } \left\{1, \ldots, B\right\} \\
#'  \end{aligned}
#' }
#' The value of \eqn{\tau} can be set based on the realized adjustment factor matrix \eqn{\mathbf{A}} or by choosing \eqn{\tau} prior to generating the adjustment factor matrix \eqn{\mathbf{A}} so that \eqn{\tau} is likely to be large enough to prevent negative bootstrap weights.
#'
#' If the adjustment factors are rescaled in this manner, it is important to adjust the scale factor used in estimating the variance with the bootstrap replicates, which becomes \eqn{\frac{\tau^2}{B}} instead of \eqn{\frac{1}{B}}.
#' \deqn{
#'  \begin{aligned}
#'  \textbf{Prior to rescaling: } v_B\left(\hat{T}_y\right) &= \frac{1}{B}\sum_{b=1}^B\left(\hat{T}_y^{*(b)}-\hat{T}_y\right)^2 \\
#'  \textbf{After rescaling: } v_B\left(\hat{T}_y\right) &= \frac{\tau^2}{B}\sum_{b=1}^B\left(\hat{T}_y^{S*(b)}-\hat{T}_y\right)^2 \\
#'  \end{aligned}
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
#' library(survey)
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
#'
make_gen_boot_factors <- function(Sigma, num_replicates, tau = "auto") {

  n <- nrow(Sigma)

  if (!isSymmetric.matrix(Sigma)) {
    stop("`Sigma` must be a symmetric matrix.")
  }

  # Generate replicate factors by simulating from a multivariate normal distribution
  replicate_factors <- t(
    MASS::mvrnorm(n = num_replicates,
                  mu = rep(1, times = n),
                  Sigma = Sigma,
                  empirical = FALSE)
  )

  # (Potentially) rescale to avoid negative weights
  if (length(tau) != 1 || is.na(tau) || (tau != "auto" & !is.numeric(tau)) || (is.numeric(tau) & tau < 0)) {
    stop("`tau` must be either 'auto' or a single positive number.")
  }

  if (any(replicate_factors < 0)) {
    if (tau == "auto") {
      rescaling_constant <- max(1 - replicate_factors) + 0.01
    } else {
      rescaling_constant <- tau
    }
    rescaled_replicate_factors <- (replicate_factors + (rescaling_constant-1))/rescaling_constant
  } else {
    rescaling_constant <- 1
    rescaled_replicate_factors <- replicate_factors
  }
  attr(rescaled_replicate_factors, 'tau') <- rescaling_constant
  attr(rescaled_replicate_factors, 'scale') <- (rescaling_constant^2)/num_replicates
  attr(rescaled_replicate_factors, 'rscales') <- rep(1, times = num_replicates)

  # Set column names
  colnames(rescaled_replicate_factors) <- sprintf("REP_%s", seq_len(num_replicates))

  # Return result
  return(rescaled_replicate_factors)
}

#' @title Convert a survey design object to a generalized boostrap replicate design
#' @description Converts a survey design object to a replicate design object
#' with replicate weights formed using the generalized bootstrap method.
#' The generalized survey bootstrap is a method for forming bootstrap replicate weights
#' from a textbook variance estimator, provided that the variance estimator
#' can be represented as a quadratic form whose matrix is positive semi-definite
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
#'   \item{\strong{"Stratified Multistage SRS"}: }{The usual stratified multistage variance estimator
#'   based on estimating the variance of cluster totals within strata at each stage.}
#'   \item{\strong{"Ultimate Cluster"}: }{The usual variance estimator based on estimating
#'   the variance of first-stage cluster totals within first-stage strata.}
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
#' @param compress This reduces the computer memory required to represent the replicate weights and has no
#' impact on estimates.
#' @param mse If \code{TRUE}, compute variances from sums of squares around the point estimate from the full-sample weights,
#' If \code{FALSE}, compute variances from sums of squares around the mean estimate from the replicate weights.
#' @return
#' A replicate design object, with class \code{svyrep.design}, which can be used with the usual functions,
#' such as \code{svymean()} or \code{svyglm()}. \cr \cr
#' Use \code{weights(..., type = 'analysis')} to extract the matrix of replicate weights.
#' Use \code{as_data_frame_with_weights()} to convert the design object to a data frame with columns
#' for the full-sample and replicate weights.
#' @export
#' @seealso Use \code{\link[svrep]{estimate_boot_reps_for_target_cv}} to help choose the number of bootstrap replicates. \cr
#' For greater customization of the method, \code{\link[svrep]{make_quad_form_matrix}} can be used to
#' represent several common variance estimators as a quadratic form's matrix,
#' which can then be used as an input to \code{\link[svrep]{make_gen_boot_factors}}.
#' @section Statistical Details:
#' Let \eqn{v( \hat{T_y})} be the textbook variance estimator for an estimated population total \eqn{\hat{T}_y} of some variable \eqn{y}.
#' The base weight for case \eqn{i} in our sample is \eqn{w_i}, and we let \eqn{\breve{y}_i} denote the weighted value \eqn{w_iy_i}.
#' Suppose we can represent our textbook variance estimator as a quadratic form: \eqn{v(\hat{T}_y) = \breve{y}\Sigma\breve{y}^T},
#' for some \eqn{n \times n} matrix \eqn{\Sigma}.
#' The only constraint on \eqn{\Sigma} is that, for our sample, it must be symmetric and positive semi-definite (in other words, it should never lead to a negative variance estimate, no matter what the value of \eqn{\breve{y}} is).
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
#'     \begin{aligned}
#'     v_B\left(\hat{T}_y\right) &=\mathbf{\breve{y}}^{\prime}\Sigma_B \mathbf{\breve{y}} \\
#'     \textit{where}& \\
#'     \boldsymbol{\Sigma}_B &= \frac{\sum_{b=1}^B\left(\mathbf{a}^{(b)}-\mathbf{1}_n\right)\left(\mathbf{a}^{(b)}-\mathbf{1}_n\right)^{\prime}}{B}
#'     \end{aligned}
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
#'  \begin{aligned}
#'  a_k^{S,(b)} &= \frac{a_k^{(b)} + \tau - 1}{\tau} \\
#'  \textit{where } \tau &\geq 1 - a_k^{(b)} \geq 1 \\
#'  &\textit{for all }k \textit{ in } \left\{ 1,\ldots,n \right\} \\
#'  &\textit{and all }b \textit{ in } \left\{1, \ldots, B\right\} \\
#'  \end{aligned}
#' }
#' The value of \eqn{\tau} can be set based on the realized adjustment factor matrix \eqn{\mathbf{A}} or by choosing \eqn{\tau} prior to generating the adjustment factor matrix \eqn{\mathbf{A}} so that \eqn{\tau} is likely to be large enough to prevent negative bootstrap weights.
#'
#' If the adjustment factors are rescaled in this manner, it is important to adjust the scale factor used in estimating the variance with the bootstrap replicates, which becomes \eqn{\frac{\tau^2}{B}} instead of \eqn{\frac{1}{B}}.
#' \deqn{
#'  \begin{aligned}
#'  \textbf{Prior to rescaling: } v_B\left(\hat{T}_y\right) &= \frac{1}{B}\sum_{b=1}^B\left(\hat{T}_y^{*(b)}-\hat{T}_y\right)^2 \\
#'  \textbf{After rescaling: } v_B\left(\hat{T}_y\right) &= \frac{\tau^2}{B}\sum_{b=1}^B\left(\hat{T}_y^{S*(b)}-\hat{T}_y\right)^2 \\
#'  \end{aligned}
#' }
#' When sharing a dataset that uses rescaled weights from a generalized survey bootstrap, the documentation for the dataset should instruct the user to use replication scale factor \eqn{\frac{\tau^2}{B}} rather than \eqn{\frac{1}{B}} when estimating sampling variances.
#' @section Variance Estimators:
#' The \strong{Horvitz-Thompson} variance estimator:
#' \deqn{
#'   v(\hat{Y}) = \sum_{i \in s}\sum_{j \in s} (1 - \frac{\pi_i \pi_j}{\pi_{ij}}) \frac{y_i}{\pi_i} \frac{y_j}{\pi_j}
#' }
#' The \strong{Yates-Grundy} variance estimator:
#' \deqn{
#'   v(\hat{Y}) = -\frac{1}{2}\sum_{i \in s}\sum_{j \in s} (1 - \frac{\pi_i \pi_j}{\pi_{ij}}) (\frac{y_i}{\pi_i} - \frac{y_j}{\pi_j})^2
#' }
#' The \strong{Stratified Multistage SRS} variance estimator is the recursive variance estimator
#' proposed by Bellhouse (1985) and used in the 'survey' package's function \link[survey]{svyrecvar}.
#' The estimator can be used for any number of sampling stages. For illustration, we describe its use
#' for two sampling stages.
#' \deqn{
#'   \begin{aligned}
#'   v(\hat{Y}) &= \hat{V}_1 + \hat{V}_2 \\
#'   \hat{V}_1 &= \sum_{h=1}^{H} (1 - \frac{n_h}{N_h})\frac{n_h}{n_h - 1} \sum_{i=1}^{n_h} (y_{hi.} - \bar{y}_{hi.})^2 \\
#'   \hat{V}_2 &= \sum_{h=1}^{H} \frac{n_h}{N_h} \sum_{i=1}^{n_h}v_{hi}(y_{hi.}) \\
#'   \end{aligned}
#' }
#' where:
#' \deqn{
#'   \begin{aligned}
#'   n_h &\textit{ is the number of sampled clusters in stratum }h \\
#'   N_h &\textit{ is the number of population clusters in stratum }h \\
#'   y_{hi.} &\textit{ is the weighted cluster total in cluster }i\textit{ of stratum }h \\
#'   \bar{y}_{hi.} &= \frac{1}{n_h}\sum_{i=1}^{n_h}y_{hi.} \\
#'   v_{hi}(y_{hi.}) &\textit{ is the estimated sampling variance of }y_{hi.}
#'   \end{aligned}
#' }
#' The \strong{Ultimate Cluster} variance estimator is simply the stratified multistage SRS
#' variance estimator, but ignoring variances from later stages of sampling.
#' \deqn{
#'   v(\hat{Y}) = \hat{V}_1
#' }
#' This is the variance estimator used in the 'survey' package when the user specifies
#' \code{option(survey.ultimate.cluster = TRUE)} or uses \code{svyrecvar(..., one.stage = TRUE)}.
#' When the first-stage sampling fractions are small, analysts often omit the finite population corrections \eqn{(1-\frac{n_h}{N_h})}
#' when using the ultimate cluster estimator.
#' \cr \cr
#' The \strong{SD1} and \strong{SD2} variance estimators are "successive difference"
#' estimators sometimes used for systematic sampling designs.
#' Ash (2014) describes each estimator as follows:
#' \deqn{
#' \begin{aligned}
#'   \hat{v}_{S D 1}(\hat{Y}) &= \left(1-\frac{n}{N}\right) \frac{n}{2(n-1)} \sum_{k=2}^n\left(\breve{y}_k-\breve{y}_{k-1}\right)^2 \\
#'   \hat{v}_{S D 2}(\hat{Y}) &= \left(1-\frac{n}{N}\right) \frac{1}{2}\left[\sum_{k=2}^n\left(\breve{y}_k-\breve{y}_{k-1}\right)^2+\left(\breve{y}_n-\breve{y}_1\right)^2\right] \\
#' \end{aligned}
#' }
#' where \eqn{\breve{y}_k = y_k/\pi_k} is the weighted value of unit \eqn{k}
#' with selection probability \eqn{\pi_k}. The SD1 estimator is recommended by Wolter (1984).
#' The SD2 estimator is the basis of the successive difference replication estimator commonly
#' used for systematic sampling designs. See Ash (2014) for details.
#' \cr \cr
#' For multistage samples, SD1 and SD2 are applied to the clusters at each stage, separately by stratum.
#' For later stages of sampling, the variance estimate from a stratum is multiplied by the product
#' of sampling fractions from earlier stages of sampling. For example, at a third stage of sampling,
#' the variance estimate from a third-stage stratum is multiplied by \eqn{\frac{n_1}{N_1}\frac{n_2}{N_2}},
#' which is the product of sampling fractions from the first-stage stratum and second-stage stratum.
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
#' @export
#'
#' @examples
#'
#'library(survey)
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
as_gen_boot_design <- function(design, variance_estimator = NULL,
                                replicates = 500, tau = "auto",
                                mse = getOption("survey.replicates.mse"),
                                compress = TRUE) {
  UseMethod("as_gen_boot_design", design)
}

#' @export
as_gen_boot_design.survey.design <- function(design, variance_estimator = NULL,
                                             replicates = 500, tau = "auto",
                                             mse = getOption("survey.replicates.mse"),
                                             compress = TRUE) {

  is_pps_design <- isTRUE(design$pps)

  if (variance_estimator %in% c("Horvitz-Thompson", "Yates-Grundy")) {
    if (!is_pps_design) {
      sprintf("For `variance_estimator='%s'`, must use a PPS design. Please see the help page for `survey::svydesign()`.") |>
        stop()
    }
    if ((variance_estimator == "Yates-Grundy") & (!design$variance %in% c("YG"))) {
      sprintf("Must specify `variance='YG'` when creating the survey design object.`") |>
        stop()
    }
    if ((variance_estimator == "Horvitz-Thompson") & (!design$variance %in% c("HT"))) {
      sprintf("Must specify `variance='HT'` when creating the survey design object.`") |>
        stop()
    }
    Sigma <- dpps_ht[['dcheck']][[1]]$dcheck |> as.matrix()
  }

  if (variance_estimator %in% c("SD1", "SD2")) {
    sprintf("For `variance_estimator='%s', assumes rows of data are sorted in the same order used in sampling.",
            variance_estimator) |> message()
    Sigma <- make_quad_form_matrix(
      variance_estimator = variance_estimator,
      cluster_ids = design$cluster,
      strata_ids = design$strata,
      strata_pop_sizes = design$fpc$popsize,
      sort_order = seq_len(nrow(design))
    )
  }
  if (variance_estimator == "Ultimate Cluster") {
    Sigma <- make_quad_form_matrix(
      variance_estimator = variance_estimator,
      cluster_ids = design$cluster,
      strata_ids = design$strata,
      strata_pop_sizes = design$fpc$popsize,
      sort_order = NULL
    )
  }

  adjustment_factors <- make_gen_boot_factors(
    Sigma = Sigma,
    num_replicates = replicates,
    tau = tau
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

  return(rep_design)
}
