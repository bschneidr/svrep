#' @title Variance Estimators
#' @description This help page describes variance estimators
#' which are commonly used for survey samples. These variance estimators
#' can be used as the basis of the generalized replication methods, implemented
#' with the functions \code{\link[svrep]{as_fays_gen_rep_design}()},
#' \code{\link[svrep]{as_gen_boot_design}()},
#' \code{\link[svrep]{make_fays_gen_rep_factors}()},
#' or \code{\link[svrep]{make_gen_boot_factors}()}
#' @rdname variance-estimators
#' @name variance-estimators
#' @section Shared Notation:
#' Let \eqn{s} denote the selected sample of size \eqn{n}, with elements \eqn{i=1,\dots,n}.
#' Element \eqn{i} in the sample had probability \eqn{\pi_i} of being included in the sample.
#' The \emph{pair} of elements \eqn{ij} was sampled with probability \eqn{\pi_{ij}}.
#'
#' The population total for a variable is denoted \eqn{Y = \sum_{i \in U}y_i},
#' and the Horvitz-Thompson estimator for \eqn{\hat{Y}} is
#' denoted \eqn{\hat{Y} = \sum_{i \in s} y_i/\pi_i}. For convenience,
#' we denote \eqn{\breve{y}_i = y_i/\pi_i}.
#'
#' The true sampling variance of \eqn{\hat{Y}} is denoted \eqn{V(\hat{Y})},
#' while an estimator of this sampling variance is denoted \eqn{v(\hat{Y})}.
#' @section Horvitz-Thompson:
#' The \strong{Horvitz-Thompson} variance estimator:
#' \deqn{
#'   v(\hat{Y}) = \sum_{i \in s}\sum_{j \in s} (1 - \frac{\pi_i \pi_j}{\pi_{ij}}) \frac{y_i}{\pi_i} \frac{y_j}{\pi_j}
#' }
#' @section Yates-Grundy:
#' The \strong{Yates-Grundy} variance estimator:
#' \deqn{
#'   v(\hat{Y}) = -\frac{1}{2}\sum_{i \in s}\sum_{j \in s} (1 - \frac{\pi_i \pi_j}{\pi_{ij}}) (\frac{y_i}{\pi_i} - \frac{y_j}{\pi_j})^2
#' }
#' @section Poisson Horvitz-Thompson:
#' The \strong{Poisson Horvitz-Thompson} variance estimator
#' is simply the Horvitz-Thompson variance estimator, but
#' where \eqn{\pi_{ij}=\pi_i \times \pi_j}, which is the case for Poisson sampling.
#' \cr \cr
#' @section Stratified Multistage SRS:
#' The \strong{Stratified Multistage SRS} variance estimator is the recursive variance estimator
#' proposed by Bellhouse (1985) and used in the 'survey' package's function \link[survey]{svyrecvar}.
#' In the case of simple random sampling without replacement (with one or more stages),
#' this estimator exactly matches the Horvitz-Thompson estimator.
#'
#' The estimator can be used for any number of sampling stages. For illustration, we describe its use
#' for two sampling stages.
#' \deqn{
#'   v(\hat{Y}) = \hat{V}_1 + \hat{V}_2
#' }
#' where
#' \deqn{
#'   \hat{V}_1 = \sum_{h=1}^{H} (1 - \frac{n_h}{N_h})\frac{n_h}{n_h - 1} \sum_{i=1}^{n_h} (y_{hi.} - \bar{y}_{hi.})^2
#' }
#' and
#' \deqn{
#'   \hat{V}_2 = \sum_{h=1}^{H} \frac{n_h}{N_h} \sum_{i=1}^{n_h}v_{hi}(y_{hi.})
#' }
#' where \eqn{n_h} is the number of sampled clusters in stratum \eqn{h},
#' \eqn{N_h} is the number of population clusters in stratum \eqn{h},
#' \eqn{y_{hi.}} is the weighted cluster total in cluster \eqn{i} of stratum \eqn{h},
#' \eqn{\bar{y}_{hi.}} is the mean weighted cluster total of stratum \eqn{h},
#' (\eqn{\bar{y}_{hi.} = \frac{1}{n_h}\sum_{i=1}^{n_h}y_{hi.}}), and
#' \eqn{v_{hi}(y_{hi.})} is the estimated sampling variance of \eqn{y_{hi.}}.
#' \cr \cr
#' @section Ultimate Cluster:
#' The \strong{Ultimate Cluster} variance estimator is simply the stratified multistage SRS
#' variance estimator, but ignoring variances from later stages of sampling.
#' \deqn{
#'   v(\hat{Y}) = \hat{V}_1
#' }
#' This is the variance estimator used in the 'survey' package when the user specifies
#' \code{option(survey.ultimate.cluster = TRUE)} or uses \code{svyrecvar(..., one.stage = TRUE)}.
#' When the first-stage sampling fractions are small, analysts often omit the finite population corrections \eqn{(1-\frac{n_h}{N_h})}
#' when using the ultimate cluster estimator.
#' @section SD1 and SD2 (Successive Difference Estimators):
#' The \strong{SD1} and \strong{SD2} variance estimators are "successive difference"
#' estimators sometimes used for systematic sampling designs.
#' Ash (2014) describes each estimator as follows:
#' \deqn{
#'   \hat{v}_{S D 1}(\hat{Y}) = \left(1-\frac{n}{N}\right) \frac{n}{2(n-1)} \sum_{k=2}^n\left(\breve{y}_k-\breve{y}_{k-1}\right)^2
#' }
#' \deqn{
#'   \hat{v}_{S D 2}(\hat{Y}) = \left(1-\frac{n}{N}\right) \frac{1}{2}\left[\sum_{k=2}^n\left(\breve{y}_k-\breve{y}_{k-1}\right)^2+\left(\breve{y}_n-\breve{y}_1\right)^2\right]
#' }
#' where \eqn{\breve{y}_k = y_k/\pi_k} is the weighted value of unit \eqn{k}
#' with selection probability \eqn{\pi_k}. The SD1 estimator is recommended by Wolter (1984).
#' The SD2 estimator is the basis of the successive difference replication estimator commonly
#' used for systematic sampling designs and is more conservative. See Ash (2014) for details.
#' \cr \cr
#' For multistage samples, SD1 and SD2 are applied to the clusters at each stage, separately by stratum.
#' For later stages of sampling, the variance estimate from a stratum is multiplied by the product
#' of sampling fractions from earlier stages of sampling. For example, at a third stage of sampling,
#' the variance estimate from a third-stage stratum is multiplied by \eqn{\frac{n_1}{N_1}\frac{n_2}{N_2}},
#' which is the product of sampling fractions from the first-stage stratum and second-stage stratum.
#' @section Deville 1 and Deville 2:
#' The \strong{"Deville-1"} and \strong{"Deville-2"} variance estimators
#' are clearly described in Matei and Tillé (2005),
#' and are intended for designs that use
#' fixed-size, unequal-probability random sampling without replacement.
#' These variance estimators have been shown to be effective
#' for designs that use a fixed sample size with a high-entropy sampling method.
#' This includes most PPSWOR sampling methods,
#' but unequal-probability systematic sampling is an important exception.
#'
#' These variance estimators take the following form:
#' \deqn{
#' \hat{v}(\hat{Y}) = \sum_{i=1}^{n} c_i (\breve{y}_i - \frac{1}{\sum_{i=k}^{n}c_k}\sum_{k=1}^{n}c_k \breve{y}_k)^2
#' }
#' where \eqn{\breve{y}_i = y_i/\pi_i} is the weighted value of the the variable of interest,
#' and \eqn{c_i} depend on the method used:
#' \itemize{
#'   \item \strong{"Deville-1"}: 
#'     \deqn{c_i=\left(1-\pi_i\right) \frac{n}{n-1}}
#'   \item \strong{"Deville-2"}: 
#'     \deqn{c_i = (1-\pi_i) \left[1 - \sum_{k=1}^{n} \left(\frac{1-\pi_k}{\sum_{k=1}^{n}(1-\pi_k)}\right)^2 \right]^{-1}}
#' }
#' In the case of simple random sampling without replacement (SRSWOR),
#' these estimators are both identical to the usual stratified multistage SRS estimator
#' (which is itself a special case of the Horvitz-Thompson estimator).
#'
#' For multistage samples, "Deville-1" and "Deville-2" are applied to the clusters at each stage, separately by stratum.
#' For later stages of sampling, the variance estimate from a stratum is multiplied by the product
#' of sampling probabilities from earlier stages of sampling. For example, at a third stage of sampling,
#' the variance estimate from a third-stage stratum is multiplied by \eqn{\pi_1 \times \pi_{(2 | 1)}},
#' where \eqn{\pi_1} is the sampling probability of the first-stage unit
#' and \eqn{\pi_{(2|1)}} is the sampling probability of the second-stage unit
#' within the first-stage unit.
#' @section Deville-Tillé:
#' See Section 6.8 of Tillé (2020) for more detail on this estimator,
#' including an explanation of its quadratic form.
#' See Deville and Tillé (2005) for the results of a simulation study
#' comparing this and other alternative estimators for balanced sampling.
#'
#' The estimator can be written as follows:
#' \deqn{
#'   v(\hat{Y})=\sum_{k \in S} \frac{c_k}{\pi_k^2}\left(y_k-\hat{y}_k^*\right)^2,
#' }
#' where
#' \deqn{
#'   \hat{y}_k^*=\mathbf{z}_k^{\top}\left(\sum_{\ell \in S} c_{\ell} \frac{\mathbf{z}_{\ell} \mathbf{z}_{\ell}^{\prime}}{\pi_{\ell}^2}\right)^{-1} \sum_{\ell \in S} c_{\ell} \frac{\mathbf{z}_{\ell} y_{\ell}}{\pi_{\ell}^2}
#' }
#' and \eqn{\mathbf{z}_k} denotes the vector of auxiliary variables for observation \eqn{k}
#' included in sample \eqn{S}, with inclusion probability \eqn{\pi_k}. The value \eqn{c_k} is set to \eqn{\frac{n}{n-q}(1-\pi_k)},
#' where \eqn{n} is the number of observations and \eqn{q} is the number of auxiliary variables.
#' @section Aronow-Samii:
#' The Aronow-Samii variance estimator is intended for use in 
#' non-measurable survey designs, where some pairs of units in the population have zero
#' pairwise inclusion probabilities (i.e., where \eqn{\pi_{ij} = 0} for some \eqn{i} and \eqn{j}).
#' An important example is when a design samples only one unit from a stratum
#' with multiple units. Unbiased variance estimation is not in general possible for such designs.
#' 
#' The estimator proposed by Aronow and Samii (2013) is conservative in the sense that
#' its bias is guaranteed to be weakly greater than zero. Aronow and Samii (2013)
#' provide a detailed discussion of this estimator's properties, concluding
#' that the estimator should be effective for designs that sample one unit
#' per stratum. However, they caution against using the estimator
#' for systematic samples, as in such cases it may be too upwardly biased to be useful.
#' 
#' This estimator requires the joint probabilities for the entire \emph{population},
#' or at least for each stratum in the population. Thus, it cannot
#' be used with the functions \code{as_fays_gen_rep_design()} or \code{as_gen_boot_design()}.
#' Instead, it can only be used with \code{make_quad_form_matrix()}.
#' 
#' Let \eqn{I_i} equal 1 if unit \eqn{i} is included in the sample and equal if not,
#' and let \eqn{I_j} equal 1 if unit \eqn{j} is included in the sample and equal 0 if not.
#' 
#' Then the Aronow-Samii variance estimator is written as follows:
#' \deqn{
#'   v(\hat{Y}) = T_1 + T_2
#' }
#' 
#' where
#' \deqn{
#'   T_1 = \sum_{i \in U}\sum_{j \in \{U: \pi_{ij} > 0\}} I_i I_j (1 - \frac{\pi_i \pi_j}{\pi_{ij}}) \frac{y_i}{\pi_i} \frac{y_j}{\pi_j}
#' }
#' 
#' and
#' 
#' \deqn{
#'   T_2 = \sum_{i \in U}\sum_{j \in \{U: \pi_{ij} = 0\}} I_i \frac{y_i^2}{2 \pi_i} + I_j \frac{y_j^2}{2 \pi_j}
#' }
#' 
#' @references
#' 
#' Aronow, P., and Samii, C. (2013). "\emph{Conservative variance estimation 
#' for sampling designs with zero pairwise inclusion probabilities}."
#' \strong{Survey Methodology}, Statistics Canada, 39(1), 231-241.
#' 
#' Ash, S. (2014). "\emph{Using successive difference replication for estimating variances}."
#' \strong{Survey Methodology}, Statistics Canada, 40(1), 47–59.
#'
#' Bellhouse, D.R. (1985). "\emph{Computing Methods for Variance Estimation in Complex Surveys}."
#' \strong{Journal of Official Statistics}, Vol.1, No.3.
#'
#' Deville, J.‐C., and Tillé, Y. (2005). "\emph{Variance approximation under balanced sampling.}"
#' \strong{Journal of Statistical Planning and Inference}, 128, 569–591.
#'
#' Tillé, Y. (2020). "\emph{Sampling and estimation from finite populations}." (I. Hekimi, Trans.). Wiley.
#'
#' Matei, Alina, and Yves Tillé. (2005).
#' “\emph{Evaluation of Variance Approximations and Estimators
#' in Maximum Entropy Sampling with Unequal Probability and Fixed Sample Size.}”
#' \strong{Journal of Official Statistics}, 21(4):543–70.
NULL
