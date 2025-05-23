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
#' @section Beaumont-Emond:
#' The \strong{"Beaumont-Emond"} variance estimator was proposed by Beaumont and Emond (2022),
#' intended for designs that use fixed-size, unequal-probability random sampling without replacement.
#' The variance estimator is simply the Horvitz-Thompson
#' variance estimator with the following approximation for the joint inclusion
#' probabilities.
#' \deqn{
#'   \pi_{kl} \approx \pi_k \pi_l \frac{n - 1}{(n-1) + \sqrt{(1-\pi_k)(1-\pi_l)}}
#' }
#' In the case of cluster sampling, this approximation is
#' applied to the clusters rather than the units within clusters,
#' with \eqn{n} denoting the number of sampled clusters. and the probabilities \eqn{\pi}
#' referring to the cluster's sampling probability. For stratified samples,
#' the joint probability for units \eqn{k} and \eqn{l} in different strata
#' is simply the product of \eqn{\pi_k} and \eqn{\pi_l}.
#' 
#' For multistage samples, this approximation is applied to the clusters at each stage, separately by stratum.
#' For later stages of sampling, the variance estimate from a stratum is multiplied by the product
#' of sampling probabilities from earlier stages of sampling. For example, at a third stage of sampling,
#' the variance estimate from a third-stage stratum is multiplied by \eqn{\pi_1 \times \pi_{(2 | 1)}},
#' where \eqn{\pi_1} is the sampling probability of the first-stage unit
#' and \eqn{\pi_{(2|1)}} is the sampling probability of the second-stage unit
#' within the first-stage unit.
#' 
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
#' @section BOSB:
#' This kernel-based variance estimator was proposed by Breidt, Opsomer, and Sanchez-Borrego (2016),
#' for use with samples selected using systematic sampling or where only a single
#' sampling unit is selected from each stratum (sometimes referred to as "fine stratification").
#' 
#' Suppose there are \eqn{n} sampled units, and
#' for each unit \eqn{i} there is a numeric population characteristic \eqn{x_i}
#' and there is a weighted total \eqn{\hat{Y}_i}, where
#' \eqn{\hat{Y}_i} is only observed in the selected sample but \eqn{x_i}
#' is known prior to sampling.
#' 
#' The variance estimator has the following form:
#' 
#' \deqn{
#'   \hat{V}_{ker}=\frac{1}{C_d} \sum_{i=1}^n (\hat{Y}_i-\sum_{j=1}^n d_j(i) \hat{Y}_j)^2
#' }
#' 
#' The terms \eqn{d_j(i)} are kernel weights given by
#' 
#' \deqn{
#'   d_j(i)=\frac{K(\frac{x_i-x_j}{h})}{\sum_{j=1}^n K(\frac{x_i-x_j}{h})}
#' }
#' 
#' where \eqn{K(\cdot)} is a symmetric, bounded kernel function
#' and \eqn{h} is a bandwidth parameter. The normalizing constant \eqn{C_d} 
#' is computed as:
#' 
#' \deqn{
#'   C_d=\frac{1}{n} \sum_{i=1}^n(1-2 d_i(i)+\sum_{j=1}^H d_j^2(i))
#' }
#' 
#' For most functions in the 'svrep' package, the kernel function
#' is the Epanechnikov kernel and the bandwidth is automatically selected
#' to yield the smallest possible nonempty kernel window, as was recommended
#' by Breidt, Opsomer, and Sanchez-Borrego (2016). That's the case for
#' the functions \code{as_fays_gen_rep_design()}, \code{as_gen_boot_design()},
#' \code{make_quad_form_matrix()}, etc. However, users can construct the quadratic
#' form matrix of this variance estimator using a different kernel and a different bandwidth
#' by directly working with the function \code{make_kernel_var_matrix()}.
#' 
#' 
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
#' @references
#' Ash, S. (2014). "\emph{Using successive difference replication for estimating variances}."
#' \strong{Survey Methodology}, Statistics Canada, 40(1), 47-59.
#' 
#' Beaumont, J.F.; Émond, N. (2022). "\emph{A Bootstrap Variance Estimation Method for Multistage Sampling and Two-Phase Sampling When Poisson Sampling Is Used at the Second Phase}."
#' \strong{Stats}, \emph{5}: 339-357.
#' https://doi.org/10.3390/stats5020019
#'
#' Bellhouse, D.R. (1985). "\emph{Computing Methods for Variance Estimation in Complex Surveys}."
#' \strong{Journal of Official Statistics}, Vol.1, No.3.
#' 
#' Breidt, F. J., Opsomer, J. D., & Sanchez-Borrego, I. (2016). 
#' "\emph{Nonparametric Variance Estimation Under Fine Stratification: An Alternative to Collapsed Strata}." 
#' \strong{Journal of the American Statistical Association}, 111(514), 822-833. https://doi.org/10.1080/01621459.2015.1058264
#'
#' Deville, J.C., and Tillé, Y. (2005). "\emph{Variance approximation under balanced sampling.}"
#' \strong{Journal of Statistical Planning and Inference}, 128, 569-591.
#'
#' Tillé, Y. (2020). "\emph{Sampling and estimation from finite populations}." (I. Hekimi, Trans.). Wiley.
#'
#' Matei, Alina, and Yves Tillé. (2005).
#' "\emph{Evaluation of Variance Approximations and Estimators
#' in Maximum Entropy Sampling with Unequal Probability and Fixed Sample Size.}"
#' \strong{Journal of Official Statistics}, 21(4):543-70.
NULL
