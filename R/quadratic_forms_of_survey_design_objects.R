#' @title Determine the quadratic form matrix of a variance estimator for a survey design object
#' @description Determines the quadratic form matrix of a specified variance estimator,
#' by parsing the information stored in a survey design object created using
#' the 'survey' package.
#' @param design A survey design object created using the 'survey' (or 'srvyr') package,
#' with class \code{'survey.design'} or \code{'svyimputationList'}. Also accepts two-phase design objects
#' with class \code{'twophase2'}; see the section below titled "Two-Phase Designs" for
#' more information about handling of two-phase designs.
#' @param variance_estimator The name of the variance estimator
#' whose quadratic form matrix should be created. \cr
#' See the section "Variance Estimators" below.
#' Options include:
#' \itemize{
#'   \item{\strong{"Yates-Grundy"}: }{The Yates-Grundy variance estimator based on
#'   first-order and second-order inclusion probabilities.}
#'   \item{\strong{"Horvitz-Thompson"}: }{The Horvitz-Thompson variance estimator based on
#'   first-order and second-order inclusion probabilities.}
#'   \item{\strong{"Poisson Horvitz-Thompson"}: }{The Horvitz-Thompson variance estimator
#'   based on assuming Poisson sampling with specified first-order inclusion probabilities.}
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
#' @return A matrix representing the quadratic form of a specified variance estimator,
#' based on extracting information about clustering, stratification,
#' and selection probabilities from the survey design object.
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
#'   \hat{v}_{S D 1}(\hat{Y}) = \left(1-\frac{n}{N}\right) \frac{n}{2(n-1)} \sum_{k=2}^n\left(\breve{y}_k-\breve{y}_{k-1}\right)^2
#' }
#' \deqn{
#'   \hat{v}_{S D 2}(\hat{Y}) = \left(1-\frac{n}{N}\right) \frac{1}{2}\left[\sum_{k=2}^n\left(\breve{y}_k-\breve{y}_{k-1}\right)^2+\left(\breve{y}_n-\breve{y}_1\right)^2\right]
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
#' @section Two-Phase Designs:
#' For a two-phase design, \code{variance_estimator} should be a list of variance estimators,
#' with two elements, such as \code{list('Ultimate Cluster', 'Poisson Horvitz-Thompson')}.
#' In two-phase designs, only the following estimators may be used for the second phase:
#' \itemize{
#'   \item "Ultimate Cluster"
#'   \item "Stratified Multistage SRS"
#'   \item "Poisson Horvitz-Thompson"
#' }
#'
#' @references
#' - Ash, S. (2014). "\emph{Using successive difference replication for estimating variances}."
#' \strong{Survey Methodology}, Statistics Canada, 40(1), 47–59.
#' \cr \cr
#' - Bellhouse, D.R. (1985). "\emph{Computing Methods for Variance Estimation in Complex Surveys}."
#' \strong{Journal of Official Statistics}, Vol.1, No.3.
#' \cr \cr
#' - Särndal, C.-E., Swensson, B., & Wretman, J. (1992). "\emph{Model Assisted Survey Sampling}." Springer New York.
#' @examples
#' \dontrun{
#' # Example 1: Bootstrap based on the successive-difference estimator ----
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
#' # Example 2: Two-phase design (second phase is nonresponse) ----
#'
#'   ## Create a survey design object
#'   twophase_design <- twophase(
#'     data = library_stsys_sample,
#'     strata = list(~ SAMPLING_STRATUM, NULL),
#'     id = list(~ 1, ~ 1),
#'     fpc = list(~ STRATUM_POP_SIZE, NULL),
#'     subset = ~ I(RESPONSE_STATUS == "Survey Respondent")
#'   )
#'
#'   ## Convert to a bootstrap design
#'   twophase_genboot_design <- as_gen_boot_design(
#'     design = twophase_design,
#'     variance_estimator = list(
#'       "SD2",
#'       "Poisson Horvitz-Thompson"
#'     ),
#'     replicates = 2000
#'   )
#' }
#' @export

get_design_quad_form <- function(design, variance_estimator) {
  UseMethod("get_design_quad_form", design)
}

#' @export
get_design_quad_form.survey.design <- function(design, variance_estimator) {

  accepted_variance_estimators <- c(
    "Yates-Grundy", "Horvitz-Thompson",
    "Poisson Horvitz-Thompson",
    "Ultimate Cluster", "Stratified Multistage SRS",
    "SD1", "SD2"
  )

  if (is.null(variance_estimator)) {
    stop("Must specify a value for `variance_estimator`.")
  }

  if (!variance_estimator %in% accepted_variance_estimators) {
    sprintf("`%s` is not a supported variance estimator, or else there is a typo.",
            variance_estimator) |> stop()
  }
  if (length(variance_estimator) > 1) {
    stop("Can only specify one estimator for `variance_estimator`.")
  }

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
    Sigma <- design[['dcheck']][[1]]$dcheck |> as.matrix()

    if (variance_estimator == "Yates-Grundy") {
      Sigma <- - Sigma
      diag(Sigma) <- diag(Sigma) - rowSums(Sigma)
      Sigma <- - Sigma
    }
  }

  if (variance_estimator %in% c("Poisson Horvitz-Thompson")) {
    Sigma <- diag(1 - design$prob)
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
  if (variance_estimator %in% c("Ultimate Cluster", "Stratified Multistage SRS")) {
    Sigma <- make_quad_form_matrix(
      variance_estimator = variance_estimator,
      cluster_ids = design$cluster,
      strata_ids = design$strata,
      strata_pop_sizes = design$fpc$popsize,
      sort_order = NULL
    )
  }
  return(Sigma)
}

#' @export
get_design_quad_form.twophase2 <- function(design, variance_estimator) {

  # Check that `variance_estimator` is correctly specified
  accepted_phase1_estimators <- c(
    "Yates-Grundy", "Horvitz-Thompson",
    "Poisson Horvitz-Thompson",
    "Ultimate Cluster", "Stratified Multistage SRS",
    "SD1", "SD2"
  )
  accepted_phase2_estimators <- c(
    "Ultimate Cluster", "Stratified Multistage SRS",
    "Poisson Horvitz-Thompson"
  )

  if (is.null(variance_estimator) || (any(sapply(variance_estimator, is.null)))) {
    stop("Must specify a value for `variance_estimator`.")
  }

  if (!is.list(variance_estimator) || (length(variance_estimator) != 2)) {
    stop("For a two-phase design, must specify `variance_estimator` as a list with two elements.")
  }

  if (!variance_estimator[[1]] %in% accepted_phase1_estimators) {
    sprintf("`%s` is not a supported variance estimator for the first phase of a two-phase design, or else there is a typo.",
            variance_estimator[[1]]) |> stop()
  }
  if (!variance_estimator[[2]] %in% accepted_phase2_estimators) {
    sprintf("`%s` is not a supported variance estimator for the second phase of a two-phase design, or else there is a typo.",
            variance_estimator[[2]]) |> stop()
  }
  if (any(sapply(variance_estimator, length) > 1)) {
    stop("Can only specify one estimator for each element of `variance_estimator`.")
  }

  # Extract quadratic form for first phase,
  # and subset to only cases selected in the second phase sample
  Sigma_phase1 <- get_design_quad_form(
    design = design$phase1$full,
    variance_estimator = variance_estimator[[1]]
  )
  Sigma_phase1 <- Sigma_phase1[design$subset, design$subset]

  # Extract the quadratic form for second phase
  # (conditional on first phase sample)
  Sigma_phase2 <- get_design_quad_form(
    design = design$phase2,
    variance_estimator = variance_estimator[[2]]
  )

  # Obtain phase 2 conditional joint inclusion probabilities
  phase2_joint_prob <- ht_matrix_to_joint_probs(Sigma_phase2)

  # Combine the quadratic forms from the two phases
  Sigma <- make_twophase_quad_form(
    sigma_1 = Sigma_phase1,
    sigma_2 = Sigma_phase2,
    phase_2_joint_probs = phase2_joint_prob
  )

  return(Sigma)
}
