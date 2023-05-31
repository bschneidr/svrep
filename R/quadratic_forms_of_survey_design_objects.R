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
#' @param ensure_psd If \code{TRUE} (the default), ensures
#' that the result is a positive semidefinite matrix. This
#' is necessary if the quadratic form is used as an input for
#' replication methods such as the generalized bootstrap.
#' For mathematical details, please see the documentation for the function \code{get_nearest_psd_matrix()}.
#' The approximation method is discussed by Beaumont and Patak (2012)
#' in the context of forming replicate weights for two-phase samples.
#' The authors argue that this approximation should
#' lead to only a small overestimation of variance.
#' @return A matrix representing the quadratic form of a specified variance estimator,
#' based on extracting information about clustering, stratification,
#' and selection probabilities from the survey design object.
#' @inheritSection make_quad_form_matrix Variance Estimators
#' @inheritSection as_gen_boot_design Two-Phase Designs
#' @references
#' - Ash, S. (2014). "\emph{Using successive difference replication for estimating variances}."
#' \strong{Survey Methodology}, Statistics Canada, 40(1), 47–59.
#' \cr \cr
#' - Beaumont, Jean-François, and Zdenek Patak. (2012). "\emph{On the Generalized Bootstrap for Sample Surveys with Special Attention to Poisson Sampling: Generalized Bootstrap for Sample Surveys.}"
#' \strong{International Statistical Review} 80 (1): 127–48.
#' \cr \cr
#' - Bellhouse, D.R. (1985). "\emph{Computing Methods for Variance Estimation in Complex Surveys}."
#' \strong{Journal of Official Statistics}, Vol.1, No.3.
#' \cr \cr
#' - Särndal, C.-E., Swensson, B., & Wretman, J. (1992). "\emph{Model Assisted Survey Sampling}." Springer New York.
#' @examples
#' \dontrun{
#' # Example 1: Quadratic form for successive-difference variance estimator ----
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
#'    ## Obtain quadratic form
#'    quad_form_matrix <- get_design_quad_form(
#'      design = design_obj,
#'      variance_estimator = "SD2"
#'    )
#'
#'    ## Estimate variance of estimated population total
#'    y <- design_obj$variables$LIBRARIA
#'    wts <- weights(design_obj, type = 'sampling')
#'    y_wtd <- as.matrix(y) * wts
#'    y_wtd[is.na(y_wtd)] <- 0
#'
#'    pop_total <- sum(y_wtd)
#'
#'    var_est <- t(y_wtd) %*% quad_form_matrix %*% y_wtd
#'    std_error <- sqrt(var_est)
#'
#'    print(pop_total); print(std_error)
#'
#'    # Compare to estimate from assuming SRS
#'    svytotal(x = ~ LIBRARIA, na.rm = TRUE,
#'             design = design_obj)
#'
#' # Example 2: Two-phase design (second phase is nonresponse) ----
#'
#'   ## Estimate response propensities, separately by stratum
#'   library_stsys_sample[['RESPONSE_PROB']] <- svyglm(
#'     design = design_obj,
#'     formula = I(RESPONSE_STATUS == "Survey Respondent") ~ SAMPLING_STRATUM,
#'     family = quasibinomial('logistic')
#'   ) |> predict(type = 'response')
#'
#'   ## Create a survey design object,
#'   ## where nonresponse is treated as a second phase of sampling
#'   twophase_design <- twophase(
#'     data = library_stsys_sample,
#'     strata = list(~ SAMPLING_STRATUM, NULL),
#'     id = list(~ 1, ~ 1),
#'     fpc = list(~ STRATUM_POP_SIZE, NULL),
#'     probs = list(NULL, ~ RESPONSE_PROB),
#'     subset = ~ I(RESPONSE_STATUS == "Survey Respondent")
#'   )
#'
#'   ## Obtain quadratic form for the two-phase variance estimator,
#'   ## where first phase variance contribution estimated
#'   ## using the successive differences estimator
#'   ## and second phase variance contribution estimated
#'   ## using the Horvitz-Thompson estimator
#'   ## (with joint probabilities based on assumption of Poisson sampling)
#'   get_design_quad_form(
#'     design = twophase_design,
#'     variance_estimator = list(
#'       "SD2",
#'       "Poisson Horvitz-Thompson"
#'     )
#'   )
#' }
#' @export

get_design_quad_form <- function(design, variance_estimator,
                                 ensure_psd = FALSE) {
  UseMethod("get_design_quad_form", design)
}

#' @export
get_design_quad_form.survey.design <- function(design, variance_estimator,
                                               ensure_psd = FALSE) {

  accepted_variance_estimators <- c(
    "Yates-Grundy", "Horvitz-Thompson",
    "Poisson Horvitz-Thompson",
    "Ultimate Cluster", "Stratified Multistage SRS",
    "SD1", "SD2", "Deville-1", "Deville-2"
  )

  if (is.null(variance_estimator)) {
    stop("Must specify a value for `variance_estimator`.")
  }

  if (length(variance_estimator) > 1) {
    stop("Can only specify one estimator for `variance_estimator`.")
  }

  if (!variance_estimator %in% accepted_variance_estimators) {
    sprintf("`%s` is not a supported variance estimator, or else there is a typo.",
            variance_estimator) |> stop()
  }

  is_pps_design <- isTRUE(design$pps)

  if (variance_estimator %in% c("Horvitz-Thompson", "Yates-Grundy")) {
    if (!is_pps_design) {
      sprintf("For `variance_estimator='%s'`, must use a PPS design. Please see the help page for `survey::svydesign()`.",
              variance_estimator) |>
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
  if (variance_estimator %in% c("Deville-1", "Deville-2")) {
    Sigma <- make_quad_form_matrix(
      variance_estimator = variance_estimator,
      cluster_ids = design$cluster,
      strata_ids = design$strata,
      probs = design$allprob,
      strata_pop_sizes = NULL,
      sort_order = NULL
    )
  }

  if (ensure_psd && !is_psd_matrix(Sigma)) {
    informative_msg <- paste0(
      "The variance estimator does not have a positive semidefinite quadratic form.",
      " Since `ensure_psd = TRUE`, the quadratic form matrix is approximated by the nearest positive semidefinite matrix."
    )
    message(informative_msg)
    Sigma <- get_nearest_psd_matrix(Sigma)
  }

  return(Sigma)
}

#' @export
get_design_quad_form.twophase2 <- function(design, variance_estimator,
                                           ensure_psd = FALSE) {

  # Check that `variance_estimator` is correctly specified
  accepted_phase1_estimators <- c(
    "Yates-Grundy", "Horvitz-Thompson",
    "Poisson Horvitz-Thompson",
    "Ultimate Cluster", "Stratified Multistage SRS",
    "SD1", "SD2", "Deville-1", "Deville-2"
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

  if (any(sapply(variance_estimator, length) > 1)) {
    stop("Can only specify one estimator for each element of `variance_estimator`.")
  }

  if (!variance_estimator[[1]] %in% accepted_phase1_estimators) {
    sprintf("`%s` is not a supported variance estimator for the first phase of a two-phase design, or else there is a typo.",
            variance_estimator[[1]]) |> stop()
  }
  if (!variance_estimator[[2]] %in% accepted_phase2_estimators) {
    sprintf("`%s` is not a supported variance estimator for the second phase of a two-phase design, or else there is a typo.",
            variance_estimator[[2]]) |> stop()
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
    phase_2_joint_probs = phase2_joint_prob,
    ensure_psd = ensure_psd
  )

  if (ensure_psd && !is_psd_matrix(Sigma)) {
    informative_msg <- paste0(
      "The combined two-phase variance estimator",
      " does not have a positive semidefinite quadratic form.",
      " Since `ensure_psd = TRUE`, the quadratic form matrix",
      " is approximated by the nearest positive semidefinite matrix."
    )
    message(informative_msg)
    Sigma <- get_nearest_psd_matrix(Sigma)
  }

  return(Sigma)
}
