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
#'   \item \strong{"Yates-Grundy"}: \cr The Yates-Grundy variance estimator based on
#'     first-order and second-order inclusion probabilities.
#'   \item \strong{"Horvitz-Thompson"}: \cr The Horvitz-Thompson variance estimator based on
#'     first-order and second-order inclusion probabilities.
#'   \item \strong{"Poisson Horvitz-Thompson"}: \cr The Horvitz-Thompson variance estimator
#'     based on assuming Poisson sampling, with first-order inclusion probabilities
#'     inferred from the sampling probabilities of the survey design object.
#'   \item \strong{"Stratified Multistage SRS"}: \cr The usual stratified multistage variance estimator
#'     based on estimating the variance of cluster totals within strata at each stage.
#'   \item \strong{"Ultimate Cluster"}: \cr The usual variance estimator based on estimating
#'     the variance of first-stage cluster totals within first-stage strata.
#'   \item \strong{"Deville-1"}: \cr A variance estimator for unequal-probability
#'     sampling without replacement, described in Matei and Tillé (2005)
#'     as "Deville 1".
#'   \item \strong{"Deville-2"}: \cr A variance estimator for unequal-probability
#'     sampling without replacement, described in Matei and Tillé (2005) as "Deville 2".
#'   \item \strong{"Deville-Tille": } \cr A variance estimator useful
#'     for balanced sampling designs, proposed by Deville and Tillé (2005).
#'   \item \strong{"SD1"}: \cr The non-circular successive-differences variance estimator described by Ash (2014),
#'     sometimes used for variance estimation for systematic sampling.
#'   \item\strong{"SD2"}:  \cr The circular successive-differences variance estimator described by Ash (2014).
#'     This estimator is the basis of the "successive-differences replication" estimator commonly used
#'     for variance estimation for systematic sampling.
#'   \item \strong{"Beaumont-Emond"}: \cr The variance estimator of Beaumont and Emond (2022)
#'     for multistage unequal-probability sampling without replacement.
#'   \item \strong{"BOSB"}: \cr The kernel-based variance estimator proposed by
#'     Breidt, Opsomer, and Sanchez-Borrego (2016) for use with systematic samples
#'     or other finely stratified designs. Uses the Epanechnikov kernel
#'     with the bandwidth automatically chosen to result in the smallest possible
#'     nonempty kernel window.
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
#' @param aux_var_names Only required if \code{variance_estimator = "Deville-Tille"}
#' or if \code{variance_estimator = "BOSB"}.
#' For the Deville-Tillé estimator, this should be a character vector of 
#' variable names for auxiliary variables to be used in the variance estimator.
#' For the BOSB estimator, this should be a string giving a single variable name
#' to use as an auxiliary variable in the kernel-based variance estimator
#' of Breidt, Opsomer, and Sanchez-Borrego (2016).
#' @return A matrix representing the quadratic form of a specified variance estimator,
#' based on extracting information about clustering, stratification,
#' and selection probabilities from the survey design object.
#' @inheritSection make_quad_form_matrix Variance Estimators
#' @inheritSection as_gen_boot_design Two-Phase Designs
#' @references
#' - Ash, S. (2014). "\emph{Using successive difference replication for estimating variances}."
#' \strong{Survey Methodology}, Statistics Canada, 40(1), 47-59.
#' \cr \cr
#' - Beaumont, Jean-François, and Zdenek Patak. (2012). "\emph{On the Generalized Bootstrap for Sample Surveys with Special Attention to Poisson Sampling: Generalized Bootstrap for Sample Surveys.}"
#' \strong{International Statistical Review} 80 (1): 127-48.
#' \cr \cr
#' - Bellhouse, D.R. (1985). "\emph{Computing Methods for Variance Estimation in Complex Surveys}."
#' \strong{Journal of Official Statistics}, Vol.1, No.3.
#' \cr \cr
#' - Breidt, F. J., Opsomer, J. D., & Sanchez-Borrego, I. (2016). 
#' "\emph{Nonparametric Variance Estimation Under Fine Stratification: An Alternative to Collapsed Strata}." 
#' \strong{Journal of the American Statistical Association}, 111(514), 822-833. https://doi.org/10.1080/01621459.2015.1058264
#' \cr \cr
#' - Deville, J.C., and Tillé, Y. (2005). "\emph{Variance approximation under balanced sampling.}"
#' \strong{Journal of Statistical Planning and Inference}, 128, 569-591.
#' \cr \cr
#' - Särndal, C.E., Swensson, B., & Wretman, J. (1992). "\emph{Model Assisted Survey Sampling}." Springer New York.
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
#' # Example 2: Kernel-based variance estimator ----
#' 
#'    Q_BOSB <- get_design_quad_form(
#'      design             = design_obj,
#'      variance_estimator = "BOSB",
#'      aux_var_names      = "SAMPLING_SORT_ORDER"
#'    )
#'    
#'    var_est <- t(y_wtd) %*% Q_BOSB %*% y_wtd
#'    std_error <- sqrt(var_est)
#'    
#'    print(pop_total); print(std_error)
#'
#' # Example 3: Two-phase design (second phase is nonresponse) ----
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
                                 ensure_psd = FALSE,
                                 aux_var_names = NULL) {
  UseMethod("get_design_quad_form", design)
}

#' @export
get_design_quad_form.survey.design <- function(design, variance_estimator,
                                               ensure_psd = FALSE,
                                               aux_var_names = NULL) {

  accepted_variance_estimators <- c(
    "Yates-Grundy", "Horvitz-Thompson",
    "Poisson Horvitz-Thompson",
    "Ultimate Cluster", "Stratified Multistage SRS",
    "SD1", "SD2", "BOSB", 
    "Deville-1", "Deville-2", "Deville-Tille", "Beaumont-Emond"
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
    Sigma <- design[['dcheck']][[1]]$dcheck |>
      as("symmetricMatrix")

    if (variance_estimator == "Yates-Grundy") {
      Sigma <- - Sigma
      diag(Sigma) <- Matrix::diag(Sigma) - Matrix::rowSums(Sigma)
      Sigma <- - Sigma
    }
  }

  if (variance_estimator %in% c("Poisson Horvitz-Thompson")) {
    Sigma <- Matrix::diag(1 - design$prob)
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
  if (variance_estimator %in% c("Deville-1", "Deville-2", "Beaumont-Emond")) {
    Sigma <- make_quad_form_matrix(
      variance_estimator = variance_estimator,
      cluster_ids = design$cluster,
      strata_ids = design$strata,
      probs = design$allprob,
      strata_pop_sizes = NULL,
      sort_order = NULL
    )
  }
  if (variance_estimator %in% c("Deville-Tille")) {

    if (is.null(aux_var_names)) {
      stop("For `variance_estimator='Deville-Tille', must supply a vector of variable names to `aux_var_names`.")
    }

    if (!all(aux_var_names %in% colnames(design$variables))) {
      stop("Some of `aux_var_names` do not show up as columns in the design object.")
    }

    aux_vars_matrix <- model.matrix(
      object = reformulate(termlabels = aux_var_names, intercept = FALSE),
      data = design$variables[,aux_var_names,drop=FALSE]
    )

    Sigma <- make_quad_form_matrix(
      variance_estimator = variance_estimator,
      cluster_ids = design$cluster,
      strata_ids = design$strata,
      probs = design$allprob,
      aux_vars = aux_vars_matrix,
      strata_pop_sizes = NULL,
      sort_order = NULL
    )
  }
  
  if (variance_estimator %in% c("BOSB")) {
    
    if (is.null(aux_var_names) || (length(aux_var_names) != 1)) {
      stop("For `variance_estimator='BOSB', must supply a single variable name to `aux_var_names`.")
    }
    
    if (!all(aux_var_names %in% colnames(design$variables))) {
      stop("Some of `aux_var_names` do not show up as columns in the design object.")
    }
    
    if (!is.numeric(design$variables[,aux_var_names,drop=TRUE])) {
      stop("The variable specified by `aux_var_names` should be a numeric variable.")
    }
    
    aux_var_matrix <- model.matrix(
      object = reformulate(termlabels = aux_var_names, intercept = FALSE),
      data = design$variables[,aux_var_names,drop=FALSE]
    )
    
    Sigma <- make_quad_form_matrix(
      variance_estimator = variance_estimator,
      cluster_ids        = design$cluster,
      strata_ids         = design$strata,
      probs              = design$allprob,
      aux_vars           = aux_var_matrix,
      strata_pop_sizes   = NULL,
      sort_order         = NULL
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
                                           ensure_psd = FALSE,
                                           aux_var_names = NULL) {

  # Check that `variance_estimator` is correctly specified
  accepted_phase1_estimators <- c(
    "Yates-Grundy", "Horvitz-Thompson",
    "Poisson Horvitz-Thompson",
    "Ultimate Cluster", "Stratified Multistage SRS",
    "SD1", "SD2", "Deville-1", "Deville-2", "Deville-Tille", "Beaumont-Emond", "BOSB"
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
    variance_estimator = variance_estimator[[1]],
    aux_var_names = aux_var_names
  )
  Sigma_phase1 <- Sigma_phase1[
    as.logical(design$subset), 
    as.logical(design$subset)
  ]

  # Extract the quadratic form for second phase
  # (conditional on first phase sample)
  Sigma_phase2 <- get_design_quad_form(
    design = design$phase2,
    variance_estimator = variance_estimator[[2]],
    aux_var_names = NULL
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

#' @title Produce a compressed representation of a survey design object
#'
#' @param design A survey design object
#' @param vars_to_keep (Optional) A character vector
#' of variables in the design to keep in the compressed design.
#' By default, none of the variables are retained.
#' @return A list with two elements. The \code{design_subset}
#' element is a a design object with only the minimal rows
#' needed to represent the survey design.
#' The \code{index} element links each row of the original design
#' to a row of \code{design_subset}, so that the design can be "uncompressed."
#' @keywords internal
compress_design <- function(design, vars_to_keep = NULL) {

  UseMethod("compress_design", design)
}

#' @export
compress_design.survey.design <- function(design, vars_to_keep = NULL) {

  if (is.null(vars_to_keep)) {
    vars_to_keep <- 0
  }

  if ((!is.null(design$pps)) && (design$pps != FALSE)) {
    compressed_design_structure <- list(
      design_subset = design,
      index = seq_len(nrow(design))
    )
  } else {
    design_structure <- cbind(design$strata, design$cluster)
    tmp <- apply(design_structure, 1, function(x) paste(x, collapse = "\r"))
    unique_elements <- !duplicated(design_structure)
    compressed_design_structure <- list(
      design_subset = design |> (\(design_obj) {
        # Reduce memory usage by dropping variables
        design_obj$variables <- design_obj$variables[,vars_to_keep,drop=FALSE]
        # Subset to only unique strata/cluster combos
        design_obj[unique_elements,]
      })(),
      index = match(tmp, tmp[unique_elements])
    )
  }

  return(compressed_design_structure)
}

#' @export
compress_design.DBIsvydesign <- function(design, vars_to_keep = NULL) {
  
  if (!is.null(vars_to_keep)) {
    design$variables <- getvars(
      formula      = vars_to_keep,
      dbconnection = design$db$connection,
      tables       = design$db$tablename,
      updates      = design$updates,
      subset       = design$subset
    )
  }

  # Produce a (potentially) compressed survey design object
  if ((!is.null(design$pps)) && (design$pps != FALSE)) {
    compressed_design_structure <- list(
      design_subset = design,
      index = seq_len(nrow(design))
    )
  } else {
    design_structure <- cbind(design$strata, design$cluster)
    tmp <- apply(design_structure, 1, function(x) paste(x, collapse = "\r"))
    unique_elements <- !duplicated(design_structure)
    compressed_design_structure <- list(
      design_subset = design |> (\(design_obj) {
        # Reduce memory usage by dropping variables
        if (!is.null(design_obj$variables)) {
          design_obj$variables <- design_obj$variables[unique_elements,vars_to_keep,drop=FALSE]
        }
        # Subset to only unique strata/cluster/weight/fpc combos
        design_obj$strata <- design_obj$strata[unique_elements,, drop = FALSE]
        design_obj$cluster <- design_obj$cluster[unique_elements,, drop = FALSE]
        if (!is.null(design_obj$allprob)) {
          design_obj$allprob <- design_obj$allprob[unique_elements,, drop = FALSE]
        }
        if (!is.null(design_obj$fpc$sampsize)) {
          design_obj$fpc$sampsize <- design_obj$fpc$sampsize[unique_elements,, drop = FALSE]
        }
        if (!is.null(design_obj$fpc$popsize)) {
          design_obj$fpc$popsize <- design_obj$fpc$popsize[unique_elements,, drop = FALSE]
        }
        design_obj$prob <- design_obj$prob[unique_elements]
        return(design_obj)
      })(),
      index = match(tmp, tmp[unique_elements])
    )
  }

  return(compressed_design_structure)
  }
