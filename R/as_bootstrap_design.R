#' @title Convert a survey design object to a bootstrap replicate design
#' @description Converts a survey design object to a replicate design object
#' with replicate weights formed using a bootstrap method. Supports stratified,
#' cluster samples with one or more stages of sampling. At each stage of sampling,
#' either simple random sampling (with or without replacement)
#' or unequal probability sampling (with or without replacement) may be used.
#' @param design A survey design object created using the 'survey' (or 'srvyr') package,
#' with class \code{'survey.design'} or \code{'svyimputationList'}.
#' @param type The type of bootstrap to use, which should be chosen based
#' on its applicability to the sampling method used for the survey.
#' The available types are the following: \cr
#' \itemize{
#'  \item \strong{"Rao-Wu-Yue-Beaumont"} (the default): \cr
#'    The bootstrap method of Beaumont and Émond (2022), which is a generalization of the Rao-Wu-Yue bootstrap,
#'    and is applicable to a wide variety of designs, including single-stage and multistage stratified designs.
#'    The design may have different sampling methods used at different stages.
#'    Each stage of sampling may potentially be PPS (i.e., use unequal probabilities), with or without replacement,
#'    and may potentially use Poisson sampling. \cr \cr
#'    For a stratum with a fixed sample size of \eqn{n} sampling units, resampling in each replicate resamples \eqn{(n-1)} sampling units with replacement.
#'  \item \strong{"Rao-Wu"}: \cr
#'    The basic Rao-Wu \eqn{(n-1)} bootstrap method, which is only applicable to single-stage designs or
#'    multistage designs where the first-stage sampling fractions are small (and can thus be ignored).
#'    Accommodates stratified designs. All sampling within a stratum must be simple random sampling with or without replacement,
#'    although the first-stage sampling is effectively treated as sampling without replacement.
#'  \item \strong{"Antal-Tille"}: \cr
#'    The doubled half bootstrap method proposed by Antal and Tillé (2014), which is only applicable to single-stage designs or
#'    multistage designs where the first-stage sampling fractions are small (and can thus be ignored).
#'    Accomodates stratified designs. Sampling within each stratum can be simple random sampling or unequal probability sampling
#'    with or without replacement.
#'  \item \strong{"Preston"}: \cr
#'    Preston's multistage rescaled bootstrap, which is applicable to single-stage designs or multistage designs
#'    with arbitrary sampling fractions. Accommodates stratified designs. All sampling within a stratum must be
#'    simple random sampling with or without replacement.
#'  \item \strong{"Canty-Davison"}: \cr
#'    The Canty-Davison bootstrap, which is only applicable to single-stage designs, with arbitrary sampling fractions.
#'    Accommodates stratified designs. All sampling with a stratum must be simple random sampling with or without replacement.
#' }
#' @param replicates Number of bootstrap replicates (should be as large as possible, given computer memory/storage limitations).
#' A commonly-recommended default is 500.
#' @param compress Use a compressed representation of the replicate weights matrix.
#' This reduces the computer memory required to represent the replicate weights and has no
#' impact on estimates.
#' @param mse If \code{TRUE}, compute variances from sums of squares around the point estimate from the full-sample weights.
#' If \code{FALSE}, compute variances from sums of squares around the mean estimate from the replicate weights.
#' @param samp_method_by_stage (Optional). By default, this function will automatically determine the sampling method used at each stage.
#' However, this argument can be used to ensure the correct sampling method is identified for each stage. \cr
#' Accepts a vector with length equal to the number of stages of sampling.
#' Each element should be one of the following: \cr
#' \itemize{
#'  \item \code{"SRSWOR"} - Simple random sampling, without replacement
#'  \item \code{"SRSWR"} - Simple random sampling, with replacement
#'  \item \code{"PPSWOR"} - Unequal probabilities of selection, without replacement
#'  \item \code{"PPSWR"} - Unequal probabilities of selection, with replacement
#'  \item \code{"Poisson"} -  Poisson sampling: each sampling unit is selected into the sample at most once, with potentially different probabilities of inclusion for each sampling unit.
#' }
#' @return
#' A replicate design object, with class \code{svyrep.design}, which can be used with the usual functions,
#' such as \code{svymean()} or \code{svyglm()}.
#'
#' Use \code{weights(..., type = 'analysis')} to extract the matrix of replicate weights. \cr
#' Use \code{as_data_frame_with_weights()} to convert the design object to a data frame with columns
#' for the full-sample and replicate weights.
#' @export
#' @seealso Use \code{\link[svrep]{estimate_boot_reps_for_target_cv}} to help choose the number of bootstrap replicates.
#'
#' The underlying function for the Rao-Wu-Yue-Beaumont bootstrap
#' is \code{\link[svrep]{make_rwyb_bootstrap_weights}},
#' and the underlying function for the Antal-Tillé bootstrap 
#' is \code{\link[svrep]{make_doubled_half_bootstrap_weights}}.
#' Other bootstrap methods are implemented using functions from the 'survey' package,
#' including: \code{\link[survey]{bootweights}} (Canty-Davison),
#' \code{\link[survey]{subbootweights}} (Rao-Wu),
#' and \code{\link[survey]{mrbweights}} (Preston).
#'
#'
#' For systematic samples, one-PSU-per-stratum designs, or other especially complex sample designs,
#' one can use the generalized survey bootstrap method. See \code{\link[svrep]{as_gen_boot_design}} or \code{\link[svrep]{make_gen_boot_factors}}.
#'
#' @references
#' Antal, E. and Tillé, Y. (2014). 
#' "A new resampling method for sampling designs without replacement: The doubled half bootstrap." 
#' \strong{Computational Statistics}, \emph{29}(5), 1345-1363. https://doi.org/10.1007/s00180-014-0495-0
#' 
#' Beaumont, J.-F.; Émond, N. (2022).
#' "A Bootstrap Variance Estimation Method for Multistage Sampling and Two-Phase Sampling When Poisson Sampling Is Used at the Second Phase."
#' \strong{Stats}, \emph{5}: 339-357.
#' https://doi.org/10.3390/stats5020019
#'
#' Canty, A.J.; Davison, A.C. (1999).
#' "Resampling-based variance estimation for labour force surveys."
#' \strong{The Statistician}, \emph{48}: 379-391.
#'
#' Preston, J. (2009).
#' "Rescaled bootstrap for stratified multistage sampling."
#' \strong{Survey Methodology}, \emph{35}(2): 227-234.
#'
#' Rao, J.N.K.; Wu, C.F.J.; Yue, K. (1992).
#' "Some recent work on resampling methods for complex surveys."
#' \strong{Survey Methodology}, \emph{18}: 209-217.
#'
#'
#' @examples
#' library(survey)
#' # Example 1: A multistage sample with two stages of SRSWOR
#'
#'   ## Load an example dataset from a multistage sample, with two stages of SRSWOR
#'   data("mu284", package = 'survey')
#'   multistage_srswor_design <- svydesign(data = mu284,
#'                                         ids = ~ id1 + id2,
#'                                         fpc = ~ n1 + n2)
#'
#'   ## Convert the survey design object to a bootstrap design
#'   set.seed(2022)
#'   bootstrap_rep_design <- as_bootstrap_design(multistage_srswor_design,
#'                                               replicates = 500)
#'
#'   ## Compare std. error estimates from bootstrap versus linearization
#'   data.frame(
#'     'Statistic' = c('total', 'mean', 'median'),
#'     'SE (bootstrap)' = c(SE(svytotal(x = ~ y1, design = bootstrap_rep_design)),
#'                          SE(svymean(x = ~ y1, design = bootstrap_rep_design)),
#'                          SE(svyquantile(x = ~ y1, quantile = 0.5,
#'                                         design = bootstrap_rep_design))),
#'     'SE (linearization)' = c(SE(svytotal(x = ~ y1, design = multistage_srswor_design)),
#'                              SE(svymean(x = ~ y1, design = multistage_srswor_design)),
#'                              SE(svyquantile(x = ~ y1, quantile = 0.5,
#'                                             design = multistage_srswor_design))),
#'     check.names = FALSE
#'   )
#'
#' # Example 2: A multistage-sample,
#' # first stage selected with unequal probabilities without replacement
#' # second stage selected with simple random sampling without replacement
#'
#'   data("library_multistage_sample", package = "svrep")
#'
#'   multistage_pps <- svydesign(data = library_multistage_sample,
#'                               ids = ~ PSU_ID + SSU_ID,
#'                               fpc = ~ PSU_SAMPLING_PROB + SSU_SAMPLING_PROB,
#'                               pps = "brewer")
#'
#'   bootstrap_rep_design <- as_bootstrap_design(
#'     multistage_pps, replicates = 500,
#'     samp_method_by_stage = c("PPSWOR", "SRSWOR")
#'   )
#'
#'   ## Compare std. error estimates from bootstrap versus linearization
#'   data.frame(
#'       'Statistic' = c('total', 'mean'),
#'       'SE (bootstrap)' = c(
#'           SE(svytotal(x = ~ TOTCIR, na.rm = TRUE,
#'                       design = bootstrap_rep_design)),
#'           SE(svymean(x = ~ TOTCIR, na.rm = TRUE,
#'                      design = bootstrap_rep_design))),
#'       'SE (linearization)' = c(
#'           SE(svytotal(x = ~ TOTCIR, na.rm = TRUE,
#'                       design = multistage_pps)),
#'           SE(svymean(x = ~ TOTCIR, na.rm = TRUE,
#'                      design = multistage_pps))),
#'       check.names = FALSE
#'   )

as_bootstrap_design <- function(design,
                                type = "Rao-Wu-Yue-Beaumont",
                                replicates = 500,
                                compress = TRUE,
                                mse = getOption("survey.replicates.mse"),
                                samp_method_by_stage = NULL) {
  UseMethod("as_bootstrap_design", design)
}

#' @export
as_bootstrap_design.survey.design <- function(design,
                                              type = "Rao-Wu-Yue-Beaumont",
                                              replicates = 500,
                                              compress = TRUE,
                                              mse = getOption("survey.replicates.mse"),
                                              samp_method_by_stage = NULL) {

  type <- tolower(type)
  permissible_types <- c("rao-wu-yue-beaumont", "preston", "rao-wu",
                         "canty-davison", "antal-tille")

  if (is.null(type) || is.na(type) || (length(type) != 1) || !all(type %in% permissible_types)) {
    stop("Invalid value of `type`. Must use either 'Rao-Wu-Yue-Beaumont', 'Preston', 'Rao-Wu', 'Antal-Tille', or 'Canty-Davison'.")
  }

  if (type == "rao-wu-yue-beaumont") {
    # Extract information from the survey design object
    is_pps_design <- isTRUE(design$pps)
    samp_unit_ids_by_stage <- design$cluster
    strata_ids_by_stage <- design$strata
    fpcs_by_stage <- design$fpc
    number_of_stages <- ncol(samp_unit_ids_by_stage)

    samp_sizes_by_stage <- as.matrix(fpcs_by_stage$sampsize)
    if (is.null(fpcs_by_stage$popsize)) {
      pop_sizes_by_stage <- matrix(Inf, nrow = nrow(samp_unit_ids_by_stage), ncol = ncol(samp_unit_ids_by_stage))
    } else {
      pop_sizes_by_stage <- as.matrix(fpcs_by_stage$popsize)
    }

    samp_unit_sel_probs_by_stage <- design[['allprob']]
    if (ncol(samp_unit_sel_probs_by_stage) < ncol(pop_sizes_by_stage)) {
      samp_unit_sel_probs_by_stage <- samp_sizes_by_stage/pop_sizes_by_stage
    }

    # Determine which stages were with-replacement
    with_replacement_stages <- apply(X = pop_sizes_by_stage, MARGIN = 2, FUN = function(N) all(is.infinite(N)))
    # For WR stages of non-PPS designs, selection probabilities will be 0
    if (!is_pps_design & all(with_replacement_stages)) {
      samp_unit_sel_probs_by_stage <- samp_sizes_by_stage / pop_sizes_by_stage
    }

    if (is.null(samp_method_by_stage)) {
      if (is_pps_design) {
        samp_method_by_stage <- rep("PPSWOR", times = number_of_stages)
      } else {
        samp_method_by_stage <- rep("SRSWOR", times = number_of_stages)
        samp_method_by_stage[with_replacement_stages] <- "SRSWR"
      }
    }

    adjustment_factors <- make_rwyb_bootstrap_weights(
      num_replicates = replicates,
      samp_unit_ids = samp_unit_ids_by_stage,
      strata_ids = strata_ids_by_stage,
      samp_unit_sel_probs = samp_unit_sel_probs_by_stage,
      samp_method_by_stage = samp_method_by_stage,
      allow_final_stage_singletons = TRUE,
      output = "factors"
    )

    rep_design <- survey::svrepdesign(
      variables = design$variables,
      weights = weights(design, type = "sampling"),
      repweights = adjustment_factors, combined.weights = FALSE,
      compress = compress, mse = mse,
      scale = 1/replicates, rscales = rep(1, times = replicates),
      type = "bootstrap"
    )

  }
  
  if (type == "antal-tille") {
    
    adjustment_factors <- make_doubled_half_bootstrap_weights(
      num_replicates       = replicates,
      samp_unit_ids        = design$cluster[, 1, drop = TRUE],
      strata_ids           = design$strata[, 1, drop = TRUE],
      samp_unit_sel_probs  = design$allprob[, 1, drop = TRUE],
      output = "factors"
    )
    
    rep_design <- survey::svrepdesign(
      variables  = design$variables,
      weights    = weights(design, type = "sampling"),
      repweights = adjustment_factors, combined.weights = FALSE,
      compress   = compress, 
      mse        = mse,
      scale      = 1/replicates, 
      rscales    = rep(1, times = replicates),
      type       = "bootstrap"
    )
    
  }

  if (!type %in% c("rao-wu-yue-beaumont", "antal-tille")) {
    type <- switch(type,
                   'preston' = "mrbbootstrap",
                   'mrbootstrap' = "mrbootstrap",
                   'rao-wu' = "subbootstrap",
                   'subbootstrap' = "subbootstrap",
                   'canty-davison' = "bootstrap",
                   'bootstrap' = "bootstrap")
    rep_design <- survey::as.svrepdesign(design = design,
                                         type = type,
                                         compress = compress,
                                         mse = mse,
                                         replicates = replicates)
  }

  if (inherits(design, 'tbl_svy') && ('package:srvyr' %in% search())) {
    rep_design <- srvyr::as_survey_rep(
      rep_design
    )
  }

  rep_design$call <- sys.call(which = -1)

  return(rep_design)
}

#' @export
as_bootstrap_design.DBIsvydesign <- function(design,
                                             type = "Rao-Wu-Yue-Beaumont",
                                             replicates = 500,
                                             compress = TRUE,
                                             mse = getOption("survey.replicates.mse"),
                                             samp_method_by_stage = NULL) {

  rep_design <- NextMethod(design)

  # Replace 'variables' with a database connection
  # and make the object have the appropriate class
  rep_design$variables <- NULL
  if (design$db$dbtype == "ODBC") {
    stop("'RODBC' no longer supported. Use the odbc package")
  } else {
    db <- DBI::dbDriver(design$db$dbtype)
    dbconn <- DBI::dbConnect(db, design$db$dbname)
  }
  rep_design$db <- list(
    dbname = design$db$dbname, tablename = design$db$tablename,
    connection = dbconn,
    dbtype = design$db$dbtype
  )
  class(rep_design) <- c(
    "DBIrepdesign", "DBIsvydesign",
    setdiff(class(rep_design), c("DBIrepdesign", "DBIsvydesign"))
  )

  rep_design$call <- sys.call(which = -1)

  return(rep_design)
}
