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
#'  \item{\strong{"Rao-Wu-Yue-Beaumont"} (the default): }{\cr
#'  The bootstrap method of Beaumont and Émond (2022), which is a generalization of the Rao-Wu-Yue bootstrap,
#'  and is applicable to a wide variety of designs, including single-stage and multistage stratified designs.
#'  The design may have different sampling methods used at different stages.
#'  Each stage of sampling may potentially be PPS (i.e., use unequal probabilities), with or without replacement,
#'  and may potentially use Poisson sampling. \cr \cr
#'  For a stratum with a fixed sample size of \eqn{n} sampling units, resampling in each replicate resamples \eqn{(n-1)} sampling units with replacement.}
#'  \item{\strong{"Rao-Wu"}: }{\cr
#'  The basic Rao-Wu \eqn{(n-1)} bootstrap method, which is only applicable to single-stage designs or
#'  multistage designs where the first-stage sampling fractions are small (and can thus be ignored).
#'  Accommodates stratified designs. All sampling within a stratum must be simple random sampling with or without replacement,
#'  although the first-stage sampling is effectively treated as sampling without replacement.}
#'  \item{\strong{"Preston"}: }{\cr
#'  Preston's multistage rescaled bootstrap, which is applicable to single-stage designs or multistage designs
#'  with arbitrary sampling fractions. Accommodates stratified designs. All sampling within a stratum must be
#'  simple random sampling with or without replacement.}
#'  \item{\strong{"Canty-Davison"}: }{\cr
#'  The Canty-Davison bootstrap, which is only applicable to single-stage designs, with arbitrary sampling fractions.
#'  Accommodates stratified designs. All sampling with a stratum must be simple random sampling with or without replacement.}
#' }
#' @param replicates Number of bootstrap replicates (should be as large as possible, given computer memory/storage limitations).
#' A commonly-recommended default is 500.
#' @param compress Use a compressed representation of the replicate weights matrix.
#' This reduces the computer memory required to represent the replicate weights and has no
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
#' @seealso Use \code{\link[svrep]{estimate_boot_reps_for_target_cv}} to help choose the number of bootstrap replicates.
#'
#'          For some complex designs, one can use \code{\link[svrep]{make_rwyb_bootstrap_weights}} to create
#'          Rao-Wu-Yue-Beaumont bootstrap weights or adjustment factors given information for each stage of sampling
#'          (the type of sampling, strata IDs, cluster IDs, selection probabilities, etc.).
#'
#'          \cr \cr
#'
#'          For systematic samples, one-PSU-per-stratum designs, or other especially complex sample designs,
#'           one can use the generalized survey bootstrap method. See \code{\link[svrep]{make_gen_boot_factors}}.
#'
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
#'   bootstrap_rep_design <- as_bootstrap_design(multistage_pps, replicates = 100)
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

as_bootstrap_design <- function(design, type = "Rao-Wu-Yue-Beaumont", replicates = 500, compress = TRUE, mse = getOption("survey.replicates.mse")) {
  UseMethod("as_bootstrap_design", design)
}

#' @export
as_bootstrap_design.survey.design <- function(design,
                                              type = "Rao-Wu-Yue-Beaumont",
                                              replicates = 500,
                                              compress = TRUE,
                                              mse = getOption("survey.replicates.mse")) {

  type <- tolower(type)

  if (type == "rao-wu-yue-beaumont") {
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
    with_replacement_stages <- apply(X = pop_sizes_by_stage, MARGIN = 2, FUN = function(N) all(is.infinite(N)))
    is_pps_design <- isTRUE(design$pps)

    samp_method_by_stage <- rep("PPSWOR", times = number_of_stages)
    if (is_pps_design) {
      samp_unit_sel_probs_by_stage <- design[['allprob']]
    } else {
      samp_unit_sel_probs_by_stage <- samp_sizes_by_stage / pop_sizes_by_stage
      samp_method_by_stage[with_replacement_stages] <- "SRSWR"
    }

    adjustment_factors <- make_rwyb_bootstrap_weights(
      num_replicates = replicates,
      samp_unit_ids = samp_unit_ids_by_stage,
      strata_ids = strata_ids_by_stage,
      samp_unit_sel_probs = samp_unit_sel_probs_by_stage,
      samp_method_by_stage = samp_method_by_stage,
      allow_singletons_at_final_stage = TRUE,
      output = "factors"
    )

    rep_design <- survey::svrepdesign(
      variables = design$variables,
      weights = weights(design, type = "sampling"),
      repweights = adjustment_factors * weights(design, type = "sampling"), combined.weights = TRUE,
      compress = compress, mse = mse,
      scale = 1/replicates, rscales = rep(1, times = replicates),
      type = "bootstrap"
    )

  }

  if (type != "rao-wu-yue-beaumont") {
    type <- switch(type,
                   'preston' = "mrbbootstrap",
                   'mrbootstrap' = "mrbootstrap",
                   'rao-wu' = "subbootstrap",
                   'subbootstrap' = "subbootstrap",
                   'canty-davison' = "bootstrap",
                   'bootstrap' = "canty-davison")
    rep_design <- survey::as.svrepdesign(design = design,
                                         type = type,
                                         compress = compress,
                                         mse = mse)
  }

  rep_design$call <- sys.call(which = -1)

  return(rep_design)
}
