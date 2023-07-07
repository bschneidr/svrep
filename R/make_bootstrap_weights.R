#' Create bootstrap replicate weights for a general survey design,
#' using the Rao-Wu-Yue-Beaumont bootstrap method
#' @description Creates bootstrap replicate weights for a multistage stratified sample design
#' using the method of Beaumont and Émond (2022), which is a generalization of the Rao-Wu-Yue bootstrap. \cr \cr
#' The design may have different sampling methods used at different stages.
#' Each stage of sampling may potentially use unequal probabilities (with or without replacement)
#' and may potentially use Poisson sampling.
#' @param num_replicates Positive integer giving the number of bootstrap replicates to create
#' @param samp_unit_ids Matrix or data frame of sampling unit IDs for each stage of sampling
#' @param strata_ids Matrix or data frame of strata IDs for each sampling unit at each stage of sampling
#' @param samp_unit_sel_probs Matrix or data frame of selection probabilities for each sampling unit
#' at each stage of sampling.
#' @param samp_method_by_stage A vector with length equal to the number of stages of sampling,
#' corresponding to the number of columns in \code{samp_unit_ids}.
#' This describes the method of sampling used at each stage.
#' Each element should be one of the following: \cr
#' \itemize{
#'  \item{"SRSWOR"}{ - Simple random sampling, without replacement}
#'  \item{"SRSWR"}{ - Simple random sampling, with replacement}
#'  \item{"PPSWOR"}{ - Unequal probabilities of selection, without replacement}
#'  \item{"PPSWR"}{ - Unequal probabilities of selection, with replacement}
#'  \item{"Poisson"}{ -  Poisson sampling: each sampling unit is selected into the sample at most once, with potentially different probabilities of inclusion for each sampling unit.}
#' }
#' @param allow_final_stage_singletons Logical value indicating whether to allow
#' non-certainty singleton strata at the final sampling stage (rather than throw an error message). \cr
#' If \code{TRUE}, the sampling unit in a non-certainty singleton stratum will have its final-stage adjustment factor
#' calculated as if it was selected with certainty at the final stage (i.e., its adjustment factor will be 1),
#' and then its final bootstrap weight will be calculated by combining this adjustment factor
#' with its final-stage selection probability.
#' @param output Either \code{"weights"} (the default) or \code{"factors"}.
#' Specifying \code{output = "factors"} returns a matrix of replicate adjustment factors which can later be multiplied by
#' the full-sample weights to produce a matrix of replicate weights. Specifying \code{output = "weights"}
#' returns the matrix of replicate weights, where the full-sample weights are inferred using \code{samp_unit_sel_probs}.
#' @details
#' Beaumont and Émond (2022) describe a general algorithm for forming bootstrap replicate weights
#' for multistage stratified samples, based on the method of Rao-Wu-Yue, with extensions
#' to sampling without replacement and use of unequal probabilities of selection
#' (i.e., sampling with probability proportional to size) as well as Poisson sampling. These methods
#' are guaranteed to produce nonnegative replicate weights and provide design-unbiased and design-consistent variance estimates for totals,
#' for designs where sampling uses one or more of the following methods:
#' \itemize{
#'  \item{"SRSWOR"}{ - Simple random sampling, without replacement}
#'  \item{"SRSWR"}{ - Simple random sampling, with replacement}
#'  \item{"PPSWR"}{ - Unequal probabilities of selection, with replacement}
#'  \item{"Poisson"}{ -  Poisson sampling: each sampling unit is selected into the sample at most once, with potentially different probabilities of inclusion for each sampling unit.}
#' }
#'
#' For designs where at least one stage's strata have sampling without replacement with unequal probabilities of selection ("PPSWOR"),
#' the bootstrap method of Beaumont and Émond (2022) is guaranteed to produce nonnegative weights, but is not design-unbiased,
#' since the method only approximates the joint selection probabilities which would be needed for unbiased estimation.
#' \cr \cr
#' Unless any stages use simple random sampling without replacement, the resulting bootstrap replicate weights
#' are guaranteed to all be strictly positive, which may be useful for calibration or analyses of domains with small sample sizes.
#' If any stages use simple random sampling without replacement, it is possible for some replicate weights to be zero.
#' \cr \cr
#' If there is survey nonresponse, it may be useful to represent the response/nonresponse as an additional
#' stage of sampling, where sampling is conducted with Poisson sampling
#' where each unit's "selection probability" at that stage is its response propensity (which typically has to be estimated).
#' @seealso If the survey design can be accurately represented using \code{\link[survey]{svydesign}},
#'  then it is easier to simply use \code{\link[svrep]{as_bootstrap_design}} with argument \code{type = "Rao-Wu-Yue-Beaumont"}.
#' @references
#' Beaumont, J.-F.; Émond, N. (2022).
#' "A Bootstrap Variance Estimation Method for Multistage Sampling and Two-Phase Sampling When Poisson Sampling Is Used at the Second Phase."
#' \strong{Stats}, \emph{5}: 339–357.
#' https://doi.org/10.3390/stats5020019
#'
#' Rao, J.N.K.; Wu, C.F.J.; Yue, K. (1992).
#' "Some recent work on resampling methods for complex surveys."
#' \strong{Surv. Methodol.}, \emph{18}: 209–217.
#'
#' @return A matrix of with the same number of rows as \code{samp_unit_ids}
#' and the number of columns equal to the value of the argument \code{num_replicates}.
#' Specifying \code{output = "factors"} returns a matrix of replicate adjustment factors which can later be multiplied by
#' the full-sample weights to produce a matrix of replicate weights.
#' Specifying \code{output = "weights"} returns the matrix of replicate weights,
#' where the full-sample weights are inferred using \code{samp_unit_sel_probs}.
#' @export
#' @seealso Use \code{\link[svrep]{estimate_boot_reps_for_target_cv}} to help choose the number of bootstrap replicates.
#' @examples
#' \dontrun{
#'  library(survey)
#'
#'  # Example 1: A multistage sample with two stages of SRSWOR
#'
#'      ## Load an example dataset from a multistage sample, with two stages of SRSWOR
#'      data("mu284", package = 'survey')
#'      multistage_srswor_design <- svydesign(data = mu284,
#'                                            ids = ~ id1 + id2,
#'                                            fpc = ~ n1 + n2)
#'
#'      ## Create bootstrap replicate weights
#'      set.seed(2022)
#'      bootstrap_replicate_weights <- make_rwyb_bootstrap_weights(
#'        num_replicates = 5000,
#'        samp_unit_ids = multistage_srswor_design$cluster,
#'        strata_ids = multistage_srswor_design$strata,
#'        samp_unit_sel_probs = multistage_srswor_design$fpc$sampsize /
#'                              multistage_srswor_design$fpc$popsize,
#'        samp_method_by_stage = c("SRSWOR", "SRSWOR")
#'      )
#'
#'      ## Create a replicate design object with the survey package
#'      bootstrap_rep_design <- svrepdesign(
#'        data = multistage_srswor_design$variables,
#'        repweights = bootstrap_replicate_weights,
#'        weights = weights(multistage_srswor_design, type = "sampling"),
#'        type = "bootstrap"
#'      )
#'
#'      ## Compare std. error estimates from bootstrap versus linearization
#'      data.frame(
#'        'Statistic' = c('total', 'mean', 'median'),
#'        'SE (bootstrap)' = c(SE(svytotal(x = ~ y1, design = bootstrap_rep_design)),
#'                             SE(svymean(x = ~ y1, design = bootstrap_rep_design)),
#'                             SE(svyquantile(x = ~ y1, quantile = 0.5,
#'                                            design = bootstrap_rep_design))),
#'        'SE (linearization)' = c(SE(svytotal(x = ~ y1, design = multistage_srswor_design)),
#'                                 SE(svymean(x = ~ y1, design = multistage_srswor_design)),
#'                                 SE(svyquantile(x = ~ y1, quantile = 0.5,
#'                                                design = multistage_srswor_design))),
#'        check.names = FALSE
#'      )
#'
#'  # Example 2: A single-stage sample selected with unequal probabilities, without replacement
#'
#'      ## Load an example dataset of U.S. counties states with 2004 Presidential vote counts
#'      data("election", package = 'survey')
#'      pps_wor_design <- svydesign(data = election_pps,
#'                                  pps = "overton",
#'                                  fpc = ~ p, # Inclusion probabilities
#'                                  ids = ~ 1)
#'
#'      ## Create bootstrap replicate weights
#'      set.seed(2022)
#'      bootstrap_replicate_weights <- make_rwyb_bootstrap_weights(
#'        num_replicates = 5000,
#'        samp_unit_ids = pps_wor_design$cluster,
#'        strata_ids = pps_wor_design$strata,
#'        samp_unit_sel_probs = pps_wor_design$prob,
#'        samp_method_by_stage = c("PPSWOR")
#'      )
#'
#'      ## Create a replicate design object with the survey package
#'      bootstrap_rep_design <- svrepdesign(
#'        data = pps_wor_design$variables,
#'        repweights = bootstrap_replicate_weights,
#'        weights = weights(pps_wor_design, type = "sampling"),
#'        type = "bootstrap"
#'      )
#'
#'      ## Compare std. error estimates from bootstrap versus linearization
#'      data.frame(
#'        'Statistic' = c('total', 'mean'),
#'        'SE (bootstrap)' = c(SE(svytotal(x = ~ Bush, design = bootstrap_rep_design)),
#'                             SE(svymean(x = ~ I(Bush/votes),
#'                                        design = bootstrap_rep_design))),
#'        'SE (Overton\'s PPS approximation)' = c(SE(svytotal(x = ~ Bush,
#'                                                            design = pps_wor_design)),
#'                                                SE(svymean(x = ~ I(Bush/votes),
#'                                                           design = pps_wor_design))),
#'        check.names = FALSE
#'      )
#' }

make_rwyb_bootstrap_weights <- function(num_replicates = 100,
                                        samp_unit_ids, strata_ids,
                                        samp_unit_sel_probs,
                                        samp_method_by_stage = rep("PPSWOR", times = ncol(samp_unit_ids)),
                                        allow_final_stage_singletons = TRUE,
                                        output = "weights") {

  number_of_stages <- ncol(samp_unit_ids)
  number_of_ultimate_units <- nrow(samp_unit_ids)

  # Check validity of inputs
  if (!is.matrix(samp_unit_sel_probs)) {
    samp_unit_sel_probs <- as.matrix(samp_unit_sel_probs)
  }
  if ((ncol(strata_ids) != number_of_stages) || (nrow(strata_ids) != number_of_ultimate_units)) {
    stop("`strata_ids` must have the same number of rows and columns as `samp_unit_ids`.")
  }
  if (any(is.na(strata_ids))) {
    stop("`strata_ids` should not have any missing values.")
  }
  if (any(is.na(samp_unit_ids))) {
    stop("`samp_unit_ids` should not have any missing values.")
  }
  if ((ncol(samp_unit_sel_probs) != number_of_stages) || (nrow(samp_unit_sel_probs) != number_of_ultimate_units)) {
    stop("`samp_unit_sel_probs` must have the same number of rows and columns as `samp_unit_ids`.")
  }
  if (length(samp_method_by_stage) != number_of_stages) {
    stop("The length of `samp_method_by_stage` should match the number of columns of `samp_unit_ids`.")
  }
  samp_method_by_stage <- toupper(samp_method_by_stage)
  if (!all(samp_method_by_stage %in% c("SRSWR", "SRSWOR", "PPSWR", "PPSWOR", "POISSON"))) {
    stop('Each element of `samp_method_by_stage` must be one of the following: "SRSWR", "SRSWOR", "PPSWR", "PPSWOR", or "Poisson"')
  }
  if (any(is.na(samp_unit_sel_probs[,samp_method_by_stage != "SRSWR"]))) {
    stop("For stages where sampling is not 'SRSWR', the corresponding column of `samp_unit_sel_probs` should not have any missing values.")
  }
  if ((length(num_replicates) != 1) || !is.numeric(num_replicates) || (!num_replicates > 0)) {
    stop("Must specify a single, positive number for the argument `num_replicates`.")
  }

  # Initialize 3-dimensional array: ultimate unit X stage X replicates
  adjustment_factors <- matrix(
    data = 1,
    nrow = number_of_ultimate_units,
    ncol = num_replicates
  )

  # Make sure each stage's sampling units are nested within strata
  # and each stage's sampling units are nested
  # within previous stage sampling units
  samp_unit_ids[,1] <- interaction(strata_ids[, 1, drop = TRUE],
                                   samp_unit_ids[, 1, drop = TRUE],
                                   sep = " | ", drop = TRUE) |>
    as.numeric()
  stage <- 2L
  while (stage <= number_of_stages) {
    strata_ids[,stage] <- interaction(
      samp_unit_ids[, stage-1L, drop=TRUE],
      strata_ids[, stage, drop=TRUE],
      sep = " | ", drop = TRUE
    ) |> as.numeric()
    samp_unit_ids[,stage] <- interaction(
      strata_ids[, stage, drop = TRUE],
      samp_unit_ids[, stage, drop = TRUE],
      sep = " | ", drop = TRUE
    ) |> as.numeric()
    stage <- stage + 1L
  }

  # Calculate replicate adjustment factors for each stage
  stage <- 1L
  while (stage <= number_of_stages) {

    samp_units <- samp_unit_ids[,stage,drop=TRUE]
    sel_probs <- samp_unit_sel_probs[,stage,drop=TRUE]
    strata <- strata_ids[,stage,drop=TRUE]

    # Determine each stratum's sample size and resample size
    distinct_strata_ids <- unique(strata)
    H <- length(distinct_strata_ids)

    for (h in seq_len(H)) {
      stratum_indices <- which(strata == distinct_strata_ids[h])
      n_h <- samp_units[stratum_indices] |> unique() |> length()
      m_h <- n_h - 1

      # Get each stratum's list of sampling units and their probabilities
      distinct_samp_units <- unique(samp_units[stratum_indices])
      distinct_sel_probs <- sel_probs[stratum_indices][!duplicated(samp_units[stratum_indices])]
      certainty_flags <- (distinct_sel_probs >= 1)

      is_noncertainty_singleton <- (n_h == 1) && !all(certainty_flags)

      if (is_noncertainty_singleton) {
        if ((stage < number_of_stages) | !allow_final_stage_singletons) {
          if (samp_method_by_stage[stage] != "POISSON") {
            error_msg <- sprintf(
              paste("Cannot form bootstrap adjustment factors for a stratum at stage %s, ",
                    "since the stratum has only one sampling unit, ",
                    "and that single sampling unit was not selected with certainty.",
                    ifelse(
                      stage == number_of_stages,
                      paste0(
                        " Setting `allow_final_stage_singletons = TRUE` will avoid this error message",
                        " by calculating that sampling unit's final-stage adjustment factor as if it was selected with certainty at the final stage,",
                        " and then calculating the bootstrap weight using this adjustment factor combined with the final-stage selection probability."
                      ),
                      ""
                    ),
                    sep = ""),
              stage
            )
            stop(error_msg)
          }
        }
      }

      # For later stages of sampling, get the previous stage selection probability of the higher-level sampling unit
      if (stage > 1) {
        row_to_select <- stratum_indices[1]
        prior_stage_probs <- samp_unit_sel_probs[row_to_select, seq_len(stage-1), drop=TRUE]
        sel_prob_of_higher_unit_from_previous_stage <- Reduce(prior_stage_probs, f = `*`)
      } else {
        sel_prob_of_higher_unit_from_previous_stage <- 1
      }

      if (samp_method_by_stage[stage] != "POISSON") {

        # Determine 'multiplicities' (number of times each PSU is resampled)
        multiplicities <- stats::rmultinom(
          n = num_replicates,
          size = m_h,
          prob = rep(n_h, times = n_h)^(-1)
        )

        # Calculate adjustment from Equation (21) of Beaumont & Emond 2022 in the cases of PPSWOR or SRSWOR,
        #                        or Equation (20) in the case of SRSWR (which is the Rao-Wu-Yue method)
        mstar_k <- multiplicities
        certainty_flags <- certainty_flags
        is_singleton_stratum <- n_h == 1

        if (is_singleton_stratum) {
          a_k <- 1
        } else if (samp_method_by_stage[stage] != "SRSWR") {
          pi_k <- distinct_sel_probs
          a_k <- 1 - sqrt((m_h*(1-pi_k))/(n_h - 1)) + sqrt((m_h*(1-pi_k))/(n_h - 1)) * (n_h/m_h) * mstar_k
        } else if (samp_method_by_stage[stage] == "SRSWR") {
          a_k <- 1 - sqrt(m_h/(n_h - 1)) + sqrt(m_h/(n_h - 1)) * (n_h/m_h) * mstar_k
        }
        a_beaumont_emond <- as.matrix(a_k)
      }

      if (samp_method_by_stage[stage] == "POISSON") {
        delta <- 1 - distinct_sel_probs
        a_beaumont_emond <- sapply(
          X = delta,
          FUN = function(scale_param) {
            stats::rgamma(n = num_replicates,
                          shape = 1/scale_param,
                          scale = scale_param)
          }
        ) |> t()
      }

      # For PPSWOR, "calibrate" adjustments to sum to actual sample size, using Equation (24) of Beaumont & Emond 2022
      if (samp_method_by_stage[stage] == "PPSWOR") {

        n_noncertainty <- sum(!certainty_flags)

        a_sum <- colSums(a_beaumont_emond[!certainty_flags,,drop=FALSE])
        adjustment_factor <- n_noncertainty / a_sum

        a_k_cal <- a_beaumont_emond
        a_k_cal[!certainty_flags,] <- t(
          apply(X = a_beaumont_emond[!certainty_flags,,drop=FALSE],
                MARGIN = 1, FUN = function(x) x*adjustment_factor)
        )

        a_beaumont_emond_cal <- a_k_cal
      } else {
        a_beaumont_emond_cal <- a_beaumont_emond
      }

      # For all stages other than the first, update the adjustment factor
      # based on the selection probability of the higher-level unit from previous stage of sampling

      if (stage > 1) {
        a_prelim <- a_beaumont_emond_cal
        Delta_prev_kk <- 1 - sel_prob_of_higher_unit_from_previous_stage
        a_updated <- 1 - sqrt((1 - Delta_prev_kk)/(1 + Delta_prev_kk)) + sqrt((1 - Delta_prev_kk)/(1 + Delta_prev_kk)) * a_prelim
        a_final <- a_updated
      } else {
        a_final <- a_beaumont_emond_cal
      }

      # For each row of the data, get the adjustment factor for the sampling unit it belongs to
      for (hi in seq_along(stratum_indices)) {
        k <- (distinct_samp_units == samp_units[stratum_indices][hi])
        i <- stratum_indices[hi]
        if (length(a_final) == 1) {
          adjustment_factors[i,] <- adjustment_factors[i,,drop=TRUE] * as.vector(a_final)
        } else {
          adjustment_factors[i,] <- adjustment_factors[i,,drop=TRUE] * as.vector(a_final[k,])
        }
      }
    }

    stage <- stage + 1L
  }

  result <- adjustment_factors

  if (output == "weights") {
    # Calculate overall sampling weight for a unit
    # as the inverse of product of selection probabilities from all stages
    overall_sampling_weights <- apply(X = samp_unit_sel_probs,
                                      MARGIN = 1,
                                      FUN = Reduce, f = `*`) ^ (-1)


    # Create replicate weights by multiplying adjustment factors by sampling weights
    replicate_weights <- apply(X = adjustment_factors,
                               MARGIN = 2,
                               FUN = function(rep_factor) rep_factor * overall_sampling_weights)
    result <- replicate_weights
  }

  return(result)
}
