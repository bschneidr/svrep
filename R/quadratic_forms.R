#' @title Represent a variance estimator as a quadratic form
#' @description Common variance estimators for estimated population totals can be represented as a quadratic form.
#' Given a choice of variance estimator and information about the sample design,
#' this function constructs the matrix of the quadratic form.
#' \cr \cr
#' In notation, let
#' \eqn{v(\hat{Y}) = \mathbf{\breve{y}}^{\prime}\mathbf{\Sigma}\mathbf{\breve{y}}},
#' where \eqn{\breve{y}} is the vector of weighted values, \eqn{y_i/\pi_i, \space i=1,\dots,n}.
#' This function constructs the \eqn{n \times n} matrix of the quadratic form, \eqn{\mathbf{\Sigma}}.
#' @param variance_estimator The name of the variance estimator
#' whose quadratic form matrix should be created. See the section "Variance Estimators" below.
#' Options include:
#' \itemize{
#'   \item \strong{"Yates-Grundy"}: The Yates-Grundy variance estimator based on
#'     first-order and second-order inclusion probabilities. If this is used,
#'     the argument \code{joint_probs} must also be used.
#'   \item \strong{"Horvitz-Thompson"}: The Horvitz-Thompson variance estimator based on
#'     first-order and second-order inclusion probabilities. If this is used,
#'     the argument \code{joint_probs} must also be used.
#'   \item \strong{"Stratified Multistage SRS"}: The usual stratified multistage variance estimator
#'     based on estimating the variance of cluster totals within strata at each stage.
#'     If this option is used, then it is necessary to also use the arguments
#'     \code{strata_ids}, \code{cluster_ids}, \code{strata_pop_sizes}, and \code{strata_pop_sizes}.
#'   \item \strong{"Ultimate Cluster"}: The usual variance estimator based on estimating
#'     the variance of first-stage cluster totals within first-stage strata.
#'     If this option is used, then it is necessary to also use the arguments
#'     \code{strata_ids}, \code{cluster_ids}, \code{strata_pop_sizes}.
#'     Optionally, to use finite population correction factors, one can also use the argument \code{strata_pop_sizes}.
#'   \item \strong{"Deville-1"}: A variance estimator for unequal-probability
#'     sampling without replacement, described in Matei and Tillé (2005)
#'     as "Deville 1". If this option is used, then it is necessary to also use the arguments
#'     \code{strata_ids}, \code{cluster_ids}, and \code{probs}.
#'   \item \strong{"Deville-2"}:  A variance estimator for unequal-probability
#'     sampling without replacement, described in Matei and Tillé (2005)
#'     as "Deville 2". If this option is used, then it is necessary to also use the arguments
#'     \code{strata_ids}, \code{cluster_ids}, and \code{probs}.
#'   \item \strong{"SD1"}: The non-circular successive-differences variance estimator described by Ash (2014),
#'     sometimes used for variance estimation for systematic sampling.
#'   \item \strong{"SD2"}: The circular successive-differences variance estimator described by Ash (2014).
#'     This estimator is the basis of the "successive-differences replication" estimator commonly used
#'     for variance estimation for systematic sampling.
#'   \item \strong{"BOSB"}: \cr The kernel-based variance estimator proposed by
#'     Breidt, Opsomer, and Sanchez-Borrego (2016) for use with systematic samples
#'     or other finely stratified designs. Uses the Epanechnikov kernel
#'     with the bandwidth automatically chosen to result in the smallest possible
#'     nonempty kernel window.
#'   \item \strong{"Deville-Tille"}: The estimator of Deville and Tillé (2005),
#'     developed for balanced sampling using the cube method.
#'   \item \strong{"Beaumont-Emond"}: The variance estimator of Beaumont and Emond (2022)
#'     for multistage unequal-probability sampling without replacement.
#' }
#' @param probs Required if \code{variance_estimator} equals \code{"Deville-1"}, \code{"Deville-2"}, or \code{"Deville-Tille"}.
#' This should be a matrix or data frame of sampling probabilities.
#' If there are multiple stages of sampling,
#' then \code{probs} can have multiple columns,
#' with one column for each level of sampling to be accounted for by the variance estimator.
#' @param joint_probs Only used if \code{variance_estimator = "Horvitz-Thompson"} or \code{variance_estimator = "Yates-Grundy"}.
#' This should be a matrix of joint inclusion probabilities.
#' Element \code{[i,i]} of the matrix is the first-order inclusion probability of unit \code{i},
#' while element \code{[i,j]} is the joint inclusion probability of units \code{i} and \code{j}.
#' @param cluster_ids Required unless \code{variance_estimator} equals \code{"Horvitz-Thompson"} or \code{"Yates-Grundy"}.
#' This should be a matrix or data frame of cluster IDs. If there are multiple stages of sampling,
#' then \code{cluster_ids} can have multiple columns,
#' with one column for each level of sampling to be accounted for by the variance estimator.
#' @param strata_ids Required if \code{variance_estimator} equals \code{"Stratified Multistage SRS"}
#' or \code{"Ultimate Cluster"}.
#' This should be a matrix or data frame of strata IDs. If there are multiple stages of sampling,
#' then \code{strata_ids} can have multiple columns,
#' with one column for each level of sampling to be accounted for by the variance estimator.
#' @param strata_pop_sizes Required if \code{variance_estimator} equals \code{"Stratified Multistage SRS"},
#' but can optionally be used if \code{variance_estimator} equals \code{"Ultimate Cluster"}, \code{"SD1"}, or \code{"SD2"}.
#' If there are multiple stages of sampling,
#' then \code{strata_pop_sizes} can have multiple columns,
#' with one column for each level of sampling to be accounted for by the variance estimator.
#' @param sort_order Required if \code{variance_estimator} equals \code{"SD1"} or \code{"SD2"}.
#' This should be a vector that orders the rows of data into the order used for sampling.
#' @param aux_vars Required if \code{variance_estimator} equals \code{"Deville-Tille"}.
#' A matrix of auxiliary variables.
#' @section Variance Estimators:
#' See \link[svrep]{variance-estimators} for a
#' description of each variance estimator.
#' @section Arguments required for each variance estimator:
#' Below are the arguments that are required or optional for each variance estimator.
#'
#' | variance_estimator       | probs      | joint_probs | cluster_ids | strata_ids  | strata_pop_sizes | sort_order | aux_vars |
#' | ------------------------ |-----------:| -----------:| -----------:| -----------:| ----------------:|-----------:|---------:|
#' | Stratified Multistage SRS|            |             | Required    | Required    | Required         |            |          |
#' | Ultimate Cluster         |            |             | Required    | Required    | Optional         |            |          |
#' | SD1                      |            |             | Required    | Optional    | Optional         | Required   |          |
#' | SD2                      |            |             | Required    | Optional    | Optional         | Required   |          |
#' | Deville-1                | Required   |             | Required    | Optional    |                  |            |          |
#' | Deville-2                | Required   |             | Required    | Optional    |                  |            |          |
#' | Beaumont-Emond           | Required   |             | Required    | Optional    |                  |            |          |
#' | Deville-Tille            | Required   |             | Required    | Optional    |                  |            | Required |
#' | BOSB                     |            |             | Required    | Optional    |                  |            | Required |
#' | Yates-Grundy             |            | Required    |             |             |                  |            |          |
#' | Horvitz-Thompson         |            | Required    |             |             |                  |            |          |
#' @md
#' @return The matrix of the quadratic form representing the variance estimator.
#' @export
#' @seealso
#' See \link[svrep]{variance-estimators} for a
#' description of each variance estimator.
#'
#' For a two-phase design, the function
#' \link[svrep]{make_twophase_quad_form} combines
#' the quadratic form matrix from each phase.
#' @examples
#' \dontrun{
#' # Example 1: The Horvitz-Thompson Estimator
#'   library(survey)
#'   data("election", package = "survey")
#'
#'   ht_quad_form_matrix <- make_quad_form_matrix(variance_estimator = "Horvitz-Thompson",
#'                                                joint_probs = election_jointprob)
#'   ##_ Produce variance estimate
#'   wtd_y <- as.matrix(election_pps$wt * election_pps$Bush)
#'   t(wtd_y) %*% ht_quad_form_matrix %*% wtd_y
#'
#'   ##_ Compare against result from 'survey' package
#'   svytotal(x = ~ Bush,
#'            design = svydesign(data=election_pps,
#'                               variance = "HT",
#'                               pps = ppsmat(election_jointprob),
#'                               ids = ~ 1, fpc = ~ p)) |> vcov()
#'
#' # Example 2: Stratified multistage Sample ----
#'
#'   data("mu284", package = 'survey')
#'   multistage_srswor_design <- svydesign(data = mu284,
#'                                         ids = ~ id1 + id2,
#'                                         fpc = ~ n1 + n2)
#'
#'   multistage_srs_quad_form <- make_quad_form_matrix(
#'     variance_estimator = "Stratified Multistage SRS",
#'     cluster_ids = mu284[,c('id1', 'id2')],
#'     strata_ids = matrix(1, nrow = nrow(mu284), ncol = 2),
#'     strata_pop_sizes = mu284[,c('n1', 'n2')]
#'   )
#'
#'   wtd_y <- as.matrix(weights(multistage_srswor_design) * mu284$y1)
#'   t(wtd_y) %*% multistage_srs_quad_form %*% wtd_y
#'
#'   ##_ Compare against result from 'survey' package
#'   svytotal(x = ~ y1, design = multistage_srswor_design) |> vcov()
#'
#' # Example 3: Successive-differences estimator ----
#'
#'   data('library_stsys_sample', package = 'svrep')
#'
#'   sd1_quad_form <- make_quad_form_matrix(
#'     variance_estimator = 'SD1',
#'     cluster_ids = library_stsys_sample[,'FSCSKEY',drop=FALSE],
#'     strata_ids = library_stsys_sample[,'SAMPLING_STRATUM',drop=FALSE],
#'     strata_pop_sizes = library_stsys_sample[,'STRATUM_POP_SIZE',drop=FALSE],
#'     sort_order = library_stsys_sample[['SAMPLING_SORT_ORDER']]
#'   )
#'
#'   wtd_y <- as.matrix(library_stsys_sample[['TOTCIR']] /
#'                       library_stsys_sample$SAMPLING_PROB)
#'   wtd_y[is.na(wtd_y)] <- 0
#'
#'   t(wtd_y) %*% sd1_quad_form %*% wtd_y
#'
#' # Example 4: Deville estimators ----
#'
#'  data('library_multistage_sample', package = 'svrep')
#'
#'  deville_quad_form <- make_quad_form_matrix(
#'      variance_estimator = 'Deville-1',
#'      cluster_ids = library_multistage_sample[,c("PSU_ID", "SSU_ID")],
#'      strata_ids = cbind(rep(1, times = nrow(library_multistage_sample)),
#'                         library_multistage_sample$PSU_ID),
#'      probs = library_multistage_sample[,c("PSU_SAMPLING_PROB",
#'                                           "SSU_SAMPLING_PROB")]
#'  )
#' }
make_quad_form_matrix <- function(variance_estimator = "Yates-Grundy",
                                  probs = NULL,
                                  joint_probs = NULL,
                                  cluster_ids = NULL,
                                  strata_ids = NULL,
                                  strata_pop_sizes = NULL,
                                  sort_order = NULL,
                                  aux_vars = NULL) {

  accepted_variance_estimators <- c(
    "Yates-Grundy", "Horvitz-Thompson",
    "Ultimate Cluster", "Stratified Multistage SRS",
    "SD1", "SD2", "BOSB",
    "Deville-1", "Deville-2", "Beaumont-Emond", "Deville-Tille"
  )

  if (length(variance_estimator) > 1) {
    stop("Can only specify one estimator for `variance_estimator`.")
  }

  if (!variance_estimator %in% accepted_variance_estimators) {
    sprintf("`%s` is not a supported variance estimator, or else there is a typo.",
            variance_estimator) |> stop()
  }

  if (variance_estimator %in% c("Yates-Grundy", "Horvitz-Thompson")) {
    if (is.null(joint_probs)) {
      sprintf("For `variance_estimator='%s'`, must supply a matrix to the argument `joint_probs`.",
              variance_estimator) |>
        stop()
    }
    if (!is.matrix(joint_probs) || any(is.na(joint_probs)) || any(joint_probs < 0)) {
      stop("`joint_probs` must be a matrix of values between 0 and 1, with no missing values.")
    }
    number_of_ultimate_units <- ncol(joint_probs)
    use_sparse_matrix <- FALSE
  }

  # Check inputs and assemble all necessary information
  # for estimators of stratified/clustered designs
  if (variance_estimator %in% c("Stratified Multistage SRS", "Ultimate Cluster", "SD1", "SD2", "Deville-1", "Deville-2", "Deville-Tille", "Beaumont-Emond", "BOSB")) {

    use_sparse_matrix <- TRUE

    # Ensure the minimal set of inputs is supplied
    if (variance_estimator %in% c("Stratified Multistage SRS", "Ultimate Cluster", "Deville-1", "Deville-2", "Deville-Tille", "Beaumont-Emond", "BOSB")) {
      if (is.null(cluster_ids) || is.null(strata_ids)) {
        sprintf(
          "For `variance_estimator='%s'`, must supply a matrix or data frame to both `strata_ids` and `cluster_ids`",
          variance_estimator
        ) |> stop()
      }
    }
    if (variance_estimator %in% c("Deville-1", "Deville-2", "Deville-Tille", "Beaumont-Emond")) {
      if (is.null(probs)) {
        sprintf(
          "For `variance_estimator='%s'`, must supply a matrix or data frame to `probs`.",
          variance_estimator
        ) |> stop()
      }
    }

    if (variance_estimator == "Stratified Multistage SRS") {
      if (is.null(strata_pop_sizes)) {
        stop("For `variance_estimator='Stratified Multistage SRS'`, must supply a matrix or data frame to `strata_pop_sizes`.")
      }
    }
    if (variance_estimator %in% c("SD1", "SD2", "Deville-1", "Deville-2", "Deville-Tille", "Beaumont-Emond", "BOSB")) {
      if (is.null(cluster_ids)) {
        sprintf(
          "For `variance_estimator='%s'`, must supply a matrix or data frame to `cluster_ids`",
          variance_estimator
        ) |> stop()
      }
    }
    if (variance_estimator %in% c("SD1", "SD2")) {
      if (is.null(sort_order)) {
        sprintf(
          "For `variance_estimator='%s'`, must supply a vector to `sort_order`",
          variance_estimator
        ) |> stop()
      }
    }
    if (variance_estimator %in% c("Deville-Tille", "BOSB")) {
      if (missing(aux_vars) || is.null(aux_vars)) {
        sprintf(
          "For `variance_estimator='%s', must supply a matrix to `aux_vars`",
          variance_estimator
        ) |> stop()
      }
      if (!is.matrix(aux_vars)) {
        stop("`aux_vars` must be a matrix.")
      }
    }

    # Determine number of stages and ultimate units,
    # and make clusters/strata explicit if they weren't supplied
    if (is.null(cluster_ids) && !is.null(strata_ids)) {
      number_of_stages <- ncol(strata_ids)
      number_of_ultimate_units <- nrow(strata_ids)
      cluster_ids <- Matrix::Matrix(data = rep(seq_len(number_of_ultimate_units),
                                       times = number_of_stages),
                            nrow = number_of_ultimate_units, ncol = number_of_stages,
                            byrow = FALSE)
    } else if (!is.null(cluster_ids) && is.null(strata_ids)) {
      number_of_stages <- ncol(cluster_ids)
      number_of_ultimate_units <- nrow(cluster_ids)
      strata_ids <- Matrix::Matrix(1,
                           nrow = number_of_ultimate_units,
                           ncol = number_of_stages)
    } else {
      number_of_stages <- ncol(strata_ids)
      number_of_ultimate_units <- nrow(strata_ids)
    }

    if ((variance_estimator %in% c("Deville-Tille", "BOSB")) && (number_of_stages > 1)) {
      sprintf(
        "For `variance_estimator = '%s', the quadratic form only takes into account the first stage of sampling.",
        variance_estimator
      ) |> message()
    }

    # Make sure each stage's sampling units are nested within strata
    # and each stage's sampling units are nested
    # within previous stage sampling units
    cluster_ids[,1] <- interaction(strata_ids[, 1, drop = TRUE],
                                   cluster_ids[, 1, drop = TRUE],
                                   sep = " | ", drop = TRUE)
    if (variance_estimator != "Ultimate Cluster") {
      stage <- 2L
      while (stage <= number_of_stages) {
        strata_ids[,stage] <- interaction(
          cluster_ids[, stage-1L, drop=TRUE],
          strata_ids[, stage, drop=TRUE],
          sep = " | ", drop = TRUE
        )
        cluster_ids[,stage] <- interaction(
          strata_ids[, stage, drop = TRUE],
          cluster_ids[, stage, drop = TRUE],
          sep = " | ", drop = TRUE
        )
        stage <- stage + 1L
      }
    }

    if (is.null(strata_pop_sizes)) {
      strata_pop_sizes <- Matrix::Matrix(data = Inf,
                                         nrow = number_of_ultimate_units,
                                         ncol = number_of_stages)
    }
  }

  if (variance_estimator == "Horvitz-Thompson") {
    n <- number_of_ultimate_units
    quad_form_matrix <- Matrix::Matrix(nrow = n, ncol = n) |> as("symmetricMatrix")
    for (i in seq_len(n)) {
      for (j in seq(from = i, to = n, by = 1)) {
        quad_form_matrix[i,j] <- 1 - (joint_probs[i,i] * joint_probs[j,j])/joint_probs[i,j]
      }
    }
    quad_form_matrix[lower.tri(quad_form_matrix)] <- t(quad_form_matrix)[lower.tri(quad_form_matrix)]
  }

  if (variance_estimator == "Yates-Grundy") {
    n <- number_of_ultimate_units
    quad_form_matrix <- -(1 - as(Matrix::tcrossprod(diag(joint_probs)), "symmetricMatrix")/as(joint_probs, "symmetricMatrix"))
    diag(quad_form_matrix) <- Matrix::diag(quad_form_matrix) - Matrix::rowSums(quad_form_matrix)
    quad_form_matrix <- -quad_form_matrix
  }

  if (variance_estimator == "Ultimate Cluster") {
    quad_form_matrix <- Matrix::Matrix(data = 0,
                                       nrow = number_of_ultimate_units,
                                       ncol = number_of_ultimate_units) |>
      as("symmetricMatrix")
    if (use_sparse_matrix) {
      quad_form_matrix <- quad_form_matrix |> as("CsparseMatrix")
    }

    # Generate quadratic form for each stratum
    for (stratum_id in unique(strata_ids[,1,drop=TRUE])) {
      stratum_indices <- which(strata_ids[,1,drop=TRUE] == stratum_id)
      stratum_pop_size <- strata_pop_sizes[,1,drop=TRUE][stratum_indices[1]]
      n_clusters <- length(unique(cluster_ids[stratum_indices, 1, drop = TRUE]))
      quad_form_matrix[stratum_indices,stratum_indices] <- distribute_matrix_across_clusters(
        cluster_level_matrix = make_srswor_matrix(n = n_clusters,
                                                  f = (n_clusters/stratum_pop_size)),
        cluster_ids = cluster_ids[stratum_indices, 1, drop = TRUE],
        rows = TRUE, cols = TRUE
      )
    }
  }

  if (variance_estimator %in% c("Stratified Multistage SRS", "Deville-1", "Deville-2", "Deville-Tille", "Beaumont-Emond", "BOSB")) {
    quad_form_matrix <- Matrix::Matrix(data = 0,
                                       nrow = number_of_ultimate_units,
                                       ncol = number_of_ultimate_units) |>
      as("symmetricMatrix")

    if (use_sparse_matrix) {
      quad_form_matrix <- quad_form_matrix |> as("CsparseMatrix")
    }

    # Obtain matrix of cluster sample sizes by stage
    count <- function(x) sum(!duplicated(x))
    sampsize <- matrix(ncol = ncol(cluster_ids), nrow = nrow(cluster_ids))
    for (i in seq_len(number_of_stages)) {
      split(sampsize[, i], strata_ids[, i]) <- lapply(split(cluster_ids[,i], strata_ids[, i]), count)
    }
    # Iterate over stages
    stage <- 1L
    while (stage <= number_of_stages) {
      # Generate quadratic form for each stratum
      for (stratum_id in unique(strata_ids[,stage,drop=TRUE])) {

        stratum_indices <- which(strata_ids[,stage,drop=TRUE] == stratum_id)

        if (variance_estimator == "Stratified Multistage SRS") {

          # Get quadratic form at current stage
          stratum_pop_size <- strata_pop_sizes[,stage,drop=TRUE][stratum_indices[1]]
          n_clusters <- sampsize[,stage,drop=TRUE][stratum_indices[1]]

          Q_current <- distribute_matrix_across_clusters(
            cluster_level_matrix = make_srswor_matrix(n = n_clusters,
                                                      f = (n_clusters/stratum_pop_size)),
            cluster_ids = cluster_ids[stratum_indices, stage, drop = TRUE],
            rows = TRUE, cols = TRUE
          )
          # Get product of sampling probabilities from previous stages
          if (stage > 1) {
            prev_n_clusters <- sampsize[stratum_indices[1],seq_len(stage-1),drop=TRUE]
            prev_stratum_pop_size <- strata_pop_sizes[stratum_indices[1],seq_len(stage-1),drop=TRUE]
            prev_samp_fraction <- (
              Reduce(f = `*`, x = prev_n_clusters) / Reduce(f = `*`, x = prev_stratum_pop_size)
            )
          } else {
            prev_samp_fraction <- 1
          }
          prev_stages_samp_prob <- prev_samp_fraction
        }

        if (variance_estimator %in% c("Deville-1", "Deville-2", "Beaumont-Emond")) {

          # Get quadratic form at current stage
          current_cluster_ids <- cluster_ids[stratum_indices, stage, drop = TRUE]
          current_probs <- probs[stratum_indices, stage, drop = TRUE]
          cluster_probs <- current_probs[!duplicated(current_cluster_ids)]

          Q_current <- distribute_matrix_across_clusters(
            cluster_level_matrix = make_ppswor_approx_matrix(
              probs = cluster_probs,
              method = variance_estimator
            ),
            cluster_ids = current_cluster_ids,
            rows = TRUE, cols = TRUE
          )
          # Get product of sampling probabilities from previous stages
          if (stage > 1) {
            probs_at_prev_stages <- probs[stratum_indices[1],seq_len(stage-1),drop=TRUE]
            prev_stages_samp_prob <- Reduce(f = `*`, x = probs_at_prev_stages)
          } else {
            prev_stages_samp_prob <- 1
          }

        }

        if ((variance_estimator %in% c("Deville-Tille")) && (stage == 1)) {
          # Get quadratic form at current stage
          current_cluster_ids <- cluster_ids[stratum_indices, stage, drop = TRUE]
          current_probs <- probs[stratum_indices, stage, drop = TRUE]
          cluster_probs <- current_probs[!duplicated(current_cluster_ids)]
          current_aux_vars <- aux_vars[stratum_indices, , drop = FALSE]
          cluster_aux_vars <- current_aux_vars[!duplicated(current_cluster_ids), , drop = FALSE]

          Q_current <- distribute_matrix_across_clusters(
            cluster_level_matrix = make_deville_tille_matrix(
              probs = cluster_probs,
              aux_vars = cluster_aux_vars
            ),
            cluster_ids = current_cluster_ids,
            rows = TRUE, cols = TRUE
          )

          prev_stages_samp_prob <- 1
        }
        
        if ((variance_estimator %in% c("BOSB")) && (stage == 1)) {
          # Get quadratic form at current stage
          current_cluster_ids <- cluster_ids[stratum_indices, stage, drop = TRUE]
          current_aux_vars <- aux_vars[stratum_indices, , drop = FALSE]
          cluster_aux_vars <- current_aux_vars[!duplicated(current_cluster_ids), , drop = FALSE]
          
          Q_current <- distribute_matrix_across_clusters(
            cluster_level_matrix = make_kernel_var_matrix(
              x         = cluster_aux_vars[, 1, drop = TRUE],
              kernel    = "Epanechnikov",
              bandwidth = "auto"
            ),
            cluster_ids = current_cluster_ids,
            rows = TRUE, cols = TRUE
          )
          
          prev_stages_samp_prob <- 1
        }

        # Add overall variance contribution from current stage/stratum sampling
        quad_form_matrix[stratum_indices,stratum_indices] <- (
          quad_form_matrix[stratum_indices,stratum_indices] +
            (Q_current * prev_stages_samp_prob)
        )
        
        if (any(is.nan(quad_form_matrix))) {}

      }
      stage <- stage + 1L
    }

  }

  if (variance_estimator %in% c("SD1", "SD2")) {
    n <- number_of_ultimate_units
    # Initialize quadratic form matrix
    quad_form_matrix <- Matrix::Matrix(data = 0, nrow = n, ncol = n) |>
      as("symmetricMatrix")
    if (use_sparse_matrix) {
      quad_form_matrix <- quad_form_matrix |> as("CsparseMatrix")
    }
    sorted_quad_form_matrix <- quad_form_matrix
    # Sort the inputs, compile into a dataframe
    sorted_df <- data.frame('Row_ID' = seq_len(n),
                            'Stratum' = strata_ids[,1,drop=TRUE],
                            'Cluster' = cluster_ids[,1,drop=TRUE],
                            'Sort_Order' = sort_order,
                            'Stratum_Pop_Size' = strata_pop_sizes[,1,drop=TRUE],
                            stringsAsFactors = FALSE)
    sorted_df <- sorted_df[order(sorted_df[['Sort_Order']],
                                 sorted_df[['Stratum']],
                                 sorted_df[['Cluster']]),]
    sorted_df[['New_Order']] <- seq_len(n)
    inverse_sort_map <- sorted_df[['New_Order']][order(sorted_df[['Row_ID']])]

    # Generate quadratic form for each stratum
    for (stratum_id in unique(sorted_df[['Stratum']])) {
      stratum_indices <- which(sorted_df[['Stratum']] == stratum_id)
      stratum_pop_size <- sorted_df[['Stratum_Pop_Size']][stratum_indices[1]]

      n_clusters <- length(
        unique(sorted_df[stratum_indices, 'Cluster', drop = TRUE])
      )

      sorted_quad_form_matrix[stratum_indices,stratum_indices] <- distribute_matrix_across_clusters(
        cluster_level_matrix = make_sd_matrix(n = n_clusters,
                                              f = (n_clusters/stratum_pop_size),
                                              type = variance_estimator),
        cluster_ids = sorted_df[stratum_indices, 'Cluster', drop = TRUE],
        rows = TRUE, cols = TRUE
      )
    }
    # Arrange matrix rows/columns to match the original order of the input data
      quad_form_matrix <- sorted_quad_form_matrix[inverse_sort_map,inverse_sort_map]
  }

  return(quad_form_matrix)
}

#' @title Create a quadratic form's matrix to represent a successive-difference variance estimator
#' @description A successive-difference variance estimator can be represented
#' as a quadratic form. This function determines the matrix of the quadratic form.
#' @param n Number of rows or columns for the matrix
#' @param f A single number between \code{0} and \code{1},
#' representing the sampling fraction. Default value is \code{0}.
#' @param type Either "SD1" or "SD2". See the "Details" section for definitions.
#' @return A matrix of dimension \code{n}
#' @details
#' Ash (2014) describes each estimator as follows:
#' \deqn{
#'   \hat{v}_{SD1}(\hat{Y}) = (1-f) \frac{n}{2(n-1)} \sum_{k=2}^n\left(\breve{y}_k-\breve{y}_{k-1}\right)^2
#' }
#' \deqn{
#'   \hat{v}_{SD2}(\hat{Y}) = \frac{1}{2}(1-f)\left[\sum_{k=2}^n\left(\breve{y}_k-\breve{y}_{k-1}\right)^2+\left(\breve{y}_n-\breve{y}_1\right)^2\right]
#' }
#' where \eqn{\breve{y}_k} is the weighted value \eqn{y_k/\pi_k} of unit \eqn{k}
#' with selection probability \eqn{\pi_k}, and \eqn{f} is the sampling fraction \eqn{\frac{n}{N}}.
#' @references
#' Ash, S. (2014). "\emph{Using successive difference replication for estimating variances}."
#' \strong{Survey Methodology}, Statistics Canada, 40(1), 47-59.
#' @keywords internal
make_sd_matrix <- function(n, f = 0, type = "SD1") {

  if (!is.numeric(n) || (length(n) != 1) || is.na(n) || (n %% 1 != 0) || (n < 1)) {
    stop("`n` must be an integer greater than or equal to 1")
  }
  if (!is.numeric(f) || (length(f) != 1) || is.na(f) || (f < 0) || (f > 1)) {
    stop("`f` must be a single number between 0 and 1")
  }

  if (n == 1) {
    C_matrix <- Matrix::Matrix(0, nrow = 1, ncol = 1)
  }
  if (type == "SD1" && (n > 1)) {
    C_matrix <- as(Matrix::Diagonal(n) * 2, "symmetricMatrix")
    C_matrix[1,1] <- 1
    C_matrix[n,n] <- 1
    for (i in seq_len(n)) {
      if (i < n) {
        C_matrix[i,i+1] <- -1
        C_matrix[i+1,i] <- -1
      }
      if (i > 1) {
        C_matrix[i,i-1] <- -1
        C_matrix[i-1,i] <- -1
      }
    }
    C_matrix <- (n/(n-1)) * C_matrix
  }
  if (type == "SD2" && (n > 1)) {
    C_matrix <- as(Matrix::Diagonal(n) * 2, "symmetricMatrix")
    for (i in seq_len(n)) {
      if (i < n) {
        C_matrix[i,i+1] <- -1
        C_matrix[i+1,i] <- -1
      }
      if (i > 1) {
        C_matrix[i,i-1] <- -1
        C_matrix[i-1,i] <- -1
      }
    }
    if (n > 2) {
      C_matrix[1,n] <- -1
      C_matrix[n,1] <- -1
    } else if (n == 2) {
      C_matrix[1,n] <- -2
      C_matrix[n,1] <- -2
    }
  }

  C_matrix <- ((1-f)/2) * C_matrix

  return(C_matrix)
}


#' @title Create a quadratic form's matrix to represent the basic variance estimator
#' for a total under simple random sampling without replacement
#' @description The usual variance estimator for simple random sampling without replacement
#' can be represented as a quadratic form.
#' This function determines the matrix of the quadratic form.
#' @param n Sample size
#' @param f A single number between \code{0} and \code{1},
#' representing the sampling fraction. Default value is \code{0}.
#' @return A symmetric matrix of dimension \code{n}
#' @details
#' The basic variance estimator of a total for simple random sampling without replacement is as follows:
#' \deqn{
#' \hat{v}(\hat{Y}) = (1 - f)\frac{n}{n - 1} \sum_{i=1}^{n} (y_i - \bar{y})^2
#' }
#' where \eqn{f} is the sampling fraction \eqn{\frac{n}{N}}. \cr \cr
#' If \eqn{f=0}, then the matrix of the quadratic form has all non-diagonal elements equal to \eqn{-(n-1)^{-1}},
#' and all diagonal elements equal to \eqn{1}. If \eqn{f > 0}, then each element
#' is multiplied by \eqn{(1-f)}. \cr \cr
#' If \eqn{n=1}, then this function returns a \eqn{1 \times 1} matrix whose sole element equals \eqn{0}
#' (essentially treating the sole sampled unit as a selection made with probability \eqn{1}).
#' @keywords internal
make_srswor_matrix <- function(n, f = 0) {
  if (!is.numeric(n) || (length(n) != 1) || is.na(n) || (n %% 1 != 0) || (n < 1)) {
    stop("`n` must be an integer greater than or equal to 1")
  }
  if (!is.numeric(f) || (length(f) != 1) || is.na(f) || (f < 0) || (f > 1)) {
    stop("`f` must be a single number between 0 and 1")
  }
  if (n == 1) {
    C_matrix <- Matrix::Matrix(0, nrow = 1, ncol = 1)
  } else {
    C_matrix <- Matrix::Matrix(-1/(n-1), nrow = n, ncol = n)
    diag(C_matrix) <- 1
    C_matrix <- (1-f) * C_matrix
  }
  return(C_matrix)
}

#' @title Create a quadratic form's matrix to represent a variance estimator
#' for PPSWOR designs, based on commonly-used approximations
#' @description Several variance estimators for designs that use
#' unequal probability sampling without replacement (i.e., PPSWOR),
#' variance estimation tends to be more accurate
#' when using an approximation estimator that uses the first-order
#' inclusion probabilities (i.e., the basic sampling weights)
#' and ignores the joint inclusion probabilities.
#' This function returns the matrix of the quadratic form used
#' to represent such variance estimators.
#' @param probs A vector of first-order inclusion probabilities
#' @param method A string specifying the approximation method to use.
#' See the "Details" section below.
#' Options include:
#' \itemize{
#'   \item "Deville-1"
#'   \item "Deville-2"
#'   \item "Beaumont-Emond"
#' }

#' @return A symmetric matrix whose dimension matches the length of \code{probs}.
#' @section Deville's Estimators: 
#' The "Deville-1" and "Deville-2" approximations have been shown to be effective
#' for designs that use a fixed sample size with a high-entropy sampling method.
#' This includes most PPSWOR sampling methods,
#' but unequal-probability systematic sampling is an important exception.
#'
#' Deville's variance estimators generally take the following form:
#' \deqn{
#' \hat{v}(\hat{Y}) = \sum_{i=1}^{n} c_i (\breve{y}_i - \frac{1}{\sum_{i=k}^{n}c_k}\sum_{k=1}^{n}c_k \breve{y}_k)^2
#' }
#' where \eqn{\breve{y}_i = y_i/\pi_i} is the weighted value of the the variable of interest,
#' and \eqn{c_i} are constants that depend on the approximation method used.  \cr \cr
#' The matrix of the quadratic form, denoted \eqn{\Sigma}, has
#' its \eqn{ij}-th entry defined as follows:
#' \deqn{
#'   \sigma_{ii} = c_i (1 - \frac{c_i}{\sum_{k=1}^{n}c_k}) \textit{ when } i = j \\
#'   \sigma_{ij}=\frac{-c_i c_j}{\sum_{k=1}^{n}c_k} \textit{ when } i \neq j \\
#' }
#' When \eqn{\pi_{i} = 1} for every unit, then \eqn{\sigma_{ij}=0} for all \eqn{i,j}.
#' If there is only one sampling unit, then \eqn{\sigma_{11}=0}; that is, the unit is treated as if it was sampled with certainty.
#'
#' The constants \eqn{c_i} are defined for each approximation method as follows,
#' with the names taken directly from Matei and Tillé (2005).
#' \itemize{
#'   \item \strong{"Deville-1"}:
#'     \deqn{c_i=\left(1-\pi_i\right) \frac{n}{n-1}}
#'   \item \strong{"Deville-2"}:
#'     \deqn{c_i = (1-\pi_i) \left[1 - \sum_{k=1}^{n} \left(\frac{1-\pi_k}{\sum_{k=1}^{n}(1-\pi_k)}\right)^2 \right]^{-1}}
#' }
#' Both of the approximations \strong{"Deville-1"} and \strong{"Deville-2"} were shown
#' in the simulation studies of Matei and Tillé (2005) to perform much better
#' in terms of MSE compared to the strictly-unbiased
#' Horvitz-Thompson and Yates-Grundy variance estimators.
#' In the case of simple random sampling without replacement (SRSWOR),
#' these estimators are identical to the usual Horvitz-Thompson variance estimator.
#' 
#' @section Beaumont-Emond Estimator: 
#' Beaumont and Emond (2022) proposed a variance estimator for unequal probability
#' sampling without replacement. This estimator is simply the Horvitz-Thompson
#' variance estimator with the following approximation for the joint inclusion
#' probabilities.
#' \deqn{
#'   \pi_{kl} \approx \pi_k \pi_l \frac{n - 1}{(n-1) + \sqrt{(1-\pi_k)(1-\pi_l)}}
#' }
#' In the case of cluster sampling, this approximation should be 
#' applied to the clusters rather than the units within clusters.
#'
#' @references
#' Matei, Alina, and Yves Tillé. 2005.
#' "Evaluation of Variance Approximations and Estimators
#' in Maximum Entropy Sampling with Unequal Probability and Fixed Sample Size."
#' Journal of Official Statistics 21(4):543-70.
#'
#' @keywords internal
make_ppswor_approx_matrix <- function(probs, method = "Deville-1") {

  n <- length(probs)
  one_minus_pi <- 1 - probs

  if (method %in% c("Deville-1", "Deville-2")) {
    if (method == "Deville-1") {
      c_k <- (1 - probs) * (n/(n-1))
    }
    if (method == "Deville-2") {
      c_k <- (1 - probs) / (
        1 - sum( (one_minus_pi/sum(one_minus_pi))^2 )
      )
    }
    
    c_sum <- sum(c_k)
    
    if ((n == 1) || (c_sum == 0)) {
      Sigma <- Matrix::Matrix(0, nrow = length(c_k), ncol = length(c_k))
    } else {
      Sigma <- outer(c_k, -c_k) / c_sum
      diag(Sigma) <- c_k*(1 - c_k/c_sum)
      Sigma <- as(Sigma, "symmetricMatrix")
    }
  }
  if (method == "Beaumont-Emond") {
    if ((n == 1)) {
      Sigma <- Matrix::Matrix(0, n, n)
    } else {
      Sigma <- - tcrossprod(sqrt(one_minus_pi)) / (n-1)
      diag(Sigma) <- one_minus_pi
      Sigma <- as(Sigma, "symmetricMatrix")
    }
  }

  return(Sigma)
}

#' @title Create a quadratic form's matrix
#' for a Deville-Tillé variance estimator for balanced samples
#' @description Creates the quadratic form matrix for a
#' variance estimator for balanced samples,
#' proposed by Deville and Tillé (2005).
#' @param probs A vector of first-order inclusion probabilities
#' @param aux_vars A matrix of auxiliary variables,
#' with the number of rows matching the number of elements of \code{probs}.
#'
#' @return A symmetric matrix whose dimension matches the length of \code{probs}.
#' @keywords internal
#' @details
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
#'
#' See Li, Chen, and Krenzke (2014) for an example of this estimator's use
#' as the basis for a generalized replication estimator. See Breidt and Chauvet (2011)
#' for a discussion of alternative simulation-based estimators for the specific
#' application of variance estimation for balanced samples selected using the cube method.
#'
#' @references
#' - Breidt, F.J. and Chauvet, G. (2011).
#' "Improved variance estimation for balanced samples drawn via the cube method."
#' Journal of Statistical Planning and Inference, 141, 411-425.
#'
#' - Deville, J.C., and Tillé, Y. (2005). "\emph{Variance approximation under balanced sampling.}"
#' \strong{Journal of Statistical Planning and Inference}, 128, 569-591.
#'
#' - Li, J., Chen, S., and Krenzke, T. (2014).
#' "Replication Variance Estimation for Balanced Sampling: An Application to the PIAAC Study."
#' Proceedings of the Survey Research Methods Section, 2014: 985-994. Alexandria, VA: American Statistical Association.
#' http://www.asasrms.org/Proceedings/papers/1984_094.pdf.
#'
#' - Tillé, Y. (2020). "\emph{Sampling and estimation from finite populations}." (I. Hekimi, Trans.). Wiley.
make_deville_tille_matrix <- function(probs, aux_vars) {

  if (all(probs == 1)) {
    return(matrix(0, length(probs), length(probs)))
  }

  n <- length(probs)

  q <- ncol(aux_vars)

  if (is.null(q) || (q < 1)) {
    stop("`aux_vars` must be a matrix with at least one column")
  }
  if (q >= n) {
    error_msg <- paste(
      "The number of columns of `aux_vars` exceeds the number of observations:",
      "the estimator is undefined in this case."
    )
    stop(error_msg)
  }

  # Calculate the constants used for each observation
  c_k <- (1 - probs) * (n/(n-q))

  # Calculate the hat matrix for creating y-star
  H <- wls_hat_matrix(X = aux_vars, w = c_k/(probs^2))

  # Subtract the hat matrix from the identity matrix
  I_minus_H <- diag(n) - H

  # Weight each entry
  wtd_I_minus_H <- I_minus_H * sqrt(c_k/(probs^2))

  Sigma <- diag(probs) %*% crossprod(wtd_I_minus_H, wtd_I_minus_H) %*% diag(probs)

  return(Sigma)
}

#' @title Make a quadratic form matrix for the kernel-based variance estimator
#' of Breidt, Opsomer, and Sanchez-Borrego (2016)
#' @description Constructs the quadratic form matrix
#' for the kernel-based variance estimator of Breidt, Opsomer, and Sanchez-Borrego (2016).
#' The bandwidth is automatically chosen to result
#' in the smallest possible nonempty kernel window.
#' @param x A numeric vector, giving the values
#' of an auxiliary variable.
#' @param kernel The name of a kernel function. 
#' Currently only "Epanechnikov" is supported.
#' @param bandwidth The bandwidth to use for the kernel.
#' The default value is `"auto"`, which means that the bandwidth
#' will be chosen automatically to produce the smallest window size
#' while ensuring that every unit has a nonempty window,
#' as suggested by Breidt, Opsomer, and Sanchez-Borrego (2016).
#' Otherwise, the user can supply their own value, which can be a 
#' single positive number.
#'
#' @return The quadratic form matrix for the variance estimator,
#' with dimension equal to the length of `x`. The resulting
#' object has an attribute `bandwidth` that can be retrieved
#' using `attr(Q, 'bandwidth')`
#' 
#' @details
#' 
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
#' If \eqn{n=2}, then the estimator is simply the estimator
#' used for simple random sampling without replacement.
#' 
#' If \eqn{n=1}, then the matrix simply has an entry equal to 0.
#' 
#' @export
#' @md
#' 
#' @references
#' Breidt, F. J., Opsomer, J. D., & Sanchez-Borrego, I. (2016). 
#' "\emph{Nonparametric Variance Estimation Under Fine Stratification: An Alternative to Collapsed Strata}." 
#' \strong{Journal of the American Statistical Association}, 111(514), 822-833. https://doi.org/10.1080/01621459.2015.1058264
#' @examples
#' # The auxiliary variable has the same value for all units
#' make_kernel_var_matrix(c(1, 1, 1))
#' 
#' # The auxiliary variable differs across units
#' make_kernel_var_matrix(c(1, 2, 3))
#' 
#' # View the bandwidth that was automatically selected
#' Q <- make_kernel_var_matrix(c(1, 2, 4))
#' attr(Q, 'bandwidth')
#' 
make_kernel_var_matrix <- function(x, kernel = "Epanechnikov", bandwidth = "auto") {
  
  if (!is.vector(x) || !is.numeric(x)) {
    stop("Only numeric vectors are supported for the argument `x`.")
  }
  
  kernel_fn <- switch(
    kernel,
    "Epanechnikov" = function(u) 0.75 * (1-pmin(abs(u), 1)^2),
    stop("The only value allowed for `kernel` is 'Epanechnikov'.")
  )
  
  # Rescale 'x' to the range [-1,1]
  max_abs_x <- max(abs(x))
  if (max_abs_x > 0) {
    x <- x/max_abs_x
  }
  
  # Construct matrix of pairwise differences
  
  H     <- length(x)
  
  diffs <- outer(x, x, FUN = `-`)

  diffs <- round(diffs, 8)
  
  # Determine bandwidth
  
  if (!is.numeric(bandwidth) && !is.character(bandwidth)) {
    stop("`bandwidth` must be either 'auto' or a single positive number.") 
  }
  if (length(bandwidth) != 1) {
    stop("`bandwidth` must be either 'auto' or a single positive number.")
  }
  if (is.na(bandwidth)) {
    stop("`bandwidth` must be either 'auto' or a single positive number.")
  }
  
  if (is.numeric(bandwidth)) {
    if (bandwidth <= 0) {
      stop("`bandwidth` must be either 'auto' or a single positive number.")
    }
    bw <- bandwidth
  }
  if (is.character(bandwidth)) {
    if (bandwidth != "auto") {
      stop("`bandwidth` must be either 'auto' or a single positive number.")
    }
  }
  
  # Choose bandwidth to use the smallest window
  # such that each unit has a nonempty window
  if (bandwidth == "auto") {
    if (H > 2) {
      bw <- 1
      min_bw <- min(abs(upper.tri(diffs)))
      for (i in seq_len(H)) {
        abs_diffs     <- abs(diffs[i,-i])
        min_abs_diff  <- min(abs_diffs)
        max_abs_diff  <- max(abs_diffs)
        
        if (bw > min_abs_diff) {
          if (max_abs_diff > min_abs_diff) {
            bw <- pmax(pmin(max_abs_diff, bw), min_bw)
          }
        }
        if (bw <= min_abs_diff) {
          candidate_bw_values <- abs_diffs[abs_diffs > bw]
          if (length(candidate_bw_values) > 0) {
            bw <- min(candidate_bw_values)
          }
        }
      }
    }
    if (H == 2) {
      bw <- Inf
    }
    if (H == 1) {
      bw <- 1
    }
  }
  K     <- kernel_fn(diffs/bw)
  
  # Construct quadratic form for the variance estimator
  
  D_ij  <- K / rowSums(K)
  
  if (H > 1) {
    C_d   <- 1 - (1/H) * (2*sum(diag(D_ij)) - sum(D_ij^2))
  }
  if (H == 1) {
    C_d   <- 1
  }
  
  Q     <- (1/C_d) * crossprod(diag(H) - D_ij)
  
  Q <- as(Q, "symmetricMatrix")
  
  # Save the bandwidth as an attribute
  attr(Q, "bandwidth") <- bw * max_abs_x
  
  return(Q)
}


#' @title Helper function to turn a cluster-level matrix into an element-level matrix
#' by duplicating rows or columns of the matrix
#' @description Turns a cluster-level matrix into an element-level matrix
#' by suitably duplicating rows or columns of the matrix.
#' @param cluster_level_matrix A square matrix, whose number of rows/columns matches the number of clusters.
#' @param cluster_ids A vector of cluster identifiers.
#' If \code{rows=TRUE}, the number of unique elements of \code{cluster_ids}
#' must match the number of rows of \code{cluster_level_matrix}.
#' If \code{cols=TRUE}, the number of unique elements of \code{cluster_ids}
#' must match the number of columns of \code{cluster_level_matrix}.
#' @param rows Whether to duplicate rows of the \code{cluster_level_matrix} for elements from the same cluster.
#' @param cols Whether to duplicate columns of the \code{cluster_level_matrix} for elements from the same cluster.
#' @return The input \code{cluster_level_matrix} has its rows/columns
#' duplicated so that the number of rows (if \code{rows=TRUE}) or columns (if \code{cols=TRUE})
#' equals the length of \code{cluster_ids}.
#' @keywords internal
distribute_matrix_across_clusters <- function(cluster_level_matrix, cluster_ids, rows = TRUE, cols = TRUE) {
  if (any(is.na(cluster_ids))) {
    stop("`cluster_ids` cannot have any missing values.")
  }
  ## Get indices for each cluster
  cluster_indices <- as.numeric(
    factor(cluster_ids,
           levels = unique(cluster_ids))
  )
  n_clusters <- max(cluster_indices)
  if ((nrow(cluster_level_matrix) != n_clusters) | (ncol(cluster_level_matrix) != n_clusters)) {

  }
  if (rows && cols) {
    if ((n_clusters != nrow(cluster_level_matrix)) || (n_clusters != ncol(cluster_level_matrix))) {
      stop("The number of rows/columns of `cluster_level_matrix` must match the number of unique elements of `cluster_ids`.")
    }
    result <- cluster_level_matrix[cluster_indices,cluster_indices,drop=FALSE]
  } else if (rows && (!cols)) {
    if ((n_clusters != nrow(cluster_level_matrix))) {
      stop("The number of rows of `cluster_level_matrix` must match the number of unique elements of `cluster_ids`.")
    }
    result <- cluster_level_matrix[cluster_indices,,drop=FALSE]
  } else if ((!rows) && cols) {
    if ((n_clusters != ncol(cluster_level_matrix))) {
      stop("The number of columns of `cluster_level_matrix` must match the number of unique elements of `cluster_ids`.")
    }
    result <- cluster_level_matrix[,cluster_indices,drop=FALSE]
  } else {
    stop("Must set `rows=TRUE` and/or `cols=TRUE`")
  }
  return(result)
}

#' @title Check whether a matrix is positive semidefinite
#' @description Check whether a matrix is positive semidefinite, based on checking for symmetric and negative eigenvalues.
#'
#' @param X A matrix with no missing or infinite values.
#' @param tolerance Tolerance for controlling whether
#' a tiny computed eigenvalue will actually be considered negative.
#' Computed negative eigenvalues will be considered
#' negative if they are less than which are less than
#' \code{-abs(tolerance * max(eigen(X)$values))}.
#' A small nonzero tolerance is recommended
#' since eigenvalues are nearly always computed with some floating-point error.
#'
#' @return A logical value. \code{TRUE} if the matrix is deemed positive semidefinite.
#' Negative otherwise (including if \code{X} is not symmetric).
#'
#' @seealso The function \code{\link[svrep]{get_nearest_psd_matrix}()}
#' can be used to approximate a symmetric matrix which is not positive semidefinite,
#' by a similar positive semidefinite matrix.
#'
#' @examples
#' X <- matrix(
#'   c(2, 5, 5,
#'     5, 2, 5,
#'     5, 5, 2),
#'   nrow = 3, byrow = TRUE
#' )
#'
#' is_psd_matrix(X)
#'
#' eigen(X)$values
#' @export
is_psd_matrix <- function(X, tolerance = sqrt(.Machine$double.eps)) {
  if (inherits(X, 'symmetricMatrix')) {
    symmetric <- TRUE
  } else {
    symmetric <- Matrix::isSymmetric(X)
  }

  if (!symmetric) {
    result <- FALSE
  } else {
    if ((sum(is.na(X)) + sum(is.infinite(X))) > 0) {
      stop("The matrix `X` should not have any missing or infinite values.")
    }
    eigenvalues <- eigen(X, only.values = TRUE)$values
    result <- all(eigenvalues >= -tolerance * abs(eigenvalues[1]))
  }
  return(result)
}

#' @title Approximates a symmetric, real matrix by the nearest positive
#' semidefinite matrix.
#'
#' @description Approximates a symmetric, real matrix by the nearest positive
#' semidefinite matrix in the Frobenius norm, using the method of Higham (1988).
#' For a real, symmetric matrix, this is equivalent to "zeroing out" negative eigenvalues.
#' See the "Details" section for more information.
#'
#' @param X A symmetric, real matrix with no missing values.
#'
#' @details
#' Let \eqn{A} denote a symmetric, real matrix which is not positive semidefinite.
#' Then we can form the spectral decomposition \eqn{A=\Gamma \Lambda \Gamma^{\prime}},
#' where \eqn{\Lambda} is the diagonal matrix
#' whose entries are eigenvalues of \eqn{A}.
#' The method of Higham (1988) is to  approximate
#' \eqn{A} with \eqn{\tilde{A} = \Gamma \Lambda_{+} \Gamma^{\prime}},
#' where the \eqn{ii}-th entry of \eqn{\Lambda_{+}} is \eqn{\max(\Lambda_{ii}, 0)}.
#'
#' @return The nearest positive semidefinite matrix
#' of the same dimension as \code{X}.
#'
#' @references
#' - Higham, N. J. (1988). "\emph{Computing a nearest symmetric positive semidefinite matrix.}" Linear Algebra and Its Applications, 103, 103-118.
#' @export
#'
#' @examples
#' X <- matrix(
#'   c(2, 5, 5,
#'     5, 2, 5,
#'     5, 5, 2),
#'   nrow = 3, byrow = TRUE
#' )
#' get_nearest_psd_matrix(X)
get_nearest_psd_matrix <- function(X) {

  eigen_decomposition <- eigen(X)
  eigen_vectors <- eigen_decomposition$vectors
  eigen_values <- eigen_decomposition$values

  updated_eigen_values <- pmax(eigen_values, 0)
  updated_eigen_values <- abs(updated_eigen_values)

  X <- eigen_vectors %*% diag(updated_eigen_values) %*% t(eigen_vectors)
  return(X)
}

#' @title Compute the matrix of joint inclusion probabilities
#' from the quadratic form of a Horvitz-Thompson variance estimator.
#'
#' @param ht_quad_form The matrix of the quadratic form
#' representing the Horvitz-Thompson variance estimator.
#' @details The quadratic form matrix of the Horvitz-Thompson variance estimator
#' has \eqn{ij}-th entry equal to \eqn{(1-\frac{\pi_i \pi_j}{\pi_{ij}})}.
#' The matrix of joint probabilties has \eqn{ij}-th entry equal to \eqn{\pi_{ij}}.
#' @return The matrix of joint inclusion probabilities
#' @keywords internal
ht_matrix_to_joint_probs <- function(ht_quad_form) {
  first_order_probs <- 1 - Matrix::diag(ht_quad_form)
  joint_probs <- ((1 - ht_quad_form)^(-1)) * outer(first_order_probs, first_order_probs)
  return(joint_probs)
}

#' @title Create the "hat matrix" for weighted least squares regression
#'
#' @param X Matrix of predictor variables, with \code{n} rows
#' @param w Vector of weights (should all be nonnegative), of length \code{n}
#'
#' @return An \eqn{n \times n} matrix. This is the "hat matrix" for a WLS regression.
#' @keywords internal
#'
wls_hat_matrix <- function(X, w) {

  # Temporarily drop entries with zero weights
  n <- length(w)
  which_zero <- (w == 0)
  X <- X[(w != 0), , drop = FALSE]
  w <- w[(w != 0)]

  # Compute the hat matrix
  rw <- sqrt(w)
  X1 <- rw * X

  QR <- qr.default(X1, LAPACK = TRUE)
  Q <- qr.qy(QR, diag(1, nrow = nrow(QR$qr), ncol = QR$rank))

  Q1 <- (1 / rw) * Q
  Q2 <- rw * Q

  H <- tcrossprod(Q1, Q2)

  # Restore entries with zero weights
  result <- matrix(0, nrow = n, ncol = n)
  result[!which_zero, !which_zero] <- H

  return(result)
}

#' @title Combine quadratic forms from each phase of a two phase design
#' @description This function combines quadratic forms from each phase of a two phase design,
#' so that the combined variance of the entire two-phase sampling design can be estimated.
#' @param sigma_1 The quadratic form for the first phase variance estimator,
#' subsetted to only include cases selected in the phase two sample.
#' @param sigma_2 The quadratic form for the second phase variance estimator,
#' conditional on the selection of the first phase sample.
#' @param phase_2_joint_probs The matrix of conditional joint
#' inclusion probabilities for the second phase, given the selected
#' first phase sample.
#' @param ensure_psd If \code{TRUE} (the default), ensures
#' that the result is a positive semidefinite matrix. This
#' is necessary if the quadratic form is used as an input for
#' replication methods such as the generalized bootstrap.
#' For details, see the help section entitled
#' "Ensuring the Result is Positive Semidefinite".
#' @return A quadratic form matrix that can be used to estimate
#' the sampling variance from a two-phase sample design.
#' @section Statistical Details:
#' The two-phase variance estimator has a quadratic form matrix \eqn{\boldsymbol{\Sigma}_{ab}} given by:
#' \deqn{
#'   \boldsymbol{\Sigma}_{ab} = {W}^{-1}_b(\boldsymbol{\Sigma}_{a^\prime} \circ D_b ){W}^{-1}_b + \boldsymbol{\Sigma}_b
#' }
#' The first term estimates the variance contribution from the first phase of sampling,
#' while the second term estimates the variance contribution from the second phase of sampling. \cr
#'
#' The full quadratic form of the variance estimator is:
#' \deqn{
#'   v(\hat{t_y}) = \breve{\breve{y^{'}}} \boldsymbol{\Sigma}_{ab} \breve{\breve{y}}
#' }
#' where the weighted variable \eqn{\breve{\breve{y}}_k = \frac{y_k}{\pi_{ak}\pi_{bk}}},
#' is formed using the first phase inclusion probability, denoted \eqn{\pi_{ak}}, and
#' the conditional second phase inclusion probability (given the selected first phase sample),
#' denoted \eqn{\pi_{bk}}. \cr
#'
#' The notation for this estimator is as follows: \cr
#' \itemize{
#'   \item \eqn{n_a} denotes the first phase sample size.
#'
#'   \item \eqn{n_b} denotes the second phase sample size.
#'
#'   \item \eqn{\boldsymbol{\Sigma}_a} denotes the matrix of dimension \eqn{n_a \times n_a}
#'   representing the quadratic form for the variance estimator
#'   used for the full first-phase design.
#'
#'   \item \eqn{\boldsymbol{\Sigma}_{a^\prime}} denotes the matrix of dimension \eqn{n_b \times n_b}
#'   formed by subsetting the rows and columns of \eqn{\boldsymbol{\Sigma}_a} to only include
#'   cases selected in the second-phase sample.
#'
#'   \item \eqn{\boldsymbol{\Sigma}_{b}} denotes
#'   the matrix of dimension \eqn{n_b \times n_b} representing the Horvitz-Thompson
#'   estimator of variance for the second-phase sample, conditional on the selected
#'   first-phase sample.
#'
#'   \item \eqn{\boldsymbol{D}_b} denotes the \eqn{n_b \times n_b} matrix of weights formed by the inverses of
#'   the second-phase joint inclusion probabilities, with element \eqn{kl} equal to \eqn{\pi_{bkl}^{-1}},
#'   where \eqn{\pi_{bkl}} is the conditional probability that units \eqn{k} and \eqn{l} are included
#'   in the second-phase sample, given the selected first-phase sample. Note that this
#'   matrix will often not be positive semidefinite, and so the two-phase variance estimator
#'   has a quadratic form which is not necessarily positive semidefinite.
#'
#'   \item \eqn{\boldsymbol{W}_b} denotes the diagonal \eqn{n_b \times n_b} matrix
#'   whose \eqn{k}-th diagonal entry is the second-phase weight \eqn{\pi_{bk}^{-1}},
#'   where \eqn{\pi_{bk}} is the conditional probability that unit \eqn{k}
#'   is included in the second-phase sample, given the selected first-phase sample.
#' }
#' @section Ensuring the Result is Positive Semidefinite:
#' Note that the matrix \eqn{(\boldsymbol{\Sigma}_{a^\prime} \circ D_b )} may not be
#' positive semidefinite, since the matrix \eqn{D_b} is not guaranteed to be positive semidefinite.
#' If \eqn{(\boldsymbol{\Sigma}_{a^\prime} \circ D_b )} is found not to be positive semidefinite,
#' then it is approximated by the nearest positive semidefinite matrix in the Frobenius norm,
#' using the method of Higham (1988). \cr \cr
#' This approximation is discussed by Beaumont and Patak (2012) in the context
#' of forming replicate weights for two-phase samples. The authors argue that
#' this approximation should lead to only a small overestimation of variance. \cr \cr
#' Since \eqn{(\boldsymbol{\Sigma}_{a^\prime} \circ D_b )}
#' is a real, symmetric matrix, this is equivalent to "zeroing out" negative eigenvalues.
#' To be more precise, denote \eqn{A=(\boldsymbol{\Sigma}_{a^\prime} \circ D_b )}.
#' Then we can form the spectral decomposition \eqn{A=\Gamma \Lambda \Gamma^{\prime}}, where \eqn{\Lambda} is the diagonal matrix
#' whose entries are eigenvalues of \eqn{A}. The method of Higham (1988)
#' is to  approximate
#' \eqn{A} with \eqn{\tilde{A} = \Gamma \Lambda_{+} \Gamma^{\prime}},
#' where the \eqn{ii}-th entry of \eqn{\Lambda_{+}} is \eqn{\max(\Lambda_{ii}, 0)}.
#'
#' @references
#' See Section 7.5 of Tillé (2020) or Section 9.3 of Särndal, Swensson, and Wretman (1992)
#' for an overview of variance estimation for two-phase sampling. In the case where
#' the Horvitz-Thompson variance estimator is used for both phases, the method used in this function
#' is equivalent to equation (9.3.8) of Särndal, Swensson, and Wretman (1992)
#' and equation (7.7) of Tillé (2020). However, this function can be used
#' for any combination of first-phase and second-phase variance estimators,
#' provided that the joint inclusion probabilities from the second-phase design
#' are available and are all nonzero.
#' \cr \cr
#' - Beaumont, Jean-François, and Zdenek Patak. (2012). "On the Generalized Bootstrap for Sample Surveys with Special Attention to Poisson Sampling: Generalized Bootstrap for Sample Surveys."
#' International Statistical Review 80 (1): 127-48.
#' \cr \cr
#' - Higham, N. J. (1988). "\emph{Computing a nearest symmetric positive semidefinite matrix.}" Linear Algebra and Its Applications, 103, 103-118.
#' \cr \cr
#' - Särndal, C.-E., Swensson, B., & Wretman, J. (1992). "\emph{Model Assisted Survey Sampling}." Springer New York.
#' \cr \cr
#' - Tillé, Y. (2020). "\emph{Sampling and estimation from finite populations}." (I. Hekimi, Trans.). Wiley.
#' @md
#' @export
#' @seealso
#' For each phase of sampling, the function
#' \link[svrep]{make_quad_form_matrix} can be used to create
#' the appropriate quadratic form matrix.
#' @examples
#' \dontrun{
#'
#' ## ---------------------- Example 1 ------------------------##
#' ## First phase is a stratified multistage sample            ##
#' ## Second phase is a simple random sample                   ##
#' ##----------------------------------------------------------##
#' data('library_multistage_sample', package = 'svrep')
#'
#' # Load first-phase sample
#'   twophase_sample <- library_multistage_sample
#'
#' # Select second-phase sample
#'   set.seed(2022)
#'
#'   twophase_sample[['SECOND_PHASE_SELECTION']] <- sampling::srswor(
#'     n = 100,
#'     N = nrow(twophase_sample)
#'   ) |> as.logical()
#'
#' # Declare survey design
#'   twophase_design <- twophase(
#'     method = "full",
#'     data = twophase_sample,
#'     # Identify the subset of first-phase elements
#'     # which were selected into the second-phase sample
#'     subset = ~ SECOND_PHASE_SELECTION,
#'     # Describe clusters, probabilities, and population sizes
#'     # at each phase of sampling
#'     id = list(~ PSU_ID + SSU_ID,
#'               ~ 1),
#'     probs = list(~ PSU_SAMPLING_PROB + SSU_SAMPLING_PROB,
#'                  NULL),
#'     fpc = list(~ PSU_POP_SIZE + SSU_POP_SIZE,
#'                NULL)
#'   )
#'
#' # Get quadratic form matrix for the first phase design
#'   first_phase_sigma <- get_design_quad_form(
#'     design = twophase_design$phase1$full,
#'     variance_estimator = "Stratified Multistage SRS"
#'   )
#'
#' # Subset to only include cases sampled in second phase
#'
#'   first_phase_sigma <- first_phase_sigma[twophase_design$subset,
#'                                          twophase_design$subset]
#'
#' # Get quadratic form matrix for the second-phase design
#'   second_phase_sigma <- get_design_quad_form(
#'     design = twophase_design$phase2,
#'     variance_estimator = "Ultimate Cluster"
#'   )
#'
#' # Get second-phase joint probabilities
#'   n <- twophase_design$phase2$fpc$sampsize[1,1]
#'   N <- twophase_design$phase2$fpc$popsize[1,1]
#'
#'   second_phase_joint_probs <- Matrix::Matrix((n/N)*((n-1)/(N-1)),
#'                                      nrow = n, ncol = n)
#'   diag(second_phase_joint_probs) <- rep(n/N, times = n)
#'
#' # Get quadratic form for entire two-phase variance estimator
#'   twophase_quad_form <- make_twophase_quad_form(
#'    sigma_1 = first_phase_sigma,
#'    sigma_2 = second_phase_sigma,
#'    phase_2_joint_probs = second_phase_joint_probs
#'  )
#'
#'  # Use for variance estimation
#'
#'    rep_factors <- make_gen_boot_factors(
#'      Sigma = twophase_quad_form,
#'      num_replicates = 500
#'    )
#'
#'    library(survey)
#'
#'    combined_weights <- 1/twophase_design$prob
#'
#'    twophase_rep_design <- svrepdesign(
#'      data = twophase_sample |>
#'        subset(SECOND_PHASE_SELECTION),
#'      type = 'other',
#'      repweights = rep_factors,
#'      weights = combined_weights,
#'      combined.weights = FALSE,
#'      scale = attr(rep_factors, 'scale'),
#'      rscales = attr(rep_factors, 'rscales')
#'    )
#'
#'    svymean(x = ~ LIBRARIA, design = twophase_rep_design)
#'
#'
#' ## ---------------------- Example 2 ------------------------##
#' ## First phase is a stratified systematic sample            ##
#' ## Second phase is nonresponse, modeled as Poisson sampling ##
#' ##----------------------------------------------------------##
#'
#' data('library_stsys_sample', package = 'svrep')
#'
#' # Determine quadratic form for full first-phase sample variance estimator
#'
#'   full_phase1_quad_form <- make_quad_form_matrix(
#'     variance_estimator = "SD2",
#'     cluster_ids = library_stsys_sample[,'FSCSKEY',drop=FALSE],
#'     strata_ids = library_stsys_sample[,'SAMPLING_STRATUM',drop=FALSE],
#'     strata_pop_sizes = library_stsys_sample[,'STRATUM_POP_SIZE',drop=FALSE],
#'     sort_order = library_stsys_sample$SAMPLING_SORT_ORDER
#'   )
#'
#' # Identify cases included in phase two sample
#' # (in this example, respondents)
#'   phase2_inclusion <- (
#'     library_stsys_sample$RESPONSE_STATUS == "Survey Respondent"
#'   )
#'   phase2_sample <- library_stsys_sample[phase2_inclusion,]
#'
#' # Estimate response propensities
#'
#'   response_propensities <- glm(
#'     data = library_stsys_sample,
#'     family = quasibinomial('logit'),
#'     formula = phase2_inclusion ~ 1,
#'     weights = 1/library_stsys_sample$SAMPLING_PROB
#'   ) |>
#'     predict(type = "response",
#'             newdata = phase2_sample)
#'
#' # Estimate conditional joint inclusion probabilities for second phase
#'
#'   phase2_joint_probs <- outer(response_propensities, response_propensities)
#'   diag(phase2_joint_probs) <- response_propensities
#'
#' # Determine quadratic form for variance estimator of second phase
#' # (Horvitz-Thompson estimator for nonresponse modeled as Poisson sampling)
#'
#'   phase2_quad_form <- make_quad_form_matrix(
#'     variance_estimator = "Horvitz-Thompson",
#'     joint_probs = phase2_joint_probs
#'   )
#'
#' # Create combined quadratic form for entire design
#'
#'  twophase_quad_form <- make_twophase_quad_form(
#'    sigma_1 = full_phase1_quad_form[phase2_inclusion, phase2_inclusion],
#'    sigma_2 = phase2_quad_form,
#'    phase_2_joint_probs = phase2_joint_probs
#'  )
#'
#'  combined_weights <- 1/(phase2_sample$SAMPLING_PROB * response_propensities)
#'
#' # Use for variance estimation
#'
#'   rep_factors <- make_gen_boot_factors(
#'     Sigma = twophase_quad_form,
#'     num_replicates = 500
#'   )
#'
#'   library(survey)
#'
#'   twophase_rep_design <- svrepdesign(
#'     data = phase2_sample,
#'     type = 'other',
#'     repweights = rep_factors,
#'     weights = combined_weights,
#'     combined.weights = FALSE,
#'     scale = attr(rep_factors, 'scale'),
#'     rscales = attr(rep_factors, 'rscales')
#'   )
#'
#'   svymean(x = ~ LIBRARIA, design = twophase_rep_design)
#' }
make_twophase_quad_form <- function(sigma_1, sigma_2, phase_2_joint_probs,
                                    ensure_psd = TRUE) {

  # Diagonal matrix whose entries are second-phase first-order probabilities
  phase2_prob_matrix <- Matrix::diag(Matrix::diag(phase_2_joint_probs))

  # Weighted version of `Sigma_1`
  wtd_sigma_1 <- sigma_1 / phase_2_joint_probs

  # If necessary, approximate `wtd_sigma_1`
  # with the nearest positive semidefinite matrix
  if (ensure_psd && !is_psd_matrix(wtd_sigma_1)) {
    paste(
      "Approximating (sigma_1/phase_2_joint_probs) with the nearest positive semidefinite matrix,",
      "since the matrix (1/phase_2_joint_probs) is not positive semidefinite.",
      "This is expected to result in a small overestimation of variance.",
      "See `help('make_twophase_quad_form', package = 'svrep')` for details."
    ) |> warning()

    wtd_sigma_1 <- get_nearest_psd_matrix(wtd_sigma_1)
  }

  # Combine the quadratic forms from the two phases
  Sigma <- `+`(
    phase2_prob_matrix %*% wtd_sigma_1 %*% phase2_prob_matrix,
    sigma_2
  )

  return(Sigma)
}
