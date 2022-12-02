suppressWarnings({
  suppressPackageStartupMessages({
    library(survey)
    library(dplyr)
    library(svrep)
    library(testthat)
  })
})

data('library_stsys_sample', package = 'svrep')
set.seed(2014)

library_stsys_sample <- library_stsys_sample |>
  mutate(
    TOTCIR = ifelse(is.na(TOTCIR), 0, TOTCIR),
    TOTSTAFF = ifelse(is.na(TOTSTAFF), 0, TOTSTAFF)
  )

# Check basic successive-difference quadratic forms ----

  successive_diffs <- function(y, type = "SD2") {
    n <- length(y)
    sum_of_squares <- sum(diff(y)^2)
    if (n == 1) {
      result <- 0
    }
    if (type == "SD1" && (n > 1)) {
      result <- sum_of_squares * (n/(2*(n-1)))
    }
    if (type == "SD2" && (n > 1)) {
      sum_of_squares <- sum_of_squares + (tail(y, 1) - head(y, 1))^2
      result <- sum_of_squares/2
    }
    return(result)
  }

  wtd_y <- library_stsys_sample[['TOTCIR']]/library_stsys_sample[['SAMPLING_PROB']]
  wtd_y <- wtd_y[order(library_stsys_sample[['SAMPLING_SORT_ORDER']])]

  test_that(
    "Helper function for successive difference quadratic forms works correctly", {
      # Basic results correct
      expect_equal(
        object = t(wtd_y) %*% svrep:::make_sd_matrix(n = length(wtd_y), type = "SD2") %*% wtd_y,
        expected = successive_diffs(wtd_y, type = "SD2") |> as.matrix()
      )
      expect_equal(
        object = t(wtd_y) %*% svrep:::make_sd_matrix(n = length(wtd_y), type = "SD1") %*% wtd_y,
        expected = successive_diffs(wtd_y, type = "SD1") |> as.matrix()
      )
      # Checks on FPC
      expect_equal(
        object = t(wtd_y) %*% svrep:::make_sd_matrix(n = length(wtd_y), type = "SD1",
                                                     f = 0.9) %*% wtd_y,
        expected = (1-0.9) * successive_diffs(wtd_y, type = "SD1") |> as.matrix()
      )
      expect_equal(
        object = t(wtd_y) %*% svrep:::make_sd_matrix(n = length(wtd_y), type = "SD2",
                                                     f = 0.9) %*% wtd_y,
        expected = (1-0.9) * successive_diffs(wtd_y, type = "SD2") |> as.matrix()
      )
      # Correct result for only a single unit
      expect_equal(
        object = svrep:::make_sd_matrix(n = 1),
        expected = matrix(0, nrow = 1, ncol = 1)
      )
      # Correct result for only two units
      expect_equal(
        object = c(99, 16) %*% svrep:::make_sd_matrix(n = 2, type = "SD1") %*% c(99,16),
        expected = as.matrix(successive_diffs(c(99,16), type = "SD1"))
      )
      expect_equal(
        object = c(99, 16) %*% svrep:::make_sd_matrix(n = 2, type = "SD2") %*% c(99,16),
        expected = as.matrix(successive_diffs(c(99,16), type = "SD2"))
      )
      # Checks on numeric arguments
      expect_error(
        object = svrep:::make_sd_matrix(n = 0),
        regexp = "must be an integer greater than"
      )
      expect_error(
        object = svrep:::make_sd_matrix(n = 2, f = 1.1),
        regexp = "must be a single number between"
      )
    })

# Check basic SRSWOR quadratic form ----

  wtd_y <- library_stsys_sample[['TOTCIR']]/library_stsys_sample[['SAMPLING_PROB']]
  n <- length(wtd_y)

  wtd_y_matrix <- as.matrix(
    library_stsys_sample[,c("TOTCIR", "TOTSTAFF")]/library_stsys_sample[['SAMPLING_PROB']]
  )

  test_that(
    "Helper function for SRSWOR quadratic forms works correctly", {
      # Basic results correct
      expect_equal(
        object = t(wtd_y) %*% svrep:::make_srswor_matrix(n = length(wtd_y)) %*% wtd_y,
        expected = n * cov(as.matrix(wtd_y))
      )
      # Matches 'survey' package
      expect_equal(
        object = t(wtd_y_matrix) %*% svrep:::make_srswor_matrix(n = length(wtd_y)) %*% wtd_y_matrix,
        expected = svyrecvar(
          x = wtd_y_matrix,
          clusters = as.matrix(seq_len(n)),
          stratas = as.matrix(rep(1, times = n)),
          fpcs = list('sampsize' = as.matrix(rep(n, times = n)),
                      'popsize' = as.matrix(rep(Inf, times = n))),
          lonely.psu = "remove", one.stage = TRUE
        )
      )
      # Checks on FPC
      expect_equal(
        object = t(wtd_y) %*% svrep:::make_srswor_matrix(n = length(wtd_y),
                                                         f = 0.9) %*% wtd_y,
        expected = (1-0.9) * n *  cov(as.matrix(wtd_y)) |> as.matrix()
      )
      # Correct result for only a single unit
      expect_equal(
        object = svrep:::make_srswor_matrix(n = 1),
        expected = matrix(0, nrow = 1, ncol = 1)
      )
      # Checks on numeric arguments
      expect_error(
        object = svrep:::make_srswor_matrix(n = 0),
        regexp = "must be an integer greater than"
      )
      expect_error(
        object = svrep:::make_srswor_matrix(n = 2, f = 1.1),
        regexp = "must be a single number between"
      )
    })

# Check Horvitz-Thompson results ----

  ## Load an example dataset that uses unequal probability sampling
  data('election', package = 'survey')
  y_wtd <- election_pps$Bush/diag(election_jointprob)

  ## Create matrix to represent the Horvitz-Thompson estimator as a quadratic form
  n <- nrow(election_pps)
  pi <- election_jointprob
  horvitz_thompson_matrix <- matrix(nrow = n, ncol = n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      horvitz_thompson_matrix[i,j] <- 1 - (pi[i,i] * pi[j,j])/pi[i,j]
    }
  }

  ht_quad_form <- make_quad_form_matrix(
    variance_estimator = "Horvitz-Thompson",
    joint_probs = election_jointprob
  )

  test_that(
    "Generates correct quadratic form for Horvitz-Thompson", {
    expect_equal(object = ht_quad_form, expected = horvitz_thompson_matrix)
  })

  svy_pkg_result <- svydesign(
    data = election_pps,
    id = ~1, fpc = ~p,
    pps = ppsmat(election_jointprob),
    variance = "HT"
  ) |> svytotal(x = ~ Bush) |> vcov()

  quad_form_result <- t(y_wtd) %*% ht_quad_form %*% y_wtd

  test_that(
    "Generates correct quadratic form for Horvitz-Thompson", {
      expect_equal(object = ht_quad_form, expected = horvitz_thompson_matrix)
      expect_equal(object = quad_form_result, expected = svy_pkg_result)
    })

# Check Sen-Yates-Grundy results ----

  y_wtd <- election_pps$Bush/diag(election_jointprob)

  syg_result <- as.matrix(0)
  for (i in seq_along(y_wtd)) {
    for (j in seq_along(y_wtd)) {
      pi_ij <- election_jointprob[i,j]
      pi_i <- election_jointprob[i,i]
      pi_j <- election_jointprob[j,j]
      contribution <- -0.5*(1 - (pi_i*pi_j)/pi_ij) * (y_wtd[i] - y_wtd[j])^2
      syg_result <- syg_result + contribution
    }
  }

  yg_quad_form_matrix <- make_quad_form_matrix(variance_estimator = "Yates-Grundy",
                                               joint_probs = election_jointprob)

  quad_form_result <- t(y_wtd) %*% yg_quad_form_matrix %*% y_wtd


  svy_pkg_result <- svydesign(
    data = election_pps,
    id = ~1, fpc = ~p,
    pps = ppsmat(election_jointprob),
    variance = "YG"
  ) |> svytotal(x = ~ Bush) |> vcov()

  test_that(
    "Generates correct quadratic form for Yates-Grundy", {
    expect_equal(object = quad_form_result, expected = syg_result)
    expect_equal(object = quad_form_result, expected = svy_pkg_result)
  })



# Check SD1 and SD2 results ----

  ##_ Basic correctness checks
  sd1_quad_form <- make_quad_form_matrix(
    variance_estimator = 'SD1',
    cluster_ids = library_stsys_sample[,'FSCSKEY',drop=FALSE],
    strata_ids = library_stsys_sample[,'SAMPLING_STRATUM',drop=FALSE],
    strata_pop_sizes = library_stsys_sample[,'STRATUM_POP_SIZE',drop=FALSE],
    sort_order = library_stsys_sample[['SAMPLING_SORT_ORDER']]
  )
  sd2_quad_form <- make_quad_form_matrix(
    variance_estimator = 'SD2',
    cluster_ids = library_stsys_sample[,'FSCSKEY',drop=FALSE],
    strata_ids = library_stsys_sample[,'SAMPLING_STRATUM',drop=FALSE],
    strata_pop_sizes = library_stsys_sample[,'STRATUM_POP_SIZE',drop=FALSE],
    sort_order = library_stsys_sample[['SAMPLING_SORT_ORDER']]
  )

  wtd_y <- as.matrix(library_stsys_sample[['TOTCIR']] /
                      library_stsys_sample$SAMPLING_PROB)

  expected_results <- library_stsys_sample |>
    arrange(SAMPLING_SORT_ORDER) |>
    group_by(SAMPLING_STRATUM) |>
    summarize(
      v_h_sd1 = (1 - unique(SAMPLING_PROB)) * successive_diffs(y = TOTCIR/SAMPLING_PROB, type = "SD1"),
      v_h_sd2 = (1 - unique(SAMPLING_PROB)) * successive_diffs(y = TOTCIR/SAMPLING_PROB, type = "SD2")
    ) |>
    summarize(
      v_sd1 = sum(v_h_sd1), v_sd2 = sum(v_h_sd2)
    ) |> unlist()

  test_that(
    "`make_quad_form_matrix()` works correctly for 'SD1' and 'SD2' with stratified single-stage samples", {
      expect_equal(
        object = as.numeric(t(wtd_y) %*% sd1_quad_form %*% wtd_y),
        expected = unname(expected_results['v_sd1'])
      )
      expect_equal(
        object = as.numeric(t(wtd_y) %*% sd2_quad_form %*% wtd_y),
        expected = unname(expected_results['v_sd2'])
      )
    })

  ##_ Check that works if rows are not in sampling sort order

  shuffled_data <- library_stsys_sample[sample(nrow(library_stsys_sample),
                                               size = nrow(library_stsys_sample)),]
  sorted_data <- shuffled_data |>
    arrange(SAMPLING_SORT_ORDER)

  test_that(
    "`make_quad_form_matrix()` with 'SD1' and 'SD2' works correctly when data frame isn't sorted", {

      unsorted_quad_form_matrix <- make_quad_form_matrix(
        variance_estimator = 'SD1',
        cluster_ids = shuffled_data[,'FSCSKEY',drop=FALSE],
        strata_ids = shuffled_data[,'SAMPLING_STRATUM',drop=FALSE],
        strata_pop_sizes = shuffled_data[,'STRATUM_POP_SIZE',drop=FALSE],
        sort_order = shuffled_data[['SAMPLING_SORT_ORDER']]
      )

      sorted_quad_form_matrix <- make_quad_form_matrix(
        variance_estimator = 'SD1',
        cluster_ids = sorted_data[,'FSCSKEY',drop=FALSE],
        strata_ids = sorted_data[,'SAMPLING_STRATUM',drop=FALSE],
        strata_pop_sizes = sorted_data[,'STRATUM_POP_SIZE',drop=FALSE],
        sort_order = sorted_data[['SAMPLING_SORT_ORDER']]
      )

      y_wtd_sorted <- sorted_data$TOTCIR/sorted_data$SAMPLING_PROB
      y_wtd_unsorted <- shuffled_data$TOTCIR/shuffled_data$SAMPLING_PROB

      sorted_result <- t(y_wtd_sorted) %*% sorted_quad_form_matrix %*% y_wtd_sorted
      unsorted_result <- t(y_wtd_unsorted) %*% unsorted_quad_form_matrix %*% y_wtd_unsorted

      expect_equal(
        object = sorted_result,
        expected = unsorted_result
      )

    })

  test_that(
    "successive-differences estimate less than ST-SRSWOR approximation, for sys sample from sorted frame", {

      expect_lt(
        object = as.numeric(t(wtd_y) %*% sd1_quad_form %*% wtd_y),
        expected = svydesign(
          data = library_stsys_sample,
          strata = ~ SAMPLING_STRATUM,
          fpc = ~ SAMPLING_PROB,
          ids = ~ 1
        ) |> svytotal(x = ~TOTCIR) |> vcov()
      )
  })

# Check "Ultimate Cluster" results ----

  test_that(
    "Correct results for ultimate cluster", {

      # Correct result for a single-stage stratified sample
      wtd_y_matrix <- as.matrix(
        library_stsys_sample[,c("TOTCIR", "TOTSTAFF")]
      )/library_stsys_sample[['SAMPLING_PROB']]
      quad_UC <- make_quad_form_matrix(
        variance_estimator = "Ultimate Cluster",
        cluster_ids = library_stsys_sample[,'FSCSKEY',drop=FALSE],
        strata_ids = library_stsys_sample[,'SAMPLING_STRATUM',drop=FALSE],
        strata_pop_sizes = library_stsys_sample[,'STRATUM_POP_SIZE',drop=FALSE]
      )
      expect_equal(
        object = t(wtd_y_matrix) %*% quad_UC %*% wtd_y_matrix,
        expected = svydesign(
          data = library_stsys_sample,
          ids = ~ FSCSKEY,
          strata = ~ SAMPLING_STRATUM,
          fpc = ~ STRATUM_POP_SIZE
        ) |> svytotal(x = ~ TOTCIR + TOTSTAFF, na.rm = TRUE) |>
          vcov()
      )

      # Correct result for a multistage cluster sample
      wtd_y_matrix <- as.matrix(
        library_multistage_sample[,c("TOTCIR", "TOTSTAFF")]
      )/library_multistage_sample[['SAMPLING_PROB']]
      wtd_y_matrix[rowSums(is.na(wtd_y_matrix)) > 0] <- 0
      quad_UC <- make_quad_form_matrix(
        variance_estimator = "Ultimate Cluster",
        cluster_ids = library_multistage_sample[,c("PSU_ID", "SSU_ID"),drop=FALSE],
        strata_ids = data.frame('PSU_STRATUM' = rep(1, times = nrow(library_multistage_sample)),
                                'SSU_STRATUM' = rep(1, times = nrow(library_multistage_sample))),
        strata_pop_sizes = library_multistage_sample[,c("PSU_POP_SIZE", "SSU_POP_SIZE"),drop=FALSE]
      )
      expect_equal(
        object = t(wtd_y_matrix) %*% quad_UC %*% wtd_y_matrix,
        expected = svydesign(
          data = library_multistage_sample |>
            mutate(SAMPLING_WEIGHT = 1/SAMPLING_PROB),
          weights = ~ SAMPLING_WEIGHT,
          ids = ~ PSU_ID,
          strata = NULL,
          fpc = ~ PSU_POP_SIZE
        ) |> svytotal(x = ~ TOTCIR + TOTSTAFF, na.rm = TRUE) |>
          vcov()
      )
  })

# Check "Stratified Multistage SRS" results ----

  test_that(
    "Correct results for stratified multistage SRS", {

      # Correct result for a single-stage stratified sample
      wtd_y_matrix <- as.matrix(
        library_stsys_sample[,c("TOTCIR", "TOTSTAFF")]
      )/library_stsys_sample[['SAMPLING_PROB']]
      quad_UC <- make_quad_form_matrix(
        variance_estimator = "Stratified Multistage SRS",
        cluster_ids = library_stsys_sample[,'FSCSKEY',drop=FALSE],
        strata_ids = library_stsys_sample[,'SAMPLING_STRATUM',drop=FALSE],
        strata_pop_sizes = library_stsys_sample[,'STRATUM_POP_SIZE',drop=FALSE]
      )
      expect_equal(
        object = t(wtd_y_matrix) %*% quad_UC %*% wtd_y_matrix,
        expected = svydesign(
          data = library_stsys_sample,
          ids = ~ FSCSKEY,
          strata = ~ SAMPLING_STRATUM,
          fpc = ~ STRATUM_POP_SIZE
        ) |> svytotal(x = ~ TOTCIR + TOTSTAFF, na.rm = TRUE) |>
          vcov()
      )

      # Correct result for a multistage cluster sample
      data('mu284', package = 'survey')
      mu284 <- mu284 |> arrange(id1, id2)
      mu284_design <- survey::svydesign(
        data = mu284,
        ids = ~ id1 + id2,
        fpc = ~ n1 + n2
      )
      wtd_y_matrix <- as.matrix(
        mu284_design$variables$y1
      ) * weights(mu284_design)
      wtd_y_matrix[rowSums(is.na(wtd_y_matrix)) > 0] <- 0
      colnames(wtd_y_matrix) <- c("y1")

      # Use function to produce quadratic form matrix
      qf_matrix <- make_quad_form_matrix(
        variance_estimator = "Stratified Multistage SRS",
        cluster_ids = mu284_design$cluster,
        strata_ids = mu284_design$strata,
        strata_pop_sizes = mu284_design$fpc$popsize
      )

      # Determine expected correct quadratic form matrix
      quad_form_UC <- make_quad_form_matrix(
        variance_estimator = "Ultimate Cluster",
        cluster_ids = mu284_design$cluster,
        strata_ids = mu284_design$strata,
        strata_pop_sizes = mu284_design$fpc$popsize
      )
      later_stage_quad_form <- make_quad_form_matrix(
        variance_estimator = "Ultimate Cluster",
        cluster_ids = mu284_design$cluster[,2,drop=FALSE],
        strata_ids = mu284_design$strata[,2,drop=FALSE],
        strata_pop_sizes = mu284_design$fpc$popsize[,2,drop=FALSE]
      )
      expected_quad_form <- quad_form_UC + (0.1 * later_stage_quad_form)

      # Compare to survey package
      expect_equal(
        object = t(wtd_y_matrix) %*% qf_matrix %*% wtd_y_matrix,
        expected = mu284_design |>
          svytotal(x = ~ y1, na.rm = TRUE) |>
          vcov()
      )
    })

# Ensure function checks inputs for issues ----

  test_that(
    "Informative errors for bad inputs", {

      # Horvitz-Thompson / Yates-Grundy
      expect_error(
        object = make_quad_form_matrix(variance_estimator = "Yates-Grundy"),
        regexp = "must supply a matrix.+joint_probs"
      )
      expect_error(
        object = make_quad_form_matrix(variance_estimator = "Yates-Grundy",
                                       joint_probs = matrix(NA, 2, 2)),
        regexp = "must be a matrix.+no missing values"
      )
      # SD1 and SD2
      expect_error(
        object = make_quad_form_matrix(variance_estimator = "SD1"),
        regexp = "must supply a matrix or data frame to `cluster_ids`"
      )
      expect_error(
        object = make_quad_form_matrix(variance_estimator = "SD1",
                                       cluster_ids = data.frame(ID = c(1,2))),
        regexp = "must supply a vector to `sort_order`"
      )

  })


