suppressWarnings({
  suppressPackageStartupMessages({
    library(survey)
    library(dplyr)
    library(labelled)
    library(svrep)
  })
})

data('library_stsys_sample', package = 'svrep')
set.seed(2014)

library_stsys_sample <- library_stsys_sample |>
  mutate(
    TOTCIR = ifelse(is.na(TOTCIR), 0, TOTCIR)
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
      expect_equal(
        object = t(wtd_y) %*% svrep:::make_sd_matrix(n = length(wtd_y), type = "SD2") %*% wtd_y,
        expected = successive_diffs(wtd_y, type = "SD2") |> as.matrix()
      )
      expect_equal(
        object = t(wtd_y) %*% svrep:::make_sd_matrix(n = length(wtd_y), type = "SD1") %*% wtd_y,
        expected = successive_diffs(wtd_y, type = "SD1") |> as.matrix()
      )
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
      expect_equal(
        object = svrep:::make_sd_matrix(n = 1),
        expected = matrix(0, nrow = 1, ncol = 1)
      )
      expect_error(
        object = svrep:::make_sd_matrix(n = 0),
        regexp = "must be an integer greater than"
      )
      expect_error(
        object = svrep:::make_sd_matrix(n = 2, f = 1.1),
        regexp = "must be a single number between"
      )
    })

# Check Horvitz-Thompson results ----

  # Load an example dataset that uses unequal probability sampling ----
  data('election', package = 'survey')

  # Create matrix to represent the Horvitz-Thompson estimator as a quadratic form ----
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

  test_that(
    "Generates correct quadratic form for Yates-Grundy", {
    expect_equal(object = quad_form_result, expected = syg_result)
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




