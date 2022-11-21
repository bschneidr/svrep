suppressPackageStartupMessages(library(survey))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(labelled))
suppressPackageStartupMessages(library(svrep))

data('library_stsys_sample', package = 'svrep')

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

  small_example <- library_stsys_sample |>
    arrange(SAMPLING_STRATUM) |>
    head(4)

  small_example <- small_example[sample(x = nrow(small_example),
                                        size = nrow(small_example),
                                        replace = FALSE),]
  sorted_small_example <- small_example |>
    arrange(SAMPLING_SORT_ORDER)

  test_that(
    "`make_quad_form_matrix()` with 'SD1' and 'SD2' works correctly when data frame isn't sorted", {

      row_ordering <- rank(small_example[['SAMPLING_SORT_ORDER']])

      unsorted_quad_form_matrix <- make_quad_form_matrix(
        variance_estimator = 'SD1',
        cluster_ids = small_example[,'FSCSKEY',drop=FALSE],
        strata_ids = small_example[,'SAMPLING_STRATUM',drop=FALSE],
        strata_pop_sizes = small_example[,'STRATUM_POP_SIZE',drop=FALSE],
        sort_order = small_example[['SAMPLING_SORT_ORDER']]
      )

      sorted_quad_form_matrix <- make_quad_form_matrix(
        variance_estimator = 'SD1',
        cluster_ids = sorted_small_example[,'FSCSKEY',drop=FALSE],
        strata_ids = sorted_small_example[,'SAMPLING_STRATUM',drop=FALSE],
        strata_pop_sizes = sorted_small_example[,'STRATUM_POP_SIZE',drop=FALSE],
        sort_order = sorted_small_example[['SAMPLING_SORT_ORDER']]
      )

      y_wtd_sorted <- sorted_small_example$TOTCIR/sorted_small_example$SAMPLING_PROB
      y_wtd_unsorted <- small_example$TOTCIR/small_example$SAMPLING_PROB

      expect_equal(
        object = unsorted_quad_form_matrix[row_ordering, row_ordering],
        expected = sorted_quad_form_matrix
      )
      expect_equal(
        object = t(y_wtd_sorted) %*% sorted_quad_form_matrix %*% y_wtd_sorted,
        expected = t(y_wtd_unsorted) %*% unsorted_quad_form_matrix %*% y_wtd_unsorted
      )
    })


