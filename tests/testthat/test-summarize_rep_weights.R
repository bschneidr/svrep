suppressPackageStartupMessages(library(survey))
suppressPackageStartupMessages(library(dplyr))

# Create example data ----
set.seed(1999)
data(api)

cluster_sample <- apiclus1
cluster_sample[['response_status']] <- sample(x = c('R', 'NR', 'UE', 'IE'),
                                              size = nrow(cluster_sample),
                                              replace = TRUE)
cluster_design <- svydesign(id=~dnum, data=cluster_sample, weights = ~ pw)
cluster_rep_design <- as.svrepdesign(cluster_design, type="JK1")
num_replicates <- ncol(cluster_rep_design$repweights)

# Test of error messages ----
test_that("Informative error messages for bad inputs", {

  expect_error(
    object = summarize_rep_weights(cluster_design),
    regexp = "must be a replicate design object"
  )

  expect_error(
    object = summarize_rep_weights(cluster_rep_design, type = "silent"),
    regexp = "must be one of"
  )
  expect_error(
    object = summarize_rep_weights(cluster_rep_design, type = NULL),
    regexp = "must be one of"
  )
})

test_that("Throws informative error for missing `by` variables", {
  expect_error(summarize_rep_weights(
    cluster_rep_design, type = 'overall', by = c("stype", "is_cromulent")
  ), regexp = "The following `by` variables are missing from the data.+cromulent",
  fixed = FALSE
  )
})

# Test the arguments work as expected ----

test_that("`type` argument works", {

  expect_equal(
    object = names(summarize_rep_weights(cluster_rep_design, type = 'both')),
    expected = c("overall", "specific")
  )

  expect_equal(
    object = names(summarize_rep_weights(cluster_rep_design, type = 'overall')),
    expected = c("nrows", "ncols", "degf_svy_pkg", "rank",
                 "avg_wgt_sum", "sd_wgt_sums",
                 "min_rep_wgt", "max_rep_wgt")
  )

  expect_equal(
    object = names(summarize_rep_weights(cluster_rep_design, type = 'specific')),
    expected = c("Rep_Column", "N", "N_NONZERO", "SUM",
                 "MEAN", "CV", "MIN", "MAX")
  )
})

test_that("`by` argument works", {

  ## Check that results have expected columns and number of rows
  grouped_result <- summarize_rep_weights(cluster_rep_design, type = 'both',
                                          by = c("stype", "sch.wide"))

  expect_equal(
    object = colnames(grouped_result$overall),
    expected = c("stype", "sch.wide",
                 "nrows", "ncols", "degf_svy_pkg", "rank",
                 "avg_wgt_sum", "sd_wgt_sums",
                 "min_rep_wgt", "max_rep_wgt")
  )
  expect_equal(
    object = colnames(grouped_result$specific),
    expected = c("stype", "sch.wide",
                 "Rep_Column", "N", "N_NONZERO", "SUM",
                 "MEAN", "CV", "MIN", "MAX")
  )

  num_groups <- nrow(unique(cluster_rep_design$variables[,c("stype", "sch.wide")]))

  expect_equal(
    object = nrow(grouped_result$overall),
    expected = num_groups
  )
  expect_equal(
    object = nrow(grouped_result$specific),
    expected = num_replicates * num_groups
  )

  ## Expected values for basic summaries
  grouped_result <- summarize_rep_weights(cluster_rep_design, type = 'both',
                                          by = c("sch.wide"))

  expected_summaries <- by(data = weights(cluster_rep_design, type = 'analysis'),
     INDICES = cluster_rep_design$variables$sch.wide,
     FUN = function(rep_wts_matrix_subset) {
       summaries <- list(
        'matrix_rank' = qr(x = rep_wts_matrix_subset, tol = 1e-05)[['rank']],
        'nrows' = nrow(rep_wts_matrix_subset),
        'n_nonzero' = apply(rep_wts_matrix_subset, MARGIN = 2,
                            FUN = function(wts) sum(wts > 0))
       )
     })

  expect_equal(
    object = grouped_result$overall$rank,
    expected = c(expected_summaries$No$matrix_rank, expected_summaries$Yes$matrix_rank)
  )

  expect_equal(
    object = grouped_result$overall$nrows,
    expected = c(expected_summaries$No$nrows, expected_summaries$Yes$nrows)
  )

  expect_equal(
    object = grouped_result$specific$N_NONZERO |> unname(),
    expected = c(expected_summaries$No$n_nonzero,
                 expected_summaries$Yes$n_nonzero) |>
      unname()
  )

})

test_that("Works correctly when there are zero rows of data", {
  expect_equal(
    object = cluster_rep_design[0,] |>  summarize_rep_weights(type = 'both'),
    expected = list(
      overall = structure(list(nrows = 0L, ncols = 15L, degf_svy_pkg = NA_integer_,
                               rank = NA_integer_, avg_wgt_sum = NA_real_, sd_wgt_sums = NA_real_,
                               min_rep_wgt = NA_real_, max_rep_wgt = NA_real_),
                          class = "data.frame", row.names = c(NA, -1L)),
      specific = structure(list(Rep_Column = 1:15, N = rep(0,times =15),
                                N_NONZERO = rep(0, times = 15),
                                SUM = rep(NA_real_, times = 15),
                                MEAN = rep(NA_real_, times = 15),
                                CV = rep(NA_real_, times = 15),
                                MIN = rep(NA_real_, times = 15),
                                MAX = rep(NA_real_, times = 15)),
                           class = "data.frame", row.names = c(NA, -15L))
    )
  )
})
