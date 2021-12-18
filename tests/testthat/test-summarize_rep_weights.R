suppressPackageStartupMessages(library(survey))

# Create example data ----
set.seed(1999)
data(api)

cluster_sample <- apiclus1
cluster_sample[['response_status']] <- sample(x = c('R', 'NR', 'UE', 'IE'),
                                              size = nrow(cluster_sample),
                                              replace = TRUE)
cluster_design <- svydesign(id=~dnum, data=cluster_sample, weights = ~ pw)
cluster_rep_design <- as.svrepdesign(cluster_design, type="JK1")

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

# Test the arguments work as expected ----

test_that("`type` argument works", {

  expect_equal(
    object = names(summarize_rep_weights(cluster_rep_design, type = 'both')),
    expected = c("overall", "specific")
  )

  expect_equal(
    object = names(summarize_rep_weights(cluster_rep_design, type = 'overall')),
    expected = c("nrows", "ncols", "degf_svy_pkg", "rank",
                 "min_rep_wgt", "max_rep_wgt", "max_CV", "SD_wgt_sums")
  )

  expect_equal(
    object = names(summarize_rep_weights(cluster_rep_design, type = 'specific')),
    expected = c("Rep_Column", "N", "N_NONZERO", "SUM",
                 "MEAN", "CV", "L", "MIN", "MAX")
  )
})
