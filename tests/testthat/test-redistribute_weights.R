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

basic_example_output <- redistribute_weights(
  design = cluster_rep_design,
  reduce_if = response_status %in% c("NR"),
  increase_if = response_status %in% c("R")
)

reduced_cases <- cluster_rep_design$variables[['response_status']] %in% "NR"
increased_cases <- cluster_rep_design$variables[['response_status']] %in% c("R")
unchanged_cases <- !reduced_cases & !increased_cases

orig_fullsample_wt_sum <- sum(cluster_rep_design$pweights)
updated_fullsample_wt_sum <- sum(basic_example_output$pweights)
orig_repwt_sums <- apply(weights(cluster_rep_design, 'analysis'), MARGIN = 2, sum)
updated_repwt_sums <- apply(weights(basic_example_output, 'analysis'), MARGIN = 2, sum)

# Test of expected results from basic example ----

test_that("Full sample weights updated correctly", {

  expect_equal(
    object = unique(unname(basic_example_output$pweights[reduced_cases])),
    expected = 0
  )

  expect_equal(
    object = basic_example_output$pweights[unchanged_cases],
    expected = cluster_rep_design$pweights[unchanged_cases]
  )

  expect_equal(object = updated_fullsample_wt_sum,
               expected = orig_fullsample_wt_sum)
})

test_that("Replicate weights updated correctly", {

  uncompressed_input_wts <- weights(cluster_rep_design, 'analysis')
  uncompressed_result_wts <- weights(basic_example_output, 'analysis')
  dimnames(uncompressed_input_wts) <- NULL
  dimnames(uncompressed_result_wts) <- NULL

  expect_equal(
    object = unique(as.vector(uncompressed_result_wts[reduced_cases,])),
    expected = 0
  )

  expect_equal(
    object = uncompressed_result_wts[unchanged_cases,],
    expected = uncompressed_input_wts[unchanged_cases,]
  )

  expect_equal(object = updated_repwt_sums,
               expected = orig_repwt_sums)
})

# Test of example with grouping variables ----

cluster_rep_design[['variables']][['stype_E_or_M']] <- ifelse(
  cluster_rep_design[['variables']][['stype']] %in% c("E", "M"), 1, 0
)

grouped_example_output <- redistribute_weights(
  design = cluster_rep_design,
  reduce_if = response_status %in% c("NR"),
  increase_if = response_status %in% c("R"),
  by = c("stype_E_or_M", "stype")
)

reduced_cases <- cluster_rep_design$variables[['response_status']] %in% "NR"
increased_cases <- cluster_rep_design$variables[['response_status']] %in% c("R")
unchanged_cases <- !reduced_cases & !increased_cases

orig_fullsample_wt_sum <- sum(cluster_rep_design$pweights)
updated_fullsample_wt_sum <- sum(grouped_example_output$pweights)
orig_repwt_sums <- apply(weights(cluster_rep_design, 'analysis'), MARGIN = 2, sum)
updated_repwt_sums <- apply(weights(grouped_example_output, 'analysis'), MARGIN = 2, sum)

test_that("With grouping variables, full sample weights updated correctly", {

  expect_equal(
    object = unique(unname(grouped_example_output$pweights[reduced_cases])),
    expected = 0
  )

  expect_equal(
    object = grouped_example_output$pweights[unchanged_cases],
    expected = cluster_rep_design$pweights[unchanged_cases]
  )

  expect_equal(object = updated_fullsample_wt_sum,
               expected = orig_fullsample_wt_sum)
})

test_that("With grouping variables, replicate weights updated correctly", {

  uncompressed_input_wts <- weights(cluster_rep_design, 'analysis')
  uncompressed_result_wts <- weights(grouped_example_output, 'analysis')
  dimnames(uncompressed_input_wts) <- NULL
  dimnames(uncompressed_result_wts) <- NULL

  expect_equal(
    object = unique(as.vector(uncompressed_result_wts[reduced_cases,])),
    expected = 0
  )

  expect_equal(
    object = uncompressed_result_wts[unchanged_cases,],
    expected = uncompressed_input_wts[unchanged_cases,]
  )

  expect_equal(object = updated_repwt_sums,
               expected = orig_repwt_sums)
})

# Test that works with replicate design created with provided replicate weights ----

standalone_wts <- `colnames<-`(weights(cluster_rep_design, 'analysis'),
                               paste0("Rep_Set_", 1:ncol(cluster_rep_design$repweights$weights)))

data_w_repweights <- cbind(cluster_rep_design$variables,
                           as.data.frame(standalone_wts))

hand_created_rep_design <- svrepdesign(
  data = data_w_repweights,
  weights = ~ pw,
  repweights = "Rep_Set_", type = "JK1",
  scale = (ncol(standalone_wts) - 1)/ncol(standalone_wts),
  combined = TRUE
)

test_that("Works with replicate design created using provided weights in data", {

  result_autocreated_wts <- weights(redistribute_weights(cluster_rep_design,
                                                         increase_if = response_status == "R",
                                                         reduce_if = response_status == "NR"),
                                    'analysis')
  result_handcreated_wts <- weights(redistribute_weights(hand_created_rep_design,
                                                         increase_if = response_status == "R",
                                                         reduce_if = response_status == "NR"),
                                    'analysis')
  dimnames(result_autocreated_wts) <- NULL
  dimnames(result_handcreated_wts) <- NULL

  expect_equal(
    expected = result_autocreated_wts,
    object = result_handcreated_wts
  )
})

# Tests of expected error messages ----
test_that("Throws error if input is not a replicate design object.", {
  expect_error(redistribute_weights(
    design = cluster_design,
    reduce_if = response_status %in% c("UE"),
    increase_if = response_status %in% c("R", "NR", "IE")
  ), regexp = "`design` must be a replicate design object.", fixed = TRUE)
})
test_that("Throws informative error for non-logical results from `reduce_if`/`increase_if`", {
  expect_error(redistribute_weights(
    design = cluster_rep_design,
    reduce_if = 2,
    increase_if = response_status %in% c("R", "NR", "IE")
  ), regexp = "The expressions supplied to `reduce_if` and `increase_if` must result in logical values of TRUE or FALSE.",
  fixed = TRUE
  )
})
test_that("Throws informative error for missing values in `reduce_if`/`increase_if`", {
  expect_error(redistribute_weights(
    design = cluster_rep_design,
    reduce_if = sample(c(TRUE, FALSE, NA),
                       size = nrow(cluster_rep_design),
                       replace = TRUE),
    increase_if = response_status %in% c("R", "NR", "IE")
  ), regexp = "The result of the expressions supplied to `reduce_if` and `increase_if` must be TRUE or FALSE, not NA.",
  fixed = TRUE
  )
})
test_that("Throws informative error for conflicting results from `reduce_if`/`increase_if`", {
  expect_error(redistribute_weights(
    design = cluster_rep_design,
    reduce_if = response_status %in% c("NR", "IE"),
    increase_if = response_status %in% c("R", "NR", "IE")
  ), regexp = "`reduce_if` and `increase_if` conflict: they imply that some cases should have weights simultaneously reduced and increased.",
  fixed = TRUE
  )
})
test_that("Throws informative error for missing `by` variables", {
  expect_error(redistribute_weights(
    design = cluster_rep_design,
    reduce_if = response_status %in% c("UE"),
    increase_if = response_status %in% c("R", "NR", "IE"),
    by = "nonexistantvar"
  ), regexp = "The following `by` variables are missing from the data.+nonexistantvar",
  fixed = FALSE
  )
})
test_that("Throws informative error for missing `reduce_if`/`increase_if` arguments", {
  expect_error(redistribute_weights(
    design = cluster_rep_design,
    increase_if = response_status %in% c("R", "NR", "IE"),
    by = "stype"
  ), regexp = "Must supply expressions",
  fixed = FALSE
  )
})
