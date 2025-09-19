
# Create example data ----
set.seed(2023)

# Create an example survey design object ----

sample_data <- data.frame(
  Y       = c(8.11, 7.42, 14.54),
  STRATUM = c(1,1,1),
  PSU     = c(1,2,3),
  FPC     = c(0.5, 0.5, 0.5),
  ALT_WGT = c(10, 15, 20)
)

survey_design <- svydesign(
  data    = sample_data,
  strata  = ~ STRATUM,
  ids     = ~ PSU,
  weights = ~ 1,
  fpc     = ~ FPC
)

orig_comp_rep_design <- as.svrepdesign(
  design = survey_design,
  type = "JKn", mse = TRUE,
  compress = TRUE
) |> compressWeights()

orig_uncomp_rep_design <- as.svrepdesign(
  design = survey_design,
  type = "JKn", mse = TRUE
)

# Check replicates added in correct place ----

test_that("Able to correctly specify location", {

  # Check correct placement for `location = "first"`
  added_first <- add_inactive_replicates(
    design = orig_comp_rep_design,
    n_to_add = 2,
    location = "first"
  )
  expect_equal(
    object = added_first |> weights(type = "analysis"),
    expected = matrix(
      c(1,1,1,1,1,1,
        0, 1.5, 1.5,
        1.5, 0, 1.5,
        1.5, 1.5, 0),
      nrow = 3, ncol = 5,
      byrow = FALSE
    )
  )
  expect_equal(object = added_first$rscales,
               expected = c(1,1,rep(1/3, times = 3)))

  # Check correct placement for `location = "last"`
  added_last <- add_inactive_replicates(
    design = orig_comp_rep_design,
    n_to_add = 2,
    location = "last"
  )
  expect_equal(
    object = added_last |> weights(type = "analysis"),
    expected = matrix(
      c(0, 1.5, 1.5,
        1.5, 0, 1.5,
        1.5, 1.5, 0,
        1,1,1,1,1,1),
      nrow = 3, ncol = 5,
      byrow = FALSE
    )
  )
  expect_equal(object = added_last$rscales,
               expected = c(rep(1/3, times = 3), 1, 1))

  # Check correct results for `location = "random"`
  expect_equal(
    object = add_inactive_replicates(
      design = orig_comp_rep_design,
      n_to_add = 2,
      location = "random"
    ) |> weights(type = "analysis") |> t() |> cov(),
    expected = matrix(
      c(0, 1.5, 1.5,
        1.5, 0, 1.5,
        1.5, 1.5, 0,
        1,1,1,1,1,1),
      nrow = 3, ncol = 5,
      byrow = FALSE
    ) |> t() |> cov()
  )

})

test_that("Correct results when `n_total` is LTE existing number of replicates", {

  expect_equal(
    object = orig_uncomp_rep_design,
    expected = orig_uncomp_rep_design |>
      add_inactive_replicates(n_total = ncol(orig_uncomp_rep_design$repweights))
  )

})

# Edge case handling ----

test_that("Emits warning when `mse=FALSE`", {
  expect_warning(
    object = as.svrepdesign(
      design = survey_design,
      type = "JKn", mse = FALSE
    ) |> add_inactive_replicates(n_to_add = 2, location = 'last'),
    regexp = "variance estimates may differ before and after adding inactive replicates"
  )
})

test_that("Correct results regardless of `combined.weights` option", {

  alt_survey_design <- svydesign(
    data    = sample_data,
    strata  = ~ STRATUM,
    ids     = ~ PSU,
    weights = ~ ALT_WGT
  )

  rep_design_separate <- alt_survey_design |> 
    as.svrepdesign(type = "JKn", mse = TRUE, compress = FALSE)

  rep_design_combined <- svrepdesign(
    data = alt_survey_design,
    weights = weights(rep_design_separate, type = 'sampling'),
    repweights = weights(rep_design_separate, type = 'analysis'),
    type = 'JKn', rscales = rep(2/3, times = 3), mse = TRUE
  )

  expect_equal(
    object = rep_design_separate |> 
      add_inactive_replicates(n_to_add = 2, location = 'last') |>
      weights(type = 'analysis'),
    expected = rep_design_combined |> 
      add_inactive_replicates(n_to_add = 2, location = 'last') |>
      weights(type = 'analysis')
  )

})
