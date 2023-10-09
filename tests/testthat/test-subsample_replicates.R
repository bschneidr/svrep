suppressPackageStartupMessages(library(survey))

# Create example data ----
set.seed(2023)

# Create an example survey design object

sample_data <- data.frame(
  Y       = c(8.11, 7.42, 14.54, 18.21, 9.18, 13.83, 7.57, 9.87),
  STRATUM = c(1,1,1,1,2,2,2,2),
  PSU     = c(1,2,3,4,5,6,7,8),
  FPC     = rep(c(0.5, 0.3), each = 4)
)

survey_design <- svydesign(
  data    = sample_data,
  strata  = ~ STRATUM,
  ids     = ~ PSU,
  weights = ~ 1,
  fpc     = ~ FPC
)

# Test subsampling of replicates

test_that("Correct number of replicates and `scale` attribute after subsampling", {
  jkn_design <- survey_design |> as.svrepdesign(type = "JKn", compress = TRUE)
  subsamp_design <- jkn_design |> subsample_replicates(n_reps = 5)

  expect_equal(
    object = ncol(subsamp_design$repweights[['weights']]),
    expected = 5
  )
  expect_equal(
    object = subsamp_design$scale,
    expected = (8/5) * jkn_design$scale
  )

  jkn_design <- survey_design |> as.svrepdesign(type = "JKn", compress = FALSE)
  subsamp_design <- jkn_design |> subsample_replicates(n_reps = 5)

  expect_equal(
    object = ncol(subsamp_design$repweights),
    expected = 5
  )
  expect_equal(
    object = subsamp_design$scale,
    expected = (8/5) * jkn_design$scale
  )
})
