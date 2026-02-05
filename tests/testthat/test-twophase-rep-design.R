suppressWarnings({
  suppressPackageStartupMessages({
    library(survey)
    library(dplyr)
    library(svrep)
    library(testthat)
  })
})

# Create example data
phase_one_data <- data.frame(
  id = c(1, 2, 3, 4, 5, 6, 7, 8),
  x  = c(109, 95, 86, 102, 106, 105, 106, 115),
  y  = c(103, 60, 92, 76, 104, 132, 127, 88),
  PHASE_ONE_PSU     = c(1, 2, 3, 4, 5, 6, 7, 8),
  PHASE_ONE_WGT     = c(95, 105, 99, 101, 85, 115, 101, 99),
  PHASE_TWO_STRATA  = rep(c(1,2), each = 4),
  PHASE_TWO_PROB    = rep(c(0.5, 0.75), each = 4),
  PHASE_TWO_SAMPLED = c(TRUE, TRUE, FALSE, FALSE, 
                        TRUE, TRUE, TRUE, FALSE)
)

# Create replicate weights for the first phase 
phase_one_design <- svydesign(
  data    = phase_one_data,
  ids     = ~ PHASE_ONE_PSU,
  weights = ~ PHASE_ONE_WGT
)

set.seed(2026)

phase_one_rep_design <- phase_one_design |>
  as.svrepdesign(type = "JK1", mse = TRUE)

# Derive a replicate design for the two-phase sample
phase_two_rep_design <- derive_twophase_rep_design(
  design               = phase_one_rep_design,
  phase_two_indicators = "PHASE_TWO_SAMPLED",
  phase_two_probs      = "PHASE_TWO_PROB",
  phase_two_strata     = "PHASE_TWO_STRATA"
)

test_that("Correct weights from `derive_twophase_rep_design()`", {

  # Check exactly correct values for full-sample weights
  full_weights <- phase_two_rep_design |> weights(type = "sampling")

  expected_full_weight <- phase_one_rep_design |> weights(type = "sampling")
  expected_full_weight <- expected_full_weight / phase_one_rep_design$variables$PHASE_TWO_PROB
  expected_full_weight <- expected_full_weight[phase_one_rep_design$variables$PHASE_TWO_SAMPLED]
  expected_full_weight[1:2] <- expected_full_weight[1:2] * (400 / sum(expected_full_weight[1:2]))
  expected_full_weight[3:5] <- expected_full_weight[3:5] * (400 / sum(expected_full_weight[3:5]))

  expect_equal(unname(full_weights), unname(expected_full_weight))

  # Replicate weights should yield identical variance estimates
  # for strata identifiers, compared to phase one sample
  expect_equal(
    svytotal(x = ~ factor(PHASE_TWO_STRATA), design = phase_two_rep_design) |> vcov(),
    svytotal(x = ~ factor(PHASE_TWO_STRATA), design = phase_one_rep_design) |> vcov()
  )
})

test_that("Correct results regardless of `combined.weights`", {

  separate_weights_design <- phase_one_design |>
    as.svrepdesign(type = "JK1", mse = TRUE)
  combined_weights_design <- separate_weights_design

  combined_weights_design$repweights <- separate_weights_design |>
    weights(type = "analysis")
  combined_weights_design$combined.weights <- TRUE

  combined_weights_result <- derive_twophase_rep_design(
    design               = phase_one_rep_design,
    phase_two_indicators = "PHASE_TWO_SAMPLED",
    phase_two_probs      = "PHASE_TWO_PROB",
    phase_two_strata     = "PHASE_TWO_STRATA"
  )
  separate_weights_result <- derive_twophase_rep_design(
    design               = phase_one_rep_design,
    phase_two_indicators = "PHASE_TWO_SAMPLED",
    phase_two_probs      = "PHASE_TWO_PROB",
    phase_two_strata     = "PHASE_TWO_STRATA"
  )

  expect_equal(
    weights(combined_weights_result, type = "analysis"), 
    weights(separate_weights_result, type = "analysis")
  )
})

test_that("Informative errors from `derive_twophase_rep_design()`", {
  expect_error(
    derive_twophase_rep_design(
      design               = phase_one_rep_design,
      phase_two_indicators = "PHASE_TWO_SAMPLED",
      phase_two_probs      = "DERP",
      phase_two_strata     = "PHASE_TWO_STRATA"
    ),
    regexp = "name of a variable"
  )
  expect_error(
    derive_twophase_rep_design(
      design               = phase_one_rep_design |> transform(STRATA = rep(NA, times = 8)),
      phase_two_indicators = "PHASE_TWO_SAMPLED",
      phase_two_probs      = "PHASE_TWO_PROB",
      phase_two_strata     = "STRATA"
    ),
    regexp = "cannot have any missing values"
  )
})
