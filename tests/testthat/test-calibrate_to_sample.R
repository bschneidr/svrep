# library(svrep)
# library(testthat)

# Prepare example data ----
  suppressPackageStartupMessages(library(survey))
  data(api, package = 'survey')
  set.seed(1999)

  primary_survey <- svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc) |>
    as.svrepdesign(type = "JK1") |> srvyr::as_survey()

  control_survey <- svydesign(id = ~ 1, fpc = ~fpc, data = apisrs) |>
    as.svrepdesign(type = "JK1")

test_that("Duplicates primary replicates when control survey has more replicates", {

  col_matching <- c(25L, 55L, 98L, 160L, 199L, 2L, 122L, 27L, 72L, 74L, 111L, 152L,
                    10L, 38L, 126L, 63L, 41L, 68L, 65L, 28L, 83L, NA, 175L, 80L,
                    106L, 43L, 97L, 150L, 61L, 139L, 90L, 129L, 168L, 173L, 67L,
                    42L, 105L, 84L, NA, 51L, 180L, 182L, 66L, 200L, 88L, 110L, 196L,
                    5L, NA, 189L, 136L, 174L, 102L, NA, 85L, 195L, 193L, 142L, 33L,
                    162L, 154L, 147L, 12L, 128L, 19L, 130L, 59L, 104L, 187L, 48L,
                    179L, 114L, 81L, 69L, 158L, 113L, 135L, 23L, 165L, 22L, 125L,
                    118L, 159L, 197L, 149L, NA, 8L, 7L, 34L, 186L, 191L, 123L, NA,
                    192L, 120L, 185L, 121L, 71L, 45L, 49L, 183L, 86L, 116L, 92L,
                    29L, 64L, 198L, 79L, 44L, 6L, 167L, 24L, 26L, 31L, 157L, 144L,
                    188L, 134L, 103L, 32L, 58L, 89L, 14L, 39L, 95L, NA, 91L, 93L,
                    21L, 146L, NA, 99L, 137L, 100L, 166L, 50L, 156L, 117L, 46L, 190L,
                    17L, 109L, 62L, 194L, 171L, 76L, 151L, 54L, 172L, 47L, 13L, 53L,
                    107L, 20L, 115L, 60L, 82L, 161L, 177L, 36L, 15L, 18L, 141L, 170L,
                    52L, 11L, 40L, 169L, 133L, 163L, 131L, 181L, 56L, 57L, 184L,
                    145L, 77L, 124L, 153L, 73L, 127L, 16L, 112L, 178L, 108L, 101L,
                    70L, NA, 148L, 1L, 140L, 138L, 143L, 9L, 75L, 176L, 164L, 87L,
                    3L, 78L, 132L, 119L, 96L, 4L, 94L, 35L, 37L, NA, 30L, 155L)

  # Check for informative message
  expect_message(
    suppressWarnings({
      calibrated_rep_design <- calibrate_to_sample(
        primary_rep_design = primary_survey,
        control_rep_design = control_survey,
        cal_formula = ~ stype + enroll,
        control_col_matches = col_matching
      )
    }),
    regexp = "primary survey has fewer replicates than the control survey",
    label = "Informative message is displayed"
  )

  # Check for correct number of columns in the result
  n_control_reps <- ncol(control_survey$repweights)
  n_primary_reps <- ncol(primary_survey$repweights)
  n_dupes <- ceiling(n_control_reps / n_primary_reps)

  expect_equal(object = ncol(calibrated_rep_design$repweights),
               expected = n_primary_reps * n_dupes,
               label = "Result has number of duplications needed")

})

test_that("Basic example gives correct results", {

  # Check for informative message
  expect_message(
    suppressWarnings({
      calibrated_rep_design <- calibrate_to_sample(
        primary_rep_design = primary_survey,
        control_rep_design = primary_survey,
        cal_formula = ~ stype + enroll,
        control_col_matches = NULL
      )
    }),
    regexp = "Matching.+will be done at random",
    label = "Informative message on random replicate matching is displayed"
  )

  epsilon_to_use <- 1e-7

  suppressMessages(
    suppressWarnings({
      calibrated_rep_design <- calibrate_to_sample(
        primary_rep_design = primary_survey,
        control_rep_design = control_survey,
        cal_formula = ~ stype + enroll,
        epsilon = epsilon_to_use
      )
    })
  )

  col_matching <- calibrated_rep_design$control_column_matches

  # Check that calibration replicates were calibrated to intended control replicates
  calibrated_replicate_estimates <- svytotal(x = ~ stype,
                                             design = calibrated_rep_design,
                                             return.replicates = TRUE) |>
    getElement("replicates")

  # Full-sample calibration targets
  control_estimates <- svytotal(x = ~ stype,
                                design = control_survey,
                                return.replicates = TRUE)
  # Replicate-specific calibration targets
  control_replicate_estimates <- getElement(control_estimates,
                                            'replicates')
  a_r <- rep(sqrt(control_survey$scale / calibrated_rep_design$scale),
             times = length(col_matching))

  calibration_targets <- sapply(seq_along(col_matching), function(i) {
    control_rep <- col_matching[i]
    if (is.na(control_rep)) {
      control_rep_est <- coef(control_estimates)
      return(control_rep_est)
    } else {
      control_rep_est <- control_replicate_estimates[control_rep,]
      coef(control_estimates) + a_r[i] * (control_rep_est - coef(control_estimates))
    }
  }) |> t()

  # Check that relative error of every estimate is below epsilon
  misfit <- abs(calibration_targets - calibrated_replicate_estimates)
  relative_error <- misfit / (1 + abs(calibration_targets))

  expect_lt(
    object = max(relative_error),
    expected = epsilon_to_use, expected.label = sprintf("specified epsilon %s",
                                                        epsilon_to_use),
    label = "Relative error of estimated calibration totals"
  )

  # Check that variance-covariance matrix of control totals is reproduced

  vcov_calibrated <- svytotal(x = ~ stype + enroll,
           design = calibrated_rep_design,
           return.replicates = FALSE) |>
    vcov() |> as.matrix() |> `attr<-`('means', NULL)

  vcov_control <- svytotal(x = ~ stype + enroll,
           design = control_survey,
           return.replicates = FALSE) |>
    vcov() |> as.matrix() |> `attr<-`('means', NULL)

  expect_equal(
    object = vcov_calibrated,
    expected = vcov_control,
    tolerance = 1e-07
  )

})

test_that("Able to manually specify column matching when control survey has more replicates", {

  epsilon_to_use <- 1e-7

  col_matching <- c(25L, 55L, 98L, 160L, 199L, 2L, 122L, 27L, 72L, 74L, 111L, 152L,
                    10L, 38L, 126L, 63L, 41L, 68L, 65L, 28L, 83L, NA, 175L, 80L,
                    106L, 43L, 97L, 150L, 61L, 139L, 90L, 129L, 168L, 173L, 67L,
                    42L, 105L, 84L, NA, 51L, 180L, 182L, 66L, 200L, 88L, 110L, 196L,
                    5L, NA, 189L, 136L, 174L, 102L, NA, 85L, 195L, 193L, 142L, 33L,
                    162L, 154L, 147L, 12L, 128L, 19L, 130L, 59L, 104L, 187L, 48L,
                    179L, 114L, 81L, 69L, 158L, 113L, 135L, 23L, 165L, 22L, 125L,
                    118L, 159L, 197L, 149L, NA, 8L, 7L, 34L, 186L, 191L, 123L, NA,
                    192L, 120L, 185L, 121L, 71L, 45L, 49L, 183L, 86L, 116L, 92L,
                    29L, 64L, 198L, 79L, 44L, 6L, 167L, 24L, 26L, 31L, 157L, 144L,
                    188L, 134L, 103L, 32L, 58L, 89L, 14L, 39L, 95L, NA, 91L, 93L,
                    21L, 146L, NA, 99L, 137L, 100L, 166L, 50L, 156L, 117L, 46L, 190L,
                    17L, 109L, 62L, 194L, 171L, 76L, 151L, 54L, 172L, 47L, 13L, 53L,
                    107L, 20L, 115L, 60L, 82L, 161L, 177L, 36L, 15L, 18L, 141L, 170L,
                    52L, 11L, 40L, 169L, 133L, 163L, 131L, 181L, 56L, 57L, 184L,
                    145L, 77L, 124L, 153L, 73L, 127L, 16L, 112L, 178L, 108L, 101L,
                    70L, NA, 148L, 1L, 140L, 138L, 143L, 9L, 75L, 176L, 164L, 87L,
                    3L, 78L, 132L, 119L, 96L, 4L, 94L, 35L, 37L, NA, 30L, 155L)

  suppressMessages(
    suppressWarnings({
      calibrated_rep_design <- calibrate_to_sample(
        primary_rep_design = primary_survey,
        control_rep_design = control_survey,
        cal_formula = ~ stype + enroll,
        control_col_matches = col_matching,
        epsilon = epsilon_to_use
      )
    })
  )

  # Check that object correctly returns the column matching
  expect_equal(object = calibrated_rep_design$control_column_matches,
               expected = col_matching,
               label = "Object correctly returns specified column matching")

  # Check that calibration replicates were calibrated to intended control replicates
  calibrated_replicate_estimates <- svytotal(x = ~ stype,
                                             design = calibrated_rep_design,
                                             return.replicates = TRUE) |>
    getElement("replicates")

  # Full-sample calibration targets
  control_estimates <- svytotal(x = ~ stype,
                                design = control_survey,
                                return.replicates = TRUE)
  # Replicate-specific calibration targets
  control_replicate_estimates <- getElement(control_estimates,
                                            'replicates')
  a_r <- rep(sqrt(control_survey$scale / calibrated_rep_design$scale),
             times = length(col_matching))

  calibration_targets <- sapply(seq_along(col_matching), function(i) {
    control_rep <- col_matching[i]
    if (is.na(control_rep)) {
      control_rep_est <- coef(control_estimates)
      return(control_rep_est)
    } else {
      control_rep_est <- control_replicate_estimates[control_rep,]
      coef(control_estimates) + a_r[i] * (control_rep_est - coef(control_estimates))
    }
  }) |> t()

  # Check that relative error of every estimate is below epsilon
  misfit <- abs(calibration_targets - calibrated_replicate_estimates)
  relative_error <- misfit / (1 + abs(calibration_targets))

  expect_lt(
    object = max(relative_error),
    expected = epsilon_to_use, expected.label = sprintf("specified epsilon %s",
                                                        epsilon_to_use),
    label = "Relative error of estimated calibration totals"
  )

})

test_that("Throws error if convergence is not achieved", {

  epsilon_to_use <- 1e-30
  max_iterations <- 2

  expect_error(
    object = {
      suppressMessages({
        suppressWarnings({
          calibrated_rep_design <- calibrate_to_sample(
            primary_rep_design = primary_survey,
            control_rep_design = control_survey,
            cal_formula = ~ stype + enroll,
            epsilon = epsilon_to_use,
            maxit = max_iterations
          )
        })
      })
    },
    regexp = "Convergence was not achieved",
    label = "Informative error message for failure to converge"
  )

})

# Check that function works for more specialized classes ----

test_that(
  desc = "Returns `tbl_svy` if any input is a `tbl_svy` and 'srvyr' is loaded", {
    library(srvyr)
    expect_true(
      suppressMessages({
        suppressWarnings({
          calibrate_to_sample(
            primary_rep_design = primary_survey |> as_survey(),
            control_rep_design = control_survey,
            cal_formula = ~ stype
          ) |> inherits('tbl_svy')
        })
      })
    )
  }
)
