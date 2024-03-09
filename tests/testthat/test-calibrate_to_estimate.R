# library(svrep)
# library(testthat)

# Prepare example data ----
  suppressPackageStartupMessages(library(survey))
  data(api, package = 'survey')
  set.seed(1999)

  primary_survey <- svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc) |>
    as.svrepdesign(type = "JK1", mse = FALSE)

  control_survey <- svydesign(id = ~ 1, fpc = ~fpc, data = apisrs) |>
    as.svrepdesign(type = "JK1", mse = FALSE)

test_that("Basic example gives correct results", {

  estimated_controls <- svytotal(x = ~ stype + enroll,
                                 design = control_survey)
  control_point_estimates <- coef(estimated_controls)
  control_vcov_estimate <- vcov(estimated_controls)

  # Check for informative message
  expect_message(
    suppressWarnings({
      calibrated_rep_design <- calibrate_to_estimate(
        rep_design = primary_survey,
        estimate = control_point_estimates,
        vcov_estimate = control_vcov_estimate,
        cal_formula = ~ stype + enroll
      )
    }),
    regexp = "Selection.+will be done at random",
    label = "Informative message on random selection of columns is displayed"
  )
  expect_warning({
    suppressMessages(
      calibrated_rep_design <- calibrate_to_estimate(
        rep_design = primary_survey,
        estimate = control_point_estimates,
        vcov_estimate = control_vcov_estimate,
        cal_formula = ~ stype + enroll
      )
    )
    },
    regexp = "centered around full-sample",
    label = "Informative message on setting `mse` to TRUE"
  )

  epsilon_to_use <- 1e-7
  col_selection <- c(1,5,10,15)

  suppressMessages(
    suppressWarnings({
      calibrated_rep_design <- calibrate_to_estimate(
        rep_design = primary_survey,
        estimate = control_point_estimates,
        vcov_estimate = control_vcov_estimate,
        cal_formula = ~ stype + enroll,
        col_selection = col_selection,
        epsilon = epsilon_to_use
      )
    })
  )

  # Check that calibration replicates were calibrated to intended control replicates
  calibrated_replicate_estimates <- svytotal(x = ~ stype + enroll,
                                             design = calibrated_rep_design,
                                             return.replicates = TRUE) |>
    getElement("replicates")


    ##_ Calculate spectral decomposition
    eigen_decomposition <- eigen(x = control_vcov_estimate,
                                 symmetric = TRUE)

    ##_ Calculate matrix of calibration targets
    v <- sapply(X = seq_along(eigen_decomposition$values),
                FUN = function(k) {
                  truncated_eigenvalue <- ifelse(eigen_decomposition$values[k] < 0,
                                                 0, eigen_decomposition$values[k])
                  sqrt(truncated_eigenvalue) * eigen_decomposition$vectors[,k]
              })

    calibration_targets <- matrix(data = control_point_estimates,
                                  nrow = ncol(primary_survey$repweights),
                                  ncol = length(control_point_estimates),
                                  byrow = TRUE)
    A_primary <- primary_survey$scale
    rscales_primary <- primary_survey$rscales

    for (i in seq_len(length(col_selection))) {
      i_star <- col_selection[i]
      calibration_targets[i_star,] <-  control_point_estimates + (sqrt(1/(A_primary * rscales_primary[i_star])) * v[,i])
    }

    ##_ Check that relative error of every estimate is below epsilon
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

test_that("Throws error if convergence is not achieved", {

  epsilon_to_use <- 1e-30
  max_iterations <- 2

  estimated_controls <- svytotal(x = ~ stype + enroll,
                                 design = control_survey)
  control_point_estimates <- coef(estimated_controls)
  control_vcov_estimate <- vcov(estimated_controls)

  expect_error(
    object = {
      suppressMessages({
        suppressWarnings({
          calibrated_rep_design <- calibrate_to_estimate(
            rep_design = primary_survey,
            estimate = control_point_estimates,
            vcov_estimate = control_vcov_estimate,
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

# Works for more specialized classes of survey design ----

test_that("Returns `tbl_svy` if the input is a `tbl_svy` and 'srvyr' is loaded", {
  estimated_controls <- svytotal(x = ~ stype,
                                 design = control_survey)
  control_point_estimates <- coef(estimated_controls)
  control_vcov_estimate <- vcov(estimated_controls)

  library(srvyr)

  suppressMessages({
    suppressWarnings({
      calibrated_rep_design <- calibrate_to_estimate(
        rep_design = primary_survey |> as_survey(),
        estimate = control_point_estimates,
        vcov_estimate = control_vcov_estimate,
        cal_formula = ~ stype
      )
    })
  })

  expect_true({
    calibrated_rep_design |> inherits(what = 'tbl_svy')
  })
})
