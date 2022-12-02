suppressWarnings({
  suppressPackageStartupMessages({
    library(survey)
    library(dplyr)
    library(svrep)
    library(testthat)
  })
})

# Create an example bootstrap survey design object ----
 library(survey)
 data('api', package = 'survey')

 boot_design <- svydesign(id=~1,strata=~stype, weights=~pw,
                          data=apistrat, fpc=~fpc) |>
   svrep::as_bootstrap_design(replicates = 100)

# Calculate estimates of interest and retain estimates from each replicate ----

 estimated_means_and_proportions <- svymean(x = ~ api00 + api99 + stype, design = boot_design,
                                            return.replicates = TRUE)
 custom_statistic <- withReplicates(design = boot_design,
                                    return.replicates = TRUE,
                                    theta = function(wts, data) {
                                      numerator <- sum(data$api00 * wts)
                                      denominator <- sum(data$api99 * wts)
                                      statistic <- numerator/denominator
                                      return(statistic)
                                    })

# Consistent results between functions ----

  current_sim_cv <- estimate_boot_sim_cv(estimated_means_and_proportions)

  test_that(
    "Consistent calculations from different functions", {
      expect_equal(
        object = estimate_boot_reps_for_target_cv(
          estimated_means_and_proportions,
          target_cv = c(current_sim_cv$SIMULATION_CV[1],
                        current_sim_cv$SIMULATION_CV[3])
        )[c(1,2),c('api00', 'stypeE')] |> as.matrix() |> diag(),
        expected = rep(current_sim_cv$N_REPLICATES[1], 2)
      )
  })


