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

  current_sim_cv <- estimate_boot_sim_cv(estimated_means_and_proportions) |>
    summary()

  test_that(
    "Consistent calculations from different functions", {
      expect_equal(
        object = estimate_boot_reps_for_target_cv(
          estimated_means_and_proportions,
          target_cv = c(current_sim_cv$SIMULATION_CV[1],
                        current_sim_cv$SIMULATION_CV[3])
        ) |> 
          summary() |>
          (\(df) df[c(1,2),c('api00', 'stypeE')])() |>
          as.matrix() |> diag(),
        expected = rep(current_sim_cv$N_REPLICATES[1], 2)
      )
  })

# Consistent outputs from main function ----

  test_that(
    "Consistent outputs from `estimate_boot_reps_for_target_cv()`", {

      output <- estimate_boot_reps_for_target_cv(
        estimated_means_and_proportions,
        target_cv = c(0.5, 0.1)
      )

      expect_equal(
        object = output[['function']](0.1) |> as.vector(),
        expected = unlist(
          output$summary[2,c('api00', 'api99', 'stypeE', 'stypeH', 'stypeM')]
        ) |> as.vector()
      )
    }
  )

# Plotting works ----

  test_that(
    "Plot method works for `estimate_boot_reps_for_target_cv()`", {

      output <- estimate_boot_reps_for_target_cv(
        estimated_means_and_proportions,
        target_cv = c(0.5, 0.1)
      )

      expect_true(
        inherits(invisible(plot(output)), 'ggplot')
      )
    }
  )

# Sanity check ----

  library(survey)
  data("election", package = "survey")

  ht_quad_form_matrix <- make_quad_form_matrix(variance_estimator = "Horvitz-Thompson",
                                               joint_probs = election_jointprob)
  ##_ Produce variance estimate
  wtd_y <- as.matrix(election_pps$wt * election_pps$Bush)
  expected_value <- as.numeric(t(wtd_y) %*% ht_quad_form_matrix %*% wtd_y)

  set.seed(2014)
  sim_results <- replicate(n = 200, expr = {
    adj_factors <- make_gen_boot_factors(Sigma = ht_quad_form_matrix,
                                         num_replicates = 50,
                                         tau = "auto")
    election_pps_bootstrap_design <- svrepdesign(
      data = election_pps,
      weights = 1 / diag(election_jointprob),
      repweights = adj_factors,
      combined.weights = FALSE,
      type = "other",
      scale = attr(adj_factors, 'scale'),
      rscales = attr(adj_factors, 'rscales')
    )

    estimate <- svytotal(x = ~ Bush, design = election_pps_bootstrap_design,
                         return.replicates = TRUE)

    sim_cv_estimate <- estimate_boot_sim_cv(estimate)[['summary']][['SIMULATION_CV']]
    var_estimate <- as.numeric(vcov(estimate))
    c('sim_cv' = sim_cv_estimate, 'var_estimate' = var_estimate)
  })

  empirical_sim_cv <- sd(sim_results['var_estimate',]) / expected_value
  mean_estimated_sim_cv <- mean(sim_results['sim_cv',])

  rel_abs_error <- abs(mean_estimated_sim_cv - empirical_sim_cv)/empirical_sim_cv

  test_that(
    "`estimate_boot_sim_cv() produces correct results", {
      expect_lt(
        object = rel_abs_error,
        expected = 0.05
      )
    }
  )
