suppressWarnings({
  suppressPackageStartupMessages({
    library(survey)
    library(dplyr)
    library(svrep)
    library(testthat)
  })
})

data('library_stsys_sample', package = 'svrep')
set.seed(2014)

library_stsys_sample <- library_stsys_sample |>
  mutate(
    TOTCIR = ifelse(is.na(TOTCIR), 0, TOTCIR),
    TOTSTAFF = ifelse(is.na(TOTSTAFF), 0, TOTSTAFF)
  )

data('mu284', package = 'survey')

mu284$phase2 <- mu284$id1 %in% c(19,31,45)
mu284 <- mu284[order(mu284$id1),]

twophase_design <- twophase(
  data = mu284,
  id = list(~id1, ~id1),
  fpc = list(~n1, NULL),
  subset = ~ phase2
)

# Same results from conversion vs. creating from scratch ----

  set.seed(2014)

  data('election', package = 'survey')

  test_that(
    "Horvitz-Thompson: Same results from conversion vs. creating from scratch", {

      # Disable test until the Matrix package adapts to change in R 4.3.2
      # (`crossprod()` becomes primitive and S3 generic in R 4.3.2)
      skip_if_not("crossprod" %in% ls(getNamespaceInfo("Matrix", "exports")))

      ## Create survey design object
      pps_design_ht <- svydesign(
        data = election_pps,
        id = ~1, fpc = ~p,
        pps = ppsmat(election_jointprob),
        variance = "HT"
      )

      ## Convert to generalized bootstrap replicate design
      set.seed(2014)
      conversion_result <- pps_design_ht |>
        as_gen_boot_design(variance_estimator = "Horvitz-Thompson",
                           replicates = 5, tau = "auto") |>
        weights(type = "replication")


      set.seed(2014)
      orig_result <- make_quad_form_matrix(variance_estimator = "Horvitz-Thompson",
                            joint_probs = election_jointprob) |>
        make_gen_boot_factors(num_replicates = 5, tau = "auto")

    expect_equal(
      object =  conversion_result,
      expected = orig_result
    )

  })

  test_that(
    "Yates-Grundy: Same results from conversion vs. creating from scratch", {

      # Disable test until the Matrix package adapts to change in R 4.3.2
      # (`crossprod()` becomes primitive and S3 generic in R 4.3.2)
      skip_if_not("crossprod" %in% ls(getNamespaceInfo("Matrix", "exports")))

      ## Create survey design object
      pps_design_yg <- svydesign(
        data = election_pps,
        id = ~1, fpc = ~p,
        pps = ppsmat(election_jointprob),
        variance = "YG"
      )

      ## Convert to generalized bootstrap replicate design
      set.seed(2014)
      conversion_result <- pps_design_yg |>
        as_gen_boot_design(variance_estimator = "Yates-Grundy",
                           replicates = 5, tau = "auto") |>
        weights(type = "replication")


      set.seed(2014)
      orig_result <- make_quad_form_matrix(variance_estimator = "Yates-Grundy",
                                           joint_probs = election_jointprob) |>
        make_gen_boot_factors(num_replicates = 5, tau = "auto")

      expect_equal(
        object = conversion_result,
        expected = orig_result
      )

    })

  test_that(
    "Ultimate Cluster: Same results from conversion vs. creating from scratch", {

      ## Create survey design object
      multistage_survey_design <- svydesign(
        data = svrep::library_multistage_sample |>
          mutate(WT = 1/SAMPLING_PROB),
        weights = ~ WT,
        ids = ~ PSU_ID + SSU_ID,
        fpc = ~ PSU_POP_SIZE + SSU_POP_SIZE
      )

      ## Convert to generalized bootstrap replicate design
      set.seed(2014)
      conversion_result <- multistage_survey_design |>
        as_gen_boot_design(variance_estimator = "Ultimate Cluster",
                           replicates = 5, tau = "auto") |>
        weights(type = "replication")


      set.seed(2014)
      orig_result <- make_quad_form_matrix(variance_estimator = "Ultimate Cluster",
                                           strata_ids = multistage_survey_design$strata,
                                           cluster_ids = multistage_survey_design$cluster,
                                           strata_pop_sizes = multistage_survey_design$fpc$popsize) |>
        make_gen_boot_factors(num_replicates = 5, tau = "auto")

      expect_equal(
        object = conversion_result,
        expected = orig_result
      )

    })

  test_that(
    "Multistage Stratified SRS: Same results from conversion vs. creating from scratch", {

      ## Create survey design object
      multistage_survey_design <- svydesign(
        data = svrep::library_multistage_sample |>
          mutate(WT = 1/SAMPLING_PROB),
        weights = ~ WT,
        ids = ~ PSU_ID + SSU_ID,
        fpc = ~ PSU_POP_SIZE + SSU_POP_SIZE
      )

      ## Convert to generalized bootstrap replicate design
      set.seed(2014)
      conversion_result <- multistage_survey_design |>
        as_gen_boot_design(variance_estimator = "Stratified Multistage SRS",
                           replicates = 5, tau = "auto") |>
        weights(type = "replication")


      set.seed(2014)
      orig_result <- make_quad_form_matrix(variance_estimator = "Stratified Multistage SRS",
                                           strata_ids = multistage_survey_design$strata,
                                           cluster_ids = multistage_survey_design$cluster,
                                           strata_pop_sizes = multistage_survey_design$fpc$popsize) |>
        make_gen_boot_factors(num_replicates = 5, tau = "auto")

      expect_equal(
        object = conversion_result,
        expected = orig_result
      )

    })

  test_that(
    "SD1: Same results from conversion vs. creating from scratch", {

      ## Create survey design object
      multistage_survey_design <- svydesign(
        data = svrep::library_multistage_sample |>
          mutate(WT = 1/SAMPLING_PROB),
        weights = ~ WT,
        ids = ~ PSU_ID + SSU_ID,
        fpc = ~ PSU_POP_SIZE + SSU_POP_SIZE
      )

      suppressMessages({
        ## Convert to generalized bootstrap replicate design
        set.seed(2014)
        conversion_result <- multistage_survey_design |>
          as_gen_boot_design(variance_estimator = "SD1",
                             replicates = 5, tau = "auto") |>
          weights(type = "replication")


        set.seed(2014)
        orig_result <- make_quad_form_matrix(variance_estimator = "SD1",
                                             strata_ids = multistage_survey_design$strata,
                                             cluster_ids = multistage_survey_design$cluster,
                                             strata_pop_sizes = multistage_survey_design$fpc$popsize,
                                             sort_order = seq_len(nrow(multistage_survey_design))) |>
          make_gen_boot_factors(num_replicates = 5, tau = "auto")
      })

      expect_equal(
        object = conversion_result,
        expected = orig_result
      )

    })

  test_that(
    "SD2: Same results from conversion vs. creating from scratch", {

      ## Create survey design object
      multistage_survey_design <- svydesign(
        data = svrep::library_multistage_sample |>
          mutate(WT = 1/SAMPLING_PROB),
        weights = ~ WT,
        ids = ~ PSU_ID + SSU_ID,
        fpc = ~ PSU_POP_SIZE + SSU_POP_SIZE
      )

      suppressMessages({
        ## Convert to generalized bootstrap replicate design
        set.seed(2014)
        conversion_result <- multistage_survey_design |>
          as_gen_boot_design(variance_estimator = "SD2",
                             replicates = 5, tau = "auto") |>
          weights(type = "replication")


        set.seed(2014)
        orig_result <- make_quad_form_matrix(variance_estimator = "SD2",
                                             strata_ids = multistage_survey_design$strata,
                                             cluster_ids = multistage_survey_design$cluster,
                                             strata_pop_sizes = multistage_survey_design$fpc$popsize,
                                             sort_order = seq_len(nrow(multistage_survey_design))) |>
          make_gen_boot_factors(num_replicates = 5, tau = "auto")
      })

      expect_equal(
        object = conversion_result,
        expected = orig_result
      )

    })

  test_that(
    "Two-phase Design: Same results from conversion vs. creating from scratch", {

      suppressMessages({
        expect_warning(
          regexp = "The sample quadratic form matrix for this design and variance estimator is not positive semidefinite.", object = {
            set.seed(2023)
            twophase_gen_boot <- twophase_design |>
              as_gen_boot_design(
                replicates = 5, tau = 1,
                variance_estimator = list('SD2', 'Ultimate Cluster'),
                psd_option = "warn"
              )
          }
        )
      })

      suppressMessages({
        suppressWarnings({
          set.seed(2023)
          gen_boot_reps <- twophase_design |>
            get_design_quad_form(
              variance_estimator = list('SD2', 'Ultimate Cluster')
            ) |>
            get_nearest_psd_matrix() |>
            make_gen_boot_factors(
              num_replicates = 5, tau = 1
            )

          expect_equal(
            object = twophase_gen_boot$repweights,
            expected = gen_boot_reps
          )
        })
      })

  })

# Rescaling functions work as expected ----

  test_that(
    "Rescaling functions work as expected", {

      suppressMessages({
        suppressWarnings({
          twophase_gen_boot <- twophase_design |>
            as_gen_boot_design(
              variance_estimator = list(
                'SD2', 'Ultimate Cluster'
              ),
              replicates = 5, tau = 1
            )
        })
      })

      rescaled_design <- twophase_gen_boot |>
        rescale_replicates(min_wgt = 0.05)
      rescaled_matrix <- twophase_gen_boot |>
        weights(type = "replication") |>
        `attr<-`('scale', twophase_gen_boot$scale) |>
          rescale_replicates(min_wgt = 0.05)

    expect_equal(
      object = rescaled_design |> weights(type = "replication"),
      expected = rescaled_matrix
    )

    expect_gte(
      object = min(rescaled_matrix), expected = 0.05
    )

    expect_equal(
      object = matrix(c(1,0.3,0.02,2,3,4), ncol = 3) |>
        `attr<-`('scale', 1) |>
        rescale_replicates(min_wgt = 0.02) |>
        `attr<-`('tau', NULL),
      expected = matrix(c(1,0.3,0.02,2,3,4), ncol = 3)
    )

  })

# Testing the `exact_vcov` option ----

  test_that(
    desc = "`exact_vcov = TRUE` gives exact variance estimate for totals", {

      # Disable test until the Matrix package adapts to change in R 4.3.2
      # (`crossprod()` becomes primitive and S3 generic in R 4.3.2)
      skip_if_not("crossprod" %in% ls(getNamespaceInfo("Matrix", "exports")))

      pps_design_yg <- svydesign(
        data = election_pps,
        id = ~1, fpc = ~p,
        pps = ppsmat(election_jointprob),
        variance = "YG"
      )

      gen_boot_est <- as_gen_boot_design(
        pps_design_yg, "Yates-Grundy",
        replicates = 40, exact_vcov = TRUE
      ) |>
        svytotal(x = ~ Bush + Kerry) |> vcov() |>
        `attr<-`('means', NULL)

      exact_est <- pps_design_yg |>
        svytotal(x = ~ Bush + Kerry) |> vcov() |>
        `attr<-`('means', NULL)

      expect_equal(object = gen_boot_est, expected = exact_est)
    }
  )

  test_that(
    desc = "`exact_vcov = TRUE` throws informative error message if there too few replicates", {

      # Disable test until the Matrix package adapts to change in R 4.3.2
      # (`crossprod()` becomes primitive and S3 generic in R 4.3.2)
      skip_if_not("crossprod" %in% ls(getNamespaceInfo("Matrix", "exports")))

      expect_error(
        regexp = "only works if.+39",
        object = {
          svydesign(
            data = election_pps,
            id = ~1, fpc = ~p,
            pps = ppsmat(election_jointprob),
            variance = "YG"
          ) |> as_gen_boot_design(
            variance_estimator = "Yates-Grundy",
            exact_vcov = TRUE,
            replicates = 39
          )
        }
      )
    }
  )


# Sanity check results ----

  test_that(
    "Bootstrap estimate with many replicates close to expected value", {

      # Disable test until the Matrix package adapts to change in R 4.3.2
      # (`crossprod()` becomes primitive and S3 generic in R 4.3.2)
      skip_if_not("crossprod" %in% ls(getNamespaceInfo("Matrix", "exports")))

      ## Create survey design object
      pps_design_yg <- svydesign(
        data = election_pps,
        id = ~1, fpc = ~p,
        pps = ppsmat(election_jointprob),
        variance = "YG"
      )

      exp_value <- svytotal(x = ~ Bush, design = pps_design_yg) |> vcov() |>
        as.numeric()

      ## Convert to generalized bootstrap replicate design
      set.seed(2014)
      conversion_result <- pps_design_yg |>
        as_gen_boot_design(variance_estimator = "Yates-Grundy",
                           replicates = 50000, tau = "auto") |>
        svytotal(x = ~ Bush) |> vcov() |> as.numeric()

      ## Check that variance estimate from bootstrap is similar to expected value
      expect_equal(
        object = conversion_result,
        expected = exp_value,
        tolerance = 0.01
      )

      ## Check that the average replicate's estimate is similar to expected value
      expect_gt(
        object =       pps_design_yg |>
          as_gen_boot_design(variance_estimator = "Yates-Grundy",
                             replicates = 50000, tau = "auto") |>
          svytotal(x = ~ Bush, return.replicates = TRUE) |>
          getElement("replicates") |>
          t.test(mu = coef(svytotal(x = ~ Bush, design = pps_design_yg))) |>
          getElement("p.value"),
        expected = 0.01
      )

    })

# Works for more specialized classes of survey designs ----

  test_that(
    desc = "Returns `tbl_svy` if the input is a `tbl_svy` and 'srvyr' is loaded", {
      library(srvyr)
      expect_true(
        twophase_design |> as_survey() |>
          as_gen_boot_design(variance_estimator = list(
            'Ultimate Cluster', 'Ultimate Cluster'
          ), replicates = 1, tau = 1) |>
          inherits('tbl_svy')
      )
    }
  )
