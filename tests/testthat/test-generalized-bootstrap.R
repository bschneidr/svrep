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

# Same results from conversion vs. creating from scratch ----

  set.seed(2014)

  data('election', package = 'survey')

  test_that(
    "Horvitz-Thompson: Same results from conversion vs. creating from scratch", {

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

# Sanity check results ----

  test_that(
    "Bootstrap estimate with many replicates close to expected value", {

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

