suppressWarnings({
  suppressPackageStartupMessages({
    library(survey)
    library(dplyr)
    library(svrep)
    library(testthat)
  })
})

set.seed(2014)

# Load example datasets ----

  data('api', package = 'survey')

  ## One-stage stratified sample (SRSWOR)
  dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)

  ## One-stage cluster sample (SRSWOR)
  dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

  ## Two-stage cluster sample (SRSWOR)
  dclus2 <- svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2)

  ## Two-stage cluster sample (SRSWR)
  dclus2wr <- svydesign(id=~dnum+snum, weights = ~pw, data=apiclus2)

  ## Library two-stage survey (first stage PPSWOR, second stage SRSWOR)
  library_multistage_survey <- svydesign(
    data = svrep::library_multistage_sample,
    fpc = ~ PSU_SAMPLING_PROB + SSU_SAMPLING_PROB,
    ids = ~ PSU_ID + SSU_ID,
    pps = "brewer"
  )

  ## Library two-stage survey with nonresponse
  multistage_survey_with_poisson_nonresponse <- svydesign(
    data = svrep::library_multistage_sample |>
      mutate(
        RESP_PROB = runif(n = nrow(svrep::library_multistage_sample),
                          min = 0.6, max = 1),
        RESP_STATUS = ifelse(
          runif(n = nrow(svrep::library_multistage_sample)) < RESP_PROB,
          1, 0)
      ),
    fpc = ~ PSU_SAMPLING_PROB + SSU_SAMPLING_PROB + RESP_PROB,
    ids = ~ PSU_ID + SSU_ID + RESP_STATUS,
    pps = "brewer"
  )

  ## One-stage PPSWOR
  data('election', package = 'survey')

  pps_design_ht <- svydesign(
    data = election_pps,
    id = ~1, fpc = ~p,
    pps = ppsmat(election_jointprob),
    variance = "HT"
  )

# Convert to RWYB bootstrap design ----

  set.seed(1999)
  dstrat_boot <- as_bootstrap_design(
    design = dstrat,
    type = "Rao-Wu-Yue-Beaumont",
    replicates = 100
  )
  set.seed(1999)
  dstrat_reps <- make_rwyb_bootstrap_weights(
    num_replicates = 100,
    samp_unit_ids = dstrat$cluster,
    strata_ids = dstrat$strata,
    samp_unit_sel_probs = dstrat$allprob,
    samp_method_by_stage = c("SRSWOR"),
    output = 'factors'
  ) |> `dimnames<-`(NULL)

  set.seed(1999)
  dclus1_boot <- as_bootstrap_design(
    design = dclus1,
    type = "Rao-Wu-Yue-Beaumont",
    replicates = 100
  )
  set.seed(1999)
  dclus1_reps <- make_rwyb_bootstrap_weights(
    num_replicates = 100,
    samp_unit_ids = dclus1$cluster,
    strata_ids = dclus1$strata,
    samp_unit_sel_probs = dclus1$allprob,
    samp_method_by_stage = c("SRSWOR"),
    output = 'factors'
  ) |> `dimnames<-`(NULL)

  set.seed(1999)
  dclus2_boot <- as_bootstrap_design(
    design = dclus2,
    type = "Rao-Wu-Yue-Beaumont",
    replicates = 100
  )
  set.seed(1999)
  dclus2_reps <- make_rwyb_bootstrap_weights(
    num_replicates = 100,
    samp_unit_ids = dclus2$cluster,
    strata_ids = dclus2$strata,
    samp_unit_sel_probs = dclus2$allprob,
    samp_method_by_stage = c("SRSWOR", "SRSWOR"),
    output = 'factors'
  ) |> `dimnames<-`(NULL)

  set.seed(1999)
  dclus2wr_boot <- as_bootstrap_design(
    design = dclus2wr,
    type = "Rao-Wu-Yue-Beaumont",
    replicates = 100
  )
  set.seed(1999)
  dclus2wr_reps <- make_rwyb_bootstrap_weights(
    num_replicates = 100,
    samp_unit_ids = dclus2$cluster,
    strata_ids = dclus2$strata,
    samp_unit_sel_probs = matrix(0,
                                 nrow = nrow(dclus2wr$fpc$sampsize),
                                 ncol = ncol(dclus2wr$fpc$sampsize)),
    samp_method_by_stage = c("SRSWR", "SRSWR"),
    output = 'factors'
  ) |> `dimnames<-`(NULL)

  set.seed(1999)
  library_multistage_boot <- as_bootstrap_design(
    design = library_multistage_survey,
    type = "Rao-Wu-Yue-Beaumont",
    replicates = 100,
    samp_method_by_stage = c("PPSWOR", "SRSWOR")
  )
  set.seed(1999)
  library_multistage_reps <- make_rwyb_bootstrap_weights(
    num_replicates = 100,
    samp_unit_ids = library_multistage_survey$cluster,
    strata_ids = library_multistage_survey$strata,
    samp_unit_sel_probs = library_multistage_survey$allprob,
    samp_method_by_stage = c("PPSWOR", "SRSWOR"),
    output = 'factors'
  ) |> `dimnames<-`(NULL)

  set.seed(1999)
  ms_w_poisson_nr_boot <- as_bootstrap_design(
    design = multistage_survey_with_poisson_nonresponse,
    type = "Rao-Wu-Yue-Beaumont",
    replicates = 100,
    samp_method_by_stage = c("PPSWOR", "SRSWOR", "Poisson")
  )
  set.seed(1999)
  ms_w_poisson_nr_reps <- make_rwyb_bootstrap_weights(
    num_replicates = 100,
    samp_unit_ids = multistage_survey_with_poisson_nonresponse$cluster,
    strata_ids = multistage_survey_with_poisson_nonresponse$strata,
    samp_unit_sel_probs = multistage_survey_with_poisson_nonresponse$allprob,
    samp_method_by_stage = c("PPSWOR", "SRSWOR", "Poisson"),
    output = 'factors'
  ) |> `dimnames<-`(NULL)

  set.seed(1999)
  election_pps_boot <- as_bootstrap_design(
    design = pps_design_ht,
    type = "Rao-Wu-Yue-Beaumont",
    replicates = 100
  )
  set.seed(1999)
  election_pps_reps <- make_rwyb_bootstrap_weights(
    num_replicates = 100,
    samp_unit_ids = pps_design_ht$cluster,
    strata_ids = pps_design_ht$strata,
    samp_unit_sel_probs = pps_design_ht$allprob,
    samp_method_by_stage = c("PPSWOR"),
    output = 'factors'
  ) |> `dimnames<-`(NULL)

  dclus1_boot |> summarize_rep_weights(type = "overall")
  dclus2_boot |> summarize_rep_weights(type = "overall")
  library_multistage_boot |> summarize_rep_weights(type = "overall")
  election_pps_boot |> summarize_rep_weights(type = "overall")

# Same results from conversion and creating from scratch ----

  test_that(
    "Same results from conversion versus creating from scratch: dstrat", {
      expect_equal(
        object = dstrat_boot$repweights,
        expected = dstrat_reps
      )
    })

  test_that(
    "Same results from conversion versus creating from scratch: dclus1", {
      expect_equal(
        object = dclus1_boot$repweights,
        expected = dclus1_reps
      )
    })

  test_that(
    "Same results from conversion versus creating from scratch: dclus2", {
      expect_equal(
        object = dclus2_boot$repweights,
        expected = dclus2_reps
      )
    })

  test_that(
    "Same results from conversion versus creating from scratch: dclus2wr", {
      expect_equal(
        object = dclus2wr_boot$repweights,
        expected = dclus2wr_reps
      )
    })

  test_that(
    "Same results from conversion versus creating from scratch: library_multistage", {
      expect_equal(
        object = library_multistage_boot$repweights,
        expected = library_multistage_reps
      )
    })

  test_that(
    "Same results from conversion versus creating from scratch: multistage_survey_with_poisson_nonresponse", {
      expect_equal(
        object = ms_w_poisson_nr_boot$repweights,
        expected = ms_w_poisson_nr_reps
      )
    })

  test_that(
    "Same results from conversion versus creating from scratch: election_pps", {
      expect_equal(
        object = election_pps_boot$repweights,
        expected = election_pps_reps
      )
    })

# Outputs either factors or weights ----

  test_that(
    "Able to correctly output either adjustment factors or weights", {
      expect_equal(
        object = withr::with_seed(2014, {
          as_bootstrap_design(dclus1, type = "Rao-Wu-Yue-Beaumont") |>
            weights(type = "analysis")
        }),
        expected = withr::with_seed(2014, {
          as_bootstrap_design(dclus1, type = "Rao-Wu-Yue-Beaumont") |>
            weights(type = "replication") * weights(dclus1)
        })
      )
      expect_equal(
        object = withr::with_seed(2014, {
          as_bootstrap_design(dclus1, type = "Rao-Wu-Yue-Beaumont",
                              replicates = 100) |>
            weights(type = "analysis")
        }),
        expected = withr::with_seed(2014, {
          make_rwyb_bootstrap_weights(
            num_replicates = 100,
            samp_unit_ids = dclus1$cluster,
            strata_ids = dclus1$strata,
            samp_unit_sel_probs = dclus1$allprob,
            samp_method_by_stage = c("SRSWOR"),
            output = 'weights'
          ) |> `dimnames<-`(NULL)
        })
      )
    })

# Expected warnings and errors ----

  data('mu284', package = 'survey')
  test_that(
    desc = "Expected error for noncertainty singleton strata", {
      expect_error(
        object =   svydesign(
          data = mu284 |>
            mutate(Stratum = case_when(
              id1 == "19" ~ "1",
              TRUE ~ "2"
            )),
          ids = ~ id1 + id2,
          strata = ~ Stratum,
          fpc = ~ n1 + n2
        ) |> as_bootstrap_design(type = "Rao-Wu-Yue-Beaumont"),
        regexp = "Cannot form bootstrap adjustment factors for a stratum at stage 1"
      )
  })

# Check for correctness ----

  ##_ SRSWOR ----
  test_that(
    "For SRSWOR, sum of weights in replicates always matches the population size", {
      set.seed(2014)
      expect_equal(
        object = svrep::library_census |>
          group_by(STABR) |>
          mutate(N_PSUS = n()) |>
          arrange(sample(x = n())) |>
          slice_head(n = 2) |>
          ungroup() %>%
          svydesign(data = ., ids = ~ 1, strata = ~ STABR, fpc = ~ N_PSUS) |>
          as_bootstrap_design(replicates = 100) |>
          summarize_rep_weights(type = "overall") |>
          select(avg_wgt_sum, sd_wgt_sums) |>
          as.numeric(),
        expected = c(9245, 0)
      )
    })

  ##_ PPSWOR ----

  test_that(
    "For PPSWOR+SRSWOR two-stage design, RWYB bootstrap has correct weights sums, and estimate close to Brewer's approximation", {
      set.seed(2014)

      # Declare a multistage design
      # where first-stage probabilities are PPSWOR sampling
      # and second-stage probabilities are based on SRSWOR
      multistage_design <- svydesign(
        data = library_multistage_sample,
        ids = ~ PSU_ID + SSU_ID,
        probs = ~ PSU_SAMPLING_PROB + SSU_SAMPLING_PROB,
        pps = "brewer"
      )

      # Convert to a bootstrap replicate design
      boot_design <- as_bootstrap_design(
        design = multistage_design,
        type = "Rao-Wu-Yue-Beaumont",
        samp_method_by_stage = c("PPSWOR", "SRSWOR"),
        replicates = 500
      )

      # Compare variance estimates
      brewer_estimate <- svytotal(x = ~ TOTCIR, na.rm = TRUE, design = multistage_design) |>
        vcov() |> as.numeric()
      boot_estimate <- svytotal(x = ~ TOTCIR, na.rm = TRUE, design = boot_design) |>
        vcov() |> as.numeric()

      expect_lt(
        object = abs(boot_estimate - brewer_estimate)/brewer_estimate,
        expected = 0.1
      )
    })




  ##_ Poisson ----

  test_that(
    "For single-stage Poisson design, RWYB bootstrap close to Horvitz-Thompson estimate", {
      set.seed(2014)

      ### Represent Poisson design as PPS design with Horvitz-Thompson estimator
      poisson_design <- multistage_survey_with_poisson_nonresponse$variables |>
        filter(RESP_STATUS == 1) %>%
        svydesign(data = .,
                  probs = ~ RESP_PROB,
                  ids = ~ 1,
                  pps = ppsmat(
                    (.$RESP_PROB %*% t(.$RESP_PROB)) |>
                      `diag<-`(.$RESP_PROB)
                  ),
                  variance = "HT")

      ### Convert to RWYB bootstrap design
      poisson_boot <- poisson_design |>
        as_bootstrap_design(samp_method_by_stage = 'Poisson',
                            replicates = 500)

      ### Calculate Horvitz-Thompson estimate
      poisson_quad_matrix <- make_quad_form_matrix(
        variance_estimator = "Horvitz-Thompson",
        joint_probs = (poisson_design$prob %*% t(poisson_design$prob)) |>
          `diag<-`(poisson_design$prob)
      )
      wtd_y <- poisson_design$variables$TOTCIR / poisson_design$prob
      wtd_y[is.na(wtd_y)] <- 0

      expected_value <- as.numeric(t(wtd_y) %*% poisson_quad_matrix %*% wtd_y)

      ### Compare bootstrap estimate to Horvitz-Thompson estimate
      rwyb_boot_est <- svytotal(x = ~ TOTCIR, design = poisson_boot, na.rm = TRUE) |>
        vcov() |> as.numeric()

      expect_lt(
        object = abs(rwyb_boot_est - expected_value)/expected_value,
        expected = 0.10
      )
    })
