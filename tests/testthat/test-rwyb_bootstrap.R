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

  ## Library two-stage survey (first stage PPSWOR, second stage SRSWOR)
  library_multistage_survey <- svydesign(
    data = svrep::library_multistage_sample,
    fpc = ~ PSU_SAMPLING_PROB + SSU_SAMPLING_PROB,
    ids = ~ PSU_ID + SSU_ID,
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
    replicates = 10
  )
  set.seed(1999)
  dstrat_reps <- make_rwyb_bootstrap_weights(
    num_replicates = 10,
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
    replicates = 10
  )
  set.seed(1999)
  dclus1_reps <- make_rwyb_bootstrap_weights(
    num_replicates = 10,
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
    replicates = 10
  )
  set.seed(1999)
  dclus2_reps <- make_rwyb_bootstrap_weights(
    num_replicates = 10,
    samp_unit_ids = dclus2$cluster,
    strata_ids = dclus2$strata,
    samp_unit_sel_probs = dclus2$allprob,
    samp_method_by_stage = c("SRSWOR", "SRSWOR"),
    output = 'factors'
  ) |> `dimnames<-`(NULL)

  set.seed(1999)
  library_multistage_boot <- as_bootstrap_design(
    design = library_multistage_survey,
    type = "Rao-Wu-Yue-Beaumont",
    replicates = 10
  )
  set.seed(1999)
  library_multistage_reps <- make_rwyb_bootstrap_weights(
    num_replicates = 10,
    samp_unit_ids = library_multistage_survey$cluster,
    strata_ids = library_multistage_survey$strata,
    samp_unit_sel_probs = library_multistage_survey$allprob,
    samp_method_by_stage = c("PPSWOR", "PPSWOR"),
    output = 'factors'
  ) |> `dimnames<-`(NULL)

  set.seed(1999)
  election_pps_boot <- as_bootstrap_design(
    design = pps_design_ht,
    type = "Rao-Wu-Yue-Beaumont",
    replicates = 10
  )
  set.seed(1999)
  election_pps_reps <- make_rwyb_bootstrap_weights(
    num_replicates = 10,
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
    "Same results from conversion versus creating from scratch: library_multistage", {
      expect_equal(
        object = library_multistage_boot$repweights,
        expected = library_multistage_reps
      )
    })

  test_that(
    "Same results from conversion versus creating from scratch: election_pps", {
      expect_equal(
        object = election_pps_boot$repweights,
        expected = election_pps_reps
      )
    })
