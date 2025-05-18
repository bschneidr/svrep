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

  ## One-stage PPSWOR
  data('election', package = 'survey')

  pps_design_ht <- svydesign(
    data = election_pps,
    id = ~1, fpc = ~p,
    pps = ppsmat(election_jointprob),
    variance = "HT"
  )

# Convert to doubled half bootstrap design ----

  set.seed(1999)
  dstrat_boot <- as_bootstrap_design(
    design = dstrat,
    type = "Antal-Tille",
    replicates = 100
  )
  set.seed(1999)
  dstrat_reps <- make_doubled_half_bootstrap_weights(
    num_replicates      = 100,
    samp_unit_ids       = dstrat$cluster[,1],
    strata_ids          = dstrat$strata[,1],
    samp_unit_sel_probs = dstrat$allprob[,1],
    output = 'factors'
  ) |> `dimnames<-`(NULL)

  set.seed(1999)
  dclus1_boot <- as_bootstrap_design(
    design = dclus1,
    type = "Antal-Tille",
    replicates = 100
  )
  set.seed(1999)
  dclus1_reps <- make_doubled_half_bootstrap_weights(
    num_replicates      = 100,
    samp_unit_ids       = dclus1$cluster[,1],
    strata_ids          = dclus1$strata[,1],
    samp_unit_sel_probs = dclus1$allprob[,1],
    output = 'factors'
  ) |> `dimnames<-`(NULL)

  set.seed(1999)
  dclus2_boot <- as_bootstrap_design(
    design = dclus2,
    type = "Antal-Tille",
    replicates = 100
  )
  set.seed(1999)
  dclus2_reps <- make_doubled_half_bootstrap_weights(
    num_replicates      = 100,
    samp_unit_ids       = dclus2$cluster[,1],
    strata_ids          = dclus2$strata[,1],
    samp_unit_sel_probs = dclus2$allprob[,1],
    output = 'factors'
  ) |> `dimnames<-`(NULL)

  set.seed(1999)
  dclus2wr_boot <- as_bootstrap_design(
    design = dclus2wr,
    type = "Antal-Tille",
    replicates = 100
  )
  set.seed(1999)
  dclus2wr_reps <- make_doubled_half_bootstrap_weights(
    num_replicates      = 100,
    samp_unit_ids       = dclus2wr$cluster[,1],
    strata_ids          = dclus2wr$strata[,1],
    samp_unit_sel_probs = dclus2wr$allprob[,1],
    output = 'factors'
  ) |> `dimnames<-`(NULL)

  set.seed(1999)
  library_multistage_boot <- as_bootstrap_design(
    design = library_multistage_survey,
    type = "Antal-Tille",
    replicates = 100
  )
  set.seed(1999)
  library_multistage_reps <- make_doubled_half_bootstrap_weights(
    num_replicates = 100,
    samp_unit_ids = library_multistage_survey$cluster[,1],
    strata_ids = library_multistage_survey$strata[,1],
    samp_unit_sel_probs = library_multistage_survey$allprob[,1],
    output = 'factors'
  ) |> `dimnames<-`(NULL)

  set.seed(1999)
  election_pps_boot <- as_bootstrap_design(
    design = pps_design_ht,
    type = "Antal-Tille",
    replicates = 100
  )
  set.seed(1999)
  election_pps_reps <- make_doubled_half_bootstrap_weights(
    num_replicates = 100,
    samp_unit_ids = pps_design_ht$cluster[,1],
    strata_ids = pps_design_ht$strata[,1],
    samp_unit_sel_probs = pps_design_ht$allprob[,1],
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
        object   = dclus2wr_boot$repweights,
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
          as_bootstrap_design(dclus1, type = "Antal-Tille") |>
            weights(type = "analysis")
        }),
        expected = withr::with_seed(2014, {
          as_bootstrap_design(dclus1, type = "Antal-Tille") |>
            weights(type = "replication") * weights(dclus1)
        })
      )
      expect_equal(
        object = withr::with_seed(2014, {
          as_bootstrap_design(dclus1, type = "Antal-Tille",
                              replicates = 100) |>
            weights(type = "analysis")
        }),
        expected = withr::with_seed(2014, {
          make_doubled_half_bootstrap_weights(
            num_replicates = 100,
            samp_unit_ids = dclus1$cluster[,1],
            strata_ids = dclus1$strata[,1],
            samp_unit_sel_probs = dclus1$allprob[,1],
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
        ) |> as_bootstrap_design(type = "Antal-Tille"),
        regexp = "Cannot form bootstrap adjustment factors for a stratum"
      )
  })

# Check for correctness ----
  
  ##_ All replicate factors are 0, 1, 2, 3 ----
  
  test_that(
    "Replicate factors take on correct set of values", {
      expect_in(dclus2wr_reps,                0:3)
      expect_in(dclus1_reps,                  0:3)
      expect_in(election_pps_boot$repweights, 0:3)
    }
  )
  
  ##_ Factors have correct sums and covariances ----
  
  test_that(
    "PPSWOR: Replicate factors have expected sums and covariances", {
      
      expect_equal(
        object   = unique(election_pps_boot$repweights |> colSums()),
        expected = 40
      )
      
      expect_true(
        object = all(1e-8 > abs(
          cov(t(election_pps_boot$repweights)) |> rowSums()
        ))
      )
      
    }
  )
  
  test_that(
    "SRSWOR: Replicate factors have expected sums and covariances", {
      
      dclus1_psu_reps <- dclus1_reps[!duplicated(dclus1$cluster$dnum),]
      
      expect_equal(
        object   = unique(dclus1_psu_reps |> colSums()),
        expected = nrow(dclus1_psu_reps)
      )
      
      expect_true(
        object = all(1e-8 > abs(
          cov(t(dclus1_psu_reps)) |> rowSums()
        ))
      )
      
    }
  )
  
  test_that(
    "SRSWR: Replicate factors have expected sums and covariances", {
      
      dclus2wr_psu_reps <- dclus2wr_reps[!duplicated(dclus2wr$cluster$dnum),]
      
      expect_equal(
        object   = unique(dclus2wr_psu_reps |> colSums()),
        expected = nrow(dclus2wr_psu_reps)
      )
      
      expect_true(
        object = all(1e-8 > abs(
          cov(t(dclus2wr_psu_reps)) |> rowSums()
        ))
      )
      
    }
  )
