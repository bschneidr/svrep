suppressWarnings({
  suppressPackageStartupMessages({
    library(survey)
    library(dplyr)
    library(svrep)
    library(testthat)
    library(srvyr)
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

# Test correct conversion for bootstrap types other than RWYB ----

  test_that(
    "`as_bootstrap_design()` works with types other than RWYB", {
      expect_equal(
        object = withr::with_seed(2014, {
          as_bootstrap_design(dclus1, type = "Canty-Davison",
                              replicates = 10) |>
            weights(type = "analysis")
        }),
        expected = withr::with_seed(2014, {
          as.svrepdesign(dclus1, type = "bootstrap",
                         replicates = 10) |>
            weights(type = "analysis")
        })
      )
      expect_equal(
        object = withr::with_seed(2014, {
          as_bootstrap_design(dclus1, type = "Rao-Wu",
                              replicates = 10) |>
            weights(type = "analysis")
        }),
        expected = withr::with_seed(2014, {
          as.svrepdesign(dclus1, type = "subbootstrap",
                         replicates = 10) |>
            weights(type = "analysis")
        })
      )
      expect_equal(
        object = withr::with_seed(2014, {
          as_bootstrap_design(dclus1, type = "Preston",
                              replicates = 10) |>
            weights(type = "analysis")
        }),
        expected = withr::with_seed(2014, {
          as.svrepdesign(dclus1, type = "mrbbootstrap",
                         replicates = 10) |>
            weights(type = "analysis")
        })
      )
  })

# Test correct conversion whether user supplies weights or not ----

  data('mu284', package = 'survey')
  data('library_multistage_sample', package = 'svrep')

  # Manually specify weights
  design_1_srs <- svydesign(
    data = mu284 |> transform(wt = (n1/5)*(n2/3)),
    ids = ~ id1 + id2,
    fpc = ~ n1  +  n2,
    weights = ~ wt
  )
  design_1_pps <- svydesign(
    data = library_multistage_sample |> transform(wt = 1/SAMPLING_PROB),
    ids = ~ PSU_ID + SSU_ID,
    fpc = ~ PSU_SAMPLING_PROB  +  SSU_SAMPLING_PROB,
    pps = "brewer",
    weights = ~ wt
  )

  # Let 'survey' package automatically figure out weights
  design_2_srs <- svydesign(
    data = mu284,
    ids = ~ id1 + id2,
    fpc = ~ n1  +  n2
  )
  design_2_pps <- svydesign(
    data = library_multistage_sample,
    ids = ~ PSU_ID + SSU_ID,
    fpc = ~ PSU_SAMPLING_PROB  +  SSU_SAMPLING_PROB,
    pps = "brewer"
  )

  test_that(
    desc = "Check that stage-specific probs can be obtained with `svydesign(..., weights)`", {
      # For multistage SRSWOR
      set.seed(2014)
      boot_1_srs <- as_bootstrap_design(design_1_srs, replicates = 10)
      set.seed(2014)
      boot_2_srs <- as_bootstrap_design(design_2_srs, replicates = 10)
      expect_equal(object = weights(boot_1_srs, type = "analysis"),
                   expected = weights(boot_2_srs, type = "analysis"))
      # For two-stage PPSWOR+SSRWOR
      set.seed(2014)
      boot_1_pps <- as_bootstrap_design(design_1_pps, replicates = 10)
      set.seed(2014)
      boot_2_pps <- as_bootstrap_design(design_2_pps, replicates = 10)
      expect_equal(object = weights(boot_1_pps, type = "analysis"),
                   expected = weights(boot_2_pps, type = "analysis"))
    }
  )

# Works for survey design objects with special classes ----

  test_that(
    desc = "Returns `tbl_svy` if the input is a `tbl_svy` and 'srvyr' is loaded", {
      expect_true(
        dstrat |> as_survey() |>
          as_bootstrap_design(replicates = 1) |>
          inherits(what = "tbl_svy")
      )
    }
  )

# Informative error message ----

  test_that(
    desc = "Informative error message if user misspecifies `type` argument.", {

      expect_error(regexp = "Must use either",
        dstrat |> as_bootstrap_design(replicates = 1,
                                      type = "bootsy-collins")
      )

    }
  )
