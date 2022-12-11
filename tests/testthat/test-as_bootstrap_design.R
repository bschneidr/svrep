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

# Test correct conversion for bootstrap types other than RWYB ----

  test_that(
    "`as_bootstrap_design()` works with types other than RWYB", {
      expect_equal(
        object = withr::with_seed(2014, {
          as_bootstrap_design(dclus1, type = "Canty-Davison") |>
            weights(type = "analysis")
        }),
        expected = withr::with_seed(2014, {
          as.svrepdesign(dclus1, type = "bootstrap") |>
            weights(type = "analysis")
        })
      )
      expect_equal(
        object = withr::with_seed(2014, {
          as_bootstrap_design(dclus1, type = "Rao-Wu") |>
            weights(type = "analysis")
        }),
        expected = withr::with_seed(2014, {
          as.svrepdesign(dclus1, type = "subbootstrap") |>
            weights(type = "analysis")
        })
      )
      expect_equal(
        object = withr::with_seed(2014, {
          as_bootstrap_design(dclus1, type = "Preston") |>
            weights(type = "analysis")
        }),
        expected = withr::with_seed(2014, {
          as.svrepdesign(dclus1, type = "mrbbootstrap") |>
            weights(type = "analysis")
        })
      )
  })
