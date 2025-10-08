suppressWarnings({
  suppressPackageStartupMessages({
    library(survey)
    library(dplyr)
    library(svrep)
    library(testthat)
    library(srvyr)
  })
})

test_that("as_bootstrap_design runs without error on strata with *mostly* weight 1", {
    
  skip_if_not_installed("survey")

  set.seed(123)
  df <- data.frame(
    strata_id = rep(1:3, each = 10),
    wgt = c(
      # Stratum 1: mostly 1s and a couple 4s
      rep(1, 8), 4, 4,
      # Stratum 2: 40-60% non-1 weights
      sample(c(1, 2, 3), size = 10, replace = TRUE, prob = c(0.6, 0.2, 0.2)),
      # Stratum 3: similar ratio but random
      sample(c(1, 2, 5), size = 10, replace = TRUE, prob = c(0.5, 0.3, 0.2))
    ),
    y = rnorm(30)
  )
  
  pps_wor_design <- survey::svydesign(
    data = df,
    strata = ~strata_id,
    weights = ~wgt,
    ids = ~1
  )
  
  expect_error(
    as_bootstrap_design(
      design     = pps_wor_design,
      type       = "Antal-Tille",
      replicates = 10
    ),
    regexp = NA  # expects no error
  )
})