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

# Create survey design objects ----
  dstrat_fpc <- svydesign(
    data = apistrat,
    id = ~1,
    strata = ~stype,
    weights = ~pw,
    fpc = ~fpc
  )

  dstrat_nofpc <- svydesign(
    data = apistrat,
    id = ~1,
    strata = ~stype,
    weights = ~pw
  )

# Test for expected formation of random groups ----

  test_that(desc = "Random groups (without var_strat) formed correctly", {

    ## Correct sizes for random groups
    jk_design <- as_random_group_jackknife_design(
      design = dstrat_nofpc,
      replicates = 13
    )

    expect_equal(
      object = jk_design$variables |>
        count(.random_group) |>
        pull("n"),
      expected = c(rep(16, times = 5),
                   rep(15, times = 8))
    )

    ## Each PSU is mapped to one and only one group
    jk_dclus1 <- as_random_group_jackknife_design(
      design = svydesign(id=~dnum, weights=~pw, data=apiclus1),
      replicates = 5
    )

    expect_equal(
      object = jk_dclus1 |> getElement("variables") |>
        distinct(dnum, .random_group) |>
        nrow(),
      expected = 15
    )

  })

  test_that(desc = "Random groups (with var_strat) formed correctly", {

    ## Correct sizes for random groups
    jk_design <- as_random_group_jackknife_design(
      design = dstrat_nofpc,
      var_strat = "stype",
      replicates = 5
    )

    expect_equal(
      object = jk_design$variables |>
        count(stype, .random_group) |>
        pull("n"),
      expected = c(rep(c(20, 10, 10), each = 5))
    )

  })

# Given the random groups, replicates formed correctly ----

  test_that(desc = "Replicates formed correctly, given random groups", {

    jk_design <- as_random_group_jackknife_design(
      design = dstrat_nofpc, replicates = 10,
      adj_method = "variance-units",
      scale_method = "variance-units"
    )

    apistrat_grouped <- apistrat
    apistrat_grouped[['GRP']] <- jk_design |>
      getElement("variables") |>
      getElement(".random_group")

    expect_equal(
      object = jk_design$repweights,
      expected = svydesign(data = apistrat_grouped,
                           ids = ~ GRP,
                           weights = ~ 1) |>
        as.svrepdesign(type = "JK1") |>
        getElement("repweights")
    )

  })

# Test correct handling of FPCs ----

  test_that(desc = "Correct handling of FPCs", {

    expect_equal(
      object = as_random_group_jackknife_design(
        design = dstrat_nofpc, var_strat_frac = 0.75,
        replicates = 5
      ) |> getElement("scale"),
      expected = ((5-1)/5) * (1 - 0.75)
    )

    expect_warning(object = {
      as_random_group_jackknife_design(
        design = dstrat_fpc, replicates = 5
      )
    }, regexp = "Ignoring finite"
    )

    expect_error(object = {
      as_random_group_jackknife_design(
        design = dstrat_nofpc,
        replicates = 5,
        var_strat_frac = 2
      )
    }, regexp = "between 0 and 1"
    )

    expect_error(object = {
      as_random_group_jackknife_design(
        design = dstrat_nofpc,
        replicates = 5,
        var_strat_frac = c(0.5, 0.2)
      )
    }, regexp = "must be either a single number"
    )

  })

# Checks for valid number of replicates ----

  test_that(desc = "Alerts user if invalid number of replicates", {
    expect_error(
      object = {as_random_group_jackknife_design(
        design = dstrat_nofpc, replicates = 51
      )},
      regexp = "fewer PSUs than"
    )
  })

# Saves random groups variable

  test_that(desc = "Saves random groups variable correctly", {
    expect_equal(object = {
      as_random_group_jackknife_design(
        design = dstrat_nofpc, replicates = 5,
        group_var_name = "VAR_UNIT"
      ) |> colnames() |> intersect("VAR_UNIT")
    }, expected = "VAR_UNIT")

    expect_warning(object = {
      as_random_group_jackknife_design(
        design = dstrat_nofpc, replicates = 5,
        group_var_name = "stype"
      )
    }, regexp = "Overwriting.+stype")
  })

# Works for more specialized classes of survey designs ----

  test_that(
    desc = "Returns `tbl_svy` if the input is a `tbl_svy` and 'srvyr' is loaded", {
      suppressPackageStartupMessages({library(srvyr)})
      expect_true(
        suppressMessages({
          dstrat_nofpc |> as_survey() |>
            as_random_group_jackknife_design(
              replicates = 10
            ) |>
            inherits('tbl_svy')
        })
      )
    }
  )
