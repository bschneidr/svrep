suppressPackageStartupMessages({
  suppressWarnings({
    library(survey)
    library(RSQLite)
  })
})

# Load example data
data('api', package = 'survey')

suppressWarnings({
  dbclus1 <- svydesign(
    id = ~dnum, weights = ~pw, fpc = ~fpc,
    data = "apiclus1",
    dbtype = "SQLite",
    dbname = system.file("api.db",package="survey")
  )
})

dclus1 <- svydesign(
  data = apiclus1,
  id = ~ dnum, weights = ~ pw, fpc = ~ fpc
)

# Test the data frame conversion functionality ----

test_that(
  "`as_data_frame_with_weights()` works for database-backed designs", {

    expect_equal(
      object = suppressWarnings(
        dbclus1 |> as_data_frame_with_weights(
          vars_to_keep = c("dnum", "snum"),
          full_wgt_name = "FULL_SAMPLE_WGT"
        )
      ),
      expected = dclus1 |> as_data_frame_with_weights(
        vars_to_keep = c("dnum", "snum"),
        full_wgt_name = "FULL_SAMPLE_WGT"
      )
    )

  }
)


# Test the generalized bootstrap functionality ----
test_that(
  desc = "Generalized bootstrap works for database-backed designs", {

    # Create generalized bootstrap replicates for database-backed design

    set.seed(2023)
    db_result <- dbclus1 |> as_gen_boot_design(
      replicates = 5, variance_estimator = "Ultimate Cluster"
    )

    # Create replicates for design in local memory
    set.seed(2023)
    non_db_result <- dclus1 |> as_gen_boot_design(
      replicates = 5, variance_estimator = "Ultimate Cluster"
    )

    # Compare two sets of replicate weights
    expect_equal(db_result$repweights, non_db_result$repweights)
    # Compare estimates
    expect_equal(
      svytotal(x = ~ api00, db_result) |> SE(),
      svytotal(x = ~ api00, non_db_result) |> SE()
    )
  }
)

# Test Fay's generalized replication functionality ----
test_that(
  desc = "Generalized replication works for database-backed designs", {

    # Averts ATLAS/MKL tests (not supported)
    skip_if(grepl(x = La_library(), pattern = "atlas|mkl",
                  ignore.case = TRUE))
    skip_if(grepl(x = extSoftVersion()[['BLAS']], pattern = "atlas|mkl",
                  ignore.case = TRUE))

    # Create Fay's generalized replication replicates for database-backed design

    set.seed(2023)
    db_result <- dbclus1 |> as_fays_gen_rep_design(
      max_replicates = 16, variance_estimator = "Ultimate Cluster"
    )

    # Create replicates for design in local memory
    set.seed(2023)
    non_db_result <- dclus1 |> as_fays_gen_rep_design(
      max_replicates = 16, variance_estimator = "Ultimate Cluster"
    )

    # Compare two sets of replicate weights
    expect_equal(db_result$repweights, non_db_result$repweights)
    # Compare estimates
    expect_equal(
      svytotal(x = ~ api00, db_result) |> SE(),
      svytotal(x = ~ api00, non_db_result) |> SE()
    )
  }
)

# Test the basic bootstrap functionality ----

test_that(
  desc = "`as_bootstrap_design()` works for database-backed designs", {

    # Create bootstrap replicates for database-backed design
    set.seed(2023)
    db_rwyb_result <- dbclus1 |> as_bootstrap_design(
      replicates = 5,
      type = "Rao-Wu-Yue-Beaumont"
    )
    set.seed(2023)
    db_survey_pkg_boot_result <- dbclus1 |> as_bootstrap_design(
      replicates = 5, type = "Preston"
    )

    # Create replicates for design in local memory
    set.seed(2023)
    nondb_rwyb_result <- dclus1 |> as_bootstrap_design(
      replicates = 5,
      type = "Rao-Wu-Yue-Beaumont"
    )
    set.seed(2023)
    nondb_survey_pkg_boot_result <- dclus1 |> as_bootstrap_design(
      replicates = 5, type = "Preston"
    )

    # Compare two sets of replicate weights
    expect_equal(
      db_rwyb_result$repweights,
      nondb_rwyb_result$repweights
    )
    expect_equal(
      db_survey_pkg_boot_result$repweights,
      nondb_survey_pkg_boot_result$repweights
    )
    # Compare estimates
    expect_equal(
      svytotal(x = ~ api00, db_rwyb_result) |> SE(),
      svytotal(x = ~ api00, nondb_rwyb_result) |> SE()
    )
  }
)

# Test jackknife functionality ----

test_that(
  desc = "`as_random_group_jackknife_design()` works for database-backed designs", {

    suppressWarnings({
      # Create bootstrap replicates for database-backed design
      set.seed(2023)
      db_jk_result <- dbclus1 |>
        as_random_group_jackknife_design(
          replicates = 5
        )

      # Create replicates for design in local memory
      set.seed(2023)
      nondb_jk_result <- dclus1 |>
        as_random_group_jackknife_design(
          replicates = 5
        )
    })

    # Compare estimates
    expect_equal(
      svytotal(x = ~ api00, db_jk_result) |> SE(),
      svytotal(x = ~ api00, nondb_jk_result) |> SE()
    )
  }
)

# Test SDR method functionality ----

test_that(
  desc = "`as_sdr_design()` works for database-backed designs", {

    suppressMessages({
      suppressWarnings({
        # Create bootstrap replicates for database-backed design
        set.seed(2023)
        db_sdr_result <- dbclus1 |>
          as_sdr_design(
            sort_variable = "dnum",
            replicates = 16
          )
  
        # Create replicates for design in local memory
        set.seed(2023)
        nondb_sdr_result <- dclus1 |>
          as_sdr_design(
            sort_variable = "dnum",
            replicates = 16
          )
      })
    })

    # Compare estimates
    expect_equal(
      svytotal(x = ~ api00, db_sdr_result)    |> SE(),
      svytotal(x = ~ api00, nondb_sdr_result) |> SE()
    )
    }
)

# Test weight redistribution functionality ----

test_that(
  desc = "Weight adjustments work for database-backed designs", {

    # Create replicate designs
    set.seed(2023)
    db_rep_design <- dbclus1 |> as_gen_boot_design(
      replicates = 5, variance_estimator = "Ultimate Cluster"
    )

    set.seed(2023)
    nondb_rep_design <- dclus1 |> as_gen_boot_design(
      replicates = 5, variance_estimator = "Ultimate Cluster"
    )

    # Compare results of weight redistribution
    db_result <- db_rep_design |>
      redistribute_weights(
        reduce_if = (stype == "H"),
        increase_if = (stype == "M"),
        by = "sch_wide"
      )

    nondb_result <- nondb_rep_design |>
      redistribute_weights(
        reduce_if = (stype == "H"),
        increase_if = (stype == "M"),
        by = "sch.wide"
      )

    expect_equal(
      object = db_result |> weights(type = "analysis"),
      expected = nondb_result |> weights(type = "analysis")
    )
    expect_equal(
      object = db_result |> weights(type = "sampling"),
      expected = nondb_result |> weights(type = "sampling")
    )

    # Check that database-backed result has the expected class
    expect_true(
      object = inherits(db_result, "DBIrepdesign")
    )
  }
)

# Test sample-based calibration methods ----

test_that(
  "`calibrate_to_sample()` works for database-backed designs", {

    set.seed(2023)
    db_primary_survey <- dbclus1 |> as_bootstrap_design(
      replicates = 5, mse = TRUE
    )
    set.seed(2023)
    nondb_primary_survey <- dclus1 |> as_bootstrap_design(
      replicates = 5, mse = TRUE
    )

    # Load example data for control survey

    control_survey <- svydesign(id = ~ 1, fpc = ~fpc, data = apisrs) |>
      as_bootstrap_design(replicates = 5)

    # Calibrate DB-backed and regular designs
    suppressMessages({
      set.seed(2023)
      calibrated_db_rep_design <- calibrate_to_sample(
        primary_rep_design = db_primary_survey,
        control_rep_design = control_survey,
        cal_formula = ~ stype + enroll,
      )
      set.seed(2023)
      calibrated_nondb_rep_design <- calibrate_to_sample(
        primary_rep_design = nondb_primary_survey,
        control_rep_design = control_survey,
        cal_formula = ~ stype + enroll,
      )
    })

    # Compare calibrated-estimates
    expect_equal(
      object = calibrated_db_rep_design |>
        svytotal(x = ~ api00),
      expected = calibrated_nondb_rep_design |>
        svytotal(x = ~ api00)
    )

  }
)

test_that(
  "`calibrate_to_estimate()` works for database-backed designs", {

    set.seed(2023)
    db_primary_survey <- dbclus1 |> as_bootstrap_design(
      replicates = 10, mse = TRUE
    )
    set.seed(2023)
    nondb_primary_survey <- dclus1 |> as_bootstrap_design(
      replicates = 10, mse = TRUE
    )

    # Load example estimates from control survey
    estimated_controls <- svytotal(
      x = ~ stype,
      design = svydesign(id = ~ 1, fpc = ~fpc, data = apisrs)
    )
    control_point_estimates <- coef(estimated_controls)
    control_vcov_estimate <- vcov(estimated_controls)

    # Calibrate DB-backed and regular designs
    suppressMessages({
      set.seed(2023)
      calibrated_db_rep_design <- calibrate_to_estimate(
        rep_design = db_primary_survey,
        estimate = control_point_estimates,
        vcov_estimate = control_vcov_estimate,
        cal_formula = ~ stype
      )
      set.seed(2023)
      calibrated_nondb_rep_design <- calibrate_to_estimate(
        rep_design = nondb_primary_survey,
        estimate = control_point_estimates,
        vcov_estimate = control_vcov_estimate,
        cal_formula = ~ stype
      )
    })

    # Compare calibrated-estimates
    expect_equal(
      object = calibrated_db_rep_design |>
        svytotal(x = ~ api00),
      expected = calibrated_nondb_rep_design |>
        svytotal(x = ~ api00)
    )

  }
)

# Disconnect from the database
dbclus1$db$connection |> dbDisconnect()
rm(dbclus1)