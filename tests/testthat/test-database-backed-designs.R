library(survey)
library(RSQLite)

# Load example data
data('api', package = 'survey')

dbclus1 <- svydesign(
  id = ~dnum, weights = ~pw, fpc = ~fpc,
  data = "apiclus1",
  dbtype = "SQLite",
  dbname = system.file("api.db",package="survey")
)

dclus1 <- svydesign(
  data = apiclus1,
  id = ~ dnum, weights = ~ pw, fpc = ~ fpc
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

# Test weight adjustment functionality ----

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
