library(survey)
library(RSQLite)

# Load example data
data('api', package = 'survey')

# Test the generalized bootstrap functionality
test_that(
  desc = "Generalized bootstrap works for database-backed designs", {

    # Create generalized bootstrap replicates for database-backed design

    dbclus1 <- svydesign(
      id = ~dnum, weights = ~pw, fpc = ~fpc,
      data = "apiclus1",
      dbtype = "SQLite",
      dbname = system.file("api.db",package="survey")
    )

    set.seed(2023)
    db_result <- dbclus1 |> as_gen_boot_design(
      replicates = 5, variance_estimator = "Ultimate Cluster"
    )

    # Create replicates for design in local memory
    set.seed(2023)
    non_db_result <- svydesign(
      data = apiclus1,
      id = ~ dnum, weights = ~ pw, fpc = ~ fpc
    ) |> as_gen_boot_design(
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
