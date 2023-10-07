suppressPackageStartupMessages(library(survey))

# Create example data ----
set.seed(1999)

#
library(survey)
set.seed(2023)

# Create an example survey design object

sample_data <- data.frame(
  Y       = c(8.11, 7.42, 14.54, 18.21, 9.18, 13.83, 7.57, 9.87),
  STRATUM = c(1,1,1,1,2,2,2,2),
  PSU     = c(1,2,3,4,5,6,7,8),
  FPC     = rep(c(0.5, 0.3), each = 4)
)

survey_design <- svydesign(
  data    = sample_data,
  strata  = ~ STRATUM,
  ids     = ~ PSU,
  weights = ~ 1,
  fpc     = ~ FPC
)

jkn_design <- survey_design |> as.svrepdesign(type = "JKn")

test_that("Shuffling replicates still gives correct variance estimates", {

  jkn_design <- survey_design |> as.svrepdesign(type = "JKn", compress = TRUE)
  expect_equal(
    object = jkn_design |> shuffle_replicates() |> svytotal(x = ~ Y) |> vcov(),
    expected = jkn_design |> svytotal(x = ~ Y) |> vcov()
  )
  jkn_design <- survey_design |> as.svrepdesign(type = "JKn", compress = FALSE)
  expect_equal(
    object = jkn_design |> shuffle_replicates() |> svytotal(x = ~ Y) |> vcov(),
    expected = jkn_design |> svytotal(x = ~ Y) |> vcov()
  )
})

