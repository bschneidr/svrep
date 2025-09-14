suppressPackageStartupMessages(library(survey))

# Create example data ----
set.seed(2023)

svy_data <- data.frame(
  WGT  = c(10,  11,  12 ),
  REP1 = c(0.0, 1.5, 1.5),
  REP2 = c(1.5, 0.0, 1.5),
  REP3 = c(1.5, 1.5, 0.0)
)

rep_svy <- svrepdesign(
  data       = svy_data,
  weights    = ~ WGT,
  repweights = "REP[1-3]",
  combined   = FALSE,
  type       = "JKn",
  scale      = 2/3,
  rscales    = c(0.75, 0.8, 0.9)
)

test_that("`get_rep_scale_coefs()` works correctly", {
  expect_equal(
    object = rep_svy |> get_rep_scale_coefs("combined"),
    expected = (2/3) * c(0.75, 0.8, 0.9)
  )
  expect_equal(
    object = rep_svy |> get_rep_scale_coefs("specific"),
    expected = c(0.75, 0.8, 0.9)
  )
  expect_equal(
    object = rep_svy |> get_rep_scale_coefs("overall"),
    expected = (2/3)
  )
})

test_that("`get_rep_type()` works correctly", {
  expect_equal(
    object = rep_svy |> get_rep_type(),
    expected = "JKn"
  )
})
