suppressWarnings({
  suppressPackageStartupMessages({
    library(survey)
    library(testthat)
  })
})

set.seed(2014)

if (packageVersion("base") <= "4.4.0") {
  sort_by <- function(x, y, ...) {
      if (inherits(y, "formula")) 
          y <- .formula2varlist(y, x)
      if (!is.list(y)) 
          y <- list(y)
      o <- do.call(order, c(unname(y), list(...)))
      x[o, , drop = FALSE]
  }
}

# Tests for row assignment methods ----

test_that(
  "Correct row assignments when available Hadamard rows are at least n", {
    
    expect_equal(
      assign_hadamard_rows(
        n = 4, hadamard_order = 4,
        circular = TRUE
      ) |> unname(),
        expected = matrix(
          c(1,2,
            2,3,
            3,4,
            4,1), 
            byrow = TRUE, ncol = 2
        )
    )

    expect_equal(
      assign_hadamard_rows(
        n = 4, hadamard_order = 8,
        circular = FALSE
      ) |> unname(),
        expected = matrix(
          c(1,2,
            2,3,
            3,4,
            4,5), 
            byrow = TRUE, ncol = 2
        )
    )
        
  }
)

test_that(
  "Correct row assignments when available Hadamard rows are less than n", {
    
    expect_equal(
      assign_hadamard_rows(
        n = 9, hadamard_order = 4, use_first_row = TRUE
      ) |> unname(),
        expected = matrix(
          c(1,2,
            2,3,
            3,4,
            4,1,
            1,3,
            3,1,
            2,4,
            4,2,
            1,2), 
            byrow = TRUE, ncol = 2
        )
    )

    expect_equal(
      assign_hadamard_rows(
        n = 9, hadamard_order = 4, use_first_row = FALSE
      ) |> unname(),
        expected = matrix(
          c(2,3,
            3,4,
            4,2,
            2,4,
            4,3,
            3,2,
            2,3,
            3,4,
            4,2), 
            byrow = TRUE, ncol = 2
        )
    )
        
  }
)

# Expected quadratic forms ----

test_that("SDRM with enough replicates matches SD2", {

  # Checks for SD2
  suppressMessages({
    rep_factors <- make_sdr_replicate_factors(5, target_number_of_replicates = 8,
                                              use_normal_hadamard = TRUE)
  })
  
  expect_equal(
    ((4/8) * tcrossprod(rep_factors - 1)) |> round(2),
    expected = make_sd_matrix(n = 5, type = "SD2") |> as.matrix()
  )

  suppressMessages({
    rep_factors <- make_sdr_replicate_factors(5, target_number_of_replicates = 8,
                                              use_normal_hadamard = TRUE)
  })
  
  expect_equal(
    ((4/8) * tcrossprod(rep_factors - 1)) |> round(2),
    expected = make_sd_matrix(n = 5, type = "SD2") |> as.matrix()
  )

  }
)

# Expected variance estimates from `as_sdr_design()`

test_that("as_sdr_design() gives expected variance estimates in basic examples", {

  # Load example stratified systematic sample
  data('library_stsys_sample', package = 'svrep')

  # First, ensure data are sorted in same order as was used in sampling
  library_stsys_sample <- library_stsys_sample |> sort_by(~ SAMPLING_SORT_ORDER)

  # Create survey design object
  simple_design_obj <- svydesign(data = library_stsys_sample, ids = ~ 1,
                                 probs = ~ SAMPLING_PROB)

  # Compare SDR variance estimate against SD2 variance estimate
  expect_equal(
    object = suppressMessages({simple_design_obj |> as_sdr_design(
      replicates = 224, sort_variable = "SAMPLING_SORT_ORDER"
    ) |> svytotal(x = ~ LIBRARIA, na.rm = TRUE) |> SE()}),
    expected = suppressMessages({simple_design_obj |> as_fays_gen_rep_design(
      max_replicates = 224, variance_estimator = "SD2"
    ) |> svytotal(x = ~ LIBRARIA, na.rm = TRUE) |> SE()}),
    tolerance = 0.0001
  )

  # Check that first-stage strata are used as the primary sort variable
  reordered_design <- library_stsys_sample |> 
    transform(RAND_ORDER = runif(n = nrow(library_stsys_sample))) |>
    transform(ORDER_WITHIN_STRATUM = gsub(x = SAMPLING_SORT_ORDER,
                                          "\\d{3}_", "")) |>
    sort_by(~ RAND_ORDER) |>
    svydesign(data = _, ids = ~ 1, strata = ~ SAMPLING_STRATUM,
              probs = ~ SAMPLING_PROB)
  
  expect_equal(
    object = suppressMessages({reordered_design |> 
      as_sdr_design(sort_variable = "ORDER_WITHIN_STRATUM", replicates = 8) |> 
      svytotal(x = ~ LIBRARIA, na.rm = TRUE)}),
    expected = suppressMessages({reordered_design |> 
      as_sdr_design(sort_variable = "SAMPLING_SORT_ORDER", replicates = 8) |> 
      svytotal(x = ~ LIBRARIA, na.rm = TRUE)})
  )
  
  # Check that FPCs are incorporated into replicate factors
  # Load example stratified systematic sample

  # Create survey design object, with and without FPC
  design_obj_with_fpcs <- svydesign(
    data = library_stsys_sample |> 
      sort_by(~ SAMPLING_SORT_ORDER),
    ids = ~ 1,
    probs = ~ SAMPLING_PROB, 
    strata = ~ SAMPLING_STRATUM, fpc = ~ STRATUM_POP_SIZE
  )
  design_obj_without_fpcs <- svydesign(
    data = library_stsys_sample |> 
      sort_by(~ SAMPLING_SORT_ORDER),
    ids = ~ 1,
    probs = ~ SAMPLING_PROB, 
    strata = ~ SAMPLING_STRATUM
  )
  
  # Create SDR replicates with and without FPCs
  suppressMessages({
    sdr_svy_w_fpcs <- design_obj_with_fpcs |> as_sdr_design(
      sort_variable = "SAMPLING_SORT_ORDER", replicates = 256
    )
  })
  
  # Check that quadratic form based on replicates has the expected FPCs
  rep_quad_form_diag <- diag(
    ((4/256) * tcrossprod(weights(sdr_svy_w_fpcs, 'replication') - 1))
  )
  
  expected_fpc <- 1 - as.vector(design_obj_with_fpcs$fpc$sampsize/design_obj_with_fpcs$fpc$popsize)
  
  expect_equal(rep_quad_form_diag, expected_fpc)
  
  # Check that variance estimate with FPC matches the SD2 estimator with FPC, for totals
  expect_equal(
    object   = suppressMessages({library_stsys_sample |> sort_by(~SAMPLING_SORT_ORDER) |> 
        transform(POP_SIZE = 1000) |> svydesign(
          data = _, ids = ~ 1, prob = ~ SAMPLING_PROB, fpc = ~ POP_SIZE
        ) |> 
        as_sdr_design(sort_variable = "SAMPLING_SORT_ORDER", replicates = 224) |> 
        svytotal(x = ~ LIBRARIA, na.rm = TRUE) |>
        SE()
    }),
    expected = suppressMessages({library_stsys_sample |> transform(POP_SIZE = 1000) |> svydesign(
      data = _, ids = ~ 1, prob = ~ SAMPLING_PROB, fpc = ~ POP_SIZE
    ) |> 
        as_fays_gen_rep_design("SD2") |> 
        svytotal(x = ~ LIBRARIA, na.rm = TRUE) |>
        SE()
    }),
    tolerance = 0.0001
  )

})
