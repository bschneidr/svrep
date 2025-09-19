suppressWarnings({
  suppressPackageStartupMessages({
    library(survey)
    library(dplyr)
    library(svrep)
    library(testthat)
  })
})

data('library_stsys_sample', package = 'svrep')
set.seed(2014)

library_stsys_sample <- library_stsys_sample |>
  mutate(
    TOTCIR = ifelse(is.na(TOTCIR), 0, TOTCIR),
    TOTSTAFF = ifelse(is.na(TOTSTAFF), 0, TOTSTAFF)
  )

data('mu284', package = 'survey')

mu284$phase2 <- mu284$id1 %in% c(19,31,45)
mu284 <- mu284[order(mu284$id1),]

twophase_design <- twophase(
  data = mu284,
  id = list(~id1, ~id1),
  fpc = list(~n1, NULL),
  subset = ~ phase2
)

# Same results from conversion vs. creating from scratch ----

set.seed(2014)

data('election', package = 'survey')

test_that(
  "Horvitz-Thompson: Same results from conversion vs. creating from scratch", {

    # Averts ATLAS/MKL tests (not supported)
    skip_if(grepl(x = La_library(), pattern = "atlas|mkl",
                  ignore.case = TRUE))
    skip_if(grepl(x = extSoftVersion()[['BLAS']], pattern = "atlas|mkl",
                  ignore.case = TRUE))

    # Disable test until the Matrix package adapts to change in R 4.3.2
    # (`crossprod()` becomes primitive and S3 generic in R 4.3.2)
    skip_if_not("crossprod" %in% ls(getNamespaceInfo("Matrix", "exports")))

    ## Create survey design object
    pps_design_ht <- svydesign(
      data = election_pps,
      id = ~1, fpc = ~p,
      pps = ppsmat(election_jointprob),
      variance = "HT"
    )

    ## Convert to generalized replication design
    set.seed(2014)
    conversion_result <- pps_design_ht |>
      as_fays_gen_rep_design(variance_estimator = "Horvitz-Thompson",
                             max_replicates = 44) |>
      weights(type = "replication") |> cov()


    set.seed(2014)
    orig_result <- make_quad_form_matrix(
      variance_estimator = "Horvitz-Thompson",
      joint_probs = election_jointprob
    ) |>
      make_fays_gen_rep_factors(max_replicates = 44) |>
      cov()

    expect_equal(
      object =  conversion_result,
      expected = orig_result,
      tolerance = 0.001
    )

  })

test_that(
  "Yates-Grundy: Same results from conversion vs. creating from scratch", {

    # Averts ATLAS/MKL tests (not supported)
    skip_if(grepl(x = La_library(), pattern = "atlas|mkl",
                  ignore.case = TRUE))
    skip_if(grepl(x = extSoftVersion()[['BLAS']], pattern = "atlas|mkl",
                  ignore.case = TRUE))

    # Disable test until the Matrix package adapts to change in R 4.3.2
    # (`crossprod()` becomes primitive and S3 generic in R 4.3.2)
    skip_if_not("crossprod" %in% ls(getNamespaceInfo("Matrix", "exports")))

    ## Create survey design object
    pps_design_yg <- svydesign(
      data = election_pps,
      id = ~1, fpc = ~p,
      pps = ppsmat(election_jointprob),
      variance = "YG"
    )

    ## Convert to generalized replication design
    set.seed(2014)
    conversion_result <- pps_design_yg |>
      as_fays_gen_rep_design(variance_estimator = "Yates-Grundy",
                             max_replicates = 44) |>
      weights(type = "replication") |> cov()


    set.seed(2014)
    orig_result <- make_quad_form_matrix(variance_estimator = "Yates-Grundy",
                                         joint_probs = election_jointprob) |>
      make_fays_gen_rep_factors(max_replicates = 44) |> cov()

    expect_equal(
      object = conversion_result,
      expected = orig_result,
      tolerance = 0.001
    )

  })

test_that(
  "Ultimate Cluster: Same results from conversion vs. creating from scratch", {

    # Averts ATLAS/MKL tests (not supported)
    skip_if(grepl(x = La_library(), pattern = "atlas|mkl",
                  ignore.case = TRUE))
    skip_if(grepl(x = extSoftVersion()[['BLAS']], pattern = "atlas|mkl",
                  ignore.case = TRUE))

    ## Create survey design object
    multistage_survey_design <- svydesign(
      data = svrep::library_multistage_sample |>
        mutate(WT = 1/SAMPLING_PROB),
      weights = ~ WT,
      ids = ~ PSU_ID + SSU_ID,
      fpc = ~ PSU_POP_SIZE + SSU_POP_SIZE
    )

    ## Convert to generalized replication design
    set.seed(2014)
    conversion_result <- multistage_survey_design |>
      as_fays_gen_rep_design(variance_estimator = "Ultimate Cluster",
                             max_replicates = 104) |>
      weights(type = "replication")


    set.seed(2014)
    orig_result <- make_quad_form_matrix(variance_estimator = "Ultimate Cluster",
                                         strata_ids = multistage_survey_design$strata,
                                         cluster_ids = multistage_survey_design$cluster,
                                         strata_pop_sizes = multistage_survey_design$fpc$popsize) |>
      make_fays_gen_rep_factors(max_replicates = 104)

    expect_equal(
      object = conversion_result,
      expected = orig_result
    )

  })

test_that(
  "Multistage Stratified SRS: Same results from conversion vs. creating from scratch", {

    # Averts ATLAS/MKL tests (not supported)
    skip_if(grepl(x = La_library(), pattern = "atlas|mkl",
                  ignore.case = TRUE))
    skip_if(grepl(x = extSoftVersion()[['BLAS']], pattern = "atlas|mkl",
                  ignore.case = TRUE))

    ## Create survey design object
    multistage_survey_design <- svydesign(
      data = svrep::library_multistage_sample |>
        mutate(WT = 1/SAMPLING_PROB),
      weights = ~ WT,
      ids = ~ PSU_ID + SSU_ID,
      fpc = ~ PSU_POP_SIZE + SSU_POP_SIZE
    )

    ## Convert to generalized replication design
    suppressMessages({
      set.seed(2014)
      conversion_result <- multistage_survey_design |>
        as_fays_gen_rep_design(variance_estimator = "Stratified Multistage SRS",
                               max_replicates = 5) |>
        weights(type = "replication")


      set.seed(2014)
      orig_result <- make_quad_form_matrix(variance_estimator = "Stratified Multistage SRS",
                                           strata_ids = multistage_survey_design$strata,
                                           cluster_ids = multistage_survey_design$cluster,
                                           strata_pop_sizes = multistage_survey_design$fpc$popsize) |>
        make_fays_gen_rep_factors(max_replicates = 5)
    })

    expect_equal(
      object = conversion_result,
      expected = orig_result
    )

  })

test_that(
  "SD1: Same results from conversion vs. creating from scratch", {

    # Averts ATLAS/MKL tests (not supported)
    skip_if(grepl(x = La_library(), pattern = "atlas|mkl",
                  ignore.case = TRUE))
    skip_if(grepl(x = extSoftVersion()[['BLAS']], pattern = "atlas|mkl",
                  ignore.case = TRUE))

    ## Create survey design object
    multistage_survey_design <- svydesign(
      data = svrep::library_multistage_sample |>
        mutate(WT = 1/SAMPLING_PROB),
      weights = ~ WT,
      ids = ~ PSU_ID + SSU_ID,
      fpc = ~ PSU_POP_SIZE + SSU_POP_SIZE
    )

    suppressMessages({
      ## Convert to generalized replication design
      set.seed(2014)
      conversion_result <- multistage_survey_design |>
        as_fays_gen_rep_design(variance_estimator = "SD1",
                               max_replicates = 5) |>
        weights(type = "replication")


      set.seed(2014)
      orig_result <- make_quad_form_matrix(variance_estimator = "SD1",
                                           strata_ids = multistage_survey_design$strata,
                                           cluster_ids = multistage_survey_design$cluster,
                                           strata_pop_sizes = multistage_survey_design$fpc$popsize,
                                           sort_order = seq_len(nrow(multistage_survey_design))) |>
        make_fays_gen_rep_factors(max_replicates = 5)
    })

    expect_equal(
      object = conversion_result,
      expected = orig_result
    )

  })

test_that(
  "SD2: Same results from conversion vs. creating from scratch", {

    # Averts ATLAS/MKL tests (not supported)
    skip_if(grepl(x = La_library(), pattern = "atlas|mkl",
                  ignore.case = TRUE))
    skip_if(grepl(x = extSoftVersion()[['BLAS']], pattern = "atlas|mkl",
                  ignore.case = TRUE))

    ## Create survey design object
    multistage_survey_design <- svydesign(
      data = svrep::library_multistage_sample |>
        mutate(WT = 1/SAMPLING_PROB),
      weights = ~ WT,
      ids = ~ PSU_ID + SSU_ID,
      fpc = ~ PSU_POP_SIZE + SSU_POP_SIZE
    )

    suppressMessages({
      ## Convert to generalized replication design
      set.seed(2014)
      conversion_result <- multistage_survey_design |>
        as_fays_gen_rep_design(variance_estimator = "SD2",
                               max_replicates = 5) |>
        weights(type = "replication")


      set.seed(2014)
      orig_result <- make_quad_form_matrix(variance_estimator = "SD2",
                                           strata_ids = multistage_survey_design$strata,
                                           cluster_ids = multistage_survey_design$cluster,
                                           strata_pop_sizes = multistage_survey_design$fpc$popsize,
                                           sort_order = seq_len(nrow(multistage_survey_design))) |>
        make_fays_gen_rep_factors(max_replicates = 5)
    })

    expect_equal(
      object = conversion_result,
      expected = orig_result
    )

  })

test_that(
  "Deville-Tille: Same results from conversion vs. creating from scratch", {

    # Averts ATLAS/MKL tests (not supported)
    skip_if(grepl(x = La_library(), pattern = "atlas|mkl",
                  ignore.case = TRUE))
    skip_if(grepl(x = extSoftVersion()[['BLAS']], pattern = "atlas|mkl",
                  ignore.case = TRUE))

    data('api', package = 'survey')

    aux_var_columns <- model.matrix(object = ~ -1 + stype, data = apistrat)

    apistrat <- cbind(apistrat, aux_var_columns)

    ## Create survey design object
    apistrat_survey_design <- svydesign(
      data = apistrat,
      weights = ~ pw,
      ids = ~ 1
    )

    suppressMessages({
      ## Convert to generalized replication design
      set.seed(2014)
      conversion_result <- apistrat_survey_design |>
        as_fays_gen_rep_design(variance_estimator = "Deville-Tille",
                               aux_var_names = c("stypeE", "stypeM", "stypeH"),
                               max_replicates = 5) |>
        weights(type = "replication")


      set.seed(2014)
      orig_result <- make_quad_form_matrix(variance_estimator = "Deville-Tille",
                                           strata_ids = apistrat_survey_design$strata,
                                           cluster_ids = apistrat_survey_design$cluster,
                                           probs = matrix(apistrat_survey_design$variables$pw^(-1),
                                                          nrow = nrow(apistrat_survey_design)),
                                           aux_vars = aux_var_columns) |>
        make_fays_gen_rep_factors(max_replicates = 5)
    })

    expect_equal(
      object = conversion_result,
      expected = orig_result
    )

  })

test_that(
  "Two-phase Design: Same results from conversion vs. creating from scratch", {

    # Averts ATLAS/MKL tests (not supported)
    skip_if(grepl(x = La_library(), pattern = "atlas|mkl",
                  ignore.case = TRUE))
    skip_if(grepl(x = extSoftVersion()[['BLAS']], pattern = "atlas|mkl",
                  ignore.case = TRUE))

    suppressMessages({
      expect_warning(
        regexp = "The sample quadratic form matrix for this design and variance estimator is not positive semidefinite.", object = {
          set.seed(2023)
          twophase_gen_repl <- twophase_design |>
            as_fays_gen_rep_design(
              max_replicates = 5,
              variance_estimator = list('SD2', 'Ultimate Cluster'),
              psd_option = "warn"
            )
        }
      )
    })

    suppressMessages({
      suppressWarnings({
        set.seed(2023)
        gen_repl_reps <- twophase_design |>
          get_design_quad_form(
            variance_estimator = list('SD2', 'Ultimate Cluster')
          ) |>
          get_nearest_psd_matrix() |>
          make_fays_gen_rep_factors(
            max_replicates = 5
          )

        expect_equal(
          object = twophase_gen_repl$repweights,
          expected = gen_repl_reps
        )
      })
    })

  })

# Rescaling functions work as expected ----

test_that(
  "Rescaling functions work as expected", {

    suppressMessages({
      suppressWarnings({
        twophase_gen_repl <- twophase_design |>
          as_fays_gen_rep_design(
            variance_estimator = list(
              'SD2', 'Ultimate Cluster'
            ),
            max_replicates = 5
          )
      })
    })

    rescaled_design <- twophase_gen_repl |>
      rescale_replicates(min_wgt = 0.05)
    rescaled_matrix <- twophase_gen_repl |>
      weights(type = "replication") |>
      `attr<-`('scale', twophase_gen_repl$scale) |>
      rescale_replicates(min_wgt = 0.05)

    expect_equal(
      object = rescaled_design |> weights(type = "replication"),
      expected = rescaled_matrix
    )

    expect_gte(
      object = min(rescaled_matrix), expected = 0.05
    )

    expect_equal(
      object = matrix(c(1,0.3,0.02,2,3,4), ncol = 3) |>
        `attr<-`('scale', 1) |>
        rescale_replicates(min_wgt = 0.02) |>
        `attr<-`('tau', NULL),
      expected = matrix(c(1,0.3,0.02,2,3,4), ncol = 3)
    )

  })

# Variance estimates give expected result ----

test_that(
  desc = "Using full number of necessary replicates gives exact variance estimate for totals", {

    # Disable test until the Matrix package adapts to change in R 4.3.2
    # (`crossprod()` becomes primitive and S3 generic in R 4.3.2)
    skip_if_not("crossprod" %in% ls(getNamespaceInfo("Matrix", "exports")))

    set.seed(2023)
    pps_design_yg <- svydesign(
      data = election_pps,
      id = ~1, fpc = ~p,
      pps = ppsmat(election_jointprob),
      variance = "YG"
    )

    Sigma <- get_design_quad_form(pps_design_yg, "Yates-Grundy")

    gen_rep_est <- as_fays_gen_rep_design(
      pps_design_yg, "Yates-Grundy",
      max_replicates = 40
    ) |>
      svytotal(x = ~ Bush + Kerry) |> vcov() |>
      `attr<-`('means', NULL)

    exact_est <- election_pps[,c("Bush", "Kerry")] |>
      apply(MARGIN = 2, function(x) x/election_pps$p) |>
      (\(X) t(X) %*% Sigma %*% X)() |>
      as.matrix()

    expect_equal(object = gen_rep_est, expected = exact_est,
                 tolerance = 0.1)
  }
)

test_that(
  desc = "Helpful message if `max_replicates` is too low.", {

    # Disable test until the Matrix package adapts to change in R 4.3.2
    # (`crossprod()` becomes primitive and S3 generic in R 4.3.2)
    skip_if_not("crossprod" %in% ls(getNamespaceInfo("Matrix", "exports")))

    expect_message(
      regexp = "number of replicates needed",
      object = {
        svydesign(
          data = election_pps,
          id = ~1, fpc = ~p,
          pps = ppsmat(election_jointprob),
          variance = "YG"
        ) |> as_fays_gen_rep_design(
          variance_estimator = "Yates-Grundy",
          max_replicates = 5
        )
      }
    )


    expect_warning(
      regexp = "not recommended",
      object = {
        suppressMessages({
          svydesign(
            data = election_pps,
            id = ~1, fpc = ~p,
            pps = ppsmat(election_jointprob),
            variance = "YG"
          ) |> as_fays_gen_rep_design(
            variance_estimator = "Yates-Grundy",
            max_replicates = 5,
            balanced = FALSE
          )
        })
      }
    )

  }
)

# Warnings for ill-advised choices -----

test_that(
  desc = "Warning when `mse = FALSE`", {
    expect_warning({
    twophase_design$phase1$full |> 
      as_fays_gen_rep_design("Ultimate Cluster", mse = FALSE) |>
      svytotal(x = ~ y1)
    }, regexp = "may produce large underestimates")
  }
)

# Works for more specialized classes of survey designs ----

test_that(
  desc = "Returns `tbl_svy` if the input is a `tbl_svy` and 'srvyr' is loaded", {
    library(srvyr)
    expect_true(
      suppressMessages({
        twophase_design |> as_survey() |>
          as_fays_gen_rep_design(variance_estimator = list(
            'Ultimate Cluster', 'Ultimate Cluster'
          ), max_replicates = 1) |>
          inherits('tbl_svy')
      })
    )
  }
)

# Test use of 'torch' for matrix decomposition ----

test_that(
  desc = "Can use 'torch' for Fay's generalized replication", {
    skip_if_not_installed("torch")
    skip_if_not(torch::torch_is_installed())
    skip_on_cran()
    
    non_torch_result <- svydesign(
      data = election_pps,
      id = ~1, fpc = ~p,
      pps = ppsmat(election_jointprob),
      variance = "YG"
    ) |> svytotal(x = ~ Bush + Kerry) |> SE()
    
    torch_result <- withr::with_options(list(svrep.torch_device = "cpu"), {
      svydesign(
        data = election_pps,
        id = ~1, fpc = ~p,
        pps = ppsmat(election_jointprob),
        variance = "YG"
      ) |> as_fays_gen_rep_design(
        variance_estimator = "Yates-Grundy",
        max_replicates = 100
      ) |> svytotal(x = ~ Bush + Kerry) |> SE()
    })
    
    expect_equal(unname(torch_result), unname(non_torch_result),
                 tolerance = 0.1)
    
  }
)