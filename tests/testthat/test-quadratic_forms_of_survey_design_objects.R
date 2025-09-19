suppressWarnings({
  suppressPackageStartupMessages({
    library(survey)
    library(dplyr)
    library(svrep)
    library(testthat)
  })
})

# Tests for two-phase designs ----
  data('mu284', package = 'survey')

  mu284$phase2 <- mu284$id1 %in% c(19,31,45)
  mu284 <- mu284[order(mu284$id1),]

  twophase_design <- twophase(
    data = mu284,
    id = list(~id1, ~id1),
    fpc = list(~n1, NULL),
    subset = ~ phase2
  )

  test_that("Expected result for two-phase design", {

    suppressMessages({
      expect_warning(
        regexp = "Approximating \\(sigma_1\\/phase_2_joint_probs\\) ", object = {
          design_quad_form <- twophase_design |> get_design_quad_form(
            variance_estimator = list('SD2', 'Ultimate Cluster'),
            ensure_psd = TRUE
          )
        }
      )
    })

    expect_silent({
      suppressMessages({
        twophase_design |> get_design_quad_form(
          variance_estimator = list('SD2', 'Ultimate Cluster'),
          ensure_psd = FALSE
        )
      })
    })

    expect_warning(
      regexp = "Approximating \\(sigma_1\\/phase_2_joint_probs\\) ", object = {
        twophase_quad_form <- make_twophase_quad_form(
          sigma_1 = make_quad_form_matrix(
            variance_estimator = "SD2",
            cluster_ids = mu284$id1 |> as.matrix(),
            strata_ids = rep(1, times = nrow(mu284)) |> as.matrix(),
            strata_pop_sizes = mu284$n1 |> as.matrix(),
            sort_order = 1:15
          )[mu284$phase2, mu284$phase2],
          sigma_2 = make_quad_form_matrix(
            variance_estimator = "Ultimate Cluster",
            cluster_ids = c(1,1,1,2,2,2,3,3,3) |> as.matrix(),
            strata_ids = rep(1, times = 9) |> as.matrix(),
            strata_pop_sizes = rep(5, times = 9) |> as.matrix()
          ),
          phase_2_joint_probs = svrep:::make_srswor_matrix(n = 3, f = 3/5) |>
            svrep:::ht_matrix_to_joint_probs() |>
            svrep:::distribute_matrix_across_clusters(
              cluster_ids = c(1,1,1,2,2,2,3,3,3)
            )
        )
      })

    expect_equal(
      object = design_quad_form,
      expected = twophase_quad_form
    )

  })

# Tests for PPS designs ----

  data('election', package = 'survey')
  yg_design <- svydesign(
    id = ~ 1,
    fpc = ~ p,
    data = election_pps,
    pps = ppsmat(election_jointprob),
    variance = "YG"
  )
  ht_design <- svydesign(
    id = ~ 1,
    fpc = ~ p,
    data = election_pps,
    pps = ppsmat(election_jointprob),
    variance = "HT"
  )

  test_that(desc = "Expected results for PPS designs", {

    expect_error(
      regexp = "must use a PPS design",
      object = {
        get_design_quad_form(
          twophase_design$phase1$sample,
          variance_estimator = "Horvitz-Thompson"
        )
      }
    )

    expect_error(
      regexp = "Must specify `variance='HT'",
      object = {
        get_design_quad_form(
          yg_design,
          variance_estimator = "Horvitz-Thompson"
        )
      }
    )
    expect_error(
      regexp = "Must specify `variance='YG'",
      object = {
        get_design_quad_form(
          ht_design,
          variance_estimator = "Yates-Grundy"
        )
      }
    )


    expect_equal(
      object = get_design_quad_form(
        yg_design,
        variance_estimator = "Yates-Grundy"
      ) |> as.matrix() |> `dimnames<-`(NULL),
      expected = make_quad_form_matrix(
        variance_estimator = "Yates-Grundy",
        joint_probs = election_jointprob
      ) |> as.matrix()
    )
    expect_equal(
      object = get_design_quad_form(
        ht_design,
        variance_estimator = "Horvitz-Thompson"
      ) |> as.matrix() |> `dimnames<-`(NULL),
      expected = make_quad_form_matrix(
        variance_estimator = "Horvitz-Thompson",
        joint_probs = election_jointprob
      ) |> as.matrix()
    )
    expect_equal(
      object = get_design_quad_form(
        ht_design,
        variance_estimator = "Poisson Horvitz-Thompson"
      ) |> as.matrix() |> `dimnames<-`(NULL),
      expected = make_quad_form_matrix(
        variance_estimator = "Horvitz-Thompson",
        joint_probs = outer(diag(election_jointprob),
                            diag(election_jointprob)) |>
          `diag<-`(diag(election_jointprob))
      ) |> as.matrix()
    )
  })

# Test for Deville-Tille estimator ----

  test_that("Correct quadratic form for Deville-Tille estimator", {

    data('api', package = 'survey')

    aux_var_columns <- model.matrix(object = ~ -1 + stype, data = apistrat)

    apistrat <- cbind(apistrat, aux_var_columns)

    ## Create survey design object
    multistage_survey_design <- svydesign(
      data = apistrat,
      weights = ~ pw,
      ids = ~ 1
    )

    expect_equal(
      object = get_design_quad_form(
        design = multistage_survey_design,
        variance_estimator = "Deville-Tille",
        aux_var_names = "stype"
      ) |> as.matrix(),
      expected = make_deville_tille_matrix(
        probs = multistage_survey_design$variables$pw^(-1),
        aux_vars = aux_var_columns
      )
    )

  })

# Test for one of the Deville variance estimators ----

  test_that("Correct quadratic form for a Deville variance estimator", {
      dev1_quad_form <- make_quad_form_matrix(
        variance_estimator = "Deville-1",
        cluster_ids = c(3,3,1,1,2,2) |> as.matrix(),
        strata_ids = c(1,1,1,1,1,1) |> as.matrix(),
        probs = c(0.258064516129032, 0.258064516129032, 0.129032258064516,
                  0.129032258064516, 0.193548387096774, 0.193548387096774) |>
          as.matrix()
      )
    
      expect_equal(
        object = svydesign(
          data = data.frame(
            psu = c(3,3,1,1,2,2),
            stratum = c(1,1,1,1,1,1),
            samp_prob = c(0.258064516129032, 0.258064516129032, 0.129032258064516,
                          0.129032258064516, 0.193548387096774, 0.193548387096774)
          ),
          ids = ~ psu,
          strata = ~ stratum,
          probs = ~ samp_prob
        ) |> get_design_quad_form("Deville-1"),
        expected = dev1_quad_form
      )
  })

# Test the BOSB variance estimator ----

  test_that("Correct quadratic form for the BOSB variance estimator", {
    
    stsys_strata <- c(rep(1, times = 100),
                      rep(2, times = 119))

    library_stsys_sample[['strata']] <- stsys_strata
    library_stsys_sample[['weights']] <- rep(1, times = nrow(library_stsys_sample))
    
    bosb_quad_form <- make_quad_form_matrix(
      variance_estimator = 'BOSB',
      cluster_ids = library_stsys_sample[,'FSCSKEY',drop=FALSE],
      strata_ids = matrix(stsys_strata,
                          nrow = nrow(library_stsys_sample),
                          ncol = 1),
      aux_vars = library_stsys_sample[['SAMPLING_SORT_ORDER']] |>
        factor() |> as.numeric() |> matrix(nrow = nrow(library_stsys_sample))
    )
    expect_equal(
      object = svydesign(
        data   = library_stsys_sample |>
          transform(
            AUX_VAR = SAMPLING_SORT_ORDER |>
              factor() |> as.numeric()
          ),
        ids     = ~ FSCSKEY,
        strata  = ~ strata,
        weights = ~ weights
      ) |> get_design_quad_form("BOSB", aux_var_names = "AUX_VAR"),
      expected = bosb_quad_form
    )
  })

# Informative messages for bad inputs ----

  test_that(
    "Informative error messages for bad inputs", {
      expect_error(
        object = {get_design_quad_form(ht_design, variance_estimator = NULL)},
        regexp = "Must specify a value"
      )
      expect_error(
        object = {get_design_quad_form(ht_design,
                                       variance_estimator = list(
                                         'SD1', 'SD2'
                                       ))},
        regexp = "Can only specify one"
      )
      expect_error(
        object = {get_design_quad_form(ht_design, "made-up")},
        regexp = "`made-up` is not a supported variance estimator"
      )
  })
