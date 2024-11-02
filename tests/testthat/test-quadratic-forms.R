suppressWarnings({
  suppressPackageStartupMessages({
    library(survey)
    library(dplyr)
    library(svrep)
    library(testthat)
  })
})

options(survey.lonely.psu = 'certainty')

data('library_stsys_sample', package = 'svrep')
set.seed(2014)

library_stsys_sample <- library_stsys_sample |>
  mutate(
    TOTCIR = ifelse(is.na(TOTCIR), 0, TOTCIR),
    TOTSTAFF = ifelse(is.na(TOTSTAFF), 0, TOTSTAFF)
  )

# Check basic successive-difference quadratic forms ----

  successive_diffs <- function(y, type = "SD2") {
    n <- length(y)
    sum_of_squares <- sum(diff(y)^2)
    if (n == 1) {
      result <- 0
    }
    if (type == "SD1" && (n > 1)) {
      result <- sum_of_squares * (n/(2*(n-1)))
    }
    if (type == "SD2" && (n > 1)) {
      sum_of_squares <- sum_of_squares + (tail(y, 1) - head(y, 1))^2
      result <- sum_of_squares/2
    }
    return(result)
  }

  wtd_y <- library_stsys_sample[['TOTCIR']]/library_stsys_sample[['SAMPLING_PROB']]
  wtd_y <- wtd_y[order(library_stsys_sample[['SAMPLING_SORT_ORDER']])]

  test_that(
    "Helper function for successive difference quadratic forms works correctly", {
      # Basic results correct
      expect_equal(
        object = as.matrix(t(wtd_y) %*% svrep:::make_sd_matrix(n = length(wtd_y), type = "SD2") %*% wtd_y),
        expected = successive_diffs(wtd_y, type = "SD2") |> as.matrix()
      )
      expect_equal(
        object = as.matrix(t(wtd_y) %*% svrep:::make_sd_matrix(n = length(wtd_y), type = "SD1") %*% wtd_y),
        expected = successive_diffs(wtd_y, type = "SD1") |> as.matrix()
      )
      # Checks on FPC
      expect_equal(
        object = as.matrix(t(wtd_y) %*% svrep:::make_sd_matrix(n = length(wtd_y), type = "SD1",
                                                     f = 0.9) %*% wtd_y),
        expected = (1-0.9) * successive_diffs(wtd_y, type = "SD1") |> as.matrix()
      )
      expect_equal(
        object = as.matrix(t(wtd_y) %*% svrep:::make_sd_matrix(n = length(wtd_y), type = "SD2",
                                                     f = 0.9) %*% wtd_y),
        expected = (1-0.9) * successive_diffs(wtd_y, type = "SD2") |> as.matrix()
      )
      # Correct result for only a single unit
      expect_equal(
        object = svrep:::make_sd_matrix(n = 1) |> as.matrix(),
        expected = matrix(0, nrow = 1, ncol = 1)
      )
      # Correct result for only two units
      expect_equal(
        object = as.matrix(c(99, 16) %*% svrep:::make_sd_matrix(n = 2, type = "SD1") %*% c(99,16)),
        expected = as.matrix(successive_diffs(c(99,16), type = "SD1"))
      )
      expect_equal(
        object = as.matrix(c(99, 16) %*% svrep:::make_sd_matrix(n = 2, type = "SD2") %*% c(99,16)),
        expected = as.matrix(successive_diffs(c(99,16), type = "SD2"))
      )
      # Checks on numeric arguments
      expect_error(
        object = svrep:::make_sd_matrix(n = 0),
        regexp = "must be an integer greater than"
      )
      expect_error(
        object = svrep:::make_sd_matrix(n = 2, f = 1.1),
        regexp = "must be a single number between"
      )
    })

# Check basic SRSWOR quadratic form ----

  wtd_y <- library_stsys_sample[['TOTCIR']]/library_stsys_sample[['SAMPLING_PROB']]
  n <- length(wtd_y)

  wtd_y_matrix <- as.matrix(
    library_stsys_sample[,c("TOTCIR", "TOTSTAFF")]/library_stsys_sample[['SAMPLING_PROB']]
  )

  test_that(
    "Helper function for SRSWOR quadratic forms works correctly", {
      # Basic results correct
      expect_equal(
        object = as.matrix(t(wtd_y) %*% svrep:::make_srswor_matrix(n = length(wtd_y)) %*% wtd_y),
        expected = n * cov(as.matrix(wtd_y))
      )
      # Matches 'survey' package
      expect_equal(
        object = as.matrix(t(wtd_y_matrix) %*% svrep:::make_srswor_matrix(n = length(wtd_y)) %*% wtd_y_matrix),
        expected = svyrecvar(
          x = wtd_y_matrix,
          clusters = as.matrix(seq_len(n)),
          stratas = as.matrix(rep(1, times = n)),
          fpcs = list('sampsize' = as.matrix(rep(n, times = n)),
                      'popsize' = as.matrix(rep(Inf, times = n))),
          lonely.psu = "remove", one.stage = TRUE
        )
      )
      # Checks on FPC
      expect_equal(
        object = as.matrix(t(wtd_y) %*% svrep:::make_srswor_matrix(n = length(wtd_y),
                                                         f = 0.9) %*% wtd_y),
        expected = (1-0.9) * n *  cov(as.matrix(wtd_y)) |> as.matrix()
      )
      # Correct result for only a single unit
      expect_equal(
        object = svrep:::make_srswor_matrix(n = 1) |> as.matrix(),
        expected = matrix(0, nrow = 1, ncol = 1)
      )
      # Checks on numeric arguments
      expect_error(
        object = svrep:::make_srswor_matrix(n = 0),
        regexp = "must be an integer greater than"
      )
      expect_error(
        object = svrep:::make_srswor_matrix(n = 2, f = 1.1),
        regexp = "must be a single number between"
      )
    })

# Check Horvitz-Thompson results ----

  ## Load an example dataset that uses unequal probability sampling
  data('election', package = 'survey')
  y_wtd <- election_pps$Bush/diag(election_jointprob)

  ## Create matrix to represent the Horvitz-Thompson estimator as a quadratic form
  ## This is the correct result that we expect the 'svrep' package to match
  n <- nrow(election_pps)
  pi <- election_jointprob
  horvitz_thompson_matrix <- matrix(nrow = n, ncol = n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      horvitz_thompson_matrix[i,j] <- 1 - (pi[i,i] * pi[j,j])/pi[i,j]
    }
  }

  ## Get the quadratic form using the 'svrep' pacakge
  ht_quad_form <- make_quad_form_matrix(
    variance_estimator = "Horvitz-Thompson",
    joint_probs = election_jointprob
  )

  ## Compare quadratic form matrix from 'svrep' to the correct result
  test_that(
    "Generates correct quadratic form for Horvitz-Thompson", {
     expect_equal(object = as.matrix(ht_quad_form), expected = horvitz_thompson_matrix)
  })

  ## Use the 'svrep' package to obtain the Horvitz-Thompson estimate
  quad_form_result <- t(y_wtd) %*% ht_quad_form %*% y_wtd

  test_that(
    "Horvitz-Thompson estimate matches the 'survey' package", {
      # Disable test until the Matrix package adapts to change in R 4.3.2
      # (`crossprod()` becomes primitive and S3 generic in R 4.3.2)
      skip_if_not("crossprod" %in% ls(getNamespaceInfo("Matrix", "exports")))

      ## Use the 'survey' package to obtain the Horvitz-Thompson estimate
      svy_pkg_result <- svydesign(
        data = election_pps,
        id = ~1, fpc = ~p,
        pps = ppsmat(election_jointprob),
        variance = "HT"
      ) |> svytotal(x = ~ Bush) |> vcov()

    expect_equal(object = as.matrix(quad_form_result), expected = svy_pkg_result)
  })

# Check Sen-Yates-Grundy results ----

  y_wtd <- election_pps$Bush/diag(election_jointprob)

  syg_result <- as.matrix(0)
  for (i in seq_along(y_wtd)) {
    for (j in seq_along(y_wtd)) {
      pi_ij <- election_jointprob[i,j]
      pi_i <- election_jointprob[i,i]
      pi_j <- election_jointprob[j,j]
      contribution <- -0.5*(1 - (pi_i*pi_j)/pi_ij) * (y_wtd[i] - y_wtd[j])^2
      syg_result <- syg_result + contribution
    }
  }

  yg_quad_form_matrix <- make_quad_form_matrix(variance_estimator = "Yates-Grundy",
                                               joint_probs = election_jointprob)

  quad_form_result <-  as.matrix(t(y_wtd) %*% yg_quad_form_matrix %*% y_wtd)

  test_that(
    "Generates correct quadratic form for Yates-Grundy", {
    # Disable test until the Matrix package adapts to change in R 4.3.2
    # (`crossprod()` becomes primitive and S3 generic in R 4.3.2)
    skip_if_not("crossprod" %in% ls(getNamespaceInfo("Matrix", "exports")))

    svy_pkg_result <- svydesign(
      data = election_pps,
      id = ~1, fpc = ~p,
      pps = ppsmat(election_jointprob),
      variance = "YG"
    ) |> svytotal(x = ~ Bush) |> vcov()

    expect_equal(object = quad_form_result, expected = syg_result)
    expect_equal(object = quad_form_result, expected = svy_pkg_result)
  })



# Check SD1 and SD2 results ----

  ##_ Basic correctness checks
  sd1_quad_form <- make_quad_form_matrix(
    variance_estimator = 'SD1',
    cluster_ids = library_stsys_sample[,'FSCSKEY',drop=FALSE],
    strata_ids = library_stsys_sample[,'SAMPLING_STRATUM',drop=FALSE],
    strata_pop_sizes = library_stsys_sample[,'STRATUM_POP_SIZE',drop=FALSE],
    sort_order = library_stsys_sample[['SAMPLING_SORT_ORDER']]
  )
  sd2_quad_form <- make_quad_form_matrix(
    variance_estimator = 'SD2',
    cluster_ids = library_stsys_sample[,'FSCSKEY',drop=FALSE],
    strata_ids = library_stsys_sample[,'SAMPLING_STRATUM',drop=FALSE],
    strata_pop_sizes = library_stsys_sample[,'STRATUM_POP_SIZE',drop=FALSE],
    sort_order = library_stsys_sample[['SAMPLING_SORT_ORDER']]
  )

  wtd_y <- as.matrix(library_stsys_sample[['TOTCIR']] /
                      library_stsys_sample$SAMPLING_PROB)

  expected_results <- library_stsys_sample |>
    arrange(SAMPLING_SORT_ORDER) |>
    group_by(SAMPLING_STRATUM) |>
    summarize(
      v_h_sd1 = (1 - unique(SAMPLING_PROB)) * successive_diffs(y = TOTCIR/SAMPLING_PROB, type = "SD1"),
      v_h_sd2 = (1 - unique(SAMPLING_PROB)) * successive_diffs(y = TOTCIR/SAMPLING_PROB, type = "SD2")
    ) |>
    summarize(
      v_sd1 = sum(v_h_sd1), v_sd2 = sum(v_h_sd2)
    ) |> unlist()

  test_that(
    "`make_quad_form_matrix()` works correctly for 'SD1' and 'SD2' with stratified single-stage samples", {
      expect_equal(
        object = as.numeric(t(wtd_y) %*% sd1_quad_form %*% wtd_y),
        expected = unname(expected_results['v_sd1'])
      )
      expect_equal(
        object = as.numeric(t(wtd_y) %*% sd2_quad_form %*% wtd_y),
        expected = unname(expected_results['v_sd2'])
      )
    })

  ##_ Check that works if rows are not in sampling sort order

  shuffled_data <- library_stsys_sample[sample(nrow(library_stsys_sample),
                                               size = nrow(library_stsys_sample)),]
  sorted_data <- shuffled_data |>
    arrange(SAMPLING_SORT_ORDER)

  test_that(
    "`make_quad_form_matrix()` with 'SD1' and 'SD2' works correctly when data frame isn't sorted", {

      unsorted_quad_form_matrix <- make_quad_form_matrix(
        variance_estimator = 'SD1',
        cluster_ids = shuffled_data[,'FSCSKEY',drop=FALSE],
        strata_ids = shuffled_data[,'SAMPLING_STRATUM',drop=FALSE],
        strata_pop_sizes = shuffled_data[,'STRATUM_POP_SIZE',drop=FALSE],
        sort_order = shuffled_data[['SAMPLING_SORT_ORDER']]
      )

      sorted_quad_form_matrix <- make_quad_form_matrix(
        variance_estimator = 'SD1',
        cluster_ids = sorted_data[,'FSCSKEY',drop=FALSE],
        strata_ids = sorted_data[,'SAMPLING_STRATUM',drop=FALSE],
        strata_pop_sizes = sorted_data[,'STRATUM_POP_SIZE',drop=FALSE],
        sort_order = sorted_data[['SAMPLING_SORT_ORDER']]
      )

      y_wtd_sorted <- sorted_data$TOTCIR/sorted_data$SAMPLING_PROB
      y_wtd_unsorted <- shuffled_data$TOTCIR/shuffled_data$SAMPLING_PROB

      sorted_result <- t(y_wtd_sorted) %*% sorted_quad_form_matrix %*% y_wtd_sorted
      unsorted_result <- t(y_wtd_unsorted) %*% unsorted_quad_form_matrix %*% y_wtd_unsorted

      expect_equal(
        object = sorted_result,
        expected = unsorted_result
      )

    })

  test_that(
    "successive-differences estimate less than ST-SRSWOR approximation, for sys sample from sorted frame", {

      expect_lt(
        object = as.numeric(t(wtd_y) %*% sd1_quad_form %*% wtd_y),
        expected = svydesign(
          data = library_stsys_sample,
          strata = ~ SAMPLING_STRATUM,
          fpc = ~ SAMPLING_PROB,
          ids = ~ 1
        ) |> svytotal(x = ~TOTCIR) |> vcov()
      )
  })
  
# Check BOSB results ----
  
  ##_ Basic correctness checks
  bosb_quad_form <- make_quad_form_matrix(
    variance_estimator = 'BOSB',
    cluster_ids = library_stsys_sample[,'FSCSKEY',drop=FALSE],
    strata_ids = matrix(c(rep(1, times = 100),
                          rep(2, times = 119)),
                        nrow = nrow(library_stsys_sample), ncol = 1),
    aux_vars = library_stsys_sample[['SAMPLING_SORT_ORDER']] |>
      factor() |> as.numeric() |> matrix(nrow = nrow(library_stsys_sample))
  )
  
  stratum_1_quad_form <- make_kernel_var_matrix(
    x = library_stsys_sample[['SAMPLING_SORT_ORDER']][1:100] |>
      factor() |> as.numeric(),
    kernel = "Epanechnikov", bandwidth = "auto"
  )
  stratum_2_quad_form <- make_kernel_var_matrix(
    x = library_stsys_sample[['SAMPLING_SORT_ORDER']][101:219] |>
      factor() |> as.numeric(),
    kernel = "Epanechnikov", bandwidth = "auto"
  )
  
  wtd_y <- as.matrix(library_stsys_sample[['TOTCIR']] /
                       library_stsys_sample$SAMPLING_PROB)
  
  test_that(
    "`make_quad_form_matrix()` works correctly for 'BOSB' with stratified single-stage samples", {
      
      actual_estimate <- as.numeric(t(wtd_y) %*% bosb_quad_form %*% wtd_y)
      
      expected_estimate <- sum(c(
        as.numeric(t(wtd_y[1:100])   %*% stratum_1_quad_form %*% wtd_y[1:100]),
        as.numeric(t(wtd_y[101:219]) %*% stratum_2_quad_form %*% wtd_y[101:219])
      ))
      
      expect_equal(
        object   = actual_estimate,
        expected = expected_estimate
      )
    })
  
  test_that(
    "Expected errors and warnings for 'BOSB'", {
      expect_error(
        object = {
          make_quad_form_matrix(
            variance_estimator = 'BOSB',
            strata_ids = matrix(c(1,2), 2, 1)
          )
        },
        regexp = "must supply.+strata.+cluster"
      )
      
      expect_error(
        object = {
          make_quad_form_matrix(
            variance_estimator = 'BOSB',
            strata_ids = matrix(c(1,1), 2, 1),
            cluster_ids = matrix(c(1,2), 2, 1)
          )
        },
        regexp = "must supply.+aux_vars"
      )
    }
  )
  
  test_that(
    "Correct quadratic forms for 'BOSB' in specific examples", {
      
      expect_equal(
        object   = make_kernel_var_matrix(c(1,1,100,100)) |> as.matrix(),
        expected = matrix(c( 1,-1, 0, 0,
                            -1, 1, 0, 0,
                             0, 0, 1,-1,
                             0, 0,-1, 1),
                          nrow = 4, ncol = 4,
                          byrow = TRUE)
      )
      
      expect_equal(
        object = make_kernel_var_matrix(100) |> as.matrix(),
        expected = matrix(0, 1, 1)
      )
      
      expect_equal(
        object   = make_kernel_var_matrix(c(1,100)) |> as.matrix(),
        expected = matrix(c( 1,-1, 
                            -1, 1),
                          nrow = 2, ncol = 2,
                          byrow = TRUE)
      )
      
      expect_equal(
        object = make_kernel_var_matrix(c(1,1,100,80))[1:2,3:4] |> as.matrix(),
        expected = matrix(0, nrow = 2, ncol = 2)
      )
      
      expect_equal(
        object   = make_kernel_var_matrix(c(1,1,2,2,3)/3) |> attr('bandwidth'),
        expected = 2/3
      )
      
      expect_equal(
        object   = make_kernel_var_matrix(c(1,1,2,2,4)) |> attr('bandwidth'),
        expected = 3/4 * 4
      )
      
      expect_equal(
        object = make_kernel_var_matrix(c(102, 103, 104)/104) |> is.nan() |> sum(),
        expected = 0
      )
      
    }
  )

# Check "Ultimate Cluster" results ----

  test_that(
    "Correct results for ultimate cluster", {

      # Correct result for a single-stage stratified sample
      wtd_y_matrix <- as.matrix(
        library_stsys_sample[,c("TOTCIR", "TOTSTAFF")]
      )/library_stsys_sample[['SAMPLING_PROB']]
      quad_UC <- make_quad_form_matrix(
        variance_estimator = "Ultimate Cluster",
        cluster_ids = library_stsys_sample[,'FSCSKEY',drop=FALSE],
        strata_ids = library_stsys_sample[,'SAMPLING_STRATUM',drop=FALSE],
        strata_pop_sizes = library_stsys_sample[,'STRATUM_POP_SIZE',drop=FALSE]
      )
      expect_equal(
        object = as.matrix(t(wtd_y_matrix) %*% quad_UC %*% wtd_y_matrix),
        expected = svydesign(
          data = library_stsys_sample,
          ids = ~ FSCSKEY,
          strata = ~ SAMPLING_STRATUM,
          fpc = ~ STRATUM_POP_SIZE
        ) |> svytotal(x = ~ TOTCIR + TOTSTAFF, na.rm = TRUE) |>
          vcov()
      )

      # Correct result for a multistage cluster sample
      wtd_y_matrix <- as.matrix(
        library_multistage_sample[,c("TOTCIR", "TOTSTAFF")]
      )/library_multistage_sample[['SAMPLING_PROB']]
      wtd_y_matrix[rowSums(is.na(wtd_y_matrix)) > 0] <- 0
      quad_UC <- make_quad_form_matrix(
        variance_estimator = "Ultimate Cluster",
        cluster_ids = library_multistage_sample[,c("PSU_ID", "SSU_ID"),drop=FALSE],
        strata_ids = data.frame('PSU_STRATUM' = rep(1, times = nrow(library_multistage_sample)),
                                'SSU_STRATUM' = rep(1, times = nrow(library_multistage_sample))),
        strata_pop_sizes = library_multistage_sample[,c("PSU_POP_SIZE", "SSU_POP_SIZE"),drop=FALSE]
      )
      expect_equal(
        object = as.matrix(t(wtd_y_matrix) %*% quad_UC %*% wtd_y_matrix),
        expected = svydesign(
          data = library_multistage_sample |>
            mutate(SAMPLING_WEIGHT = 1/SAMPLING_PROB),
          weights = ~ SAMPLING_WEIGHT,
          ids = ~ PSU_ID,
          strata = NULL,
          fpc = ~ PSU_POP_SIZE
        ) |> svytotal(x = ~ TOTCIR + TOTSTAFF, na.rm = TRUE) |>
          vcov()
      )
  })

  test_that(
    "Fills in implicit cluster or strata IDs", {
      wtd_y_matrix <- as.matrix(
        library_stsys_sample[,c("TOTCIR", "TOTSTAFF")]
      )/library_stsys_sample[['SAMPLING_PROB']]
      qf <- make_quad_form_matrix(
        variance_estimator = "SD1",
        cluster_ids = library_stsys_sample[,'FSCSKEY',drop=FALSE],
        #strata_ids = library_stsys_sample[,'SAMPLING_STRATUM',drop=FALSE],
        strata_pop_sizes = library_stsys_sample |> transmute(N = 50000),
        sort_order = library_stsys_sample$SAMPLING_SORT_ORDER
      )
      expect_equal(object = dim(qf),
                   expected = c(219, 219))
  })

# Check "Stratified Multistage SRS" results ----

  test_that(
    "Correct results for stratified multistage SRS", {

      # Correct result for a single-stage stratified sample
      wtd_y_matrix <- as.matrix(
        library_stsys_sample[,c("TOTCIR", "TOTSTAFF")]
      )/library_stsys_sample[['SAMPLING_PROB']]
      quad_UC <- make_quad_form_matrix(
        variance_estimator = "Stratified Multistage SRS",
        cluster_ids = library_stsys_sample[,'FSCSKEY',drop=FALSE],
        strata_ids = library_stsys_sample[,'SAMPLING_STRATUM',drop=FALSE],
        strata_pop_sizes = library_stsys_sample[,'STRATUM_POP_SIZE',drop=FALSE]
      )
      expect_equal(
        object = as.matrix(t(wtd_y_matrix) %*% quad_UC %*% wtd_y_matrix),
        expected = svydesign(
          data = library_stsys_sample,
          ids = ~ FSCSKEY,
          strata = ~ SAMPLING_STRATUM,
          fpc = ~ STRATUM_POP_SIZE
        ) |> svytotal(x = ~ TOTCIR + TOTSTAFF, na.rm = TRUE) |>
          vcov()
      )

      # Correct result for a multistage cluster sample
      data('mu284', package = 'survey')
      mu284 <- mu284 |> arrange(id1, id2)
      mu284_design <- survey::svydesign(
        data = mu284,
        ids = ~ id1 + id2,
        fpc = ~ n1 + n2
      )
      wtd_y_matrix <- as.matrix(
        mu284_design$variables$y1
      ) * weights(mu284_design)
      wtd_y_matrix[rowSums(is.na(wtd_y_matrix)) > 0] <- 0
      colnames(wtd_y_matrix) <- c("y1")

      # Use function to produce quadratic form matrix
      qf_matrix <- make_quad_form_matrix(
        variance_estimator = "Stratified Multistage SRS",
        cluster_ids = mu284_design$cluster,
        strata_ids = mu284_design$strata,
        strata_pop_sizes = mu284_design$fpc$popsize
      )

      # Determine expected correct quadratic form matrix
      quad_form_UC <- make_quad_form_matrix(
        variance_estimator = "Ultimate Cluster",
        cluster_ids = mu284_design$cluster,
        strata_ids = mu284_design$strata,
        strata_pop_sizes = mu284_design$fpc$popsize
      )
      later_stage_quad_form <- make_quad_form_matrix(
        variance_estimator = "Ultimate Cluster",
        cluster_ids = mu284_design$cluster[,2,drop=FALSE],
        strata_ids = mu284_design$strata[,2,drop=FALSE],
        strata_pop_sizes = mu284_design$fpc$popsize[,2,drop=FALSE]
      )
      expected_quad_form <- quad_form_UC + (0.1 * later_stage_quad_form)

      # Compare to survey package
      expect_equal(
        object = as.matrix(t(wtd_y_matrix) %*% qf_matrix %*% wtd_y_matrix),
        expected = mu284_design |>
          svytotal(x = ~ y1, na.rm = TRUE) |>
          vcov()
      )
    })
  
# Check "Beaumont-Emond" results ----
  
  test_that(
    "Helper function for Beaumont-Emond quadratic forms works correctly", {
      
      data('api', package = 'survey')
      
      # For a single-stage stratified SRSWOR sample,
      # results match usual variance estimator
      wtd_y_matrix <- as.matrix(
        library_stsys_sample[,c("TOTCIR", "TOTSTAFF")]
      )/library_stsys_sample[['SAMPLING_PROB']]
      quad_BE <- make_quad_form_matrix(
        variance_estimator = "Beaumont-Emond",
        cluster_ids = library_stsys_sample[,'FSCSKEY',drop=FALSE],
        strata_ids = library_stsys_sample[,'SAMPLING_STRATUM',drop=FALSE],
        probs = library_stsys_sample[,'SAMPLING_PROB',drop=FALSE]
      )
      expect_equal(
        object = as.matrix(t(wtd_y_matrix) %*% quad_BE %*% wtd_y_matrix),
        expected = svydesign(
          data = library_stsys_sample,
          ids = ~ FSCSKEY,
          strata = ~ SAMPLING_STRATUM,
          fpc = ~ STRATUM_POP_SIZE
        ) |> svytotal(x = ~ TOTCIR + TOTSTAFF, na.rm = TRUE) |>
          vcov()
      )
      
      # For a stratified multistage SRSWOR sample,
      # results match the usual variance estimator
      wtd_y_matrix <- as.matrix(apiclus2[,c('api00', 'api99')] * apiclus2$pw)
      apiclus2_design <- svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2)
      quad_BE <- make_quad_form_matrix(
        variance_estimator = "Beaumont-Emond",
        cluster_ids = apiclus2[,c("dnum", "snum"),drop=FALSE],
        strata_ids = data.frame('PSU_STRATUM' = rep(1, times = nrow(apiclus2)),
                                'SSU_STRATUM' = rep(1, times = nrow(apiclus2))),
        probs = apiclus2_design$allprob
      )
      expect_equal(
        object = as.matrix(t(wtd_y_matrix) %*% quad_BE %*% wtd_y_matrix),
        expected = svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2) |> 
          svytotal(x = ~ api00 + api99, na.rm = TRUE) |>
          vcov()
      )

      # Correct result for only a single unit
      expect_equal(
        object = svrep:::make_ppswor_approx_matrix(probs = 0.5) |> as.matrix(),
        expected = matrix(0, nrow = 1, ncol = 1)
      )
    })

# Check "Deville-1" and "Deville-2" results ----

  test_that(
    "Correct results for Deville-1 and Deville-2", {

      # Basic correct form for Deville-1

      dev1_quad_form <- svrep:::make_ppswor_approx_matrix(
        probs = c(1/3, 1/2, 1/4),
        method = "Deville-1"
      )

      c_i <- (1 - c(1/3, 1/2, 1/4)) * (3/2)
      exp_dev1_quad_form <- outer(
        -c_i, c_i
      ) / sum(c_i)
      diag(exp_dev1_quad_form) <- c_i * (1 - c_i/sum(c_i))

      expect_equal(object = as.matrix(dev1_quad_form), expected = exp_dev1_quad_form)

      # Basic correct form for Deville-2

      dev2_quad_form <- svrep:::make_ppswor_approx_matrix(
        probs = c(1/3, 1/2, 1/4),
        method = "Deville-2"
      )
      pi_i <- c(1/3, 1/2, 1/4)
      c_i <- (1 - pi_i) / (
        1 - sum(((1-pi_i)/sum(1-pi_i))^2)
      )
      exp_dev2_quad_form <- outer(
        -c_i, c_i
      ) / sum(c_i)
      diag(exp_dev2_quad_form) <- c_i * (1 - c_i/sum(c_i))

      expect_equal(object = as.matrix(dev2_quad_form), expected = exp_dev2_quad_form)

      # Works as expected for clusters
      dev1_quad_form <- make_quad_form_matrix(
        variance_estimator = "Deville-1",
        cluster_ids = c(3,3,1,1,2,2) |> as.matrix(),
        strata_ids = c(1,1,1,1,1,1) |> as.matrix(),
        probs = c(0.258064516129032, 0.258064516129032, 0.129032258064516,
                  0.129032258064516, 0.193548387096774, 0.193548387096774) |>
          as.matrix()
      )
      expect_equal(
        object = as.matrix(dev1_quad_form),
        expected = svrep:::make_ppswor_approx_matrix(
          probs = c(0.258064516129032, 0.129032258064516, 0.193548387096774),
          method = "Deville-1"
        ) |> svrep:::distribute_matrix_across_clusters(
          cluster_ids = c(3,3,1,1,2,2)
        ) |> as.matrix()
      )

  })

# Check Deville-Tille results ----

  test_that("Deville-Tille quadratic form gives expected estimate", {

    ## Create input data
    y <- rnorm(n = 5, mean = 5)
    aux_vars <- matrix(rnorm(n = 15), nrow = 5, ncol = 3)
    probs <- c(1/3, 1/2, 1, 1/2, 1/4)

    ## Calculate the estimator "by hand"
    c_values <- ((length(probs))/(length(probs) - ncol(aux_vars))) * (1-probs)

    denom_matrix <- lapply(seq_along(probs), function(k) {
      c_values[k] * (aux_vars[k,] %*% t(aux_vars[k,])) / (probs[k]^2)
    }) |> Reduce(f = `+`)

    inv_denom_matrix <- solve(denom_matrix)

    right_vector <- lapply(seq_along(probs), function(l) {
      c_l <- c_values[l]
      z_l <- aux_vars[l,]
      y_l <- y[l]
      pi_l <- probs[l]
      c_l * (z_l * y_l) / (pi_l^2)
    }) |> Reduce(f = `+`)

    hat_matrix <- inv_denom_matrix %*% right_vector

    y_star <- as.vector(aux_vars %*% hat_matrix)

    var_est <- sapply(seq_along(probs), function(k) {
      (c_values[k]/(probs[k]^2)) * (y[k] - y_star[k])^2
    }) |> sum()

    ## Compare the "by-hand" result to the quadratic form result
    Q <- make_deville_tille_matrix(
      probs = probs, aux_vars = aux_vars
    )

    Q_var_est <- as.vector(t(y/probs) %*% Q %*% (y/probs))

    expect_equal(
      object = Q_var_est,
      expected = var_est
    )

    ## If every unit selected with certainty, variance estimate should be 0
    expect_equal(
      object = make_deville_tille_matrix(probs = c(1,1),
                                         aux_vars = matrix(c(1,1), 2, 1)),
      expected = matrix(c(0,0,0,0), nrow = 2, ncol = 2)
    )

  })

# Ensure function checks inputs for issues ----

  test_that(
    "Informative errors for bad inputs", {

      expect_error(
        object = make_quad_form_matrix("made-up"),
        regexp = "`made-up` is not a supported variance estimator"
      )
      expect_error(
        object = make_quad_form_matrix(c("Yates-Grundy", "Horvitz-Thompson")),
        regexp = "Can only specify one"
      )

      # Horvitz-Thompson / Yates-Grundy
      expect_error(
        object = make_quad_form_matrix(variance_estimator = "Yates-Grundy"),
        regexp = "must supply a matrix.+joint_probs"
      )
      expect_error(
        object = make_quad_form_matrix(variance_estimator = "Yates-Grundy",
                                       joint_probs = matrix(NA, 2, 2)),
        regexp = "must be a matrix.+no missing values"
      )
      # SD1 and SD2
      expect_error(
        object = make_quad_form_matrix(variance_estimator = "SD1"),
        regexp = "must supply a matrix or data frame to `cluster_ids`"
      )
      expect_error(
        object = make_quad_form_matrix(variance_estimator = "SD1",
                                       cluster_ids = data.frame(ID = c(1,2))),
        regexp = "must supply a vector to `sort_order`"
      )
      # Stratified multistage SRS
      expect_error(
        object = make_quad_form_matrix(variance_estimator = "Stratified Multistage SRS",
                                       cluster_ids = data.frame(ID = c(1,2))),
        regexp = "matrix or data frame to both"
      )
      expect_error(
        object = make_quad_form_matrix(variance_estimator = "Stratified Multistage SRS",
                                       cluster_ids = data.frame(ID = c(1,2)),
                                       strata_ids = data.frame(STRATUM = c(1,1))),
        regexp = "matrix or data frame to `strata_pop_sizes`"
      )
      # Deville-1 and Deville-2
      expect_error(
        object = make_quad_form_matrix(variance_estimator = "Deville-1",
                                       cluster_ids = data.frame(ID = c(1,2)),
                                       strata_ids = data.frame(STRATUM = c(1,1))),
        regexp = "must supply a matrix or data frame to `probs`"
      )
  })

# Helper functions work ----

  test_that("`distribute_matrix_across_clusters() works", {
    expect_equal(
      object =   svrep:::distribute_matrix_across_clusters(
        cluster_level_matrix = matrix(c(1,2,3,4,
                                        5,6), nrow = 3, ncol = 2,
                                      byrow = TRUE),
        cluster_ids = c(1,2,2,3),
        cols = FALSE
      ),
      expected = structure(c(1, 3, 3, 5, 2, 4, 4, 6), dim = c(4L, 2L))
    )
    expect_equal(
      object =   svrep:::distribute_matrix_across_clusters(
        cluster_level_matrix = matrix(c(1,2,3,4,
                                        5,6), nrow = 2, ncol = 3,
                                      byrow = TRUE),
        cluster_ids = c(1,2,2,3),
        cols = TRUE, rows = FALSE
      ),
      expected = structure(c(1, 4, 2, 5, 2, 5, 3, 6), dim = c(2L, 4L))
    )
    expect_error(
      object = svrep:::distribute_matrix_across_clusters(
        cluster_level_matrix = matrix(c(1,2,3,4,
                                        5,6), nrow = 2, ncol = 3,
                                      byrow = TRUE),
        cluster_ids = c(1,2,2,3),
        cols = FALSE, rows = FALSE
      ),
      regexp = "Must set `rows=TRUE` and/or `cols=TRUE`"
    )
  })


