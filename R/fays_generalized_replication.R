make_fay_gen_rep_factors <- function(Sigma, max_replicates) {

  n <- nrow(Sigma)

  # Calculate spectral decomposition ----
  eigen_decomposition <- eigen(x = Sigma, symmetric = TRUE)
  matrix_rank <- Matrix::rankMatrix(Sigma, method = "qr")

  # Obtain eigenvectors scaled by square roots of eigenvalues ----
  v <- sapply(X = seq_along(eigen_decomposition$values),
              FUN = function(k) {
                truncated_eigenvalue <- ifelse(eigen_decomposition$values[k] < 0,
                                               0, eigen_decomposition$values[k])
                sqrt(truncated_eigenvalue) * eigen_decomposition$vectors[,k]
              })
  v <- v[, seq_len(matrix_rank), drop=FALSE]

  # Generate Hadamard matrix
  k <- matrix_rank
  H <- (2*survey::hadamard(k) - 1)
  k_prime <- ncol(H)
  shuffle_order <- sample(x = k_prime, size = k_prime)
  H <- H[shuffle_order, shuffle_order]

  # Construct replicate factors
  replicate_factors <- (v[,seq_len(k),drop=FALSE] %*% H[seq_len(k),,drop=FALSE])

  max_flipped_value <- max(-replicate_factors)
  if (max_flipped_value > 0) {
    scale_factor <- 1/max_flipped_value
  } else {
    scale_factor <- 1
  }

  replicate_factors <- 1 + (scale_factor * replicate_factors)

  scale <- (k_prime * (scale_factor^2))^(-1)

  # Potentially subsample replicate factors
  if (max_replicates >= k_prime) {
    num_replicates <- k_prime
  }
  if (max_replicates < k_prime) {
    num_replicates <- max_replicates
    scale <- scale * (k_prime/num_replicates)
    replicate_factors <- replicate_factors[,sample(num_replicates),drop=FALSE]
  }

  replicate_factors |> rowMeans()

  attr(replicate_factors, 'scale') <- scale
  attr(replicate_factors, 'rscales') <- rep(1, times = num_replicates)

  # Set column names
  colnames(replicate_factors) <- sprintf("REP_%s", seq_len(num_replicates))

  # Return result
  return(replicate_factors)
}

#' @export
as_fays_gen_rep_design <- function(design, variance_estimator = NULL,
                               max_replicates = 500,
                               psd_option = "warn",
                               mse = getOption("survey.replicates.mse"),
                               compress = TRUE) {
  UseMethod("as_fays_gen_rep_design", design)
}

#' @export
as_fays_gen_rep_design.twophase2 <- function(design, variance_estimator = NULL,
                                         max_replicates = 500,
                                         psd_option = "warn",
                                         mse = getOption("survey.replicates.mse"),
                                         compress = TRUE) {

  Sigma <- get_design_quad_form(design, variance_estimator)

  if (!is_psd_matrix(Sigma)) {
    problem_msg <- paste0(
      "The sample quadratic form matrix",
      " for this design and variance estimator",
      " is not positive semidefinite."
    )
    if (psd_option == "warn") {

      warning_msg <- paste0(
        problem_msg,
        " It will be approximated by the nearest",
        " positive semidefinite matrix."
      )
      warning(warning_msg)
      Sigma <- get_nearest_psd_matrix(Sigma)

    } else {
      error_msg <- paste0(
        problem_msg,
        " This can be handled using the `psd_option` argument."
      )
      stop(error_msg)
    }
  }

  adjustment_factors <- make_fay_gen_rep_factors(
    Sigma = Sigma,
    max_replicates = max_replicates
  )

  rep_design <- survey::svrepdesign(
    variables = design$phase1$full$variables[design$subset,,drop=FALSE],
    weights = stats::weights(design, type = "sampling"),
    repweights = as.matrix(adjustment_factors),
    combined.weights = FALSE,
    compress = compress, mse = mse,
    scale = attr(adjustment_factors, 'scale'),
    rscales = attr(adjustment_factors, 'rscales'),
    type = "other"
  )

  if (inherits(design, 'tbl_svy') && ('package:srvyr' %in% search())) {
    rep_design <- srvyr::as_survey_rep(
      rep_design
    )
  }

  rep_design$call <- sys.call(which = -1)

  return(rep_design)
}

#' @export
as_fays_gen_rep_design.survey.design <- function(design, variance_estimator = NULL,
                                             max_replicates = 500,
                                             psd_option = 'warn',
                                             mse = getOption("survey.replicates.mse"),
                                             compress = TRUE) {

  # Produce a (potentially) compressed survey design object
  if ((!is.null(design$pps)) && (design$pps != FALSE)) {
    compressed_design_structure <- list(
      design_subset = design,
      index = seq_len(nrow(design))
    )
  } else {
    design_structure <- cbind(design$strata, design$cluster)
    tmp <- apply(design_structure, 1, function(x) paste(x, collapse = "\r"))
    unique_elements <- !duplicated(design_structure)
    compressed_design_structure <- list(
      design_subset = design |> (\(design_obj) {
        # Reduce memory usage by dropping variables
        design_obj$variables <- design_obj$variables[,0,drop=FALSE]
        # Subset to only unique strata/cluster combos
        design_obj[unique_elements,]
      })(),
      index = match(tmp, tmp[unique_elements])
    )
  }

  # Get the quadratic form of the variance estimator,
  # for the compressed design object
  Sigma <- get_design_quad_form(
    compressed_design_structure$design_subset,
    variance_estimator
  )

  # Check that the matrix is positive semidefinite
  if (!is_psd_matrix(Sigma)) {
    problem_msg <- "The sample quadratic form matrix for this design and variance estimator is not positive semidefinite."
    if (psd_option == "warn") {

      warning_msg <- paste0(
        problem_msg,
        " It will be approximated by the nearest positive semidefinite matrix."
      )
      warning(warning_msg)
      Sigma <- get_nearest_psd_matrix(Sigma)

    } else {
      error_msg <- paste0(
        problem_msg, " This can be handled using the `psd_option` argument."
      )
      stop(error_msg)
    }
  }

  # Generate adjustment factors for the compressed design object
  adjustment_factors <- make_fay_gen_rep_factors(
    Sigma = Sigma,
    max_replicates = max_replicates
  )

  scale <- attr(adjustment_factors, 'scale')
  rscales <- attr(adjustment_factors, 'rscales')

  # Uncompress the adjustment factors
  adjustment_factors <- distribute_matrix_across_clusters(
    cluster_level_matrix = adjustment_factors,
    cluster_ids = compressed_design_structure$index,
    rows = TRUE, cols = FALSE
  )

  attr(adjustment_factors, 'scale') <- scale
  attr(adjustment_factors, 'rscales') <- rscales

  # Create the survey design object
  rep_design <- survey::svrepdesign(
    variables = design$variables,
    weights = stats::weights(design, type = "sampling"),
    repweights = as.matrix(adjustment_factors),
    combined.weights = FALSE,
    compress = compress, mse = mse,
    scale = scale,
    rscales = rscales,
    type = "other"
  )

  if (inherits(design, 'tbl_svy') && ('package:srvyr' %in% search())) {
    rep_design <- srvyr::as_survey_rep(
      rep_design
    )
  }

  rep_design$call <- sys.call(which = -1)

  return(rep_design)
}


# Examples ----

##_ Systematic sample ----
data('library_stsys_sample', package = 'svrep')

## First, ensure data are sorted in same order as was used in sampling
library_stsys_sample <- library_stsys_sample[
  order(library_stsys_sample$SAMPLING_SORT_ORDER),
]

## Create a survey design object
design_obj <- svydesign(
  data = library_stsys_sample,
  strata = ~ SAMPLING_STRATUM,
  ids = ~ 1,
  fpc = ~ STRATUM_POP_SIZE
)

## Convert to generalized replicate design
gen_boot_design_sd2 <- as_gen_boot_design(
  design = design_obj,
  variance_estimator = "SD2",
  replicates = 250, exact_vcov = TRUE
)
gen_rep_design_sd2 <- as_fays_gen_rep_design(
  design = design_obj,
  variance_estimator = "SD2",
  max_replicates = 250,
  mse = TRUE
)

svytotal(x = ~ TOTSTAFF, na.rm = TRUE, design = gen_boot_design_sd2)
svytotal(x = ~ TOTSTAFF, na.rm = TRUE, design = gen_rep_design_sd2)

svyquantile(x = ~ LIBRARIA, quantiles = 0.5, na.rm = TRUE,
            design = gen_boot_design_sd2, interval.type = "quantile")
svyquantile(x = ~ LIBRARIA, quantiles = 0.5, na.rm = TRUE,
            design = gen_rep_design_sd2, interval.type = "quantile")

##_ "Large" sample ----

design_obj <- svydesign(data = svrep::lou_vax_survey,
                        ids = ~ 1)

### Random groups jackknife
  group_randomly <- function(n, n_groups) {
    sample(as.numeric(cut(seq_len(n), breaks = n_groups)),
           size = n, replace = FALSE)
  }
  ran_grp_jkn_design <- svrep::lou_vax_survey |>
    transform(random_group = group_randomly(n = 1000,
                                            n_groups = 50)) |>
    svydesign(data = _, ids = ~ random_group) |>
    as.svrepdesign(type = "JK1", mse = TRUE)

### Fay's generalized replication method
  gen_rep_design <- as_fays_gen_rep_design(
    design = design_obj,
    variance_estimator = "Ultimate Cluster",
    mse = TRUE,
    max_replicates = 50
  )

  svytotal(x = ~ VAX_STATUS, na.rm = TRUE,
           design = design_obj) |> SE()

  gen_boot_design <- as_gen_boot_design(
    design = design_obj,
    variance_estimator = "Ultimate Cluster",
    mse = TRUE, exact_vcov = FALSE,
    replicates = 50
  )

  svytotal(x = ~ VAX_STATUS, na.rm = TRUE,
           design = gen_boot_design) |> SE()


  rep_results <- replicate(n = 500, expr = {
    ran_grp_jkn_design <- svrep::lou_vax_survey |>
      transform(random_group = group_randomly(n = 1000,
                                              n_groups = 20)) |>
      svydesign(data = _, ids = ~ random_group) |>
      as.svrepdesign(type = "JK1", mse = TRUE)
    gen_rep_design <- as_fay_gen_rep_design(
      design = design_obj,
      variance_estimator = "Ultimate Cluster",
      mse = TRUE,
      max_replicates = 20
    )
    list(
      'jk' = svymean(x = ~ I(as.numeric(VAX_STATUS == "Vaccinated")), na.rm = TRUE,
                      design = ran_grp_jkn_design),
      'genrep' = svymean(x = ~ I(as.numeric(VAX_STATUS == "Vaccinated")), na.rm = TRUE,
                          design = gen_rep_design)
    ) |> sapply(SE)
  })

exp_value <- (svymean(x = ~ I(as.numeric(VAX_STATUS == "Vaccinated")), na.rm = TRUE,
                      design = design_obj) |> SE())

((rep_results - (rep(exp_value, 2))^2)) |>
  rowMeans() |> sqrt()


## Simulation ----

rep_results <- replicate(
  n = 2000, expr = {
    suppressMessages({
      gen_rep_design_sd2 <- as_fay_gen_rep_design(
        design = design_obj,
        variance_estimator = "SD2",
        max_replicates = 50,
        mse = TRUE
      )
    })
    svytotal(x = ~ TOTSTAFF,
             na.rm = TRUE,
             design = gen_rep_design_sd2) |> SE()
  }
)

mean_result <- rep_results |> mean()
sd_results <- sd(rep_results)
sd_results/mean_result
