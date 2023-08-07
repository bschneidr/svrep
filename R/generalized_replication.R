as_gen_rep_design.survey.design <- function(design, variance_estimator = NULL,
                                            replicates = 50,
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

  n <- nrow(Sigma)
  decomposition <- svd(Sigma, nu = 0, nv = n)
  deltas <- decomposition$d
  lambdas <- decomposition$v

  Delta <- Matrix::Matrix(0, n, n) |> as("symmetricMatrix")
  for (k in seq_len(n)) {
    Delta <- Delta + (deltas[k] * outer(lambdas[,k], lambdas[,k]))
  }

  sqrt_D <- Matrix::Diagonal(x = sqrt(pmax(deltas, 0))) |> as("symmetricMatrix")

  H <- survey::hadamard(n = n)*(-2) + 1
  k_prime <- ncol(H)
  A <- sapply(seq_len(k_prime), function(r) {
    1 + colSums(Matrix::Diagonal(x = H[1:n,r]  * sqrt(pmax(deltas, 0))) %*% lambdas)
  })

  A <- Matrix::Matrix(1, n, k_prime)
  for (r in seq_len(k_prime)) {
    for (m in seq_len(n)) {
      A[,r] <- A[,r] + H[m,r] * sqrt(pmax(deltas[m], 0)) * lambdas[,m]
    }
  }

  (k_prime)^(-1) * (A %*% t(A))

  D <- Matrix::Diagonal(x = decomposition$d) |> as("symmetricMatrix")
  U <- decomposition$v
  V_prime <- Matrix::t(U)

  reconstructed <- U %*% D %*% V_prime

  # Generate adjustment factors for the compressed design object
  adjustment_factors <- make_gen_boot_factors(
    Sigma = Sigma,
    num_replicates = replicates,
    tau = tau,
    exact_vcov = exact_vcov
  )

  tau <- attr(adjustment_factors, 'tau')
  scale <- attr(adjustment_factors, 'scale')
  rscales <- attr(adjustment_factors, 'rscales')

  # Uncompress the adjustment factors
  adjustment_factors <- distribute_matrix_across_clusters(
    cluster_level_matrix = adjustment_factors,
    cluster_ids = compressed_design_structure$index,
    rows = TRUE, cols = FALSE
  )

  attr(adjustment_factors, 'tau') <- tau
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

  rep_design$tau <- tau

  if (inherits(design, 'tbl_svy') && ('package:srvyr' %in% search())) {
    rep_design <- srvyr::as_survey_rep(
      rep_design
    )
  }

  rep_design$call <- sys.call(which = -1)

  return(rep_design)
}
