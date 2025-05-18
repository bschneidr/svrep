#' @title Form replication factors using Fay's generalized replication method
#' @description Generate a matrix of replication factors
#' using Fay's generalized replication method.
#' This method yields a fully efficient variance estimator
#' if a sufficient number of replicates is used.
#' @param Sigma A quadratic form matrix corresponding to
#' a target variance estimator. Must be positive semidefinite.
#' @param max_replicates The maximum number of replicates to allow.
#' The function will attempt to create the minimum number of replicates
#' needed to produce a fully-efficient variance estimator.
#' If more replicates are needed than \code{max_replicates}, then the full number of replicates
#' needed will be created, but only a random subsample will be retained.
#' @param balanced If \code{balanced=TRUE}, the replicates
#' will all contribute equally to variance estimates, but
#' the number of replicates needed may slightly increase.
#' @return A matrix of replicate factors,
#' with the number of rows matching the number of rows of \code{Sigma}
#' and the number of columns less than or equal to \code{max_replicates}.
#' To calculate variance estimates using these factors,
#' use the overall scale factor given by calling
#' \code{attr(x, "scale")} on the result.
#' @section Statistical Details:
#' See Fay (1989) for a full explanation of Fay's generalized replication method.
#' This documentation provides a brief overview.
#'
#' Let \eqn{\boldsymbol{\Sigma}} be the quadratic form matrix for a target variance estimator,
#' which is assumed to be positive semidefinite.
#' Suppose the rank of \eqn{\boldsymbol{\Sigma}} is \eqn{k},
#' and so \eqn{\boldsymbol{\Sigma}} can be represented by the spectral decomposition
#' of \eqn{k} eigenvectors and eigenvalues, where the \eqn{r}-th eigenvector and eigenvalue
#' are denoted \eqn{\mathbf{v}_{(r)}} and \eqn{\lambda_r}, respectively.
#' \deqn{
#' \boldsymbol{\Sigma} = \sum_{r=1}^k \lambda_r \mathbf{v}_{(r)} \mathbf{v^{\prime}}_{(r)}
#' }
#' If \code{balanced = FALSE}, then we let \eqn{\mathbf{H}} denote an identity matrix
#' with \eqn{k' = k} rows/columns. If \code{balanced = TRUE}, then we let \eqn{\mathbf{H}} be a Hadamard matrix (with all entries equal to \eqn{1} or \eqn{-1}),
#' of order \eqn{k^{\prime} \geq k}. Let \eqn{\mathbf{H}_{mr}} denote the entry in row
#' \eqn{m} and column \eqn{r} of \eqn{\mathbf{H}}.
#'
#' Then \eqn{k^{\prime}} replicates are formed as follows.
#' Let \eqn{r} denote a given replicate, with \eqn{r = 1, ..., k^{\prime}},
#' and let \eqn{c} denote some positive constant (yet to be specified).
#'
#' The \eqn{r}-th replicate adjustment factor \eqn{\mathbf{f}_{r}} is formed as:
#' \deqn{
#'   \mathbf{f}_{r} = 1 + c \sum_{m=1}^k H_{m r} \lambda_{(m)}^{\frac{1}{2}} \mathbf{v}_{(m)}
#' }
#'
#' If \code{balanced = FALSE}, then \eqn{c = 1}. If \code{balanced = TRUE},
#' then \eqn{c = \frac{1}{\sqrt{k^{\prime}}}}.
#'
#' If any of the replicates
#' are negative, you can use \code{\link[svrep]{rescale_reps}},
#' which recalculates the replicate factors with a smaller value of \eqn{c}.
#'
#' If all \eqn{k^{\prime}} replicates are used, then variance estimates are calculated as:
#' \deqn{
#'   v_{rep}\left(\hat{T}_y\right) = \sum_{r=1}^{k^{\prime}}\left(\hat{T}_y^{*(r)}-\hat{T}_y\right)^2
#' }
#' For population totals, this replication variance estimator
#' will \emph{exactly} match the target variance estimator
#' if the number of replicates \eqn{k^{\prime}} matches the rank of \eqn{\Sigma}.
#'
#' @section The Number of Replicates:
#'
#' If \code{balanced=TRUE}, the number of replicates created
#' may need to increase slightly.
#' This is due to the fact that a Hadamard matrix
#' of order \eqn{k^{\prime} \geq k} is used to balance the replicates,
#' and it may be necessary to use order \eqn{k^{\prime} > k}.
#'
#' If the number of replicates \eqn{k^{\prime}} is too large for practical purposes,
#' then one can simply retain only a random subset of \eqn{R} of the \eqn{k^{\prime}} replicates.
#' In this case, variances are calculated as follows:
#' \deqn{
#'   v_{rep}\left(\hat{T}_y\right) = \frac{k^{\prime}}{R} \sum_{r=1}^{R}\left(\hat{T}_y^{*(r)}-\hat{T}_y\right)^2
#' }
#' This is what happens if \code{max_replicates} is less than the
#' matrix rank of \code{Sigma}: only a random subset
#' of the created replicates will be retained.
#'
#' Subsampling replicates is only recommended when
#' using \code{balanced=TRUE}, since in this case every replicate
#' contributes equally to variance estimates. If \code{balanced=FALSE},
#' then randomly subsampling replicates is valid but may
#' produce large variation in variance estimates since replicates
#' in that case may vary greatly in their contribution to variance
#' estimates.
#'
#' @section Reproducibility:
#'
#' If \code{balanced=TRUE}, a Hadamard matrix
#' is used as described above. The Hadamard matrix is
#' deterministically created using the function
#' \code{\link[survey]{hadamard}()} from the 'survey' package.
#' However, the order of rows/columns is randomly permuted
#' before forming replicates.
#'
#' In general, column-ordering of the replicate weights is random.
#' To ensure exact reproducibility, it is recommended to call
#' \code{\link[base]{set.seed}()} before using this function.
#'
#' @references
#'
#' Fay, Robert. 1989.
#' "Theory And Application Of Replicate Weighting For Variance Calculations."
#' In, 495-500. Alexandria, VA: American Statistical Association.
#' http://www.asasrms.org/Proceedings/papers/1989_033.pdf
#'
#' @seealso Use \code{\link[svrep]{rescale_reps}} to eliminate negative adjustment factors.
#' @export
#' @examples
#' \dontrun{
#'   library(survey)
#'
#' # Load an example dataset that uses unequal probability sampling ----
#'   data('election', package = 'survey')
#'
#' # Create matrix to represent the Horvitz-Thompson estimator as a quadratic form ----
#'   n <- nrow(election_pps)
#'   pi <- election_jointprob
#'   horvitz_thompson_matrix <- matrix(nrow = n, ncol = n)
#'   for (i in seq_len(n)) {
#'     for (j in seq_len(n)) {
#'       horvitz_thompson_matrix[i,j] <- 1 - (pi[i,i] * pi[j,j])/pi[i,j]
#'     }
#'   }
#'
#'   ## Equivalently:
#'
#'   horvitz_thompson_matrix <- make_quad_form_matrix(
#'     variance_estimator = "Horvitz-Thompson",
#'     joint_probs = election_jointprob
#'   )
#'
#' # Make generalized replication adjustment factors ----
#'
#'   adjustment_factors <- make_fays_gen_rep_factors(
#'     Sigma = horvitz_thompson_matrix,
#'     max_replicates = 50
#'   )
#'   attr(adjustment_factors, 'scale')
#'
#' # Compute the Horvitz-Thompson estimate and the replication estimate
#'
#' ht_estimate <- svydesign(data = election_pps, ids = ~ 1,
#'                          prob = diag(election_jointprob),
#'                          pps = ppsmat(election_jointprob)) |>
#'   svytotal(x = ~ Kerry)
#'
#' rep_estimate <- svrepdesign(
#'   data = election_pps,
#'   weights = ~ wt,
#'   repweights = adjustment_factors,
#'   combined.weights = FALSE,
#'   scale = attr(adjustment_factors, 'scale'),
#'   rscales = rep(1, times = ncol(adjustment_factors)),
#'   type = "other",
#'   mse = TRUE
#' ) |>
#'   svytotal(x = ~ Kerry)
#'
#' SE(rep_estimate)
#' SE(ht_estimate)
#' SE(rep_estimate) / SE(ht_estimate)
#' }

make_fays_gen_rep_factors <- function(
    Sigma,
    max_replicates = Inf,
    balanced = TRUE
) {

  n <- nrow(Sigma)

  # Calculate spectral decomposition ----
  eigen_decomposition <- compute_eigen_decomposition(Sigma)
  
  Sigma_rank <- Matrix::rankMatrix(Sigma, method = "qr")

  # Obtain eigenvectors scaled by square roots of eigenvalues ----
  v <- sapply(X = seq_along(eigen_decomposition$values),
              FUN = function(k) {
                if (eigen_decomposition$values[k] < 0) {
                  rep(0, times = nrow(eigen_decomposition$vectors))
                } else {
                  sqrt(eigen_decomposition$values[k]) * eigen_decomposition$vectors[,k]
                }
              })
  v <- v[, seq_len(Sigma_rank), drop=FALSE]

  # Create replication factors
  if (!balanced) {
    k_prime <- Sigma_rank
    # Shuffle the column order
    shuffle_order <- sample(x = k_prime, size = k_prime, replace = FALSE)
    replicate_factors <- v[,shuffle_order,drop=FALSE]
  }

  if (balanced) {
    # Generate Hadamard matrix
    H <- (2*survey::hadamard(Sigma_rank) - 1)
    k_prime <- ncol(H)
    shuffle_order <- sample(x = k_prime, size = k_prime, replace = FALSE)
    H <- H[shuffle_order, shuffle_order]

    # Construct replicate factors
    replicate_factors <- (
      v[,seq_len(Sigma_rank),drop=FALSE] %*% H[seq_len(Sigma_rank),,drop=FALSE]
    )
    replicate_factors <- replicate_factors / sqrt(k_prime)
  }

  replicate_factors <- 1 + replicate_factors

  # Set overall scale factor used in estimation
  scale <- 1

  # Potentially subsample replicate factors
  # (and adjust the scale factor if doing any subsampling)
  if (max_replicates >= k_prime) {
    num_replicates <- k_prime
  }
  if (max_replicates < k_prime) {
    num_replicates <- max_replicates
    scale <- scale * (k_prime/num_replicates)
    replicate_factors <- replicate_factors[,sample(num_replicates),drop=FALSE]

    if (max_replicates < Sigma_rank) {
      msg <- sprintf(
        "The number of replicates needed for fully efficient replication is %s, but `max_replicates` is set to %s.",
        k_prime, max_replicates
      )
      msg <- paste(msg, "Only a random sample of replicates will be retained.")
      message(msg)
      if (!balanced) {
        warning("Random subsampling of replicates is not recommended when `balanced=FALSE`. See the help page `?make_fays_gen_rep_factors` for details.")
      }
    }
  }

  attr(replicate_factors, 'scale') <- scale
  attr(replicate_factors, 'rscales') <- rep(1, times = num_replicates)

  # Set column names
  colnames(replicate_factors) <- sprintf("REP_%s", seq_len(num_replicates))

  # Return result
  return(replicate_factors)
}

#' @title Convert a survey design object to a replication design
#' using Fay's generalized replication method
#' @description Converts a survey design object to a replicate design object
#' with replicate weights formed using the generalized replication method of Fay (1989).
#' The generalized replication method forms replicate weights
#' from a textbook variance estimator, provided that the variance estimator
#' can be represented as a quadratic form whose matrix is positive semidefinite
#' (this covers a large class of variance estimators).
#' @param design A survey design object created using the 'survey' (or 'srvyr') package,
#' with class \code{'survey.design'} or \code{'svyimputationList'}.
#' @param variance_estimator The name of the variance estimator
#' whose quadratic form matrix should be created.
#' See \link[svrep]{variance-estimators} for a
#' detailed description of each variance estimator.
#' Options include:
#' \itemize{
#'   \item \strong{"Yates-Grundy"}: \cr The Yates-Grundy variance estimator based on
#'     first-order and second-order inclusion probabilities.
#'   \item \strong{"Horvitz-Thompson"}: \cr The Horvitz-Thompson variance estimator based on
#'     first-order and second-order inclusion probabilities.
#'   \item \strong{"Poisson Horvitz-Thompson"}: \cr The Horvitz-Thompson variance estimator
#'     based on assuming Poisson sampling, with first-order inclusion probabilities
#'     inferred from the sampling probabilities of the survey design object.
#'   \item \strong{"Stratified Multistage SRS"}: \cr The usual stratified multistage variance estimator
#'     based on estimating the variance of cluster totals within strata at each stage.
#'   \item \strong{"Ultimate Cluster"}: \cr The usual variance estimator based on estimating
#'     the variance of first-stage cluster totals within first-stage strata.
#'   \item \strong{"Deville-1"}: \cr A variance estimator for unequal-probability
#'     sampling without replacement, described in Matei and Tillé (2005)
#'     as "Deville 1".
#'   \item \strong{"Deville-2"}: \cr A variance estimator for unequal-probability
#'     sampling without replacement, described in Matei and Tillé (2005) as "Deville 2".
#'   \item \strong{"Deville-Tille": } \cr A variance estimator useful
#'     for balanced sampling designs, proposed by Deville and Tillé (2005).
#'   \item \strong{"SD1"}: \cr The non-circular successive-differences variance estimator described by Ash (2014),
#'     sometimes used for variance estimation for systematic sampling.
#'   \item \strong{"SD2"}:  \cr The circular successive-differences variance estimator described by Ash (2014).
#'     This estimator is the basis of the "successive-differences replication" estimator commonly used
#'     for variance estimation for systematic sampling.
#'   \item \strong{"BOSB"}: \cr The kernel-based variance estimator proposed by
#'     Breidt, Opsomer, and Sanchez-Borrego (2016) for use with systematic samples
#'     or other finely stratified designs. Uses the Epanechnikov kernel
#'     with the bandwidth automatically chosen to result in the smallest possible
#'     nonempty kernel window.
#'   \item\strong{"Beaumont-Emond"}: \cr The variance estimator of Beaumont and Emond (2022)
#'     for multistage unequal-probability sampling without replacement.
#' }
#' @param aux_var_names (Only used if \code{variance_estimator = "Deville-Tille")}.
#' A vector of the names of auxiliary variables used in sampling.
#' @param max_replicates The maximum number of replicates to allow (should be as large as possible, given computer memory/storage limitations).
#' A commonly-recommended default is 500. If the number of replicates needed
#' for a balanced, fully-efficient estimator is less than \code{max_replicates},
#' then only the number of replicates needed will be created.
#' If more replicates are needed than \code{max_replicates}, then the full number of replicates
#' needed will be created, but only a random subsample will be retained.
#' @param balanced If \code{balanced=TRUE}, the replicates
#' will all contribute equally to variance estimates, but
#' the number of replicates needed may slightly increase.
#' @param psd_option Either \code{"warn"} (the default) or \code{"error"}.
#' This option specifies what will happen if the target variance estimator
#' has a quadratic form matrix which is not positive semidefinite. This
#' can occasionally happen, particularly for two-phase designs. \cr
#' If \code{psd_option="error"}, then an error message will be displayed. \cr
#' If \code{psd_option="warn"}, then a warning message will be displayed,
#' and the quadratic form matrix will be approximated by the most similar
#' positive semidefinite matrix.
#' This approximation was suggested by Beaumont and Patak (2012),
#' who note that this is conservative in the sense of producing
#' overestimates of variance.
#' Beaumont and Patak (2012) argue that this overestimation is expected to be
#' small in magnitude. See \code{\link[svrep]{get_nearest_psd_matrix}}
#' for details of the approximation.
#' @param compress This reduces the computer memory required to represent the replicate weights and has no
#' impact on estimates.
#' @param mse If \code{TRUE} (the default), compute variances from sums of squares around the point estimate from the full-sample weights.
#' If \code{FALSE}, compute variances from sums of squares around the mean estimate from the replicate weights.
#' For Fay's generalized replication method, setting \code{mse = FALSE} can potentially
#' lead to large underestimates of variance.
#' @return
#' A replicate design object, with class \code{svyrep.design}, which can be used with the usual functions,
#' such as \code{svymean()} or \code{svyglm()}.
#'
#' Use \code{weights(..., type = 'analysis')} to extract the matrix of replicate weights.
#'
#' Use \code{as_data_frame_with_weights()} to convert the design object to a data frame with columns
#' for the full-sample and replicate weights.
#' @export

#' @section Statistical Details:
#' See Fay (1989) for a full description of this replication method,
#' or see the documentation in \link[svrep]{make_fays_gen_rep_factors} for implementation details.
#'
#' See \link[svrep]{variance-estimators} for a
#' description of each variance estimator available for use with
#' this function.
#'
#' Use \code{\link[svrep]{rescale_reps}} to eliminate negative adjustment factors.
#'
#' @seealso
#' For greater customization of the method, \code{\link[svrep]{make_quad_form_matrix}} can be used to
#' represent several common variance estimators as a quadratic form's matrix,
#' which can then be used as an input to \code{\link[svrep]{make_fays_gen_rep_factors}}.
#'
#' @section Two-Phase Designs:
#' For a two-phase design, \code{variance_estimator} should be a list of variance estimators' names,
#' with two elements, such as \code{list('Ultimate Cluster', 'Poisson Horvitz-Thompson')}.
#' In two-phase designs, only the following estimators may be used for the second phase:
#' \itemize{
#'   \item "Ultimate Cluster"
#'   \item "Stratified Multistage SRS"
#'   \item "Poisson Horvitz-Thompson"
#' }
#' For statistical details on the handling of two-phase designs,
#' see the documentation for \link[svrep]{make_twophase_quad_form}.
#' @references
#' The generalized replication method was first proposed in
#' Fay (1984). Fay (1989) refined the generalized replication method
#' to produce "balanced" replicates, in the sense that
#' each replicate contributes equally to variance estimates.
#' The advantage of balanced replicates is that one can
#' still obtain a reasonable variance estimate
#' by using only a random subset of the replicates.
#'
#' - Ash, S. (2014). "\emph{Using successive difference replication for estimating variances}."
#' \strong{Survey Methodology}, Statistics Canada, 40(1), 47-59.
#' \cr \cr
#' - Beaumont, J.-F.; Émond, N. (2022).
#' "\emph{A Bootstrap Variance Estimation Method for Multistage Sampling and Two-Phase Sampling When Poisson Sampling Is Used at the Second Phase}."
#' \strong{Stats}, \emph{5}: 339-357.
#' https://doi.org/10.3390/stats5020019
#' \cr \cr
#' - Breidt, F. J., Opsomer, J. D., & Sanchez-Borrego, I. (2016). 
#' "\emph{Nonparametric Variance Estimation Under Fine Stratification: An Alternative to Collapsed Strata}." 
#' \strong{Journal of the American Statistical Association}, 111(514), 822-833. https://doi.org/10.1080/01621459.2015.1058264
#' \cr \cr
#' - Deville, J. C., and Tillé, Y. (2005). "\emph{Variance approximation under balanced sampling.}"
#' \strong{Journal of Statistical Planning and Inference}, 128, 569-591.
#' \cr \cr
#' - Dippo, Cathryn, Robert Fay, and David Morganstein. 1984. "Computing Variances from Complex Samples with Replicate Weights." In, 489-94. Alexandria, VA: American Statistical Association. http://www.asasrms.org/Proceedings/papers/1984_094.pdf.
#' \cr \cr
#' - Fay, Robert. 1984. "Some Properties of Estimates of Variance Based on Replication Methods." In, 495-500. Alexandria, VA: American Statistical Association. http://www.asasrms.org/Proceedings/papers/1984_095.pdf.
#' \cr \cr
#' - Fay, Robert. 1989. "Theory And Application Of Replicate Weighting For Variance Calculations." In, 495-500. Alexandria, VA: American Statistical Association. http://www.asasrms.org/Proceedings/papers/1989_033.pdf
#' \cr \cr
#' - Matei, Alina, and Yves Tillé. (2005).
#' "\emph{Evaluation of Variance Approximations and Estimators
#' in Maximum Entropy Sampling with Unequal Probability and Fixed Sample Size.}"
#' \strong{Journal of Official Statistics}, 21(4):543-70.
#' @examples
#' if (FALSE) {
#'
#'   library(survey)
#'
#'   ## Load an example systematic sample ----
#'   data('library_stsys_sample', package = 'svrep')
#'
#'   ## First, ensure data are sorted in same order as was used in sampling
#'   library_stsys_sample <- library_stsys_sample |>
#'     sort_by(~ SAMPLING_SORT_ORDER)
#'
#'   ## Create a survey design object
#'   design_obj <- svydesign(
#'     data = library_stsys_sample,
#'     strata = ~ SAMPLING_STRATUM,
#'     ids = ~ 1,
#'     fpc = ~ STRATUM_POP_SIZE
#'   )
#'
#'   ## Convert to generalized replicate design
#'
#'   gen_rep_design_sd2 <- as_fays_gen_rep_design(
#'     design = design_obj,
#'     variance_estimator = "SD2",
#'     max_replicates = 250,
#'     mse = TRUE
#'   )
#'
#'   svytotal(x = ~ TOTSTAFF, na.rm = TRUE, design = gen_rep_design_sd2)
#' }
#' @export
as_fays_gen_rep_design <- function(design, variance_estimator = NULL,
                                   aux_var_names = NULL,
                                   max_replicates = Inf,
                                   balanced = TRUE,
                                   psd_option = "warn",
                                   mse = TRUE,
                                   compress = TRUE) {
  
  if (!mse) {
    warning("When `balanced = FALSE`, setting `mse = FALSE` may produce large underestimates of variance.")
  }
  
  UseMethod("as_fays_gen_rep_design", design)
}

#' @export
as_fays_gen_rep_design.twophase2 <- function(design, variance_estimator = NULL,
                                             aux_var_names = NULL,
                                             max_replicates = Inf,
                                             balanced = TRUE,
                                             psd_option = "warn",
                                             mse = TRUE,
                                             compress = TRUE) {

  Sigma <- get_design_quad_form(design, variance_estimator, aux_var_names = aux_var_names)

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

  adjustment_factors <- make_fays_gen_rep_factors(
    Sigma = Sigma,
    max_replicates = max_replicates,
    balanced = balanced
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
                                                 aux_var_names = NULL,
                                                 max_replicates = Inf,
                                                 balanced = TRUE,
                                                 psd_option = 'warn',
                                                 mse = TRUE,
                                                 compress = TRUE) {
  

  # Produce a (potentially) compressed survey design object
  compressed_design_structure <- compress_design(design, vars_to_keep = aux_var_names)

  # Get the quadratic form of the variance estimator,
  # for the compressed design object
  Sigma <- get_design_quad_form(
    compressed_design_structure$design_subset,
    variance_estimator,
    aux_var_names = aux_var_names
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
  adjustment_factors <- make_fays_gen_rep_factors(
    Sigma = Sigma,
    max_replicates = max_replicates,
    balanced = balanced
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

#' @export
as_fays_gen_rep_design.DBIsvydesign <- function(design, variance_estimator = NULL,
                                                aux_var_names = NULL,
                                                max_replicates = Inf,
                                                balanced = TRUE,
                                                psd_option = 'warn',
                                                mse = TRUE,
                                                compress = TRUE) {

  # Produce a (potentially) compressed survey design object
  compressed_design_structure <- compress_design(design, vars_to_keep = aux_var_names)

  # Get the quadratic form of the variance estimator,
  # for the compressed design object
  Sigma <- get_design_quad_form(
    compressed_design_structure$design_subset,
    variance_estimator,
    aux_var_names = aux_var_names
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
  adjustment_factors <- make_fays_gen_rep_factors(
    Sigma = Sigma,
    max_replicates = max_replicates,
    balanced = balanced
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
    data = data.frame(DUMMY = seq_len(nrow(adjustment_factors))),
    weights = stats::weights(design, type = "sampling"),
    repweights = as.matrix(adjustment_factors),
    combined.weights = FALSE,
    compress = compress, mse = mse,
    scale = scale,
    rscales = rscales,
    type = "other"
  )

  # Replace 'variables' with a database connection
  # and make the object have the appropriate class
  rep_design$variables <- NULL
  if (design$db$dbtype == "ODBC") {
    stop("'RODBC' no longer supported. Use the odbc package")
  } else {
    db <- DBI::dbDriver(design$db$dbtype)
    dbconn <- DBI::dbConnect(db, design$db$dbname)
  }
  rep_design$db <- list(
    dbname = design$db$dbname, tablename = design$db$tablename,
    connection = dbconn,
    dbtype = design$db$dbtype
  )
  class(rep_design) <- c("DBIrepdesign", "DBIsvydesign", class(rep_design))

  if (inherits(design, 'tbl_svy') && ('package:srvyr' %in% search())) {
    rep_design <- srvyr::as_survey_rep(
      rep_design
    )
  }

  rep_design$call <- sys.call(which = -1)

  return(rep_design)

}
