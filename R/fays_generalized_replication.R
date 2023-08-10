#' @title Form replication factors using Fay's generalized replication method
#' @description Generate a matrix of replication factors
#' using Fay's generalized replication method. The replication
#' factors are guaranteed to be nonnegative. This method yields
#' a fully efficient variance estimator if a sufficient number of
#' replicates is used.
#' @param Sigma A quadratic form matrix corresponding to
#' a target variance estimator. Must be positive semidefinite.
#' @param max_replicates The maximum number of replicates to allow.
#' The function will attempt to create the minimum number of replicates
#' needed to produce a fully-efficient variance estimator.
#' If more replicates are needed than \code{max_replicates}, then the full number of replicates
#' needed will be created, but only a random subsample will be retained.
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
#' Let \eqn{\mathbf{H}} be a Hadamard matrix (with all entries equal to \eqn{1} or \eqn{-1}),
#' of order \eqn{k^{\prime} \geq k}. Let \eqn{\mathbf{H}_{mr}} denote the entry in row
#' \eqn{m} and column \eqn{r} of \eqn{\mathbf{H}}.
#'
#' Then \eqn{k^{\prime}} replicates are formed as follows.
#' Let \eqn{r} denote a given replicate, with \eqn{r = 1, ..., k^{\prime}},
#' and let \eqn{c} denote some positive constant (yet to be specified).
#'
#' For now, let \eqn{c=1} and create the \eqn{r}-th replicate's adjustment factor \eqn{\mathbf{f}_{r}} as:
#' \deqn{
#'   \mathbf{f}_{r} = 1 + c \sum_{m=1}^k H_{m r} \lambda^{\frac{1}{2}} v_{(m)}
#' }
#' It's possible that one or more of the replicates will have a negative adjustment factor.
#' This can easily be eliminated by an appropriate choice of the constant \eqn{c}.
#' If there are any negative adjustment factors,
#' then set \eqn{c^{-1}} equal to some value just larger than the
#' absolute value of the minimum replicate adjustment factor (from across all the replicates).
#' Then recalculate \eqn{\mathbf{f}_{r}} using the updated value of \eqn{c}.
#'
#' If all \eqn{k^{\prime}} replicates are used, then variance estimates are calculated as:
#' \deqn{
#'   v_{rep}\left(\hat{T}_y\right) = \frac{c^2}{k^{\prime}} \sum_{r=1}^{k^{\prime}}\left(\hat{T}_y^{*(r)}-\hat{T}_y\right)^2
#' }
#' For population totals, this replication variance estimator
#' will \emph{exactly} match the target variance estimator.
#'
#' If the number of replicates \eqn{k^{\prime}} is too large for practical purposes,
#' then one can simply retain only a random subset of \eqn{R} of the \eqn{k^{\prime}} replicates.
#' In this case, variances are calculated as follows:
#' \deqn{
#'   v_{rep}\left(\hat{T}_y\right) = \frac{c^2}{R} \sum_{r=1}^{R}\left(\hat{T}_y^{*(r)}-\hat{T}_y\right)^2
#' }
#' This is what happens if \code{max_replicates} is less than the
#' matrix rank of \code{Sigma}: only a random subset
#' of the created replicates will be retained.
#'
#' @export
#' @examples
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
#'   adjustment_factors <- make_fay_gen_rep_factors(
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
#'   type = "other"
#' ) |>
#'   svytotal(x = ~ Kerry)
#'
#' SE(rep_estimate)
#' SE(ht_estimate)
#' SE(rep_estimate) / SE(ht_estimate)

make_fay_gen_rep_factors <- function(Sigma, max_replicates) {

  n <- nrow(Sigma)

  # Calculate spectral decomposition ----
  eigen_decomposition <- eigen(x = Sigma, symmetric = TRUE)
  Sigma_rank <- Matrix::rankMatrix(Sigma, method = "qr")

  # Obtain eigenvectors scaled by square roots of eigenvalues ----
  v <- sapply(X = seq_along(eigen_decomposition$values),
              FUN = function(k) {
                truncated_eigenvalue <- ifelse(eigen_decomposition$values[k] < 0,
                                               0, eigen_decomposition$values[k])
                sqrt(truncated_eigenvalue) * eigen_decomposition$vectors[,k]
              })
  v <- v[, seq_len(Sigma_rank), drop=FALSE]

  # Generate Hadamard matrix
  H <- (2*survey::hadamard(Sigma_rank) - 1)
  k_prime <- ncol(H)
  shuffle_order <- sample(x = k_prime, size = k_prime, replace = FALSE)
  H <- H[shuffle_order, shuffle_order]

  # Construct replicate factors
  replicate_factors <- (
    v[,seq_len(Sigma_rank),drop=FALSE] %*% H[seq_len(Sigma_rank),,drop=FALSE]
  )

  max_flipped_value <- max(-replicate_factors)
  if (max_flipped_value > 0) {
    scale_factor <- 1/(max_flipped_value*1.001)
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

    if (max_replicates < Sigma_rank) {
      sprintf(
        "The number of replicates needed for balanced, fully-efficient replication is %s, but `max_replicates` is set to %s",
        k_prime, max_replicates
      ) |> message()
    }
  }

  replicate_factors |> rowMeans()

  attr(replicate_factors, 'scale') <- scale
  attr(replicate_factors, 'rscales') <- rep(1, times = num_replicates)

  # Set column names
  colnames(replicate_factors) <- sprintf("REP_%s", seq_len(num_replicates))

  # Return result
  return(replicate_factors)
}

#' @title Convert a survey design object to a generalized replication replicate design
#' @description Converts a survey design object to a replicate design object
#' with replicate weights formed using the generalized replication method of Fay (1984).
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
#'   \item{\strong{"Yates-Grundy"}: }{The Yates-Grundy variance estimator based on
#'   first-order and second-order inclusion probabilities.}
#'   \item{\strong{"Horvitz-Thompson"}: }{The Horvitz-Thompson variance estimator based on
#'   first-order and second-order inclusion probabilities.}
#'   \item{\strong{"Poisson Horvitz-Thompson"}: }{The Horvitz-Thompson variance estimator
#'   based on assuming Poisson sampling, with first-order inclusion probabilities
#'   inferred from the sampling probabilities of the survey design object.}
#'   \item{\strong{"Stratified Multistage SRS"}: }{The usual stratified multistage variance estimator
#'   based on estimating the variance of cluster totals within strata at each stage.}
#'   \item{\strong{"Ultimate Cluster"}: }{The usual variance estimator based on estimating
#'   the variance of first-stage cluster totals within first-stage strata.}
#'   \item{\strong{"Deville-1"}: }{A variance estimator for unequal-probability
#'   sampling without replacement, described in Matei and Tillé (2005)
#'   as "Deville 1".}
#'   \item{\strong{"Deville-2"}: }{A variance estimator for unequal-probability
#'   sampling without replacement, described in Matei and Tillé (2005)
#'   as "Deville 2".}
#'   \item{\strong{"SD1"}: }{The non-circular successive-differences variance estimator described by Ash (2014),
#'   sometimes used for variance estimation for systematic sampling.}
#'   \item{\strong{"SD2"}: }{The circular successive-differences variance estimator described by Ash (2014).
#'   This estimator is the basis of the "successive-differences replication" estimator commonly used
#'   for variance estimation for systematic sampling.}
#' }
#' @param max_replicates The maximum number of replicates to allow (should be as large as possible, given computer memory/storage limitations).
#' A commonly-recommended default is 500. If the number of replicates needed
#' for a balanced, fully-efficient estimator is less than \code{max_replicates},
#' then only the number of replicates needed will be created.
#' If more replicates are needed than \code{max_replicates}, then the full number of replicates
#' needed will be created, but only a random subsample will be retained.
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
#' @param mse If \code{TRUE}, compute variances from sums of squares around the point estimate from the full-sample weights,
#' If \code{FALSE}, compute variances from sums of squares around the mean estimate from the replicate weights.
#' @return
#' A replicate design object, with class \code{svyrep.design}, which can be used with the usual functions,
#' such as \code{svymean()} or \code{svyglm()}.
#'
#' Use \code{weights(..., type = 'analysis')} to extract the matrix of replicate weights.
#'
#' Use \code{as_data_frame_with_weights()} to convert the design object to a data frame with columns
#' for the full-sample and replicate weights.
#' @export
#' @seealso
#' For greater customization of the method, \code{\link[svrep]{make_quad_form_matrix}} can be used to
#' represent several common variance estimators as a quadratic form's matrix,
#' which can then be used as an input to \code{\link[svrep]{make_fay_gen_rep_factors}}.
#'
#' See \link[svrep]{variance-estimators} for a
#' description of each variance estimator.
#' @section Statistical Details:
#' See Fay (1989) for a full description of this replication method,
#' or see the documentation in \link[svrep]{make_fay_gen_rep_factors} for implementation details.
#'
#' Let \eqn{v( \hat{T_y})} be the textbook variance estimator for an estimated population total \eqn{\hat{T}_y} of some variable \eqn{y}.
#' The base weight for case \eqn{i} in our sample is \eqn{w_i}, and we let \eqn{\breve{y}_i} denote the weighted value \eqn{w_iy_i}.
#' Suppose we can represent our textbook variance estimator as a quadratic form: \eqn{v(\hat{T}_y) = \breve{y}\Sigma\breve{y}^T},
#' for some \eqn{n \times n} matrix \eqn{\Sigma}.
#' The only constraint on \eqn{\Sigma} is that, for our sample, it must be symmetric and positive semidefinite.
#'
#' The replication process creates \eqn{B} sets of replicate weights, where the \eqn{b}-th set of replicate weights is a vector of length \eqn{n} denoted \eqn{\mathbf{a}^{(b)}}, whose \eqn{k}-th value is denoted \eqn{a_k^{(b)}}.
#' This yields \eqn{B} replicate estimates of the population total, \eqn{\hat{T}_y^{*(b)}=\sum_{k \in s} a_k^{(b)} \breve{y}_k}, for \eqn{b=1, \ldots B}, which can be used to estimate sampling variance.
#'
#' \deqn{
#'   v_B\left(\hat{T}_y\right)=\frac{\sum_{b=1}^B\left(\hat{T}_y^{*(b)}-\hat{T}_y\right)^2}{B}
#' }
#'
#' This replicate variance estimator can be written as a quadratic form:
#'
#'   \deqn{
#'     v_B\left(\hat{T}_y\right) =\mathbf{\breve{y}}^{\prime}\Sigma_B \mathbf{\breve{y}}
#'   }
#'   where
#'   \deqn{
#'     \boldsymbol{\Sigma}_B = \frac{\sum_{b=1}^B\left(\mathbf{a}^{(b)}-\mathbf{1}_n\right)\left(\mathbf{a}^{(b)}-\mathbf{1}_n\right)^{\prime}}{B}
#'   }
#'
#' Note that if the vector of adjustment factors \eqn{\mathbf{a}^{(b)}} has expectation \eqn{\mathbf{1}_n} and variance-covariance matrix \eqn{\boldsymbol{\Sigma}},
#' then we have the across-replicate expectation \eqn{E_{*}\left( \boldsymbol{\Sigma}_B \right) = \boldsymbol{\Sigma}}. Since the bootstrap process takes the sample values \eqn{\breve{y}} as fixed, the bootstrap expectation of the variance estimator is \eqn{E_{*} \left( \mathbf{\breve{y}}^{\prime}\Sigma_B \mathbf{\breve{y}}\right)= \mathbf{\breve{y}}^{\prime}\Sigma \mathbf{\breve{y}}}.
#' Thus, we can produce a replication variance estimator with the same expectation as the textbook variance estimator simply by randomly generating \eqn{\mathbf{a}^{(b)}} from a distribution with the following two conditions:
#' \cr
#'     \strong{Condition 1}: \eqn{\quad \mathbf{E}_*(\mathbf{a})=\mathbf{1}_n}
#' \cr
#'     \strong{Condition 2}: \eqn{\quad \mathbf{E}_*\left(\mathbf{a}-\mathbf{1}_n\right)\left(\mathbf{a}-\mathbf{1}_n\right)^{\prime}=\mathbf{\Sigma}}
#' \cr \cr
#' With this function, the matrix \eqn{\boldsymbol{\Sigma}} undergoes a spectral decomposition,
#' and the eigenvectors and eigenvalues are combined using a Hadamard matrix (see Fay (1984)),
#' which produces 'balanced' replicates. See
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
#' \strong{Survey Methodology}, Statistics Canada, 40(1), 47–59.
#' \cr \cr
#' - Beaumont, Jean-François, and Zdenek Patak. 2012. “On the Generalized Bootstrap for Sample Surveys with Special Attention to Poisson Sampling: Generalized Bootstrap for Sample Surveys.” International Statistical Review 80 (1): 127–48. https://doi.org/10.1111/j.1751-5823.2011.00166.x.
#' \cr \cr
#' - Bertail, and Combris. 1997. “Bootstrap Généralisé d’un Sondage.” Annales d’Économie Et de Statistique, no. 46: 49. https://doi.org/10.2307/20076068.
#' \cr \cr
#' - Dippo, Cathryn, Robert Fay, and David Morganstein. 1984. “Computing Variances from Complex Samples with Replicate Weights.” In, 489–94. Alexandria, VA: American Statistical Association. http://www.asasrms.org/Proceedings/papers/1984_094.pdf.
#' \cr \cr
#' - Fay, Robert. 1984. “Some Properties of Estimates of Variance Based on Replication Methods.” In, 495–500. Alexandria, VA: American Statistical Association. http://www.asasrms.org/Proceedings/papers/1984_095.pdf.
#' \cr \cr
#' - Fay, Robert. 1989. “Theory And Application Of Replicate Weighting For Variance Calculations.” In, 495–500. Alexandria, VA: American Statistical Association. http://www.asasrms.org/Proceedings/papers/1989_033.pdf
#' \cr \cr
#' - Matei, Alina, and Yves Tillé. (2005).
#' “\emph{Evaluation of Variance Approximations and Estimators
#' in Maximum Entropy Sampling with Unequal Probability and Fixed Sample Size.}”
#' \strong{Journal of Official Statistics}, 21(4):543–70.
#' @examples
#' library(survey)
#'
#' ## Load an example systematic sample ----
#' data('library_stsys_sample', package = 'svrep')
#'
#' ## First, ensure data are sorted in same order as was used in sampling
#' library_stsys_sample <- library_stsys_sample[
#'   order(library_stsys_sample$SAMPLING_SORT_ORDER),
#' ]
#'
#' ## Create a survey design object
#' design_obj <- svydesign(
#'   data = library_stsys_sample,
#'   strata = ~ SAMPLING_STRATUM,
#'   ids = ~ 1,
#'   fpc = ~ STRATUM_POP_SIZE
#' )
#'
#' ## Convert to generalized replicate design
#' gen_boot_design_sd2 <- as_gen_boot_design(
#'   design = design_obj,
#'   variance_estimator = "SD2",
#'   replicates = 250, exact_vcov = TRUE
#' )
#' gen_rep_design_sd2 <- as_fays_gen_rep_design(
#'   design = design_obj,
#'   variance_estimator = "SD2",
#'   max_replicates = 250,
#'   mse = TRUE
#' )
#'
#' svytotal(x = ~ TOTSTAFF, na.rm = TRUE, design = gen_boot_design_sd2)
#' svytotal(x = ~ TOTSTAFF, na.rm = TRUE, design = gen_rep_design_sd2)
#'
#' svyquantile(x = ~ LIBRARIA, quantiles = 0.5, na.rm = TRUE,
#'             design = gen_boot_design_sd2, interval.type = "quantile")
#' svyquantile(x = ~ LIBRARIA, quantiles = 0.5, na.rm = TRUE,
#'             design = gen_rep_design_sd2)
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

# ##_ "Large" sample ----
#
# design_obj <- svydesign(data = svrep::lou_vax_survey,
#                         ids = ~ 1)
#
# ### Random groups jackknife
#   group_randomly <- function(n, n_groups) {
#     sample(as.numeric(cut(seq_len(n), breaks = n_groups)),
#            size = n, replace = FALSE)
#   }
#   ran_grp_jkn_design <- svrep::lou_vax_survey |>
#     transform(random_group = group_randomly(n = 1000,
#                                             n_groups = 50)) |>
#     svydesign(data = _, ids = ~ random_group) |>
#     as.svrepdesign(type = "JK1", mse = TRUE)
#
# ### Fay's generalized replication method
#   gen_rep_design <- as_fays_gen_rep_design(
#     design = design_obj,
#     variance_estimator = "Ultimate Cluster",
#     mse = TRUE,
#     max_replicates = 50
#   )
#
#   svytotal(x = ~ VAX_STATUS, na.rm = TRUE,
#            design = design_obj) |> SE()
#
#   gen_boot_design <- as_gen_boot_design(
#     design = design_obj,
#     variance_estimator = "Ultimate Cluster",
#     mse = TRUE, exact_vcov = FALSE,
#     replicates = 50
#   )
#
#   svytotal(x = ~ VAX_STATUS, na.rm = TRUE,
#            design = gen_boot_design) |> SE()
#
#
#   rep_results <- replicate(n = 500, expr = {
#     ran_grp_jkn_design <- svrep::lou_vax_survey |>
#       transform(random_group = group_randomly(n = 1000,
#                                               n_groups = 20)) |>
#       svydesign(data = _, ids = ~ random_group) |>
#       as.svrepdesign(type = "JK1", mse = TRUE)
#     gen_rep_design <- as_fay_gen_rep_design(
#       design = design_obj,
#       variance_estimator = "Ultimate Cluster",
#       mse = TRUE,
#       max_replicates = 20
#     )
#     list(
#       'jk' = svymean(x = ~ I(as.numeric(VAX_STATUS == "Vaccinated")), na.rm = TRUE,
#                       design = ran_grp_jkn_design),
#       'genrep' = svymean(x = ~ I(as.numeric(VAX_STATUS == "Vaccinated")), na.rm = TRUE,
#                           design = gen_rep_design)
#     ) |> sapply(SE)
#   })
#
# exp_value <- (svymean(x = ~ I(as.numeric(VAX_STATUS == "Vaccinated")), na.rm = TRUE,
#                       design = design_obj) |> SE())
#
# ((rep_results - (rep(exp_value, 2))^2)) |>
#   rowMeans() |> sqrt()
#
#
# ## Simulation ----
#
# rep_results <- replicate(
#   n = 2000, expr = {
#     suppressMessages({
#       gen_rep_design_sd2 <- as_fay_gen_rep_design(
#         design = design_obj,
#         variance_estimator = "SD2",
#         max_replicates = 50,
#         mse = TRUE
#       )
#     })
#     svytotal(x = ~ TOTSTAFF,
#              na.rm = TRUE,
#              design = gen_rep_design_sd2) |> SE()
#   }
# )
#
# mean_result <- rep_results |> mean()
# sd_results <- sd(rep_results)
# sd_results/mean_result
