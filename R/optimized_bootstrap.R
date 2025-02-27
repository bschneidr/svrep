#' @title Optimization-based Bootstrap Weights
#' @description Create bootstrap replicate weights
#' using an optimization algorithm.
#' @param Sigma The matrix of the quadratic form used to represent the variance estimator. Must be positive semidefinite.
#' @param num_replicates The number of bootstrap replicates to create.
#' @param max_iter The maximum iterations to allow for the optimization algorithm.
#' @param max_loss The maximum loss to allow for the optimization algorithm.
#' @param torch_optimizer An optimization function from the 'torch' package,
#' such as \code{optim_ignite_adam}.
#' @param ... Optional arguments to pass to the optimization function
#' supplied by \code{torch_optimizer}.
#' @param .verbose Whether to periodically show the value of the loss function
#' every several iterations of the optimization algorithm.
#'
#' @return
#' A matrix with the same number of rows as \code{Sigma},
#' and the number of columns equal to \code{num_replicates}.
#' @export
#'
#' @examples
#' library(survey)
#'
#' ## Load example data
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
#' ## Obtain quadratic form
#' quad_form_matrix <- get_design_quad_form(
#'   design = design_obj,
#'   variance_estimator = "SD2"
#' )
#'
#' ## Make optimization-based bootstrap replicates
#' replication_factors <- make_optim_boot_factors(
#'   Sigma = quad_form_matrix,
#'   max_iter = 15000,
#'   max_loss = 0.005
#' )
#'
#' ## Use the replicate factors for variance estimation
#' optim_boot_design <- svrepdesign(
#'   data = library_stsys_sample,
#'   repweights = replication_factors,
#'   type = "other",
#'   scale = 1/ncol(replication_factors),
#'   rscales = rep(1, times = ncol(replication_factors)),
#'   combined.weights = FALSE,
#'   weights = ~ I(1/SAMPLING_PROB)
#' )
#'
#' optim_boot_design |> svymean(x = ~ LIBRARIA, na.rm = TRUE)

make_optim_boot_factors <- function(
    Sigma,
    num_replicates = Matrix::rankMatrix(Sigma) + 1,
    max_iter = 20000,
    max_loss = 0.0001,
    torch_optimizer = \(params, ...) torch::optim_ignite_adamw(params, lr = 0.0005, ...),
    ...,
    .verbose = TRUE
) {

  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("The `torch` package must be installed to use this function.")
  }

  n <- nrow(Sigma)

  Sigma_rank <- Matrix::rankMatrix(Sigma)

  if (num_replicates <= Sigma_rank) {
    err_msg <- paste(
      "`make_optim_boot_factors` only works if `num_replicates` exceeds the rank of `Sigma`,",
      sprintf("which is %s", Sigma_rank)
    )
    stop(err_msg)
  }

  Sigma_tensor = Sigma |> as.matrix() |> torch::torch_tensor()

  constrain_sums <- function(X) {
    X = (X / X$mean(1))
    X = (X$t() / X$mean(2))$t()
    X = (X / X$mean(1))
    X = (X$t() / X$mean(2))$t()
    X
  }

  # log_A The logarithms of the matrix A we are trying to optimize
  # target_sigma The target quadratic form matrix
  loss_fn <- function(log_A, target_sigma) {
    A = log_A$clamp(0) |> constrain_sums()
    # Loss consists of:
    # (1) Squared Frobenius distance
    #     between weights' covariance matrix
    #     and the target covariance matrix
    cov_loss = (
      A$cov(correction = 0) - target_sigma
    )$square()$sum()

    # (2) Squared Euclidean distance from mean across replicates
    #     and '1'
    mean_loss = (A$t()$mean(1) - torch::torch_ones(n))$square()$sum()

    # Total loss is the sum of (1) and (2)
    sse = cov_loss + mean_loss
    return(sse)
  }

  # Create an initial solution
  initial_solution <- make_gen_boot_factors(
    Sigma = Sigma,
    num_replicates = num_replicates,
    tau = 1, exact_vcov = TRUE
  )
  initial_solution[sign(initial_solution) < 0] <- 0
  #initial_solution[initial_solution < 0] <- 1e-10
  
  log_A = torch::torch_tensor(initial_solution,
                              requires_grad = TRUE)

  # Create an optimizer
  optimizer = torch_optimizer(params = log_A, ...)

  # Iteratively update the matrix and evaluate loss

  # Initialize loss and iteration index
  iteration <- 1L
  loss <- loss_fn(log_A, Sigma_tensor)

  if (inherits(optimizer, "optim_lbfgs")) {
    calc_loss <- function() {
      optimizer$zero_grad()
      loss_value <- loss_fn(log_A, Sigma_tensor)
      loss_value$backward()
      loss_value
    }
  }

  # Conduct the loop
  while ( (as.numeric(loss) > max_loss) & (iteration < max_iter) ) {

    # Print information for every 10-th iteration
    if (.verbose) {
      if ((iteration %% 500) == 0) {
        cat("Iteration: ", iteration, "   Loss: ", loss$item(), "\n")
      }
    }

    if (inherits(optimizer, "optim_lbfgs")) {
      optimizer$step(calc_loss)
    } else {
      # Set the gradient to zero
      optimizer$zero_grad()

      # Determine loss at current iteration
      loss = loss_fn(log_A, Sigma_tensor)

      # Calculate gradients
      loss$backward()

      # Update `log_A`
      optimizer$step()
    }

    iteration <- iteration + 1L
  }

  # Check convergence
  converged <- as.numeric(loss) < max_loss
  if (!converged) {
    warning_msg <- sprintf(
      paste0("After %s iterations, convergence was not achieved.",
             "The final value of the loss function is %s,",
             "but the specified convergence criterion is `loss < %s`"),
      max_iter, loss$item(), max_loss
    )
    warning(warning_msg)
  }

  # Create the matrix of replicates by exponentiating log_A
  A <- log_A$clamp(0) |> constrain_sums() |> as.matrix()
  attr(A, 'converged') <- converged

  return(A)
}

#' @title Convert a survey design object to a optimization-based bootstrap replicate design
#' @description Converts a survey design object to a replicate design object
#' with replicate weights formed using the optimization-based bootstrap method.
#' The generalized survey bootstrap is a method for forming bootstrap replicate weights
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
#' @param replicates Number of bootstrap replicates (should be as large as possible, given computer memory/storage limitations).
#' A commonly-recommended default is 500.

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
#' @param max_iter The maximum iterations to allow for the optimization algorithm.
#' @param max_loss The maximum loss to allow for the optimization algorithm.
#' @param torch_optimizer An optimization function from the 'torch' package,
#' such as \code{optim_ignite_adam}.
#' @param ... Optional arguments to pass to the optimization function
#' supplied by \code{torch_optimizer}.
#' @param compress This reduces the computer memory required to represent the replicate weights and has no
#' impact on estimates.
#' @param mse If \code{TRUE}, compute variances from sums of squares around the point estimate from the full-sample weights,
#' If \code{FALSE}, compute variances from sums of squares around the mean estimate from the replicate weights.
#' @param .verbose Whether to periodically show the value of the loss function
#' every several iterations of the optimization algorithm.
#' @return
#' A replicate design object, with class \code{svyrep.design}, which can be used with the usual functions,
#' such as \code{svymean()} or \code{svyglm()}.
#'
#' Use \code{weights(..., type = 'analysis')} to extract the matrix of replicate weights.
#'
#' Use \code{as_data_frame_with_weights()} to convert the design object to a data frame with columns
#' for the full-sample and replicate weights.
#' @export
#' @seealso Use \code{\link[svrep]{estimate_boot_reps_for_target_cv}} to help choose the number of bootstrap replicates. \cr
#'
#' For greater customization of the method, \code{\link[svrep]{make_quad_form_matrix}} can be used to
#' represent several common variance estimators as a quadratic form's matrix,
#' which can then be used as an input to \code{\link[svrep]{make_optim_boot_factors}}.
#'
#' See \link[svrep]{variance-estimators} for a
#' description of each variance estimator.
#' @section Statistical Details:
#' Let \eqn{v( \hat{T_y})} be the textbook variance estimator for an estimated population total \eqn{\hat{T}_y} of some variable \eqn{y}.
#' The base weight for case \eqn{i} in our sample is \eqn{w_i}, and we let \eqn{\breve{y}_i} denote the weighted value \eqn{w_iy_i}.
#' Suppose we can represent our textbook variance estimator as a quadratic form: \eqn{v(\hat{T}_y) = \breve{y}\Sigma\breve{y}^T},
#' for some \eqn{n \times n} matrix \eqn{\Sigma}.
#' The only constraint on \eqn{\Sigma} is that, for our sample, it must be symmetric and positive semidefinite.
#'
#' The bootstrapping process creates \eqn{B} sets of replicate weights, where the \eqn{b}-th set of replicate weights is a vector of length \eqn{n} denoted \eqn{\mathbf{a}^{(b)}}, whose \eqn{k}-th value is denoted \eqn{a_k^{(b)}}.
#' This yields \eqn{B} replicate estimates of the population total, \eqn{\hat{T}_y^{*(b)}=\sum_{k \in s} a_k^{(b)} \breve{y}_k}, for \eqn{b=1, \ldots B}, which can be used to estimate sampling variance.
#'
#' \deqn{
#'   v_B\left(\hat{T}_y\right)=\frac{\sum_{b=1}^B\left(\hat{T}_y^{*(b)}-\hat{T}_y\right)^2}{B}
#' }
#'
#' This bootstrap variance estimator can be written as a quadratic form:
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
#' then we have the bootstrap expectation \eqn{E_{*}\left( \boldsymbol{\Sigma}_B \right) = \boldsymbol{\Sigma}}. Since the bootstrap process takes the sample values \eqn{\breve{y}} as fixed, the bootstrap expectation of the variance estimator is \eqn{E_{*} \left( \mathbf{\breve{y}}^{\prime}\Sigma_B \mathbf{\breve{y}}\right)= \mathbf{\breve{y}}^{\prime}\Sigma \mathbf{\breve{y}}}.
#' Thus, we can produce a bootstrap variance estimator with the same expectation as the textbook variance estimator simply by randomly generating \eqn{\mathbf{a}^{(b)}} from a distribution with the following two conditions:
#' \cr
#'     \strong{Condition 1}: \eqn{\quad \mathbf{E}_*(\mathbf{a})=\mathbf{1}_n}
#' \cr
#'     \strong{Condition 2}: \eqn{\quad \mathbf{E}_*\left(\mathbf{a}-\mathbf{1}_n\right)\left(\mathbf{a}-\mathbf{1}_n\right)^{\prime}=\mathbf{\Sigma}}
#' \cr \cr
#' While there are multiple ways to generate adjustment factors satisfying these conditions,
#' the simplest general method is to simulate from a multivariate normal distribution: \eqn{\mathbf{a} \sim MVN(\mathbf{1}_n, \boldsymbol{\Sigma})}.
#' This is the method used by this function.

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
#' The generalized survey bootstrap was first proposed by Bertail and Combris (1997).
#' See Beaumont and Patak (2012) for a clear overview of the generalized survey bootstrap.
#' The generalized survey bootstrap represents one strategy for forming replication variance estimators
#' in the general framework proposed by Fay (1984) and Dippo, Fay, and Morganstein (1984).
#'
#'
#' - Beaumont, Jean-François, and Zdenek Patak. 2012. “On the generalized bootstrap for Sample Surveys with Special Attention to Poisson Sampling: optimization-based bootstrap for Sample Surveys.” International Statistical Review 80 (1): 127–48. https://doi.org/10.1111/j.1751-5823.2011.00166.x.
#' \cr \cr
#' - Bertail, and Combris. 1997. “Bootstrap Généralisé d’un Sondage.” Annales d’Économie Et de Statistique, no. 46: 49. https://doi.org/10.2307/20076068.
#' \cr \cr
#' - Dippo, Cathryn, Robert Fay, and David Morganstein. 1984. “Computing Variances from Complex Samples with Replicate Weights.” In, 489–94. Alexandria, VA: American Statistical Association. http://www.asasrms.org/Proceedings/papers/1984_094.pdf.
#' \cr \cr
#' - Fay, Robert. 1984. “Some Properties of Estimates of Variance Based on Replication Methods.” In, 495–500. Alexandria, VA: American Statistical Association. http://www.asasrms.org/Proceedings/papers/1984_095.pdf.
#' @export
#'
#' @examples
#' \dontrun{
#' library(survey)
#'
#'# Example 1: Bootstrap based on the Yates-Grundy estimator ----
#'    set.seed(2014)
#'
#'    data('election', package = 'survey')
#'
#'    ## Create survey design object
#'    pps_design_yg <- svydesign(
#'      data = election_pps,
#'      id = ~1, fpc = ~p,
#'      pps = ppsmat(election_jointprob),
#'      variance = "YG"
#'    )
#'
#'    ## Convert to optimization-based bootstrap replicate design
#'    optim_boot_design_yg <- pps_design_yg |>
#'      as_optim_boot_design(
#'        variance_estimator = "Yates-Grundy",
#'        replicates = 50,
#'        max_loss = 1e-7,
#'        torch_optimizer = \(params, ...) {
#'          torch::optim_ignite_adamw(params, lr = 0.01)
#'        }
#'      )
#'
#'    svytotal(x = ~ Bush + Kerry, design = pps_design_yg)
#'    svytotal(x = ~ Bush + Kerry, design = optim_boot_design_yg)
#'
#'# Example 2: Bootstrap based on the successive-difference estimator ----
#'
#'    data('library_stsys_sample', package = 'svrep')
#'
#'    ## First, ensure data are sorted in same order as was used in sampling
#'    library_stsys_sample <- library_stsys_sample[
#'      order(library_stsys_sample$SAMPLING_SORT_ORDER),
#'    ]
#'
#'    ## Create a survey design object
#'    design_obj <- svydesign(
#'      data = library_stsys_sample,
#'      strata = ~ SAMPLING_STRATUM,
#'      ids = ~ 1,
#'      fpc = ~ STRATUM_POP_SIZE
#'    )
#'
#'    ## Convert to optimization-based bootstrap replicate design
#'    optim_boot_design_sd2 <- as_optim_boot_design(
#'      design = design_obj,
#'      variance_estimator = "SD2",
#'      replicates = 2000
#'    )
#'
#'    ## Estimate sampling variances
#'    svytotal(x = ~ TOTSTAFF, na.rm = TRUE, design = optim_boot_design_sd2)
#'    svytotal(x = ~ TOTSTAFF, na.rm = TRUE, design = design_obj)
#'
#' # Example 3: Two-phase sample ----
#' # -- First stage is stratified systematic sampling,
#' # -- second stage is response/nonresponse modeled as Poisson sampling
#'
#'   nonresponse_model <- glm(
#'     data = library_stsys_sample,
#'     family = quasibinomial('logit'),
#'     formula = I(RESPONSE_STATUS == "Survey Respondent") ~ 1,
#'     weights = 1/library_stsys_sample$SAMPLING_PROB
#'   )
#'
#'   library_stsys_sample[['RESPONSE_PROPENSITY']] <- predict(
#'     nonresponse_model,
#'     newdata = library_stsys_sample,
#'     type = "response"
#'   )
#'
#'   twophase_design <- twophase(
#'     data = library_stsys_sample,
#'     # Identify cases included in second phase sample
#'     subset = ~ I(RESPONSE_STATUS == "Survey Respondent"),
#'     strata = list(~ SAMPLING_STRATUM, NULL),
#'     id = list(~ 1, ~ 1),
#'     probs = list(NULL, ~ RESPONSE_PROPENSITY),
#'     fpc = list(~ STRATUM_POP_SIZE, NULL),
#'   )
#'
#'   twophase_boot_design <- as_optim_boot_design(
#'     design = twophase_design,
#'     variance_estimator = list(
#'       "SD2", "Poisson Horvitz-Thompson"
#'     )
#'   )
#'
#'   svytotal(x = ~ LIBRARIA, design = twophase_boot_design)
#'
#' }
as_optim_boot_design <- function(design, variance_estimator = NULL,
                                 replicates = 500,
                                 psd_option = "warn",
                                 max_iter = 20000,
                                 max_loss = 0.0001,
                                 torch_optimizer = \(params, ...) torch::optim_ignite_adamw(params, lr = 0.0005, ...),
                                 ...,
                                 mse = getOption("survey.replicates.mse"),
                                 compress = TRUE,
                                 .verbose = TRUE) {
  UseMethod("as_optim_boot_design", design)
}

#' @export
as_optim_boot_design.twophase2 <- function(design, variance_estimator = NULL,
                                           replicates = 500,
                                           psd_option = "warn",
                                           max_iter = 20000,
                                           max_loss = 0.0001,
                                           torch_optimizer = \(params, ...) torch::optim_ignite_adamw(params, lr = 0.0005, ...),
                                           ...,
                                           mse = getOption("survey.replicates.mse"),
                                           compress = TRUE,
                                           .verbose = TRUE) {

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

  adjustment_factors <- make_optim_boot_factors(
    Sigma = Sigma,
    num_replicates = replicates,
    max_iter = max_iter,
    max_loss = max_loss,
    torch_optimizer = torch_optimizer,
    .verbose = .verbose
  )

  scale <- 1/ncol(adjustment_factors)
  rscales <- rep(1, times = ncol(adjustment_factors))

  rep_design <- survey::svrepdesign(
    variables = design$phase1$full$variables[design$subset,,drop=FALSE],
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
as_optim_boot_design.survey.design <- function(design, variance_estimator = NULL,
                                               replicates = 500,
                                               psd_option = "warn",
                                               max_iter = 20000,
                                               max_loss = 0.0001,
                                               torch_optimizer = \(params, ...) torch::optim_iginite_adamw(params, lr = 0.0005, ...),
                                               ...,
                                               mse = getOption("survey.replicates.mse"),
                                               compress = TRUE,
                                               .verbose = TRUE) {

    # Produce a (potentially) compressed survey design object
  compressed_design_structure <- compress_design(design, vars_to_keep = aux_var_names)

    # Get the quadratic form of the variance estimator,
    # for the compressed design object
  Sigma <- get_design_quad_form(
    compressed_design_structure$design_subset,
    variance_estimator,
    ...
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
  adjustment_factors <- make_optim_boot_factors(
    Sigma = Sigma,
    num_replicates = replicates,
    max_iter = max_iter,
    max_loss = max_loss,
    torch_optimizer = torch_optimizer,
    .verbose = .verbose
  )

  scale <- 1/ncol(adjustment_factors)
  rscales <- rep(1, times = ncol(adjustment_factors))

  # Uncompress the adjustment factors
  adjustment_factors <- distribute_matrix_across_clusters(
    cluster_level_matrix = adjustment_factors,
    cluster_ids = compressed_design_structure$index,
    rows = TRUE, cols = FALSE
  )

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
