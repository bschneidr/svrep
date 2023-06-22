#' @title Optimization-based Bootstrap Weights
#' @description Create bootstrap relpicate weights
#' using an optimization algorithm.
#' @param Sigma The matrix of the quadratic form used to represent the variance estimator. Must be positive semidefinite.
#' @param num_replicates The number of bootstrap replicates to create.
#' @param max_iter The maximum iterations to allow for the optimization algorithm.
#' @param max_loss The maximum loss to allow for the optimization algorithm.
#' @param torch_optimizer An optimization function from the 'torch' package,
#' such as \code{optim_adam}.
#' @param ... Optional arguments to pass to the optimization function
#' supplied by \code{torch_optimizer}.
#' @param .verbose Whether to periodically show the value of the loss function
#' every several iterations.
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
    max_iter = 15000,
    max_loss = 0.005,
    torch_optimizer = \(params, ...) torch::optim_adam(params, lr = 0.01, ...),
    ...,
    .verbose = TRUE
) {

  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("The `torch` package must be installed to use this function.")
  }

  n <- nrow(Sigma)
  Sigma_tensor = Sigma |> as.matrix() |> torch::torch_tensor()

  # log_A The logarithms of the matrix A we are trying to optimize
  # target_sigma The target quadratic form matrix
  loss_fn <- function(log_A, target_sigma) {
    A = torch::torch_exp(log_A)
    # Loss consists of:
    # (1) Frobenius distance
    #     between weights' covariance matrix
    #     and the target covariance matrix
    cov_loss = (
      A$cov(correction = 0) - target_sigma
    )$square()$sum()

    # (2) Euclidean distance from mean across replicates
    #     and '1'
    mean_loss = (A$t()$mean(1) - torch::torch_ones(n))$square()$sum()

    # Total loss is the sum of (1) and (2)
    sse = cov_loss + (n^(-1))*mean_loss
    return(sse)
  }

  # Create an initial solution
  initial_solution <- make_gen_boot_factors(
    Sigma = Sigma,
    num_replicates = num_replicates,
    tau = 'auto', exact_vcov = TRUE
  )

  log_A = torch::torch_tensor(log(initial_solution),
                              requires_grad = TRUE)

  # Create an optimizer
  optimizer = torch_optimizer(params = log_A, ...)

  # Iteratively update the matrix and evaluate loss

  # Initialize loss and iteration index
  iteration <- 1L
  loss <- loss_fn(log_A, Sigma_tensor)

  # Conduct the loop
  while ( (as.numeric(loss) > max_loss) & (iteration < max_iter) ) {

    # Print information for every 10-th iteration
    if (.verbose) {
      if ((iteration %% 500) == 0) {
        cat("Iteration: ", iteration, "   Loss: ", loss$item(), "\n")
      }
    }

    # Set the gradient to zero
    optimizer$zero_grad()

    # Determine loss at current iteration
    loss = loss_fn(log_A, Sigma_tensor)

    # Calculate gradients
    loss$backward()

    # Update `log_A`
    optimizer$step()

    iteration <- iteration + 1L
  }

  # Check convergence
  converged <- as.numeric(loss) < max_loss
  if (!converged) {
    warning_msg <- sprintf(
      paste0("After %s iterations, convergence was not achieved.",
             "The final value of the loss function is %s,",
             "but the specified convergence criterion is `loss < %s`"),
      max_iter, loss, max_loss
    )
    warning(warning_msg)
  }

  # Create the matrix of replicates by exponentiating log_A
  A <- log_A |> torch::torch_exp() |> as.matrix()
  attr(A, 'converged') <- converged

  return(A)
}
