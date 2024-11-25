#' @title Create a row-assignment matrix for Successive Difference Replication-Method
#' @description Creates a row assignment matrix:
#' for each row of a dataset of size \eqn{n}, assigns two rows
#' of a Hadamard matrix.
#' @param n The sample size of the data.
#' @param hadamard_order The order of the Hadamard matrix (i.e., the number of rows/columns)
#' @param number_of_cycles The number of cycles to use in the row assignment.
#' Must be at least as large as \code{n/hadamard_order}. Only
#' applies when \code{n} exceeds the number of available rows
#' in the Hadamard matrix. The number of available rows
#' is \code{hadamard_order} when \code{use_first_row = TRUE},
#' and \code{hadamard_order - 1} when \code{use_first_row = FALSE}.
#' @param use_first_row Whether to use the first row of the Hadamard matrix.
#' The first row of a Hadamard matrix is often all 1's, and so using the first
#' row to create replicate factors leads to the creation of a replicate
#' whose weights exactly match the full-sample weights. Thus, using the first row
#' of the Hadamard matrix may be undesirable for practical purposes,
#' even if it is valid for the purpose of variance estimation.
#' @param circular \code{TRUE} or \code{FALSE}. Only applies
#' when the number of available rows in the Hadamard matrix 
#' is at least as large as \code{n}.
#' Whether to make a circular row assignment,
#' so that the resulting successive-difference replication variance estimator
#' is equivalent to the SD2 variance estimator rather than the SD1 variance
#' estimator (see Ash 2014). The number of available rows
#' is \code{hadamard_order} when \code{use_first_row = TRUE},
#' and \code{hadamard_order - 1} when \code{use_first_row = FALSE}.
#' @details
#' Implements row-assignment methods described in Ash (2014)
#' and in Fay and Train (1995). The row-assignment method
#' depends on the number of available rows of the Hadamard matrix used.
#' The number of available rows
#' is \code{hadamard_order} when \code{use_first_row = TRUE},
#' and \code{hadamard_order - 1} when \code{use_first_row = FALSE}.
#' 
#' When the number of available Hadamard rows is at least as large as \code{n},
#' then the row assignment is as follows. Let \eqn{i=1,\dots,n} be the index
#' for the data that will receive row assignments, 
#' and let \eqn{a_j} denote row \eqn{j} of a Hadamard matrix. 
#' For \eqn{i < n}, the assignments
#' are entries \eqn{a_i} and \eqn{a_{i+1}} when \code{use_first_row = TRUE},
#' and when \code{use_first_row = FALSE}, the assignments are entries
#' \eqn{a_{i+1}} and \eqn{a_{i+2}}. The assignment for \eqn{i=n} depends on 
#' whether \code{circular = TRUE}. If \eqn{circular=TRUE}, then the assignment
#' for \eqn{i=n} is entries \eqn{a_i} and \eqn{a_1} if \code{use_first_row = TRUE},
#' and entries \eqn{a_{i+1}} and \eqn{a_2} if \code{use_first_row = FALSE}.
#' 
#' When the number of available Hadamard rows is less than \code{n},
#' then the row assignment method is the method denoted as RA1 in Ash (2014).
#' This method uses the argument \code{number_of_cycles} and 
#' does \emph{not} use the argument \code{circular}.
#' 
#' @return A matrix with \code{n} rows and two columns. Each row
#' gives the assignment of two rows of a Hadamard matrix to the row of data.
#' @references
#' Ash, S. (2014). "\emph{Using successive difference replication for estimating variances}."
#' \strong{Survey Methodology}, Statistics Canada, 40(1), 47–59.
#' @examples
#' # Assign rows of a 4x4 Hadamard matrix
#' # to 9 observations, using two cycles
#' assign_hadamard_rows(
#'   n = 9,
#'   hadamard_order = 4,
#'   number_of_cycles = 2,
#'   circular = TRUE
#' )
#'
assign_hadamard_rows <- function(
  n, hadamard_order, 
  number_of_cycles = ceiling(n/hadamard_order), 
  use_first_row    = TRUE,
  circular         = TRUE
) {

  if ((hadamard_order %% 4) > 0) {
    stop("`hadamard_order` must be a multiple of 4.")
  }

  k_A <- hadamard_order

  if (!use_first_row) {
    k_A <- k_A - 1
  }
  if (number_of_cycles < ceiling(n/k_A)) {
    error_msg <- sprintf("Given the values of `n`, `hadamard_order`, and `use_first_row`, `number_of_cycles` must be at least %s",
                         ceiling(n/k_A))
  }

  # Determine range of step sizes to use to form row-assignment matrix.
  # Cannot be larger than (order of Hadamard matrix minus one).
  if (((k_A - 1) * number_of_cycles) %% 2) {
    max_step_size <- number_of_cycles + 1L
    if (max_step_size >= (k_A - 1)) {
      max_step_size <- number_of_cycles - 1L
    }
    max_step_size <- pmin(max_step_size, k_A - 1)
  } else {
    max_step_size <- pmin(number_of_cycles, k_A - 1)
  }

  # Set initial values before looping
  row_assignments <- NULL

  start_point <- 1L
  i <- start_point
  unused_i <- seq(from = start_point,
                  to = k_A, by = 1)

  # Loop over step sizes, forming all cyclical row assignments
  # available at a given step size
  for (step_size in 1:max_step_size) {
    while (length(unused_i) > 0) {
      j <- i + step_size
      if (j > k_A) {
        if (start_point == 1) {
          j <- (j - k_A)
        } else {
          j <- (j - k_A) + 1L
        }
      }
      if (i %in% unused_i) {
        row_assignments <- c(row_assignments, c(i,j))
        unused_i <- setdiff(unused_i, i)
        i <- j
      } else {
        i <- i + 1L
      }
    }
    unused_i <- seq(from = start_point,
                    to = k_A, by = 1)
    i <- start_point
  }

  if (!use_first_row) {
    row_assignments <- row_assignments + 1L
  }

  # Compile the row assignments into a matrix
  row_assignment_matrix <- matrix(row_assignments,
                                  ncol = 2, byrow = TRUE,
                                  dimnames = list(NULL, c("row_1", "row_2")))

  # If the row assignment matrix doesn't cover all of the sample,
  # repeat rows from the assignment matrix as necessary
  if (nrow(row_assignment_matrix) < n) {
    row_assignment_matrix <- rep(as.vector(t(row_assignment_matrix)),
                                 times = ceiling( n/nrow(row_assignment_matrix))) |>
      matrix(ncol = 2, byrow = TRUE,
             dimnames = list(NULL, c("row_1", "row_2")))
    row_assignment_matrix <- row_assignment_matrix[seq_len(n),]
  }
  # If the row assignment matrix has too many rows,
  # only use the number of rows necessary
  if (nrow(row_assignment_matrix) > n) {
    row_assignment_matrix <- row_assignment_matrix[seq_len(n),]
  }
  if ((n <= k_A) && circular) {
    row_assignment_matrix[n,2] <- row_assignment_matrix[1,1] 
  }

  return(row_assignment_matrix)
}

#' Create matrix of replicate factors to use for Successive Difference Replication Method
#'
#' @param n The number of sampling units
#' @param target_number_of_replicates The target number of replicates to create.
#' This will determine the order of the Hadamard matrix to use when
#' creating replicate factors. The actual number of replicates will
#' be a multiple of 4.
#' If \code{use_normal_hadamard = FALSE}, then the actual number of replicates
#' will be a power of 4, which means that the actual number of
#' replicates might be much larger than the target.
#' @param use_normal_hadamard Whether to use a normal Hadamard matrix:
#' that is, a matrix whose first row and first column only have entries
#' equal to 1.
#' @return
#' A matrix of replicate factors, with \code{n} rows
#' and the number of columns corresponding to the order of the Hadamard matrix used
#' (which is greater than or equal to \code{target_number_of_replicates}).
#' @export
#'
#' @examples
#'
#' # Note that one of the replicates has every factor equal to 1
#' # Also note that this matches Table 1 in Ash (2014)
#' make_sdr_replicate_factors(
#'   n = 4,
#'   target_number_of_replicates = 4,
#'   use_normal_hadamard = TRUE
#' )
#'
#' # Note the difference when using a non-normal Hadamard matrix
#' rep_factors <- make_sdr_replicate_factors(
#'   n = 4,
#'   target_number_of_replicates = 4,
#'   use_normal_hadamard = FALSE
#' )
#'
#' # These replicate factors are equivalent
#' # to the SD2 variance estimator
#' tcrossprod(rep_factors - 1)
#'
#' # Compare to the quadratic form of the SD2 estimator
#' sd1_quad_form <- make_quad_form_matrix(
#'   variance_estimator = "SD1",
#'   cluster_ids        = matrix(1:4, ncol = 1),
#'   sort_order         = matrix(1:4, ncol = 1)
#' )
#'
make_sdr_replicate_factors <- function(n, target_number_of_replicates, use_normal_hadamard = FALSE) {

  # Create Hadamard matrix to use for replicate factors
  if (!use_normal_hadamard) {
    H_4 <- matrix(
      c(1,-1,1,1,
        -1,-1,-1,1,
        1,-1,-1,-1,
        1,1,-1,1),
      nrow = 4, ncol = 4,
      byrow = TRUE
    )

    H_A <- H_4
    while (ncol(H_A) < target_number_of_replicates) {
      H_A <- rbind(cbind(H_A, H_A),
                   cbind(H_A, -H_A))
    }
  }
  if (use_normal_hadamard) {
    H_A <- survey::hadamard(target_number_of_replicates - 1)
    H_A <- 2*H_A - 1 # Convert from 1/0 format to 1/-1 format
  }
  hadamard_order <- ncol(H_A)
  smallest_hadamard_order <- find_minimum_hadamard_order(target_number_of_replicates)
  sprintf("Using Hadamard matrix of order %s. If `use_normal_hadamard=TRUE`, the smallest possible order is %s.",
          hadamard_order, smallest_hadamard_order) |>
    message()


  # Assign rows of the Hadamard matrix to each observation

  row_assignment_matrix <- assign_hadamard_rows(
    n                = n,
    hadamard_order   = hadamard_order,
    number_of_cycles = ceiling(n / hadamard_order),
    use_first_row    = TRUE,
    circular         = TRUE
  )

  # Create replicate factors based on the Hadamard matrix and the row assignments
  replicate_factors <- sapply(X = 1L:hadamard_order, FUN = function(r) {
    hadamard_entries <- cbind('h_1r' = H_A[row_assignment_matrix[,'row_1'], r],
                              'h_2r' = H_A[row_assignment_matrix[,'row_2'], r])

    1 + (hadamard_entries %*% c(1,-1) * (2^(-3/2)))
  }, simplify = TRUE)

  return(replicate_factors)
}

#' @title Convert a survey design object to a successive differences replicate design
#' @description
#' Converts a survey design object to a replicate design object
#' with replicate weights formed using the successive differences replication (SDR) method.
#' The SDR method is suitable for designs that use
#' systematic sampling or finely-stratified sampling designs.
#' The rows of data in the survey design should be sorted
#' in the same order that was used for sampling.
#' @param design A survey design object created using the 'survey' (or 'srvyr') package,
#' with class \code{'survey.design'} or \code{'svyimputationList'}.
#' @param replicates The target number of replicates to create.
#' This will determine the order of the Hadamard matrix to use when
#' creating replicate factors.
#' If \code{use_normal_hadamard = TRUE}, then the actual number of replicates will be the
#' smallest \emph{multiple} of 4 that is greater or equal to the specified value of \code{replicates}.
#' If \code{use_normal_hadamard = FALSE}, then the actual number of replicates will be the
#' smallest \emph{power} of 4 that is greater or equal to the specified value of \code{replicates}.
#' @param sort_variable To create SDR replicates, the data must
#' be sorted in the order used in sampling. The name of a sorting variable
#' can be supplied as a character string. If no variable name is supplied,
#' then this function assumes that the data are already sorted into the correct order.
#' @param use_normal_hadamard Whether to use a normal Hadamard matrix:
#' that is, a matrix whose first row and first column only have entries
#' equal to 1.
#' @param compress Use a compressed representation of the replicate weights matrix.
#' This reduces the computer memory required to represent the replicate weights and has no
#' impact on estimates.
#' @param mse If \code{TRUE}, compute variances from sums of squares around the point estimate from the full-sample weights,
#' If \code{FALSE}, compute variances from sums of squares around the mean estimate from the replicate weights.
#' @references
#' Ash, S. (2014). "\emph{Using successive difference replication for estimating variances}."
#' \strong{Survey Methodology}, Statistics Canada, 40(1), 47–59.
#'
#' Fay, R.E. and Train, G.F. (1995). "\emph{Aspects of Survey and Model-Based Postcensal Estimation of
#' Income and Poverty Characteristics for States and Counties}." Joint Statistical Meetings,
#' Proceedings of the Section on Government Statistics, 154-159.
#' @return
#' A replicate design object, with class \code{svyrep.design}, which can be used with the usual functions,
#' such as \code{svymean()} or \code{svyglm()}.
#'
#' Use \code{weights(..., type = 'analysis')} to extract the matrix of replicate weights. \cr
#' Use \code{as_data_frame_with_weights()} to convert the design object to a data frame with columns
#' for the full-sample and replicate weights.
#' @export
#' @export
#'
#' @examples
#'
#' # Load example stratified systematic sample
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
#' ## Convert to SDR replicate design
#' sdr_design <- as_sdr_design(
#'   design = design_obj,
#'   replicates = 180,
#'   sort_variable = "SAMPLING_SORT_ORDER",
#'   use_normal_hadamard = TRUE
#' )
#'
#' ## Compare to generalized bootstrap
#' ## based on the SD2 estimator that SDR approximates
#' gen_boot_design <- as_gen_boot_design(
#'   design = design_obj,
#'   variance_estimator = "SD2",
#'   replicates = 180,
#'   exact_vcov = TRUE
#' )
#'
#' ## Estimate sampling variances
#' svytotal(x = ~ TOTSTAFF, na.rm = TRUE, design = sdr_design)
#' svytotal(x = ~ TOTSTAFF, na.rm = TRUE, design = gen_boot_design)

as_sdr_design <- function(
  design,
  replicates,
  sort_variable = NULL,
  use_normal_hadamard = FALSE,
  compress = TRUE,
  mse = getOption("survey.replicates.mse")
) {
  UseMethod("as_sdr_design", design)
}

#' @export
as_sdr_design.survey.design <- function(
    design,
    replicates,
    sort_variable = NULL,
    use_normal_hadamard = FALSE,
    compress = TRUE,
    mse = getOption("survey.replicates.mse")
) {

  if (is.null(sort_variable)) {
    message(
      "Since `sort_variable = NULL`, assuming rows of data are sorted in the same order used in sampling."
    )
  }

  # Produce a (potentially) compressed survey design object
  compressed_design_structure <- compress_design(design, vars_to_keep = sort_variable)


  # Create adjustment factors for the compressed design structure
  adjustment_factors <- make_sdr_replicate_factors(
    n = nrow(compressed_design_structure$design_subset),
    target_number_of_replicates = replicates,
    use_normal_hadamard = use_normal_hadamard
  )

  if (!is.null(sort_variable)) {
    sort_order <- rank(design$variables[[sort_variable]], ties.method = 'first')

    sort_data <- data.frame(
      ORIG_ORDER    = seq_len(nrow(compressed_design_structure$design)),
      PSU_STRATUM   = compressed_design_structure$design$strata[,1,drop=TRUE], 
      SORT_VARIABLE = compressed_design_structure$design$variables[[sort_variable]]
    )

    sort_data <- sort_data |> sort_by(~ PSU_STRATUM + SORT_VARIABLE)
    sort_data[['SORT_ORDER']] <- seq_len(nrow(sort_data))
    sort_data <- sort_data |> sort_by(~ ORIG_ORDER)
    adjustment_factors <- adjustment_factors[sort_data[['SORT_ORDER']], , drop = FALSE]
  }

  # Uncompress the adjustment factors
  adjustment_factors <- distribute_matrix_across_clusters(
    cluster_level_matrix = adjustment_factors,
    cluster_ids = compressed_design_structure$index,
    rows = TRUE, cols = FALSE
  )

  # Create a replicate survey design object
  rep_design <- survey::svrepdesign(
    variables = design$variables,
    weights = weights(design, type = "sampling"),
    repweights = adjustment_factors, combined.weights = FALSE,
    compress = compress, mse = mse,
    type = "successive-difference"
  )

  # Return the result
  if (inherits(design, 'tbl_svy') && ('package:srvyr' %in% search())) {
    rep_design <- srvyr::as_survey_rep(rep_design)
  }

  rep_design$call <- sys.call(which = -1)

  return(rep_design)
}
