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
#' \strong{Survey Methodology}, Statistics Canada, 40(1), 47-59.
#' @keywords internal
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
#' will be \eqn{4 \times 2^k} for some integer \eqn{k}, 
#' which means that the actual number of replicates might be much larger than the target.
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
#' print(rep_factors)
#'
#' # These replicate factors are equivalent
#' # to the SD2 variance estimator
#' tcrossprod(rep_factors - 1)
#'
#' # Compare to the quadratic form of the SD2 estimator
#' sd2_quad_form <- make_quad_form_matrix(
#'   variance_estimator = "SD2",
#'   cluster_ids        = matrix(1:4, ncol = 1),
#'   sort_order         = matrix(1:4, ncol = 1)
#' )
#' print(sd2_quad_form)
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
#' @param design A survey design object created using the 'survey' (or 'srvyr') package,
#' with class \code{'survey.design'} or \code{'svyimputationList'}.
#' @param replicates The target number of replicates to create.
#' This will determine the order of the Hadamard matrix to use when
#' creating replicate factors.
#' If \code{use_normal_hadamard = TRUE}, then the actual number of replicates will be
#' greater than or equal to \code{replicates} and determined by identifying
#' the smallest available Hadamard matrix available from the 'survey' package.
#' If \code{use_normal_hadamard = FALSE}, then the actual number of replicates will be the
#' smallest \emph{power} of 4 that is greater or equal to the specified value of \code{replicates}.
#' @param sort_variable A character string specifying the name
#' of a sorting variable. This variable should give
#' the sort order used in sampling. If the design includes strata,
#' then the replicate factors will be assigned after first sorting by the 
#' first-stage strata identifier
#' and then sorting by the value of \code{sort_variable} within each stratum.
#' @param use_normal_hadamard Whether to use a normal Hadamard matrix:
#' that is, a matrix whose first row and first column only have entries
#' equal to 1. This means that one of the replicates will be an "inactive" replicate.
#' See the "Details" section for more information.
#' @param compress Use a compressed representation of the replicate weights matrix.
#' This reduces the computer memory required to represent the replicate weights and has no
#' impact on estimates.
#' @param mse If \code{TRUE}, compute variances from sums of squares around the point estimate from the full-sample weights,
#' If \code{FALSE}, compute variances from sums of squares around the mean estimate from the replicate weights.
#' @section Statistical Overview:
#' The successive difference replication method was proposed by Fay and Train (1995)
#' as a replication method appropriate for samples selected using systematic sampling.
#' It is designed to yield variance estimates for totals that are equivalent to
#' successive difference variance estimators described in Fay and Train (1995).
#' There are different methods for forming the replicate factors depending on
#' whether the replicate variance estimator is meant to be equivalent to the 
#' SD2 variance estimator (i.e., the circular successive difference estimator) 
#' or the SD1 variance estimator (the non-circular successive difference estimator) 
#' described in Ash (2014). This function uses the approach based on the SD2 variance estimator.
#' For multistage designs, this replication method only takes into account information about the
#' first stage of sampling.
#' 
#' The scale factor to be used for variance estimation with the replicate weights
#' is \eqn{4/R}, where \eqn{R} is the number of replicates. This scale factor will 
#' be used even when there are finite population corrections; see the subsection below.
#' 
#' As an alternative to the successive difference replication estimator,
#' one can use a generalized replication method where the target variance estimator
#' is the "SD1" or "SD2" estimator. See the functions \link[svrep]{as_gen_boot_design}
#' or \link[svrep]{as_fays_gen_rep_design} for more details on generalized replication
#' and see the help section \link[svrep]{variance-estimators} for more details
#' on the "SD1" and "SD2" variance estimators.
#' 
#' @section Details on Stratification and Finite Population Corrections:
#' 
#' If the design includes strata,
#' then the replicate factors will be assigned after first sorting by the 
#' first-stage strata identifier and then sorting by the value of \code{sort_variable} 
#' within each stratum.
#' 
#' If there are finite population correction factors, then these finite population correction factors
#' will be applied to the replicate factors. This means that variance estimates with the finite population correction
#' do not require any adjustment to the overall scale factor used in variance estimation. This is the approach
#' used by the U.S. Census Bureau for the 5-year American Community Survey (ACS) replicate weights (U.S. Census Bureau, 2022, p. 12-8).
#' This approach is used regardless of whether the design has one overall finite population correction factor
#' or has different finite population correction factors for different strata.
#' 
#' @section Details on Row Assignments for Creating Replicate Factors:
#' The number of replicates must match the order of an available Hadamard matrix.
#' A Hadamard matrix can either be normal or non-normal: a normal Hadamard matrix
#' is one where the entries in the first row and in the first column are all equal to one.
#' If the user specifies \code{use_normal_hadamard = TRUE}, then there are more choices
#' of Hadamard matrix sizes available, and so greater flexibility in choosing the
#' number of replicates to create. When a normal Hadamard matrix is used, this will result
#' in the creation of an inactive replicate (sometimes referred to as a "dead" replicate),
#' which is a replicate where all the replicate factors equal one. Inactive replicates
#' are perfectly valid for variance estimation, though some users may find them
#' confusing.
#' 
#' An important part of the process of creating replicate weights is the assignment of rows of the Hadamard matrix
#' to primary sampling units. The method of Ash (2014) referred to as "RA1" is used for row assignments,
#' which means that the replication-based variance estimates for totals will
#' be equivalent to the SD2 variance estimator described by Ash (2014). The number of cycles
#' used with the "RA1" method is the smallest integer greater than \eqn{n/R}, where
#' \eqn{n} is the number of primary sample units and \eqn{R} is the number of replicates.
#' 
#' @references
#' Ash, S. (2014). "\emph{Using successive difference replication for estimating variances}."
#' \strong{Survey Methodology}, Statistics Canada, 40(1), 47-59.
#'
#' Fay, R.E. and Train, G.F. (1995). "\emph{Aspects of Survey and Model-Based Postcensal Estimation of
#' Income and Poverty Characteristics for States and Counties}." Joint Statistical Meetings,
#' Proceedings of the Section on Government Statistics, 154-159.
#' 
#' U.S. Census Bureau. (2022). "\emph{American Community Survey and Puerto Rico Community Survey Design and Methodology, Version 3.0.}" 
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
#' library(survey)
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
#'   data   = library_stsys_sample,
#'   strata = ~ SAMPLING_STRATUM,
#'   ids    = ~ 1,
#'   fpc    = ~ STRATUM_POP_SIZE
#' )
#'
#' ## Convert to SDR replicate design
#' sdr_design <- as_sdr_design(
#'   design              = design_obj,
#'   replicates          = 180,
#'   sort_variable       = "SAMPLING_SORT_ORDER",
#'   use_normal_hadamard = TRUE
#' )
#'
#' ## Compare to generalized bootstrap
#' ## based on the SD2 estimator that SDR approximates
#' gen_boot_design <- as_gen_boot_design(
#'   design             = design_obj,
#'   variance_estimator = "SD2",
#'   replicates         = 180,
#'   exact_vcov         = TRUE
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
  mse = TRUE
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
    mse = TRUE
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

  if (is.null(sort_variable)) {
    stop("Must specify a variable name for `sort_variable`.")
  }

  if (!is.null(sort_variable)) {

    if (!is.character(sort_variable) || length(sort_variable) != 1) {
      stop("`sort_variable` must be a single string.")
    }
    if (!sort_variable %in% colnames(compressed_design_structure$design_subset$variables)) {
      sprintf("The variable `%s` does  not appear in the data.", sort_variable) |> stop()
    }
    if (any(is.na(compressed_design_structure$design_subset$variables[[sort_variable]]))) {
      stop("`sort_variable` cannot have any missing values in the data.")
    }

    sort_data <- data.frame(
      ORIG_ORDER    = seq_len(nrow(compressed_design_structure$design)),
      PSU_STRATUM   = compressed_design_structure$design$strata[,1,drop=TRUE], 
      SORT_VARIABLE = compressed_design_structure$design$variables[[sort_variable]]
    )

    if (packageVersion("base") <= "4.4.0") {
      sort_by <- function(x, y, ...) {
          if (inherits(y, "formula")) 
              y <- .formula2varlist(y, x)
          if (!is.list(y)) 
              y <- list(y)
          o <- do.call(order, c(unname(y), list(...)))
          x[o, , drop = FALSE]
      }
    }

    sort_data <- sort_data |> sort_by(~ PSU_STRATUM + SORT_VARIABLE)
    sort_data[['SORT_ORDER']] <- seq_len(nrow(sort_data))
    sort_data <- sort_data |> sort_by(~ ORIG_ORDER)
    adjustment_factors <- adjustment_factors[sort_data[['SORT_ORDER']], , drop = FALSE]
  }

  # Apply first-stage FPCs
  if (!is.null(compressed_design_structure$design$fpc$popsize)) {
    fpc_factors <- sqrt(1 - (
      compressed_design_structure$design$fpc$sampsize[,1,drop = TRUE] /
      compressed_design_structure$design$fpc$popsize[,1,drop = TRUE] 
    ))
    stopifnot(nrow(adjustment_factors) == length(fpc_factors))
    adjustment_factors <- 1 + (fpc_factors * (adjustment_factors - 1))
    message("Finite population corrections are incorporated into the replicate factors.")
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

#' @export
as_sdr_design.DBIsvydesign <- function(
  design,
  replicates,
  sort_variable = NULL,
  use_normal_hadamard = FALSE,
  compress = TRUE,
  mse = TRUE
) {

rep_design <- NextMethod(design)

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
class(rep_design) <- c(
  "DBIrepdesign", "DBIsvydesign",
  setdiff(class(rep_design), c("DBIrepdesign", "DBIsvydesign"))
)

rep_design$call <- sys.call(which = -1)

return(rep_design)
}
