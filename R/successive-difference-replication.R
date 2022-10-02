#' Create a row-assignment matrix for Successive Difference Replication-Method
#'
#' @param n The sample size of the data.
#' @param hadamard_order The order of the Hadamard matrix
#' @param number_of_cycles The number of cycles to use in the row assignment.
#' Must be at least as large as \code{n/hadamard_order}.
#' @param use_first_row Whether to use the first row of the Hadamard matrix.
#' The first row of a Hadamard matrix is often all 1's, and so using the first
#' row to create replicate factors leads to the creation of a replicate
#' whose weights exactly match the full-sample weights. Thus, using the first row
#' of the Hadamard matrix may be undesirable for practical purposes,
#' even if it is valid for the purpose of variance estimation.
#' @return A matrix with \code{n} rows and two columns. Each row
#' gives the assignment of two rows of a Hadamard matrix to the row of data.
#' @export
#'
#' @examples
assign_hadamard_rows <- function(n, hadamard_order, number_of_cycles = ceiling(n/hadamard_order), use_first_row = TRUE) {

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

  return(row_assignment_matrix)
}

#' Create matrix of replicate factors to use for Successive Difference Replication Method
#'
#' @param n The sample size of the data.
#' @param target_number_of_replicates The target number of replicates to create.
#' This will determine the order of the Hadamard matrix to use when
#' creating replicate factors. The actual number of replicates will
#' be a power of four.
#'
#' @return
#' @export
#'
#' @examples
create_sdr_replicate_factors <- function(n, target_number_of_replicates, use_normal_hadamard = FALSE) {

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
    H_A <- H_A <- 2*H_A - 1 # Convert from 1/0 format to 1/-1 format
  }
  hadamard_order <- ncol(H_A)
  sprintf("Using Hadamard matrix of order %s", hadamard_order) |>
    message()

  # Assign rows of the Hadamard matrix to each observation

  row_assignment_matrix <- assign_hadamard_rows(
    n = n,
    hadamard_order = hadamard_order,
    use_first_row = TRUE
  )

  # Create replicate factors based on the Hadamard matrix and the row assignments
  replicate_factors <- sapply(X = 1L:hadamard_order, FUN = function(r) {
    hadamard_entries <- cbind('h_1r' = H_A[row_assignment_matrix[,'row_1'], r],
                              'h_2r' = H_A[row_assignment_matrix[,'row_2'], r])

    1 + (hadamard_entries %*% c(1,-1) * (2^(-3/2)))
  }, simplify = TRUE)

  return(replicate_factors)
}

sd_circular_estimator <- function(y, N) {

  # Finite population correction
    n <- length(y)
    f <- n/N

  # Contrast matrix described by Ash (2011)

    C_matrix <- diag(n) * 2
    for (i in seq_len(n)) {
      if (i < n) {
        C_matrix[i,i+1] <- -1
        C_matrix[i+1,i] <- -1
      }
      if (i > 0) {
        C_matrix[i,i-1] <- -1
        C_matrix[i-1,i] <- -1
      }
    }
    C_matrix[1,n] <- -1
    C_matrix[n,1] <- -1

  # Estimate sampling variance
    wtd_y <- y * (N/n)

    variance_estimate <- 0.5 * (1-f) * (t(wtd_y) %*% C_matrix %*% wtd_y)
    return(variance_estimate)
}

sd_noncircular_estimator <- function(y, N) {

  # Finite population correction
  n <- length(y)
  f <- n/N

  # Contrast matrix described by Ash (2011)

  C_matrix <- diag(n) * 2
  C_matrix[1,1] <- 1
  C_matrix[n,n] <- 1
  for (i in seq_len(n)) {
    if (i < n) {
      C_matrix[i,i+1] <- -1
      C_matrix[i+1,i] <- -1
    }
    if (i > 0) {
      C_matrix[i,i-1] <- -1
      C_matrix[i-1,i] <- -1
    }
  }

  # Estimate sampling variance
  wtd_y <- y * (N/n)

  variance_estimate <- 0.5 * (1-f) * (n/(n-1)) * (t(wtd_y) %*% C_matrix %*% wtd_y)
  return(variance_estimate)
}
