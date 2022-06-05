#' Create a row-assignment matrix for Successive-Difference-Replication-Method
#'
#' @param n The sample size of the data.
#' @param hadamard_order The order of the Hadamard matrix
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
assign_hadamard_rows <- function(n, hadamard_order, use_first_row = TRUE) {

  if ((hadamard_order %% 4) > 0) {
    stop("`hadamard_order` must be a multiple of 4.")
  }

  k_A <- hadamard_order

  if (!use_first_row) {
    k_A <- k_A - 1
  }

  # Determine number of cycles
  n_cycles <- ceiling(n / k_A)

  # Determine range of step sizes to use to form row-assignment matrix.
  # Cannot be larger than (order of Hadamard matrix minus one).
  if (((k_A - 1) * n_cycles) %% 2) {
    max_step_size <- n_cycles + 1L
    if (max_step_size >= (k_A - 1)) {
      max_step_size <- n_cycles - 1L
    }
  } else {
    max_step_size <- pmin(n_cycles, k_A - 1)
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
