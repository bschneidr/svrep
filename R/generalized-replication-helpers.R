#' @title Compute the eigendecomposition of a matrix
#' @details Computes the eigendecomposition of a symmetric matrix,
#' either using the base R function `eigen()` or instead using the 'torch' 
#' package if the user has 'torch' installed and
#' has set `options(svrep.torch_device = "cpu")` 
#' or `options(svrep.torch_device = "cuda")`.
#' @param Sigma A symmetric matrix
#'
#' @return A list object with two named components:
#' - `values`: A vector of eigenvalues.
#' - `vectors`: A matrix with eigenvectors corresponding to the eigenvalues.
#' @keywords internal
#' @md
compute_eigen_decomposition <- function(Sigma) {
  
  # Determine whether, and how, the torch package should be used
  use_torch <- FALSE
  
  svrep_torch_device <- getOption("svrep.torch_device")
  
  if (!is.null(svrep_torch_device) && svrep_torch_device %in% c("cpu", "cuda")) {
    if (requireNamespace("torch", quietly = TRUE)) {
      use_torch <- torch::torch_is_installed()
      if (svrep_torch_device == "cuda" && !torch::cuda_is_available()) {
        svrep_torch_device <- "cpu"
        message("Using CPU for torch, since `torch::cuda_is_available()` returned FALSE.")
      }
    }
  }
  
  
  # Compute the eigen decomposition with or without torch,
  # ensure that eigenvalues/vectors are returned in descending order
  if (!use_torch) {
    eigen_decomposition <- eigen(x = Sigma, symmetric = TRUE)
  }
  
  if (use_torch) {
    
    # Note that with 'torch', eigenvalues/vectors are in *ascending* order
    reversal_indices <- seq(from = ncol(Sigma), to = 1, by = -1)
    
    eigen_decomposition <- Sigma |> as.matrix() |> 
      torch::torch_tensor(device = torch::torch_device(
        type = svrep_torch_device
      )) |>
      torch::linalg_eigh() |>
      (\(torch_eigh_output) {
        list('values' = as.numeric(torch_eigh_output[[1]])[reversal_indices],
             'vectors' = as.matrix(torch_eigh_output[[2]])[,reversal_indices,drop=FALSE])
      })()
  }
  
  return(eigen_decomposition)
}