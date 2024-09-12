#' @title Package-level options for svrep
#' @description This help page describes the overall
#' options that can be set for your R session, 
#' using the function `options()`.
#' @rdname svrep-package-options
#' @name svrep-package-options
#' @md
#' @section Options for using the 'torch' package to speed up certain operations: 
#' The 'torch' package provides access to fast linear algebra routines
#' and is a particularly convenient approach to working with GPUs or
#' conducting multithreaded linear algebra operations.
#' 
#' Setting the following options will allow functions in 'svrep' to
#' use the 'torch' package to speed up certain computationally intensive
#' operations that occur when creating replicate weights, particularly
#' for Fay's generalized replication method or generalized bootstrap methods.
#' 
#' The option `svrep.torch_device` accepts the following options:
#' 
#' - `options(svrep.torch_device = 'none')`: The 'torch' package will not be used.
#' 
#' - `options(svrep.torch_device = 'cpu')`: The 'torch' package will be used
#'   with all operations done on the CPU.
#'   
#' - `options(svrep.torch_device = 'cuda')`: The 'torch' package will be used,
#'   with operations conducted on the GPU if possible. This requires the user's
#'   computer to have a CUDA-enabled GPU.
#' 
#' Note that precise values for matrix decompositions can vary between
#' different linear algebra libraries (including among common BLAS/LAPACK),
#' and so the replicate weights created with 'torch' may not exactly
#' match those created without 'torch'; differences will generally be small.
#' 
#' The following function from 'torch' will control
#' the number of threads used for computations, which can have a large 
#' impact on speed.
#' 
#' - `torch::set_num_threads()`: Sets the number of threads that 'torch' can use.
#' 
#' @section Relevant options from the 'survey' package: 
#' 
#' The 'survey' package has the following options which are of particular relevance
#' to users of 'svrep'.
#' 
#' - `options(survey.replicates.mse = TRUE/FALSE)`: The default value for this option
#'   is `FALSE`. This option controls the default value used for the `mse` argument in
#'   the functions `svrepdesign()` and `as.svrepdesign()`. 
#'  
#' Call `help('survey.replicates.mse', package = 'survey')` for more details.
#' 
#' In nearly all cases, it is safer to use 
#' `options(survey.replicates.mse = TRUE)`, or--better yet--
#' to always specify `svrepdesign(..., mse = TRUE)` or 
#' `as.svrepdesign(..., mse = TRUE)` when using functions with an `mse` argument.
#' 
#' For replicate weights created using Fay's generalized replication method or 
#' the generalized bootstrap, using `mse = FALSE` can result in badly biased
#' variance estimates.
#' 
#' - `options(survey.multicore = TRUE/FALSE)`: The default value for this option
#'   is `FALSE`. Setting this option to `TRUE` means that multiple processors
#'   will be used for certain variance calculations involving replicate weights,
#'   such as in `svyglm()`. 
#'   
#' This can potentially speed up calculations but is not guaranteed to do so. 
#'   
#' Call `help('survey.multicore', package = 'survey')` for more details.
#' 
NULL