#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom DBI dbGetQuery
#' @importFrom DBI dbIsValid
#' @importFrom Matrix diag
#' @importFrom Matrix isSymmetric
#' @importFrom Matrix rowSums
#' @importFrom Matrix t
#' @importFrom methods as
#' @importFrom methods is
#' @importFrom sampling inclusionprobabilities
#' @importFrom sampling UPmaxentropy
#' @importFrom sampling UPpoisson
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats na.pass
#' @importFrom stats reformulate
#' @importFrom stats setNames
#' @importFrom stats terms
#' @importFrom stats weights
#' @importFrom utils packageVersion
## usethis namespace: end
NULL

.onLoad <- function(...) {
  if (is.null(getOption("svrep.torch_device"))) {
    options(svrep.torch_device = 'none')
  }
}
