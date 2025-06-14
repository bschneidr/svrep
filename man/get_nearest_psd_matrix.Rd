% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quadratic_forms.R
\name{get_nearest_psd_matrix}
\alias{get_nearest_psd_matrix}
\title{Positive Semidefinite Matrix Approximation}
\usage{
get_nearest_psd_matrix(X)
}
\arguments{
\item{X}{A symmetric, real matrix with no missing values.}
}
\value{
The nearest positive semidefinite matrix
of the same dimension as \code{X}.
}
\description{
Approximates a symmetric, real matrix by the nearest positive
semidefinite matrix in the Frobenius norm, using the method of Higham (1988).
For a real, symmetric matrix, this is equivalent to "zeroing out" negative eigenvalues.
See the "Details" section for more information.
}
\details{
Let \eqn{A} denote a symmetric, real matrix which is not positive semidefinite.
Then we can form the spectral decomposition \eqn{A=\Gamma \Lambda \Gamma^{\prime}},
where \eqn{\Lambda} is the diagonal matrix
whose entries are eigenvalues of \eqn{A}.
The method of Higham (1988) is to  approximate
\eqn{A} with \eqn{\tilde{A} = \Gamma \Lambda_{+} \Gamma^{\prime}},
where the \eqn{ii}-th entry of \eqn{\Lambda_{+}} is \eqn{\max(\Lambda_{ii}, 0)}.
}
\examples{
X <- matrix(
  c(2, 5, 5,
    5, 2, 5,
    5, 5, 2),
  nrow = 3, byrow = TRUE
)
get_nearest_psd_matrix(X)
}
\references{
- Higham, N. J. (1988). "\emph{Computing a nearest symmetric positive semidefinite matrix.}" Linear Algebra and Its Applications, 103, 103-118.
}
