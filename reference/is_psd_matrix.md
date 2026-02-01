# Check if Matrix is Positive Semidefinite

Check whether a matrix is positive semidefinite, based on checking for
symmetric and negative eigenvalues.

## Usage

``` r
is_psd_matrix(X, tolerance = sqrt(.Machine$double.eps))
```

## Arguments

- X:

  A matrix with no missing or infinite values.

- tolerance:

  Tolerance for controlling whether a tiny computed eigenvalue will
  actually be considered negative. Computed negative eigenvalues will be
  considered negative if they are less than which are less than
  `-abs(tolerance * max(eigen(X)$values))`. A small nonzero tolerance is
  recommended since eigenvalues are nearly always computed with some
  floating-point error.

## Value

A logical value. `TRUE` if the matrix is deemed positive semidefinite.
Negative otherwise (including if `X` is not symmetric).

## See also

The function
[`get_nearest_psd_matrix()`](https://bschneidr.github.io/svrep/reference/get_nearest_psd_matrix.md)
can be used to approximate a symmetric matrix which is not positive
semidefinite, by a similar positive semidefinite matrix.

## Examples

``` r
X <- matrix(
  c(2, 5, 5,
    5, 2, 5,
    5, 5, 2),
  nrow = 3, byrow = TRUE
)

is_psd_matrix(X)
#> [1] FALSE

eigen(X)$values
#> [1] 12 -3 -3
```
