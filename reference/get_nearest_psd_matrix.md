# Positive Semidefinite Matrix Approximation

Approximates a symmetric, real matrix by the nearest positive
semidefinite matrix in the Frobenius norm, using the method of Higham
(1988). For a real, symmetric matrix, this is equivalent to "zeroing
out" negative eigenvalues. See the "Details" section for more
information.

## Usage

``` r
get_nearest_psd_matrix(X)
```

## Arguments

- X:

  A symmetric, real matrix with no missing values.

## Value

The nearest positive semidefinite matrix of the same dimension as `X`.

## Details

Let \\A\\ denote a symmetric, real matrix which is not positive
semidefinite. Then we can form the spectral decomposition \\A=\Gamma
\Lambda \Gamma^{\prime}\\, where \\\Lambda\\ is the diagonal matrix
whose entries are eigenvalues of \\A\\. The method of Higham (1988) is
to approximate \\A\\ with \\\tilde{A} = \Gamma \Lambda\_{+}
\Gamma^{\prime}\\, where the \\ii\\-th entry of \\\Lambda\_{+}\\ is
\\\max(\Lambda\_{ii}, 0)\\.

## References

\- Higham, N. J. (1988). "*Computing a nearest symmetric positive
semidefinite matrix.*" Linear Algebra and Its Applications, 103,
103-118.

## Examples

``` r
X <- matrix(
  c(2, 5, 5,
    5, 2, 5,
    5, 5, 2),
  nrow = 3, byrow = TRUE
)
get_nearest_psd_matrix(X)
#>      [,1] [,2] [,3]
#> [1,]    4    4    4
#> [2,]    4    4    4
#> [3,]    4    4    4
```
