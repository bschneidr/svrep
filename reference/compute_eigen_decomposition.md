# Compute the eigendecomposition of a matrix

Compute the eigendecomposition of a matrix

## Usage

``` r
compute_eigen_decomposition(Sigma)
```

## Arguments

- Sigma:

  A symmetric matrix

## Value

A list object with two named components:

- `values`: A vector of eigenvalues.

- `vectors`: A matrix with eigenvectors corresponding to the
  eigenvalues.

## Details

Computes the eigendecomposition of a symmetric matrix, either using the
base R function [`eigen()`](https://rdrr.io/r/base/eigen.html) or
instead using the 'torch' package if the user has 'torch' installed and
has set `options(svrep.torch_device = "cpu")` or
`options(svrep.torch_device = "cuda")`.
