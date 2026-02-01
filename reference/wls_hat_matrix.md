# Create the "hat matrix" for weighted least squares regression

Create the "hat matrix" for weighted least squares regression

## Usage

``` r
wls_hat_matrix(X, w)
```

## Arguments

- X:

  Matrix of predictor variables, with `n` rows

- w:

  Vector of weights (should all be nonnegative), of length `n`

## Value

An \\n \times n\\ matrix. This is the "hat matrix" for a WLS regression.
