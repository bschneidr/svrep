# Compute the matrix of joint inclusion probabilities from the quadratic form of a Horvitz-Thompson variance estimator.

Compute the matrix of joint inclusion probabilities from the quadratic
form of a Horvitz-Thompson variance estimator.

## Usage

``` r
ht_matrix_to_joint_probs(ht_quad_form)
```

## Arguments

- ht_quad_form:

  The matrix of the quadratic form representing the Horvitz-Thompson
  variance estimator.

## Value

The matrix of joint inclusion probabilities

## Details

The quadratic form matrix of the Horvitz-Thompson variance estimator has
\\ij\\-th entry equal to \\(1-\frac{\pi_i \pi_j}{\pi\_{ij}})\\. The
matrix of joint probabilities has \\ij\\-th entry equal to
\\\pi\_{ij}\\.
