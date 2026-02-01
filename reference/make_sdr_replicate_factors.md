# Factors for the Successive Difference Replication Method

Create matrix of replicate factors to use for Successive Difference
Replication Method

## Usage

``` r
make_sdr_replicate_factors(
  n,
  target_number_of_replicates,
  use_normal_hadamard = FALSE
)
```

## Arguments

- n:

  The number of sampling units

- target_number_of_replicates:

  The target number of replicates to create. This will determine the
  order of the Hadamard matrix to use when creating replicate factors.
  The actual number of replicates will be a multiple of 4. If
  `use_normal_hadamard = FALSE`, then the actual number of replicates
  will be \\4 \times 2^k\\ for some integer \\k\\, which means that the
  actual number of replicates might be much larger than the target.

- use_normal_hadamard:

  Whether to use a normal Hadamard matrix: that is, a matrix whose first
  row and first column only have entries equal to 1.

## Value

A matrix of replicate factors, with `n` rows and the number of columns
corresponding to the order of the Hadamard matrix used (which is greater
than or equal to `target_number_of_replicates`).

## Examples

``` r
# Note that one of the replicates has every factor equal to 1
# Also note that this matches Table 1 in Ash (2014)
make_sdr_replicate_factors(
  n = 4,
  target_number_of_replicates = 4,
  use_normal_hadamard = TRUE
)
#> Using Hadamard matrix of order 4. If `use_normal_hadamard=TRUE`, the smallest possible order is 4.
#>      [,1]      [,2]      [,3]      [,4]
#> [1,]    1 1.7071068 1.0000000 1.7071068
#> [2,]    1 0.2928932 1.7071068 1.0000000
#> [3,]    1 1.7071068 1.0000000 0.2928932
#> [4,]    1 0.2928932 0.2928932 1.0000000

# Note the difference when using a non-normal Hadamard matrix
rep_factors <- make_sdr_replicate_factors(
  n = 4,
  target_number_of_replicates = 4,
  use_normal_hadamard = FALSE
)
#> Using Hadamard matrix of order 4. If `use_normal_hadamard=TRUE`, the smallest possible order is 4.
print(rep_factors)
#>           [,1]      [,2]      [,3]      [,4]
#> [1,] 1.7071068 1.0000000 1.7071068 1.0000000
#> [2,] 0.2928932 1.0000000 1.0000000 1.7071068
#> [3,] 1.0000000 0.2928932 1.0000000 0.2928932
#> [4,] 1.0000000 1.7071068 0.2928932 1.0000000

# These replicate factors are equivalent
# to the SD2 variance estimator
tcrossprod(rep_factors - 1)
#>      [,1] [,2] [,3] [,4]
#> [1,]  1.0 -0.5  0.0 -0.5
#> [2,] -0.5  1.0 -0.5  0.0
#> [3,]  0.0 -0.5  1.0 -0.5
#> [4,] -0.5  0.0 -0.5  1.0

# Compare to the quadratic form of the SD2 estimator
sd2_quad_form <- make_quad_form_matrix(
  variance_estimator = "SD2",
  cluster_ids        = matrix(1:4, ncol = 1),
  sort_order         = matrix(1:4, ncol = 1)
)
print(sd2_quad_form)
#> 4 x 4 sparse Matrix of class "dsCMatrix"
#>                         
#> [1,]  1.0 -0.5  .   -0.5
#> [2,] -0.5  1.0 -0.5  .  
#> [3,]  .   -0.5  1.0 -0.5
#> [4,] -0.5  .   -0.5  1.0
```
