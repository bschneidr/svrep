# Add Inactive Replicates to a Survey Design Object

Adds inactive replicates to a survey design object. An inactive
replicate is a replicate that does not contribute to variance estimates
but adds to the matrix of replicate weights so that the matrix has the
desired number of columns. The new replicates' values are simply equal
to the full-sample weights.

## Usage

``` r
add_inactive_replicates(design, n_total, n_to_add, location = "last")
```

## Arguments

- design:

  A survey design object, created with either the `survey` or `srvyr`
  packages.

- n_total:

  The total number of replicates that the result should contain. If the
  design already contains `n_total` replicates (or more), then no update
  is made.

- n_to_add:

  The number of additional replicates to add. Can only use the `n_total`
  argument OR the `n_to_add` argument, not both.

- location:

  Either `"first"`, `"last"` (the default), or `"random"`. Specifies
  where the columns of new replicates should be located in the matrix of
  replicate weights. Use `"first"` to place new replicates first (i.e.,
  in the leftmost part of the matrix), `"last"` to place the new
  replicates last (i.e., in the rightmost part of the matrix). Use
  `"random"` to intersperse the new replicates in random column
  locations of the matrix; the original replicates will still be in
  their original order.

## Value

An updated survey design object, where the number of columns of
replicate weights has potentially increased. The increase only happens
if the user specifies the `n_to_add` argument instead of `n_total`, of
if the user specifies `n_total` and `n_total` is less than the number of
columns of replicate weights that the design already had.

## Statistical Details

Inactive replicates are also sometimes referred to as "dead replicates",
for example in Ash (2014). The purpose of adding inactive replicates is
to increase the number of columns of replicate weights without impacting
variance estimates. This can be useful, for example, when combining data
from a survey across multiple years, where different years use different
number of replicates, but a consistent number of replicates is desired
in the combined data file.

Suppose the initial replicate design has \\L\\ replicates, with
respective constants \\c_k\\ for \\k=1,\dots,L\\ used to estimate
variance with the formula \$\$v\_{R} = \sum\_{k=1}^L
c_k\left(\hat{T}\_y^{(k)}-\hat{T}\_y\right)^2\$\$ where \\\hat{T}\_y\\
is the estimate produced using the full-sample weights and
\\\hat{T}\_y^{(k)}\\ is the estimate from replicate \\k\\.

Inactive replicates are simply replicates that are exactly equal to the
full sample: that is, the replicate \\k\\ is called "inactive" if its
vector of replicate weights exactly equals the full-sample weights. In
this case, when using the formula above to estimate variances, these
replicates contribute nothing to the variance estimate.

If the analyst uses the variant of the formula above where the
full-sample estimate \\\hat{T}\_y\\ is replaced by the average replicate
estimate (i.e., \\L^{-1}\sum\_{k=1}^{L}\hat{T}\_y^{(k)}\\), then
variance estimates will differ before vs. after adding the inactive
replicates. For this reason, it is strongly recommend to explicitly
specify `mse=TRUE` when creating a replicate design object in R with
functions such as
[`svrepdesign()`](https://rdrr.io/pkg/survey/man/svrepdesign.html),
[`as_bootstrap_design()`](https://bschneidr.github.io/svrep/reference/as_bootstrap_design.md),
etc. If working with an already existing replicate design, you can
update the `mse` option to `TRUE` simply by using code such as
`my_design$mse <- TRUE`.

## References

Ash, S. (2014). "*Using successive difference replication for estimating
variances*." **Survey Methodology**, Statistics Canada, 40(1), 47-59.

## Examples

``` r
set.seed(2023)

# Create an example survey design object

  sample_data <- data.frame(
    PSU     = c(1,2,3)
  )

  survey_design <- svydesign(
    data = sample_data,
    ids = ~ PSU,
    weights = ~ 1
  )

  rep_design <- survey_design |>
    as.svrepdesign(type = "JK1", mse = TRUE)

# Inspect replicate weights

  rep_design |> weights(type = "analysis")
#>      [,1] [,2] [,3]
#> [1,]  0.0  1.5  1.5
#> [2,]  1.5  0.0  1.5
#> [3,]  1.5  1.5  0.0

# Inspect replicates after adding inactive replicates

  rep_design |>
    add_inactive_replicates(n_total = 5, location = "first") |>
    weights(type = "analysis")
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    1  0.0  1.5  1.5
#> [2,]    1    1  1.5  0.0  1.5
#> [3,]    1    1  1.5  1.5  0.0

  rep_design |>
    add_inactive_replicates(n_to_add = 2, location = "last") |>
    weights(type = "analysis")
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]  0.0  1.5  1.5    1    1
#> [2,]  1.5  0.0  1.5    1    1
#> [3,]  1.5  1.5  0.0    1    1

  rep_design |>
    add_inactive_replicates(n_to_add = 5, location = "random") |>
    weights(type = "analysis")
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,]    1    1    1  0.0    1  1.5    1  1.5
#> [2,]    1    1    1  1.5    1  0.0    1  1.5
#> [3,]    1    1    1  1.5    1  1.5    1  0.0
```
