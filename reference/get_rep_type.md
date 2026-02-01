# Access Type of Replication Method

Identify the type of replication method used for the replicates in a
replicate design object.

## Usage

``` r
get_rep_type(rep_design)
```

## Arguments

- rep_design:

  A replicate design object

## Value

A character string giving the type of replication method.

## Examples

``` r
data('scd', package = 'survey')

scd_design <- svydesign(
  data   = scd, 
  id     = ~ ambulance, 
  prob   = ~ 1,
  strata = ~ ESA,
  nest   = TRUE
)

scd_design |>
  as_bootstrap_design(replicates = 5) |>
  get_rep_type()
#> [1] "bootstrap"

scd_design |>
  as_fays_gen_rep_design(
    variance_estimator = "Ultimate Cluster"
  ) |>
  get_rep_type()
#> [1] "other"
```
