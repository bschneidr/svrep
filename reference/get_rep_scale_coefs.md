# Access Replication Scale Coefficients

Get the scale coefficients used for variance estimation in a replicate
design object.

## Usage

``` r
get_rep_scale_coefs(rep_design, type = "combined")
```

## Arguments

- rep_design:

  A replicate design object

- type:

  Either `'overall'`, `'specific'`, or `'combined'`. See the details
  section below. The result for `'overall'` is the single overall scale
  coefficient. The result for `'specific'` is the vector of
  replicate-specific coefficients. The result for `'combined'` is the
  product of the overall and replicate-specific coefficients.

## Value

If `type = 'overall'`, the result is a single number. Otherwise, the
result is a vector with length matching the number of replicates.

## Details

For a statistic \\\hat{\theta}\\, replication methods estimate the
sampling variance using \\R\\ replicate estimates, with the estimate for
the \\r\\-th replicate denoted \\\hat{\theta}\_r\\.

The formula for the variance estimate is the following: \$\$
v(\hat{\theta}) = C \sum\_{r=1}^{R} c_r (\hat{\theta_r} -
\hat{\theta})^2 \$\$

The terms \\C\\ and \\c_r, r=1,\dots,R\\ are scale coefficients. \\C\\
is the overall coefficient, and \\c_r, r=1,\dots,R\\ are
replicate-specific coefficients.

Specifying `get_rep_scale_coefs(type='overall')` returns the overall
coefficient \\C\\. Specifying `type='specific'` returns the
replicate-specific coefficients \\c_r, r=1,\dots,R\\.

Specifying `type='combined'` returns a vector with \\R\\ elements, where
the \\r\\-th element is \\C \times c_r\\.

## Examples

``` r
data('api', package = 'survey')

api_design <- svydesign(
  data    = apistrat, 
  id      = ~ 1, 
  strata  = ~ stype,
  weights = ~ pw,
  nest    = TRUE
)

jk_design <- api_design |>
  as_random_group_jackknife_design(
    replicates = 12
  )

jk_design |>
  get_rep_scale_coefs('overall')
#> [1] 1

jk_design |>
  get_rep_scale_coefs('specific')
#>  [1] 0.915 0.915 0.915 0.915 0.915 0.915 0.915 0.915 0.920 0.920 0.920 0.920

jk_design |>
  get_rep_scale_coefs('combined')
#>  [1] 0.915 0.915 0.915 0.915 0.915 0.915 0.915 0.915 0.920 0.920 0.920 0.920
```
