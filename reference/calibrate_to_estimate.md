# Sample-based Calibration to An Estimate

Calibrate the weights of a primary survey to match estimated totals from
a control survey, using adjustments to the replicate weights to account
for the variance of the estimated control totals. The adjustments to
replicate weights are conducted using the method proposed by Fuller
(1998). This method can be used to implement general calibration as well
as post-stratification or raking specifically (see the details for the
`calfun` parameter).

## Usage

``` r
calibrate_to_estimate(
  rep_design,
  estimate,
  vcov_estimate,
  cal_formula,
  calfun = survey::cal.linear,
  bounds = list(lower = -Inf, upper = Inf),
  verbose = FALSE,
  maxit = 50,
  epsilon = 1e-07,
  variance = NULL,
  col_selection = NULL
)
```

## Arguments

- rep_design:

  A replicate design object for the primary survey, created with either
  the `survey` or `srvyr` packages.

- estimate:

  A vector of estimated control totals. The names of entries must match
  the names from calling
  `svytotal(x = cal_formula, design = rep_design)`.

- vcov_estimate:

  A variance-covariance matrix for the estimated control totals. The
  column names and row names must match the names of `estimate`.

- cal_formula:

  A formula listing the variables to use for calibration. All of these
  variables must be included in `rep_design`.

- calfun:

  A calibration function from the `survey` package, such as
  [cal.linear](https://rdrr.io/pkg/survey/man/make.calfun.html),
  [cal.raking](https://rdrr.io/pkg/survey/man/make.calfun.html), or
  [cal.logit](https://rdrr.io/pkg/survey/man/make.calfun.html). Use
  `cal.linear` for ordinary post-stratification, and `cal.raking` for
  raking. See [calibrate](https://rdrr.io/pkg/survey/man/calibrate.html)
  for additional details.

- bounds:

  Parameter passed to
  [grake](https://rdrr.io/pkg/survey/man/calibrate.html) for
  calibration. See
  [calibrate](https://rdrr.io/pkg/survey/man/calibrate.html) for
  details.

- verbose:

  Parameter passed to
  [grake](https://rdrr.io/pkg/survey/man/calibrate.html) for
  calibration. See
  [calibrate](https://rdrr.io/pkg/survey/man/calibrate.html) for
  details.

- maxit:

  Parameter passed to
  [grake](https://rdrr.io/pkg/survey/man/calibrate.html) for
  calibration. See
  [calibrate](https://rdrr.io/pkg/survey/man/calibrate.html) for
  details.

- epsilon:

  Parameter passed to
  [grake](https://rdrr.io/pkg/survey/man/calibrate.html) for
  calibration.  
  After calibration, the absolute difference between each calibration
  target and the calibrated estimate will be no larger than `epsilon`
  times (1 plus the absolute value of the target). See
  [calibrate](https://rdrr.io/pkg/survey/man/calibrate.html) for
  details.

- variance:

  Parameter passed to
  [grake](https://rdrr.io/pkg/survey/man/calibrate.html) for
  calibration. See
  [calibrate](https://rdrr.io/pkg/survey/man/calibrate.html) for
  details.

- col_selection:

  Optional parameter to determine which replicate columns will have
  their control totals perturbed. If supplied, `col_selection` must be
  an integer vector with length equal to the length of `estimate`.

## Value

A replicate design object, with full-sample weights calibrated to totals
from `estimate`, and replicate weights adjusted to account for variance
of the control totals. The element `col_selection` indicates, for each
replicate column of the calibrated primary survey, which column of
replicate weights it was matched to from the control survey.

## Details

With the Fuller method, each of `k` randomly-selected replicate columns
from the primary survey are calibrated to control totals formed by
perturbing the `k`-dimensional vector of estimated control totals using
a spectral decomposition of the variance-covariance matrix of the
estimated control totals. Other replicate columns are simply calibrated
to the unperturbed control totals.  

Because the set of replicate columns whose control totals are perturbed
should be random, there are multiple ways to ensure that this matching
is reproducible. The user can either call
[set.seed](https://rdrr.io/r/base/Random.html) before using the
function, or supply a vector of randomly-selected column indices to the
argument `col_selection`.

## Syntax for Common Types of Calibration

For ratio estimation with an auxiliary variable `X`, use the following
options:  
- `cal_formula = ~ -1 + X`  
- `variance = 1`,  
- `cal.fun = survey::cal.linear`

For post-stratification, use the following option:

\- `cal.fun = survey::cal.linear`

For raking, use the following option:

\- `cal.fun = survey::cal.raking`

## References

Fuller, W.A. (1998). "Replication variance estimation for two-phase
samples." **Statistica Sinica**, *8*: 1153-1164.

Opsomer, J.D. and A. Erciulescu (2021). "Replication variance estimation
after sample-based calibration." **Survey Methodology**, *47*: 265-277.

## See also

\[calibrate_to_sample()\] as an alternative when the control survey has
replicate weights.

## Examples

``` r
# \donttest{

# Load example data for primary survey ----

  data(api)

  primary_survey <- svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc) |>
    as.svrepdesign(type = "JK1")

# Load example data for control survey ----

  control_survey <- svydesign(id = ~ 1, fpc = ~fpc, data = apisrs) |>
    as.svrepdesign(type = "JK1")

# Estimate control totals ----

  estimated_controls <- svytotal(x = ~ stype + enroll,
                                 design = control_survey)
  control_point_estimates <- coef(estimated_controls)
  control_vcov_estimate <- vcov(estimated_controls)

# Calibrate totals for one categorical variable and one numeric ----

  calibrated_rep_design <- calibrate_to_estimate(
    rep_design = primary_survey,
    estimate = control_point_estimates,
    vcov_estimate = control_vcov_estimate,
    cal_formula = ~ stype + enroll
  )
#> Selection of replicate columns whose control totals will be perturbed will be done at random.
#> For tips on reproducible selection, see `help('calibrate_to_estimate')`
#> Warning: Setting `mse` to TRUE; variance estimates will be centered around full-sample estimate, not mean of replicates.

# Inspect estimates before and after calibration ----

  ##_ For the calibration variables, estimates and standard errors
  ##_ from calibrated design will match those of the control survey

    svytotal(x = ~ stype + enroll, design = primary_survey)
#>             total        SE
#> stypeE    4873.97   1333.32
#> stypeH     473.86    158.70
#> stypeM     846.17    167.55
#> enroll 3404940.13 932235.03
    svytotal(x = ~ stype + enroll, design = control_survey)
#>             total        SE
#> stypeE    4397.74    196.00
#> stypeH     774.25    142.85
#> stypeM    1022.01    160.33
#> enroll 3621074.34 169519.65
    svytotal(x = ~ stype + enroll, design = calibrated_rep_design)
#>             total        SE
#> stypeE    4397.74    196.00
#> stypeH     774.25    142.85
#> stypeM    1022.01    160.33
#> enroll 3621074.34 169519.65

  ##_ Estimates from other variables will be changed as well

    svymean(x = ~ api00 + api99, design = primary_survey)
#>         mean     SE
#> api00 644.17 26.329
#> api99 606.98 26.998
    svymean(x = ~ api00 + api99, design = control_survey)
#>         mean     SE
#> api00 656.58 9.2497
#> api99 624.68 9.5003
    svymean(x = ~ api00 + api99, design = calibrated_rep_design)
#>         mean     SE
#> api00 642.69 27.336
#> api99 606.91 28.236

# Inspect weights before and after calibration ----

  summarize_rep_weights(primary_survey, type = 'overall')
#>   nrows ncols degf_svy_pkg rank avg_wgt_sum sd_wgt_sums min_rep_wgt max_rep_wgt
#> 1   183    15           14   15        6194    403.1741           0    36.26464
  summarize_rep_weights(calibrated_rep_design, type = 'overall')
#>   nrows ncols degf_svy_pkg rank avg_wgt_sum sd_wgt_sums min_rep_wgt max_rep_wgt
#> 1   183    15           14   15        6194 4.84015e-09           0    119.1069

# For reproducibility, specify which columns are randomly selected for Fuller method ----

  column_selection <- calibrated_rep_design$col_selection
  print(column_selection)
#> [1]  4 11  7  8

  calibrated_rep_design <- calibrate_to_estimate(
    rep_design = primary_survey,
    estimate = control_point_estimates,
    vcov_estimate = control_vcov_estimate,
    cal_formula = ~ stype + enroll,
    col_selection = column_selection
  )
#> Warning: Setting `mse` to TRUE; variance estimates will be centered around full-sample estimate, not mean of replicates.
# }
```
