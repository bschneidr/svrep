# Sample-based Calibration with Replicates

Calibrate the weights of a primary survey to match estimated totals from
a control survey, using adjustments to the replicate weights to account
for the variance of the estimated control totals. Both surveys must have
replicate weights. The adjustments to replicate weights are conducted
using the method proposed by Opsomer and Erciulescu (2021). This method
can be used to implement general calibration as well as
post-stratification or raking specifically (see the details for the
`calfun` parameter).

## Usage

``` r
calibrate_to_sample(
  primary_rep_design,
  control_rep_design,
  cal_formula,
  calfun = survey::cal.linear,
  bounds = list(lower = -Inf, upper = Inf),
  verbose = FALSE,
  maxit = 50,
  epsilon = 1e-07,
  variance = NULL,
  control_col_matches = NULL
)
```

## Arguments

- primary_rep_design:

  A replicate design object for the primary survey, created with either
  the `survey` or `srvyr` packages.

- control_rep_design:

  A replicate design object for the control survey.

- cal_formula:

  A formula listing the variables to use for calibration. All of these
  variables must be included in both `primary_rep_design` and
  `control_rep_design`.

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

- control_col_matches:

  Optional parameter to specify which control survey replicate is
  matched to each primary survey replicate. If the \\i-th\\ entry of
  `control_col_matches` equals \\k\\, then replicate \\i\\ in
  `primary_rep_design` is matched to replicate \\k\\ in
  `control_rep_design.` Entries of `NA` denote a primary survey
  replicate not matched to any control survey replicate. If this
  parameter is not used, matching is done at random.

## Value

A replicate design object, with full-sample weights calibrated to totals
from `control_rep_design`, and replicate weights adjusted to account for
variance of the control totals. If `primary_rep_design` had fewer
columns of replicate weights than `control_rep_design`, then the number
of replicate columns and the length of `rscales` will be increased by a
multiple `k`, and the `scale` will be updated by dividing by `k`.  
  
The element `control_column_matches` indicates, for each replicate
column of the calibrated primary survey, which column of replicate
weights it was matched to from the control survey. Columns which were
not matched to control survey replicate column are indicated by `NA`.  
  
The element `degf` will be set to match that of the primary survey to
ensure that the degrees of freedom are not erroneously inflated by
potential increases in the number of columns of replicate weights.

## Details

With the Opsomer-Erciulescu method, each column of replicate weights
from the control survey is randomly matched to a column of replicate
weights from the primary survey, and then the column from the primary
survey is calibrated to control totals estimated by perturbing the
control sample's full-sample estimates using the estimates from the
matched column of replicate weights from the control survey.  
  
If there are fewer columns of replicate weights in the control survey
than in the primary survey, then not all primary replicate columns will
be matched to a replicate column from the control survey.  

If there are more columns of replicate weights in the control survey
than in the primary survey, then the columns of replicate weights in the
primary survey will be duplicated `k` times, where `k` is the smallest
positive integer such that the resulting number of columns of replicate
weights for the primary survey is greater than or equal to the number of
columns of replicate weights in the control survey.  

Because replicate columns of the control survey are matched *at random*
to primary survey replicate columns, there are multiple ways to ensure
that this matching is reproducible. The user can either call
[set.seed](https://rdrr.io/r/base/Random.html) before using the
function, or supply a mapping to the argument `control_col_matches`.

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

Opsomer, J.D. and A. Erciulescu (2021). "Replication variance estimation
after sample-based calibration." **Survey Methodology**, *47*: 265-277.

## See also

\[calibrate_to_estimate()\] as an alternative when the control survey's
data are not available or it doesn't have replicate weights, but an
estimate and its variance-covariance are available.

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

# Calibrate totals for one categorical variable and one numeric ----

  calibrated_rep_design <- calibrate_to_sample(
    primary_rep_design = primary_survey,
    control_rep_design = control_survey,
    cal_formula = ~ stype + enroll,
  )
#> The primary survey has fewer replicates than the control survey, so columns in the primary survey will be duplicated 14 times, with suitable adjustments made to `scale` and `rscales`.
#> Matching between primary and control replicates will be done at random.
#> For tips on reproducible matching, see `help('calibrate_to_sample')`
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
#> api00 642.69 27.393
#> api99 606.91 28.302

# Inspect weights before and after calibration ----

  summarize_rep_weights(primary_survey, type = 'overall')
#>   nrows ncols degf_svy_pkg rank avg_wgt_sum sd_wgt_sums min_rep_wgt max_rep_wgt
#> 1   183    15           14   15        6194    403.1741           0    36.26464
  summarize_rep_weights(calibrated_rep_design, type = 'overall')
#>   nrows ncols degf_svy_pkg rank avg_wgt_sum  sd_wgt_sums min_rep_wgt
#> 1   183   210           47   48        6194 1.401234e-09           0
#>   max_rep_wgt
#> 1    124.8867

# For reproducibility, specify how to match replicates between surveys ----

  column_matching <- calibrated_rep_design$control_col_matches
  print(column_matching)
#> NULL

  calibrated_rep_design <- calibrate_to_sample(
    primary_rep_design = primary_survey,
    control_rep_design = control_survey,
    cal_formula = ~ stype + enroll,
    control_col_matches = column_matching
  )
#> The primary survey has fewer replicates than the control survey, so columns in the primary survey will be duplicated 14 times, with suitable adjustments made to `scale` and `rscales`.
#> Matching between primary and control replicates will be done at random.
#> For tips on reproducible matching, see `help('calibrate_to_sample')`
#> Warning: Setting `mse` to TRUE; variance estimates will be centered around full-sample estimate, not mean of replicates.
# }
```
