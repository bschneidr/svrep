# Calibrating to Estimated Control Totals

## Sample-based Calibration: An Introduction

Calibration weighting adjustments such as post-stratification or raking
are often helpful for reducing sampling variance or non-sampling errors
such as nonresponse bias. Typically, the benchmark data used for these
calibration adjustments are estimates published by agencies such as the
United States Census Bureau. For example, pollsters in the United States
frequently rake polling data so that estimates for variables such as age
or educational attainment match benchmark estimates from the American
Community Survey (ACS).

While benchmark data (also known as control totals) for raking and
calibration are often treated as the “true” population values, they are
usually themselves estimates with their own sampling variance or margin
of error. When we calibrate to estimated control totals rather than to
“true” population values, we may need to account for the variance of the
estimated control totals to ensure that calibrated estimates
appropriately reflect sampling error of both the primary survey of
interest and the survey from which the control totals were estimated.
This is especially important if the control totals have large margins of
error.

A handful of statistical methods have been developed for the problem of
conducting replication variance estimation after sample-based
calibration; see Opsomer and Erciulescu (2021) for a clear overview of
the literature on this topic. All of these methods apply calibration
weighting adjustment to full-sample weights and to each column of
replicate weights. The key “trick” of these methods is to adjust each
column of replicate weights to a slightly different set of control
totals, varying the control totals used across all of the replicates in
such a way that the variation across the columns is in a sense
proportionate to the sampling variance of the control totals.

These statistical methods differ in the way that they generate different
control totals for each column of replicate weights and in the type of
data they require the analyst to use. The method of Fuller (1998)
requires the analyst to have a variance-covariance matrix for the
estimated control totals, while the method of Opsomer and Erciulescu
(2021) requires the analyst to use the full dataset for the control
survey along with associated replicate weights.

## Functions for Implementing Sample-Based Calibration

The ‘svrep’ package provides two functions to implement sample-based
calibration.

With the function
[`calibrate_to_estimate()`](https://bschneidr.github.io/svrep/reference/calibrate_to_estimate.md),
adjustments to replicate weights are conducted using the method of
Fuller (1998), requiring a variance-covariance matrix for the estimated
control totals.

``` r
calibrate_to_estimate(
  rep_design = rep_design,
  estimate = vector_of_control_totals,
  vcov_estimate = variance_covariance_matrix_for_controls,
  cal_formula = ~ CALIBRATION_VARIABLE_1 + CALIBRATION_VARIABLE_2 + ...,
)
```

With the function
[`calibrate_to_sample()`](https://bschneidr.github.io/svrep/reference/calibrate_to_sample.md),
adjustments to replicate weights are conducted using the method proposed
by Opsomer and Erciulescu (2021), requiring a dataset with replicate
weights to use for estimating control totals and their sampling
variance.

``` r
calibrate_to_sample(
  primary_rep_design = primary_rep_design,
  control_rep_design = control_rep_design
  cal_formula = ~ CALIBRATION_VARIABLE_1 + CALIBRATION_VARIABLE_2 + ...,
)
```

For both functions, it is possible to use a variety of calibration
options from the `survey` package’s
[`calibrate()`](https://rdrr.io/pkg/survey/man/calibrate.html) function.
For example, the user can specify a specific calibration function to
use, such as `calfun = survey::cal.linear` to implement
post-stratification or `calfun = survey::cal.raking` to implement
raking. The `bounds` argument can be used to specify bounds for the
calibration weights, and the arguments such as `maxit` or `epsilon`
allow finer control over the Newton-Raphson algorithm used to implement
calibration.

## An Example Using a Vaccination Survey

To illustrate the different methods for conducting sample-based
calibration, we’ll use an example survey measuring Covid-19 vaccination
status and a handful of demographic variables, based on a simple random
sample of 1,000 residents of Louisville, Kentucky.

``` r
# Load the data
library(svrep)
#> Loading required package: survey
#> Loading required package: grid
#> Loading required package: Matrix
#> Loading required package: survival
#> 
#> Attaching package: 'survey'
#> The following object is masked from 'package:graphics':
#> 
#>     dotchart
data("lou_vax_survey")

# Inspect the first few rows
head(lou_vax_survey) |> knitr::kable()
```

| RESPONSE_STATUS | RACE_ETHNICITY                                          | SEX    | EDUC_ATTAINMENT       | VAX_STATUS | SAMPLING_WEIGHT |
|:----------------|:--------------------------------------------------------|:-------|:----------------------|:-----------|----------------:|
| Nonrespondent   | White alone, not Hispanic or Latino                     | Female | Less than high school | NA         |         596.702 |
| Nonrespondent   | Black or African American alone, not Hispanic or Latino | Female | High school or beyond | NA         |         596.702 |
| Respondent      | White alone, not Hispanic or Latino                     | Female | Less than high school | Vaccinated |         596.702 |
| Nonrespondent   | White alone, not Hispanic or Latino                     | Female | Less than high school | NA         |         596.702 |
| Nonrespondent   | White alone, not Hispanic or Latino                     | Female | High school or beyond | NA         |         596.702 |
| Respondent      | White alone, not Hispanic or Latino                     | Female | High school or beyond | Vaccinated |         596.702 |

For the purpose of variance estimation, we’ll create jackknife replicate
weights.

``` r
suppressPackageStartupMessages(
  library(survey)
)

lou_vax_survey_rep <- svydesign(
  data = lou_vax_survey,
  ids = ~ 1, weights = ~ SAMPLING_WEIGHT
) |> 
  as.svrepdesign(type = "JK1", mse = TRUE)
```

    #> Call: as.svrepdesign.default(svydesign(data = lou_vax_survey, ids = ~1, 
    #>     weights = ~SAMPLING_WEIGHT), type = "JK1", mse = TRUE)
    #> Unstratified cluster jacknife (JK1) with 1000 replicates and MSE variances.

Because the survey’s key outcome, vaccination status, is only measured
for respondents, we’ll do a quick nonresponse weighting adjustment to
help make reasonable estimates for this outcome.

``` r
# Conduct nonresponse weighting adjustment

nr_adjusted_design <- lou_vax_survey_rep |>
  redistribute_weights(
    reduce_if = RESPONSE_STATUS == "Nonrespondent",
    increase_if = RESPONSE_STATUS == "Respondent"
  ) |>
  subset(RESPONSE_STATUS == "Respondent")

# Inspect the result of the adjustment
rbind(
  'Original' = summarize_rep_weights(lou_vax_survey_rep, type = 'overall'),
  'NR-adjusted' = summarize_rep_weights(nr_adjusted_design, type = 'overall')
)[,c("nrows", "rank", "avg_wgt_sum", "sd_wgt_sums")]
#>             nrows rank avg_wgt_sum  sd_wgt_sums
#> Original     1000 1000      596702 0.000000e+00
#> NR-adjusted   502  502      596702 8.219437e-11
```

All of the work so far has given us the replicate design for the primary
survey, prepared for calibration. Now we need to obtain benchmark data
we can use for the calibration. We’ll use a Public-Use Microdata Sample
(PUMS) dataset from the ACS as our source for benchmark data on
race/ethnicity, sex, and educational attainment.

``` r
data("lou_pums_microdata")
```

``` r
# Inspect some of the rows/columns of data ----
tail(lou_pums_microdata, n = 5) |> 
  dplyr::select(AGE, SEX, RACE_ETHNICITY, EDUC_ATTAINMENT) |>
  knitr::kable()
```

| AGE | SEX    | RACE_ETHNICITY                     | EDUC_ATTAINMENT       |
|----:|:-------|:-----------------------------------|:----------------------|
|  20 | Female | Other Race, not Hispanic or Latino | High school or beyond |
|  50 | Male   | Hispanic or Latino                 | Less than high school |
|  57 | Female | Other Race, not Hispanic or Latino | High school or beyond |
|  44 | Male   | Hispanic or Latino                 | High school or beyond |
|  25 | Female | Hispanic or Latino                 | High school or beyond |

Next, we’ll prepare the PUMS data to use replication variance estimation
using provided replicate weights.

``` r
# Convert to a survey design object ----
  pums_rep_design <- svrepdesign(
      data = lou_pums_microdata,
      weights = ~ PWGTP,
      repweights = "PWGTP\\d{1,2}",
      type = "successive-difference",
      variables = ~ AGE + SEX + RACE_ETHNICITY + EDUC_ATTAINMENT,
      mse = TRUE
    )

  pums_rep_design
#> Call: svrepdesign.default(data = lou_pums_microdata, weights = ~PWGTP, 
#>     repweights = "PWGTP\\d{1,2}", type = "successive-difference", 
#>     variables = ~AGE + SEX + RACE_ETHNICITY + EDUC_ATTAINMENT, 
#>     mse = TRUE)
#> with 80 replicates and MSE variances.
```

When conduction calibration, we have to make sure that the data from the
control survey represent the same population as the primary survey.
Since the Louisville vaccination survey only represents adults, we need
to subset the control survey design to adults.

``` r
# Subset to only include adults
pums_rep_design <- pums_rep_design |> subset(AGE >= 18)
```

In addition, we need to ensure that the control survey design has
calibration variables that align with the variables in the primary
survey design of interest. This may require some data manipulation.

``` r
suppressPackageStartupMessages(
  library(dplyr)
)

# Check that variables match across data sources ----
  pums_rep_design$variables |>
    dplyr::distinct(RACE_ETHNICITY)
#>                                            RACE_ETHNICITY
#> 1 Black or African American alone, not Hispanic or Latino
#> 2                     White alone, not Hispanic or Latino
#> 3                                      Hispanic or Latino
#> 4                      Other Race, not Hispanic or Latino

  setdiff(lou_vax_survey_rep$variables$RACE_ETHNICITY,
          pums_rep_design$variables$RACE_ETHNICITY)
#> character(0)
  setdiff(lou_vax_survey_rep$variables$SEX,
          pums_rep_design$variables$SEX)
#> character(0)
  setdiff(lou_vax_survey_rep$variables$EDUC_ATTAINMENT,
          pums_rep_design$variables$EDUC_ATTAINMENT)
#> character(0)
```

``` r
# Estimates from the control survey (ACS)
svymean(
  design = pums_rep_design,
  x = ~ RACE_ETHNICITY + SEX + EDUC_ATTAINMENT
)
#>                                                                          mean
#> RACE_ETHNICITYBlack or African American alone, not Hispanic or Latino 0.19950
#> RACE_ETHNICITYHispanic or Latino                                      0.04525
#> RACE_ETHNICITYOther Race, not Hispanic or Latino                      0.04631
#> RACE_ETHNICITYWhite alone, not Hispanic or Latino                     0.70894
#> SEXMale                                                               0.47543
#> SEXFemale                                                             0.52457
#> EDUC_ATTAINMENTHigh school or beyond                                  0.38736
#> EDUC_ATTAINMENTLess than high school                                  0.61264
#>                                                                           SE
#> RACE_ETHNICITYBlack or African American alone, not Hispanic or Latino 0.0010
#> RACE_ETHNICITYHispanic or Latino                                      0.0002
#> RACE_ETHNICITYOther Race, not Hispanic or Latino                      0.0008
#> RACE_ETHNICITYWhite alone, not Hispanic or Latino                     0.0007
#> SEXMale                                                               0.0007
#> SEXFemale                                                             0.0007
#> EDUC_ATTAINMENTHigh school or beyond                                  0.0033
#> EDUC_ATTAINMENTLess than high school                                  0.0033

# Estimates from the primary survey (Louisville vaccination survey)
svymean(
  design = nr_adjusted_design,
  x = ~ RACE_ETHNICITY + SEX + EDUC_ATTAINMENT
)
#>                                                                           mean
#> RACE_ETHNICITYBlack or African American alone, not Hispanic or Latino 0.169323
#> RACE_ETHNICITYHispanic or Latino                                      0.033865
#> RACE_ETHNICITYOther Race, not Hispanic or Latino                      0.057769
#> RACE_ETHNICITYWhite alone, not Hispanic or Latino                     0.739044
#> SEXFemale                                                             0.535857
#> SEXMale                                                               0.464143
#> EDUC_ATTAINMENTHigh school or beyond                                  0.458167
#> EDUC_ATTAINMENTLess than high school                                  0.541833
#>                                                                           SE
#> RACE_ETHNICITYBlack or African American alone, not Hispanic or Latino 0.0168
#> RACE_ETHNICITYHispanic or Latino                                      0.0081
#> RACE_ETHNICITYOther Race, not Hispanic or Latino                      0.0104
#> RACE_ETHNICITYWhite alone, not Hispanic or Latino                     0.0196
#> SEXFemale                                                             0.0223
#> SEXMale                                                               0.0223
#> EDUC_ATTAINMENTHigh school or beyond                                  0.0223
#> EDUC_ATTAINMENTLess than high school                                  0.0223
```

### Raking to estimated control totals

We’ll start by raking to estimates from the ACS for race/ethnicity, sex,
and educational attainment, first using the
[`calibrate_to_sample()`](https://bschneidr.github.io/svrep/reference/calibrate_to_sample.md)
method and then using the
[`calibrate_to_estimate()`](https://bschneidr.github.io/svrep/reference/calibrate_to_estimate.md)
method. For the
[`calibrate_to_sample()`](https://bschneidr.github.io/svrep/reference/calibrate_to_sample.md)
method, we need to obtain a vector of point estimates for the control
totals, and an accompanying variance-covariance matrix for the
estimates.

``` r
acs_control_totals <- svytotal(
  x = ~ RACE_ETHNICITY + SEX + EDUC_ATTAINMENT,
  design = pums_rep_design
)

control_totals_for_raking <- list(
  'estimates' = coef(acs_control_totals),
  'variance-covariance' = vcov(acs_control_totals)
)

# Inspect point estimates
control_totals_for_raking$estimates
#> RACE_ETHNICITYBlack or African American alone, not Hispanic or Latino 
#>                                                                119041 
#>                                      RACE_ETHNICITYHispanic or Latino 
#>                                                                 27001 
#>                      RACE_ETHNICITYOther Race, not Hispanic or Latino 
#>                                                                 27633 
#>                     RACE_ETHNICITYWhite alone, not Hispanic or Latino 
#>                                                                423027 
#>                                                               SEXMale 
#>                                                                283688 
#>                                                             SEXFemale 
#>                                                                313014 
#>                                  EDUC_ATTAINMENTHigh school or beyond 
#>                                                                231136 
#>                                  EDUC_ATTAINMENTLess than high school 
#>                                                                365566

# Inspect a few rows of the control totals' variance-covariance matrix
control_totals_for_raking$`variance-covariance`[5:8,5:8] |>
  `colnames<-`(NULL)
#>                                           [,1]      [,2]        [,3]       [,4]
#> SEXMale                              355572.45 -29522.95   129208.95   196840.6
#> SEXFemale                            -29522.95 379494.65    81455.95   268515.8
#> EDUC_ATTAINMENTHigh school or beyond 129208.95  81455.95  4019242.10 -3808577.2
#> EDUC_ATTAINMENTLess than high school 196840.55 268515.75 -3808577.20  4273933.5
```

Crucially, we note that the vector of control totals has the same names
as the estimates produced by using
[`svytotal()`](https://rdrr.io/pkg/survey/man/surveysummary.html) with
the primary survey design object whose weights we plan to adjust.

``` r
svytotal(x = ~ RACE_ETHNICITY + SEX + EDUC_ATTAINMENT,
         design = nr_adjusted_design)
#>                                                                        total
#> RACE_ETHNICITYBlack or African American alone, not Hispanic or Latino 101035
#> RACE_ETHNICITYHispanic or Latino                                       20207
#> RACE_ETHNICITYOther Race, not Hispanic or Latino                       34471
#> RACE_ETHNICITYWhite alone, not Hispanic or Latino                     440989
#> SEXFemale                                                             319747
#> SEXMale                                                               276955
#> EDUC_ATTAINMENTHigh school or beyond                                  273389
#> EDUC_ATTAINMENTLess than high school                                  323313
#>                                                                            SE
#> RACE_ETHNICITYBlack or African American alone, not Hispanic or Latino 10003.0
#> RACE_ETHNICITYHispanic or Latino                                       4824.4
#> RACE_ETHNICITYOther Race, not Hispanic or Latino                       6222.7
#> RACE_ETHNICITYWhite alone, not Hispanic or Latino                     11713.1
#> SEXFemale                                                             13301.6
#> SEXMale                                                               13301.6
#> EDUC_ATTAINMENTHigh school or beyond                                  13289.2
#> EDUC_ATTAINMENTLess than high school                                  13289.2
```

To calibrate the design to the estimates, we supply the estimates and
the variance-covariance matrix to
[`calibrate_to_estimate()`](https://bschneidr.github.io/svrep/reference/calibrate_to_estimate.md),
and we supply the `cal_formula` argument with the same formula we would
use for
[`svytotal()`](https://rdrr.io/pkg/survey/man/surveysummary.html). To
use a raking adjustment, we specify `calfun = survey::cal.raking`.

``` r
raked_design <- calibrate_to_estimate(
  rep_design = nr_adjusted_design,
  estimate = control_totals_for_raking$estimates,
  vcov_estimate = control_totals_for_raking$`variance-covariance`,
  cal_formula = ~ RACE_ETHNICITY + SEX + EDUC_ATTAINMENT,
  calfun = survey::cal.raking, # Required for raking
  epsilon = 1e-9
)
#> Selection of replicate columns whose control totals will be perturbed will be done at random.
#> For tips on reproducible selection, see `help('calibrate_to_estimate')`
```

Now we can compare the estimated totals for the calibration variables to
the actual control totals. As we might intuitively expect, the estimated
totals from the survey now match the control totals, and the standard
errors for the estimated totals match the standard errors of the control
totals.

``` r
# Estimated totals after calibration
svytotal(x = ~ RACE_ETHNICITY + SEX + EDUC_ATTAINMENT,
         design = raked_design)
#>                                                                        total
#> RACE_ETHNICITYBlack or African American alone, not Hispanic or Latino 119041
#> RACE_ETHNICITYHispanic or Latino                                       27001
#> RACE_ETHNICITYOther Race, not Hispanic or Latino                       27633
#> RACE_ETHNICITYWhite alone, not Hispanic or Latino                     423027
#> SEXFemale                                                             313014
#> SEXMale                                                               283688
#> EDUC_ATTAINMENTHigh school or beyond                                  231136
#> EDUC_ATTAINMENTLess than high school                                  365566
#>                                                                            SE
#> RACE_ETHNICITYBlack or African American alone, not Hispanic or Latino  633.63
#> RACE_ETHNICITYHispanic or Latino                                       107.98
#> RACE_ETHNICITYOther Race, not Hispanic or Latino                       472.41
#> RACE_ETHNICITYWhite alone, not Hispanic or Latino                      594.14
#> SEXFemale                                                              616.03
#> SEXMale                                                                596.30
#> EDUC_ATTAINMENTHigh school or beyond                                  2004.80
#> EDUC_ATTAINMENTLess than high school                                  2067.35

# Matches the control totals!
cbind(
  'total' = control_totals_for_raking$estimates,
  'SE' = control_totals_for_raking$`variance-covariance` |>
    diag() |> sqrt()
)
#>                                                                        total
#> RACE_ETHNICITYBlack or African American alone, not Hispanic or Latino 119041
#> RACE_ETHNICITYHispanic or Latino                                       27001
#> RACE_ETHNICITYOther Race, not Hispanic or Latino                       27633
#> RACE_ETHNICITYWhite alone, not Hispanic or Latino                     423027
#> SEXMale                                                               283688
#> SEXFemale                                                             313014
#> EDUC_ATTAINMENTHigh school or beyond                                  231136
#> EDUC_ATTAINMENTLess than high school                                  365566
#>                                                                              SE
#> RACE_ETHNICITYBlack or African American alone, not Hispanic or Latino  633.6287
#> RACE_ETHNICITYHispanic or Latino                                       107.9829
#> RACE_ETHNICITYOther Race, not Hispanic or Latino                       472.4107
#> RACE_ETHNICITYWhite alone, not Hispanic or Latino                      594.1448
#> SEXMale                                                                596.2990
#> SEXFemale                                                              616.0314
#> EDUC_ATTAINMENTHigh school or beyond                                  2004.8048
#> EDUC_ATTAINMENTLess than high school                                  2067.3494
```

We can now see what effect the raking adjustment has had on our primary
estimate of interest, which is the overall Covid-19 vaccination rate.
The raking adjustment has reduced our estimate of the vaccination rate
by about one percentage point and results in a similar standard error
estimate.

``` r
estimates_by_design <- svyby_repwts(
  rep_designs = list(
    "NR-adjusted" = nr_adjusted_design,
    "Raked" = raked_design
  ),
  FUN = svytotal,
  formula = ~ RACE_ETHNICITY + SEX + EDUC_ATTAINMENT
)

t(estimates_by_design[,-1]) |>
  knitr::kable()
```

|                                                                       | NR-adjusted |       Raked |
|:----------------------------------------------------------------------|------------:|------------:|
| RACE_ETHNICITYBlack or African American alone, not Hispanic or Latino |  101035.199 | 119041.0000 |
| RACE_ETHNICITYHispanic or Latino                                      |   20207.040 |  27001.0000 |
| RACE_ETHNICITYOther Race, not Hispanic or Latino                      |   34470.833 |  27633.0000 |
| RACE_ETHNICITYWhite alone, not Hispanic or Latino                     |  440988.928 | 423027.0000 |
| SEXFemale                                                             |  319746.689 | 313014.0000 |
| SEXMale                                                               |  276955.311 | 283688.0000 |
| EDUC_ATTAINMENTHigh school or beyond                                  |  273389.363 | 231136.0000 |
| EDUC_ATTAINMENTLess than high school                                  |  323312.637 | 365566.0000 |
| se1                                                                   |   10002.951 |    633.6287 |
| se2                                                                   |    4824.430 |    107.9829 |
| se3                                                                   |    6222.719 |    472.4107 |
| se4                                                                   |   11713.138 |    594.1448 |
| se5                                                                   |   13301.627 |    616.0314 |
| se6                                                                   |   13301.627 |    596.2990 |
| se7                                                                   |   13289.206 |   2004.8048 |
| se8                                                                   |   13289.206 |   2067.3494 |

Instead of doing the raking using a vector of control totals and their
variance-covariance matrix, we could have instead done the raking by
simply supplying the two replicate design objects to the function
[`calibrate_to_sample()`](https://bschneidr.github.io/svrep/reference/calibrate_to_sample.md).
This uses the Opsomer-Erciulescu method of adjusting replicate weights,
in contrast to
[`calibrate_to_estimate()`](https://bschneidr.github.io/svrep/reference/calibrate_to_estimate.md),
which uses Fuller’s method of adjusting replicate weights.

``` r
raked_design_opsomer_erciulescu <- calibrate_to_sample(
  primary_rep_design = nr_adjusted_design,
  control_rep_design = pums_rep_design,
  cal_formula = ~ RACE_ETHNICITY + SEX + EDUC_ATTAINMENT,
  calfun = survey::cal.raking,
  epsilon = 1e-9
)
#> Matching between primary and control replicates will be done at random.
#> For tips on reproducible matching, see `help('calibrate_to_sample')`
```

We can see that the two methods yield identical point estimates from the
full-sample weights, and the standard errors match nearly exactly for
the calibration variables (race/ethnicity, sex, and educational
attainment). However, there are small but slightly more noticeable
differences in the standard errors for other variables, such as
`VAX_STATUS`, resulting from the fact that the two methods have
different methods of adjusting the replicate weights. Opsomer and
Erciulescu (2021) explain the differences between the two methods and
discuss why the the Opsomer-Erciulescu method used in
[`calibrate_to_sample()`](https://bschneidr.github.io/svrep/reference/calibrate_to_sample.md)
may have better statistical properties than the Fuller method used in
[`calibrate_to_estimate()`](https://bschneidr.github.io/svrep/reference/calibrate_to_estimate.md).

``` r
estimates_by_design <- svyby_repwts(
  rep_designs = list(
    "calibrate_to_estimate()" = raked_design,
    "calibrate_to_sample()" = raked_design_opsomer_erciulescu
  ),
  FUN = svytotal,
  formula = ~ VAX_STATUS + RACE_ETHNICITY + SEX + EDUC_ATTAINMENT
)

t(estimates_by_design[,-1]) |>
  knitr::kable()
```

|                                                                       | calibrate_to_estimate() | calibrate_to_sample() |
|:----------------------------------------------------------------------|------------------------:|----------------------:|
| VAX_STATUSUnvaccinated                                                |             282904.4084 |           282904.4084 |
| VAX_STATUSVaccinated                                                  |             313797.5916 |           313797.5916 |
| RACE_ETHNICITYBlack or African American alone, not Hispanic or Latino |             119041.0000 |           119041.0000 |
| RACE_ETHNICITYHispanic or Latino                                      |              27001.0000 |            27001.0000 |
| RACE_ETHNICITYOther Race, not Hispanic or Latino                      |              27633.0000 |            27633.0000 |
| RACE_ETHNICITYWhite alone, not Hispanic or Latino                     |             423027.0000 |           423027.0000 |
| SEXFemale                                                             |             313014.0000 |           313014.0000 |
| SEXMale                                                               |             283688.0000 |           283688.0000 |
| EDUC_ATTAINMENTHigh school or beyond                                  |             231136.0000 |           231136.0000 |
| EDUC_ATTAINMENTLess than high school                                  |             365566.0000 |           365566.0000 |
| se1                                                                   |              13429.6614 |            13386.1701 |
| se2                                                                   |              13392.3319 |            13397.3449 |
| se3                                                                   |                633.6287 |              633.6287 |
| se4                                                                   |                107.9829 |              107.9829 |
| se5                                                                   |                472.4107 |              472.4107 |
| se6                                                                   |                594.1448 |              594.1448 |
| se7                                                                   |                616.0314 |              616.0314 |
| se8                                                                   |                596.2990 |              596.2990 |
| se9                                                                   |               2004.8048 |             2004.8048 |
| se10                                                                  |               2067.3494 |             2067.3494 |

### Post-stratification

The primary difference between post-stratification and raking is that
post-stratification essentially involves only a single calibration
variable, with population benchmarks provided for each value of that
variable. In the Louisville vaccination survey, that variable is called
`POSTSTRATUM` and is based on combinations of race/ethnicity, sex, and
educational attainment.

``` r
# Create matching post-stratification variable in both datasets
  nr_adjusted_design <- nr_adjusted_design |>
    transform(POSTSTRATUM = interaction(RACE_ETHNICITY, SEX, EDUC_ATTAINMENT,
                                        sep = "|"))

  pums_rep_design <- pums_rep_design |>
    transform(POSTSTRATUM = interaction(RACE_ETHNICITY, SEX, EDUC_ATTAINMENT,
                                        sep = "|"))
  
  levels(pums_rep_design$variables$POSTSTRATUM) <- levels(
    nr_adjusted_design$variables$POSTSTRATUM
  )

# Estimate control totals
  acs_control_totals <- svytotal(
    x = ~ POSTSTRATUM,
    design = pums_rep_design
  )
  
  poststratification_totals <- list(
    'estimate' = coef(acs_control_totals),
    'variance-covariance' = vcov(acs_control_totals)
  )

# Inspect the control totals
  poststratification_totals$estimate |>
    as.data.frame() |>
    `colnames<-`('estimate') |>
    knitr::kable()
```

|                                                                                                   | estimate |
|:--------------------------------------------------------------------------------------------------|---------:|
| POSTSTRATUMBlack or African American alone, not Hispanic or Latino\|Female\|High school or beyond |    12168 |
| POSTSTRATUMHispanic or Latino\|Female\|High school or beyond                                      |     3998 |
| POSTSTRATUMOther Race, not Hispanic or Latino\|Female\|High school or beyond                      |     6190 |
| POSTSTRATUMWhite alone, not Hispanic or Latino\|Female\|High school or beyond                     |    84041 |
| POSTSTRATUMBlack or African American alone, not Hispanic or Latino\|Male\|High school or beyond   |    17648 |
| POSTSTRATUMHispanic or Latino\|Male\|High school or beyond                                        |     4132 |
| POSTSTRATUMOther Race, not Hispanic or Latino\|Male\|High school or beyond                        |     6687 |
| POSTSTRATUMWhite alone, not Hispanic or Latino\|Male\|High school or beyond                       |    96272 |
| POSTSTRATUMBlack or African American alone, not Hispanic or Latino\|Female\|Less than high school |    41944 |
| POSTSTRATUMHispanic or Latino\|Female\|Less than high school                                      |    10321 |
| POSTSTRATUMOther Race, not Hispanic or Latino\|Female\|Less than high school                      |     6753 |
| POSTSTRATUMWhite alone, not Hispanic or Latino\|Female\|Less than high school                     |   118273 |
| POSTSTRATUMBlack or African American alone, not Hispanic or Latino\|Male\|Less than high school   |    47281 |
| POSTSTRATUMHispanic or Latino\|Male\|Less than high school                                        |     8550 |
| POSTSTRATUMOther Race, not Hispanic or Latino\|Male\|Less than high school                        |     8003 |
| POSTSTRATUMWhite alone, not Hispanic or Latino\|Male\|Less than high school                       |   124441 |

To post-stratify the design, we can either supply the estimates and
their variance-covariance matrix to
[`calibrate_to_estimate()`](https://bschneidr.github.io/svrep/reference/calibrate_to_estimate.md),
or we can supply the two replicate design objects to
[`calibrate_to_sample()`](https://bschneidr.github.io/svrep/reference/calibrate_to_sample.md).
With either method, we need to supply the `cal_formula` argument with
the same formula we would use for
[`svytotal()`](https://rdrr.io/pkg/survey/man/surveysummary.html). To
use a post-stratification adjustment (rather than raking), we specify
`calfun = survey::cal.linear`.

``` r
# Post-stratify the design using the estimates
poststrat_design_fuller <- calibrate_to_estimate(
  rep_design = nr_adjusted_design,
  estimate = poststratification_totals$estimate,
  vcov_estimate = poststratification_totals$`variance-covariance`,
  cal_formula = ~ POSTSTRATUM, # Specify the post-stratification variable
  calfun = survey::cal.linear # This option is required for post-stratification
)
#> Selection of replicate columns whose control totals will be perturbed will be done at random.
#> For tips on reproducible selection, see `help('calibrate_to_estimate')`
```

``` r
# Post-stratify the design using the two samples
poststrat_design_opsomer_erciulescu <- calibrate_to_sample(
  primary_rep_design = nr_adjusted_design,
  control_rep_design = pums_rep_design,
  cal_formula = ~ POSTSTRATUM, # Specify the post-stratification variable
  calfun = survey::cal.linear # This option is required for post-stratification
)
#> Matching between primary and control replicates will be done at random.
#> For tips on reproducible matching, see `help('calibrate_to_sample')`
```

As with the raking example, we can see that the full-sample
post-stratified estimates are exactly the same for the two methods. The
standard errors for post-stratification variables are essentially
identical, while the standard errors for other variables differ
slightly.

``` r
estimates_by_design <- svyby_repwts(
  rep_designs = list(
    "calibrate_to_estimate()" = poststrat_design_fuller,
    "calibrate_to_sample()" = poststrat_design_opsomer_erciulescu
  ),
  FUN = svymean,
  formula = ~ VAX_STATUS + RACE_ETHNICITY + SEX + EDUC_ATTAINMENT
)

t(estimates_by_design[,-1]) |>
  knitr::kable()
```

|                                                                       | calibrate_to_estimate() | calibrate_to_sample() |
|:----------------------------------------------------------------------|------------------------:|----------------------:|
| VAX_STATUSUnvaccinated                                                |               0.4779776 |             0.4779776 |
| VAX_STATUSVaccinated                                                  |               0.5220224 |             0.5220224 |
| RACE_ETHNICITYBlack or African American alone, not Hispanic or Latino |               0.1994982 |             0.1994982 |
| RACE_ETHNICITYHispanic or Latino                                      |               0.0452504 |             0.0452504 |
| RACE_ETHNICITYOther Race, not Hispanic or Latino                      |               0.0463095 |             0.0463095 |
| RACE_ETHNICITYWhite alone, not Hispanic or Latino                     |               0.7089418 |             0.7089418 |
| SEXFemale                                                             |               0.4754266 |             0.4754266 |
| SEXMale                                                               |               0.5245734 |             0.5245734 |
| EDUC_ATTAINMENTHigh school or beyond                                  |               0.3873558 |             0.3873558 |
| EDUC_ATTAINMENTLess than high school                                  |               0.6126442 |             0.6126442 |
| se1                                                                   |               0.0234262 |             0.0234672 |
| se2                                                                   |               0.0234262 |             0.0234672 |
| se3                                                                   |               0.0009765 |             0.0009766 |
| se4                                                                   |               0.0001857 |             0.0001856 |
| se5                                                                   |               0.0007811 |             0.0007809 |
| se6                                                                   |               0.0007040 |             0.0007039 |
| se7                                                                   |               0.0007464 |             0.0007464 |
| se8                                                                   |               0.0007464 |             0.0007464 |
| se9                                                                   |               0.0033336 |             0.0033339 |
| se10                                                                  |               0.0033336 |             0.0033339 |

## Reproducibility

The calibration methods for
[`calibrate_to_estimate()`](https://bschneidr.github.io/svrep/reference/calibrate_to_estimate.md)
and
[`calibrate_to_sample()`](https://bschneidr.github.io/svrep/reference/calibrate_to_sample.md)
involve one element of randomization: determining which columns of
replicate weights are assigned to a given perturbation of the control
totals. In the
[`calibrate_to_sample()`](https://bschneidr.github.io/svrep/reference/calibrate_to_sample.md)
method of Fuller (1998), if the control totals are a vector of dimension
$p$, then $p$ columns of replicate weights will be calibrated to $p$
different vectors of perturbed control totals, formed using the $p$
scaled eigenvectors from a spectral decomposition of the control totals’
variance-covariance matrix (sorted in order by the largest to smallest
eigenvalues). To control which columns of replicate weights will be
calibrated to each set of perturbed control totals, we can use the
function argument `col_selection`.

``` r
# Randomly select which columns will be assigned to each set of perturbed control totals
dimension_of_control_totals <- length(poststratification_totals$estimate)

columns_to_perturb <- sample(x = 1:ncol(nr_adjusted_design$repweights),
                             size = dimension_of_control_totals)

print(columns_to_perturb)
#>  [1] 466 549 364 644 851  56 308 697 657 471 543 477 692 458  24 951

# Perform the calibration
poststratified_design <- calibrate_to_estimate(
  rep_design = nr_adjusted_design,
  estimate = poststratification_totals$estimate,
  vcov_estimate = poststratification_totals$`variance-covariance`,
  cal_formula = ~ POSTSTRATUM,
  calfun = survey::cal.linear,
  col_selection = columns_to_perturb # Specified for reproducibility
)
```

The calibrated survey design object contains an element
`perturbed_control_cols` which indicates which columns were calibrated
to the perturbed control totals; this can be useful to save and use as
an input to `col_selection` to ensure reproducibility.

``` r
poststratified_design$perturbed_control_cols
#> NULL
```

For
[`calibrate_to_sample()`](https://bschneidr.github.io/svrep/reference/calibrate_to_sample.md),
matching is done between columns of replicate weights in the primary
survey and columns of replicate weights in the control survey. The
matching is done at random unless the user specifies otherwise using the
argument `control_col_matches`. In the Louisville Vaccination Survey,
the primary survey has 1,000 replicates while the control survey has 80
columns. So we can match these 80 columns to the 1,000 replicates by
specifying 1,000 values consisting of `NA` or integers between 1 and 80.

``` r
# Randomly match the primary replicates to control replicates
set.seed(1999)

column_matching <- rep(NA, times = ncol(nr_adjusted_design$repweights))
column_matching[sample(x = 1:1000, size = 80)] <- 1:80

str(column_matching)
#>  int [1:1000] NA NA NA 34 NA NA NA 68 NA NA ...

# Perform the calibration
poststratified_design <- calibrate_to_sample(
  primary_rep_design = nr_adjusted_design,
  control_rep_design = pums_rep_design,
  cal_formula = ~ POSTSTRATUM,
  calfun = survey::cal.linear,
  control_col_matches = column_matching
)
```

The calibrated survey design object contains an element
`control_column_matches` which control survey replicate each primary
survey replicate column was matched to.

``` r
str(poststratified_design$control_column_matches)
#>  int [1:1000] NA NA NA 34 NA NA NA 68 NA NA ...
```

## References

Fuller, Wayne A. 1998. “Replication Variance Estimation for Two-Phase
Samples.” *Statistica Sinica* 8 (4): 1153–64.

Opsomer, J. D., and A. L. Erciulescu. 2021. “Replication Variance
Estimation After Sample-Based Calibration.” *Survey Methodology,
Statistics Canada* Vol. 47 (No. 2).
