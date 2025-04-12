
<!-- README.md is generated from README.Rmd. Please edit that file -->

# svrep <a href="https://bschneidr.github.io/svrep/"><img src="man/figures/logo.png" align="right" height="139" alt="svrep website" /></a>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/svrep)](https://CRAN.R-project.org/package=svrep)
[![Mentioned in Awesome Official
Statistics](https://awesome.re/mentioned-badge.svg)](https://github.com/SNStatComp/awesome-official-statistics-software)
[![Codecov test
coverage](https://codecov.io/gh/bschneidr/svrep/branch/main/graph/badge.svg)](https://app.codecov.io/gh/bschneidr/svrep?branch=main)
<!-- [![R-CMD-check](https://github.com/bschneidr/svrep/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bschneidr/svrep/actions/workflows/R-CMD-check.yaml) -->
<!-- badges: end -->

svrep extends the survey package in R by implementing flexible
replication methods such as the generalized bootstrap that can be used
to analyze data from complex sample designs. Functions from svrep can be
used to create replicate weights and use special replication-based
methods of estimating variances from complex sample surveys.

The following kinds of replication methods are currently supported:

- **Bootstrap** replication, including the Rao-Wu-Yue-Beaumont (RWYB)
  bootstrap which can handle data from stratified, multistage samples
  with arbitrary sampling fractions and unequal probability sampling.

- **Random-groups jackknife** replication, a method that can be more
  efficient than the bootstrap in certain settings while being more
  robust than the standard delete-one jackknife when estimating
  variances for quantiles and other non-smooth statistics.

- **Fay’s generalized replication method**, a highly flexible and
  efficient replication method proposed by Fay (1989). This method can
  be used for stratified multistage samples with unequal probability
  sampling, as well samples selected with more complex methods such as
  two-phase sampling or the cube method.

- **The generalized survey bootstrap**, a bootstrap form of generalized
  replication which can handle a much wider range of complex sample
  design features compared to standard bootstrap methods.

- **Successive difference replication**, a replication method that is
  especially useful for surveys with systematic sampling or fine
  stratification. This is the method used by the largest household
  survey in the United States, the American Community Survey.

The svrep package also implements sample-based calibration methods,
wherein one sample’s weights are calibrated to estimates from another
sample, for example when calibrating a second-phase sample to a
first-phase sample or when calibrating one survey using estimates from a
separate benchmark survey. svrep also provides general data manipulation
and diagnostic tools to help survey statisticians work with replicate
weights and export them to a CSV or other output file format. For
bootstrap replicates specifically, there are also special tools to help
analysts identify the number of bootstrap replicates necessary to obtain
a desired precision level.

## Installation

You can install the released version of svrep from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("svrep")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bschneidr/svrep")
```

## Citation

When using the ‘svrep’ package, please make sure to cite it in any
resulting publications. This is appreciated by the package maintainer
and helps to incentivize ongoing development, maintenance, and support.

> Schneider B. (2023). “svrep: Tools for Creating, Updating, and
> Analyzing Survey Replicate Weights”. R package version 0.7.0.
> <https://doi.org/10.32614/CRAN.package.svrep>.

When using the ‘svrep’ package, please also cite the ‘survey’ package
and R itself, since they are essential to the use of ‘svrep’. Call
`citation('svrep')`, `citation('survey')`, and `citation('base')` for
more information and to generate BibTex entries for citing these
packages as well as R.

## Example usage

### Creating replicate weights

Suppose we have data from a survey selected using a complex sampling
method such as cluster sampling. To represent the complex survey design,
we can create a survey design object using the survey package.

``` r
library(survey)
library(svrep)
data(api, package = "survey")
set.seed(2021)

# Create a survey design object for a sample
# selected using a single-stage cluster sample without replacement
dclus1 <- svydesign(
  data    = apiclus1,
  ids     = ~dnum, 
  weights = ~pw, 
  fpc     = ~fpc
)
```

To help us estimate sampling variances, we can create bootstrap
replicate weights. The function `as_bootstrap_design()` creates
bootstrap replicate weights appropriate to common complex sampling
designs, using bootstrapping methods from the ‘survey’ package as well
as additional methods such as the Rao-Wu-Yue-Beaumont method (a
generalization of the Rao-Wu bootstrap).

``` r
# Create replicate-weights survey design
orig_rep_design <- dclus1 |> as_bootstrap_design( 
  replicates = 500,
  type       = "Rao-Wu-Yue-Beaumont"
)

print(orig_rep_design)
#> Call: as_bootstrap_design(dclus1, replicates = 500, type = "Rao-Wu-Yue-Beaumont")
#> Survey bootstrap with 500 replicates.
```

For especially complex survey designs (e.g., systematic samples), the
generalized survey bootstrap can be used.

``` r
# Load example data for a stratified systematic sample
data('library_stsys_sample', package = 'svrep')

# First, ensure data are sorted in same order as was used in sampling
library_stsys_sample <- library_stsys_sample[
  order(library_stsys_sample$SAMPLING_SORT_ORDER),
]

# Create a survey design object
design_obj <- svydesign(
  data   = library_stsys_sample,
  strata = ~ SAMPLING_STRATUM,
  ids    = ~ 1,
  fpc    = ~ STRATUM_POP_SIZE
)

# Convert to generalized bootstrap replicate design
gen_boot_design_sd2 <- as_gen_boot_design(
  design             = design_obj,
  variance_estimator = "SD2",
  replicates         = 500
)
#> For `variance_estimator='SD2', assumes rows of data are sorted in the same order used in sampling.
```

For relatively simple designs, we can also use the random-groups
jackknife.

``` r
# Create random-group jackknife replicates
# for a single-stage survey with many first-stage sampling units
rand_grp_jk_design <- apisrs |>
  svydesign(data    = _, 
            ids     = ~ 1,
            weights = ~ pw) |>
  as_random_group_jackknife_design(
    replicates = 20
  )
```

### Adjusting for non-response or unknown eligibility

In social surveys, unit nonresponse is extremely common. It is also
somewhat common for respondent cases to be classified as “ineligible”
for the survey based on their response. In general, sampled cases are
typically classified as “respondents”, “nonrespondents”, “ineligible
cases”, and “unknown eligibility” cases.

``` r
# Create variable giving response status
orig_rep_design$variables[['response_status']] <- sample(
  x       = c("Respondent", "Nonrespondent",
              "Ineligible", "Unknown eligibility"),
  prob    = c(0.6, 0.2, 0.1, 0.1),
  size    = nrow(orig_rep_design),
  replace = TRUE
)

table(orig_rep_design$variables$response_status)
#> 
#>          Ineligible       Nonrespondent          Respondent Unknown eligibility 
#>                  16                  32                 119                  16
```

It is common practice to adjust weights when there is non-response or
there are sampled cases whose eligibility for the survey is unknown. The
most common form of adjustment is “weight redistribution”: for example,
weights from non-respondents are reduced to zero, and weights from
respondents are correspondingly increased so that the total weight in
the sample is unchanged. In order to account for these adjustments when
estimating variances for survey statistics, the adjustments are repeated
separately for each set of replicate weights. This process can be easily
implemented using the `redistribute_weights()` function.

``` r
# Adjust weights for unknown eligibility
ue_adjusted_design <- redistribute_weights(
  design      = orig_rep_design,
  reduce_if   = response_status %in% c("Unknown eligibility"),
  increase_if = !response_status %in% c("Unknown eligibility")
)
```

By supplying column names to the `by` argument of
`redistribute_weights()`, adjustments are conducted separately in
different groups. This can be used to conduct nonresponse weighting
class adjustments.

``` r
nr_adjusted_design <- redistribute_weights(
  design      = ue_adjusted_design,
  reduce_if   = response_status == "Nonrespondent",
  increase_if = response_status == "Respondent",
  by          = c("stype")
)
```

### Comparing estimates from different sets of weights

In order to assess whether weighting adjustments have an impact on the
estimates we care about, we want to compare the estimates from the
different sets of weights. The function `svyby_repwts()` makes it easy
to compare estimates from different sets of weights.

``` r
# Estimate overall means (and their standard errors) from each design
overall_estimates <- svyby_repwts(
  rep_designs = list('original'             = orig_rep_design,
                     'nonresponse-adjusted' = nr_adjusted_design),
  formula     = ~ api00, 
  FUN         = svymean
)
print(overall_estimates, row.names = FALSE)
#>           Design_Name    api00       se
#>  nonresponse-adjusted 641.2030 25.54368
#>              original 644.1694 23.06284

# Estimate domain means (and their standard errors) from each design
domain_estimates <- svyby_repwts(
  rep_designs = list('original'             = orig_rep_design,
                     'nonresponse-adjusted' = nr_adjusted_design),
  formula = ~ api00, 
  by      = ~ stype, 
  FUN     = svymean
)
print(domain_estimates, row.names = FALSE)
#>           Design_Name stype    api00       se
#>  nonresponse-adjusted     E 649.9188 25.56366
#>              original     E 648.8681 22.31347
#>  nonresponse-adjusted     H 603.5390 45.26079
#>              original     H 618.5714 37.39448
#>  nonresponse-adjusted     M 616.3260 36.27983
#>              original     M 631.4400 31.03957
```

We can even test for differences in estimates from the two sets of
weights and calculate confidence intervals for their difference.

``` r
estimates <- svyby_repwts(
  rep_designs = list('original'             = orig_rep_design,
                     'nonresponse-adjusted' = nr_adjusted_design),
  formula     = ~ api00, 
  FUN         = svymean
)

vcov(estimates)
#>                      nonresponse-adjusted original
#> nonresponse-adjusted             652.4793 585.5253
#> original                         585.5253 531.8947

diff_between_ests <- svycontrast(
  stat      = estimates,
  contrasts = list(
    "Original vs. Adjusted" = c(-1,1)
  )
)
print(diff_between_ests)
#>                       contrast     SE
#> Original vs. Adjusted   2.9664 3.6501
confint(diff_between_ests)
#>                           2.5 %   97.5 %
#> Original vs. Adjusted -4.187705 10.12056
```

### Diagnosing potential issues with weights

When adjusting replicate weights, there are several diagnostics which
can be used to ensure that the adjustments were carried out correctly
and that they do more good than harm. The function
`summarize_rep_weights()` helps by allowing you to quickly summarize the
replicate weights.

For example, when carrying out nonresponse adjustments, we might want to
verify that all of the weights for nonrespondents have been set to zero
in each replicate. We can use the `summarize_rep_weights()` to compare
summary statistics for each replicate, and we can use its `by` argument
to group the summaries by one or more variables.

``` r
summarize_rep_weights(
  rep_design = nr_adjusted_design,
  type       = 'specific',
  by         = "response_status"
) |> 
  subset(Rep_Column %in% 1:2)
#>          response_status Rep_Column   N N_NONZERO       SUM     MEAN        CV       MIN       MAX
#> 1             Ineligible          1  16        16  608.1360 38.00850 1.2415437 0.5632079 120.38814
#> 2             Ineligible          2  16        16  739.2634 46.20397 0.7578107 0.5422029  77.44622
#> 501        Nonrespondent          1  32         0    0.0000  0.00000       NaN 0.0000000   0.00000
#> 502        Nonrespondent          2  32         0    0.0000  0.00000       NaN 0.0000000   0.00000
#> 1001          Respondent          1 119       119 6236.0577 52.40385 1.0431318 0.6072282 151.10496
#> 1002          Respondent          2 119       119 6426.4544 54.00382 0.8345243 0.5971008 102.40567
#> 1501 Unknown eligibility          1  16         0    0.0000  0.00000       NaN 0.0000000   0.00000
#> 1502 Unknown eligibility          2  16         0    0.0000  0.00000       NaN 0.0000000   0.00000
```

At the end of the adjustment process, we can inspect the number of rows
and columns and examine the variability of the weights across all of the
replicates.

``` r
nr_adjusted_design |>
  subset(response_status == "Respondent") |>
  summarize_rep_weights(
    type = 'overall'
  )
#>   nrows ncols degf_svy_pkg rank avg_wgt_sum sd_wgt_sums min_rep_wgt max_rep_wgt
#> 1   119   500           29   30    5625.555    1257.982   0.5305136     367.826
```

### Sample-based calibration

When we rake or poststratify to estimated control totals rather than to
“true” population values, we may need to account for the variance of the
estimated control totals to ensure that calibrated estimates
appropriately reflect sampling error of both the primary survey of
interest and the survey from which the control totals were estimated.
The ‘svrep’ package provides two functions which accomplish this. The
function `calibrate_to_estimate()` requires the user to supply a vector
of control totals and its variance-covariance matrix, while the function
`calibrate_to_sample()` requires the user to supply a dataset with
replicate weights to use for estimating control totals and their
sampling variance.

As an example, suppose we have a survey measuring vaccination status of
adults in Louisville, Kentucky. For variance estimation, we use 100
bootstrap replicates.

``` r
data("lou_vax_survey")

# Load example data
lou_vax_survey <- svydesign(
  data    = lou_vax_survey,
  ids     = ~ 1, 
  weights = ~ SAMPLING_WEIGHT
) |>
  as_bootstrap_design(
    replicates = 100, 
    mse        = TRUE
  )

# Adjust for nonresponse
lou_vax_survey <- lou_vax_survey |>
  redistribute_weights(
    reduce_if   = RESPONSE_STATUS == "Nonrespondent",
    increase_if = RESPONSE_STATUS == "Respondent"
  ) |>
  subset(RESPONSE_STATUS == "Respondent")
```

To reduce nonresponse bias or coverage error for the survey, we can rake
the survey to population totals for demographic groups estimated by the
Census Bureau in the American Community Survey (ACS). To estimate the
population totals for raking purposes, we can use microdata with
replicate weights.

``` r
# Load microdata to use for estimating control totals
data("lou_pums_microdata")

acs_benchmark_survey <- survey::svrepdesign(
  data       = lou_pums_microdata,
  variables  = ~ UNIQUE_ID + AGE + SEX + RACE_ETHNICITY + EDUC_ATTAINMENT,
  weights    = ~ PWGTP, 
  repweights = "PWGTP\\d{1,2}",
  type       = "successive-difference",
  mse        = TRUE
)
```

We can see that the distribution of race/ethnicity among respondents
differs from the distribution of race/ethnicity in the ACS benchmarks.

``` r
# Compare demographic estimates from the two data sources
estimate_comparisons <- data.frame(
  'Vax_Survey'    = lou_vax_survey |>
    svymean(x = ~ RACE_ETHNICITY) |> 
    coef(),
  'ACS_Benchmark' = acs_benchmark_survey |>
    svymean(x = ~ RACE_ETHNICITY) |> 
    coef()
)
rownames(estimate_comparisons) <- gsub(x = rownames(estimate_comparisons),
                                       "RACE_ETHNICITY", "")
print(estimate_comparisons)
#>                                                         Vax_Survey ACS_Benchmark
#> Black or African American alone, not Hispanic or Latino 0.16932271    0.19949824
#> Hispanic or Latino                                      0.03386454    0.04525039
#> Other Race, not Hispanic or Latino                      0.05776892    0.04630955
#> White alone, not Hispanic or Latino                     0.73904382    0.70894182
```

There are two options for calibrating the sample to the control totals
from the benchmark survey. With the first approach, we supply point
estimates and their variance-covariance matrix to the function
`calibrate_to_estimate()`.

``` r
# Estimate control totals and their variance-covariance matrix
control_totals <- svymean(
  x      = ~ RACE_ETHNICITY + EDUC_ATTAINMENT,
  design = acs_benchmark_survey
)
point_estimates <- coef(control_totals)
vcov_estimates <- vcov(control_totals)

# Calibrate the vaccination survey to the estimated control totals
vax_survey_raked_to_estimates <- calibrate_to_estimate(
  rep_design    = lou_vax_survey,
  estimate      = point_estimates,
  vcov_estimate = vcov_estimates,
  cal_formula   = ~ RACE_ETHNICITY + EDUC_ATTAINMENT,
  calfun        = survey::cal.raking
)
```

With the second approach, we supply the control survey’s replicate
design to `calibrate_to_sample()`.

``` r
vax_survey_raked_to_acs_sample <- calibrate_to_sample(
  primary_rep_design = lou_vax_survey,
  control_rep_design = acs_benchmark_survey,
  cal_formula        = ~ RACE_ETHNICITY + EDUC_ATTAINMENT,
  calfun             = survey::cal.raking
)
```

After calibration, we can see that the estimated vaccination rate has
decreased, and the estimated standard error of the estimated vaccination
rate has increased.

``` r
# Compare the two sets of estimates
svyby_repwts(
  rep_design = list(
    'NR-adjusted'       = lou_vax_survey,
    'Raked to estimate' = vax_survey_raked_to_estimates,
    'Raked to sample'   = vax_survey_raked_to_acs_sample
  ),
  formula    = ~ VAX_STATUS,
  FUN        = svymean,
  keep.names = FALSE
)
#>         Design_Name VAX_STATUSUnvaccinated VAX_STATUSVaccinated        se1        se2
#> 1       NR-adjusted              0.4621514            0.5378486 0.01870176 0.01870176
#> 2 Raked to estimate              0.4732623            0.5267377 0.01901224 0.01901224
#> 3   Raked to sample              0.4732623            0.5267377 0.01900022 0.01900022
```

### Saving results to a data file

Once we’re satisfied with the weights, we can create a data frame with
the analysis variables and columns of final full-sample weights and
replicate weights. This format is easy to export to data files that can
be loaded into R or other software later.

``` r
data_frame_with_final_weights <- vax_survey_raked_to_estimates |>
  as_data_frame_with_weights(
    full_wgt_name  = "RAKED_WGT",
    rep_wgt_prefix = "RAKED_REP_WGT_"
  )

# Preview first 10 column names
colnames(data_frame_with_final_weights) |> head(10)
#>  [1] "RESPONSE_STATUS" "RACE_ETHNICITY"  "SEX"             "EDUC_ATTAINMENT" "VAX_STATUS"      "SAMPLING_WEIGHT"
#>  [7] "RAKED_WGT"       "RAKED_REP_WGT_1" "RAKED_REP_WGT_2" "RAKED_REP_WGT_3"
```

``` r
# Write the data to a CSV file
write.csv(
  x    = data_frame_with_final_weights,
  file = "survey-data_with-updated-weights.csv"
)
```
