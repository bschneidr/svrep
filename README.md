
<!-- README.md is generated from README.Rmd. Please edit that file -->

# svrep

<!-- badges: start -->
<!-- [![R-CMD-check](https://github.com/bschneidr/svrep/workflows/R-CMD-check/badge.svg)](https://github.com/bschneidr/svrep/actions) -->

[![Codecov test
coverage](https://codecov.io/gh/bschneidr/svrep/branch/main/graph/badge.svg)](https://app.codecov.io/gh/bschneidr/svrep?branch=main)
<!-- badges: end -->

svrep provides methods for creating, updating, and analyzing replicate
weights for surveys. Functions from svrep can be used to implement
adjustments to replicate designs (e.g. nonresponse weighting class
adjustments) and analyze their effect on the replicate weights and on
estimates of interest.

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

## Example usage

Suppose we have a replicate-weights survey design object created with
the survey package. This survey design object can include respondents,
non-respondents, and cases with unknown eligibility.

``` r
library(survey)
library(svrep)
data(api, package = "survey")
set.seed(2021)

# Create variable giving response status
apiclus1$response_status <- sample(x = c("Respondent", "Nonrespondent",
                                         "Ineligible", "Unknown eligibility"),
                                   size = nrow(apiclus1),
                                   replace = TRUE)

# Create replicate-weights survey design
dclus1 <- svydesign(data = apiclus1,
                    id = ~dnum, weights = ~pw, fpc = ~fpc)

orig_rep_design <- as.svrepdesign(dclus1)

print(orig_rep_design)
#> Call: as.svrepdesign.default(dclus1)
#> Unstratified cluster jacknife (JK1) with 15 replicates.
```

### Adjusting for non-response or unknown eligibility

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
ue_adjusted_design <- redistribute_weights(design = orig_rep_design,
                                           reduce_if = response_status %in% c("Unknown eligibility"),
                                           increase_if = !response_status %in% c("Unknown eligibility"))
```

By supplying column names to the `by` argument of
`redistribute_weights()`, adjustments are conducted separately in
different groups. This can be used to conduct nonresponse weighting
class adjustments.

``` r
nr_adjusted_design <- redistribute_weights(design = ue_adjusted_design,
                                           reduce_if = response_status == "Nonrespondent",
                                           increase_if = response_status == "Respondent",
                                           by = c("stype"))
```

### Comparing estimates from different sets of weights

In order to assess whether weighting adjustments have an impact on the
estimates we care about, we want to compare the estimates from the
different sets of weights. The function `svyby_repwts()` makes it easy
to compare estimates from different sets of weights.

``` r
# Estimate overall means (and their standard errors) from each design
overall_estimates <- svyby_repwts(
  rep_designs = list('original' = orig_rep_design,
                     'nonresponse-adjusted' = nr_adjusted_design),
  formula = ~ api00, FUN = svymean
)
print(overall_estimates, row.names = FALSE)
#>           Design_Name    api00       se
#>  nonresponse-adjusted 646.4465 30.66081
#>              original 644.1694 26.32936

# Estimate domain means (and their standard errors) from each design
domain_estimates <- svyby_repwts(
  rep_designs = list('original' = orig_rep_design,
                     'nonresponse-adjusted' = nr_adjusted_design),
  formula = ~ api00, by = ~ stype, FUN = svymean
)
print(domain_estimates, row.names = FALSE)
#>           Design_Name stype    api00       se
#>  nonresponse-adjusted     E 641.9463 34.19443
#>              original     E 648.8681 25.37430
#>  nonresponse-adjusted     H 699.5455 14.24657
#>              original     H 618.5714 46.34412
#>  nonresponse-adjusted     M 643.3429 41.47212
#>              original     M 631.4400 33.68762
```

We can even test for differences in estimates from the two sets of
weights and calculate confidence intervals for their difference.

``` r
estimates <- svyby_repwts(
  rep_designs = list('original' = orig_rep_design,
                     'nonresponse-adjusted' = nr_adjusted_design),
  formula = ~ api00, FUN = svymean
)

vcov(estimates)
#>                      nonresponse-adjusted original
#> nonresponse-adjusted             940.0856 784.1324
#> original                         784.1324 693.2352

diff_between_ests <- svycontrast(stat = estimates,
                                 contrasts = list(
                                   "Original vs. Adjusted" = c(-1,1)
                                 ))
print(diff_between_ests)
#>                       contrast     SE
#> Original vs. Adjusted  -2.2771 8.0657
confint(diff_between_ests)
#>                           2.5 %   97.5 %
#> Original vs. Adjusted -18.08562 13.53147
```

### Diagnosing potential issues with weights

When adjusting replicate weights, there are several diagnostics which
can be used to ensure that the adjustments were carried out correctly
and that they do more good than harm. The function
`summarize_rep_weights()` helps by allowing you to quickly summarize the
replicate weights.

For example, we can check to see that post-stratification correctly
removes variation in the column sums of the replicate weights.

``` r
# Post-stratify the design
post_stratified_design <- postStratify(design = orig_rep_design,
                                       strata = ~ stype + awards,
                                       population = xtabs(~stype + awards,
                                                          data = apipop))
# Compare overall summaries of the weights
list(
  'original' = summarize_rep_weights(orig_rep_design,
                                     type = 'overall'),
  'post-stratified' = summarize_rep_weights(post_stratified_design,
                                            type = 'overall')
)
#> $original
#>   nrows ncols degf_svy_pkg rank avg_wgt_sum sd_wgt_sums min_rep_wgt max_rep_wgt
#> 1   183    15           14   15        6194    403.1741           0    36.26464
#> 
#> $`post-stratified`
#>   nrows ncols degf_svy_pkg rank avg_wgt_sum  sd_wgt_sums min_rep_wgt
#> 1   183    15           14   15        6194 6.431099e-13           0
#>   max_rep_wgt
#> 1        93.4
```

When carrying out nonresponse adjustments, we might want to make sure
that column sums are the same before and after the adjustments.

``` r
# Summarize each column of replicate weights,
# before and after non-response adjustments
orig_rep_col_summaries <- summarize_rep_weights(orig_rep_design, type = 'specific')
adj_rep_col_summaries <- summarize_rep_weights(nr_adjusted_design, type = 'specific')

head(adj_rep_col_summaries)
#>   Rep_Column   N N_NONZERO      SUM     MEAN       CV MIN      MAX
#> 1          1 183        96 6237.518 34.08480 1.086404   0 179.4969
#> 2          2 183        97 6491.370 35.47197 1.066158   0 158.8797
#> 3          3 183        97 6563.900 35.86830 1.086089   0 181.0731
#> 4          4 183        94 6164.989 33.68846 1.113185   0 177.4097
#> 5          5 183        98 6563.900 35.86830 1.074456   0 181.0731
#> 6          6 183        97 6491.370 35.47197 1.080886   0 180.3158

# Compare the sums and number of nonzero entries
# before and after adjustment
weight_summaries  <- rbind(cbind(orig_rep_col_summaries, Design = 'Original'),
                           cbind(adj_rep_col_summaries, Design = 'NR-adjusted'))

weight_summaries <- weight_summaries[,c("Rep_Column", "Design", "N_NONZERO", "SUM")]
weight_summaries <- weight_summaries[order(weight_summaries$Design, decreasing = TRUE),]
weight_summaries <- weight_summaries[order(weight_summaries$Rep_Column),]
rownames(weight_summaries) <- NULL

head(weight_summaries)
#>   Rep_Column      Design N_NONZERO      SUM
#> 1          1    Original       172 6237.518
#> 2          1 NR-adjusted        96 6237.518
#> 3          2    Original       179 6491.370
#> 4          2 NR-adjusted        97 6491.370
#> 5          3    Original       181 6563.900
#> 6          3 NR-adjusted        97 6563.900
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
lou_vax_survey <- svydesign(ids = ~ 1, weights = ~ SAMPLING_WEIGHT,
                            data = lou_vax_survey) |>
  as.svrepdesign(type = "boot", replicates = 100, mse = TRUE)

# Adjust for nonresponse
lou_vax_survey <- lou_vax_survey |>
  redistribute_weights(
    reduce_if = RESPONSE_STATUS == "Nonrespondent",
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
  data = lou_pums_microdata,
  variables = ~ UNIQUE_ID + AGE + SEX + RACE_ETHNICITY + EDUC_ATTAINMENT,
  weights = ~ PWGTP, repweights = "PWGTP\\d{1,2}",
  type = "successive-difference",
  mse = TRUE
)
```

We can see that the vaccination survey seems to underrepresent
individuals who identify as Black or as Hispanic or Latino.

``` r
# Compare demographic estimates from the two data sources
estimate_comparisons <- data.frame(
  'Vax_Survey' = svymean(x = ~ RACE_ETHNICITY, design = acs_benchmark_survey) |> coef(),
  'ACS_Benchmark' = svymean(x = ~ RACE_ETHNICITY, design = lou_vax_survey) |> coef()
)
rownames(estimate_comparisons) <- gsub(x = rownames(estimate_comparisons),
                                       "RACE_ETHNICITY", "")
print(estimate_comparisons)
#>                                                         Vax_Survey
#> Black or African American alone, not Hispanic or Latino 0.19949824
#> Hispanic or Latino                                      0.04525039
#> Other Race, not Hispanic or Latino                      0.04630955
#> White alone, not Hispanic or Latino                     0.70894182
#>                                                         ACS_Benchmark
#> Black or African American alone, not Hispanic or Latino    0.16932271
#> Hispanic or Latino                                         0.03386454
#> Other Race, not Hispanic or Latino                         0.05776892
#> White alone, not Hispanic or Latino                        0.73904382
```

There are two options for calibrating the sample to the estimate
controls. With the first approach, we supply point estimates and their
variance-covariance matrix to the function `calibrate_to_estimate()`.

``` r
# Estimate control totals and their variance-covariance matrix
control_totals <- svymean(x = ~ RACE_ETHNICITY + EDUC_ATTAINMENT,
                          design = acs_benchmark_survey)
point_estimates <- coef(control_totals)
vcov_estimates <- vcov(control_totals)

# Calibrate the vaccination survey to the estimated control totals
vax_survey_raked_to_estimates <- calibrate_to_estimate(
  rep_design = lou_vax_survey,
  estimate = point_estimates,
  vcov_estimate = vcov_estimates,
  cal_formula = ~ RACE_ETHNICITY + EDUC_ATTAINMENT,
  calfun = survey::cal.raking
)
```

With the second approach, we supply the control survey’s replicate
design to `calibrate_to_sample()`.

``` r
vax_survey_raked_to_acs_sample <- calibrate_to_sample(
  primary_rep_design = lou_vax_survey,
  control_rep_design = acs_benchmark_survey,
  cal_formula = ~ RACE_ETHNICITY + EDUC_ATTAINMENT,
  calfun = survey::cal.raking
)
```

After calibration, we can see that the estimated vaccination rate has
decreased, and the estimated standard error of the estimated vaccination
rate has increased.

``` r
# Compare the two sets of estimates
svyby_repwts(
  rep_design = list(
    'NR-adjusted' = lou_vax_survey,
    'Raked to estimate' = vax_survey_raked_to_estimates,
    'Raked to sample' = vax_survey_raked_to_acs_sample
  ),
  formula = ~ VAX_STATUS,
  FUN = svymean
)
#>                         Design_Name VAX_STATUSUnvaccinated VAX_STATUSVaccinated
#> NR-adjusted             NR-adjusted              0.4621514            0.5378486
#> Raked to estimate Raked to estimate              0.4732623            0.5267377
#> Raked to sample     Raked to sample              0.4732623            0.5267377
#>                          se1        se2
#> NR-adjusted       0.02430176 0.02430176
#> Raked to estimate 0.02448676 0.02448676
#> Raked to sample   0.02446881 0.02446881
```
