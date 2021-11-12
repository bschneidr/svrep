
<!-- README.md is generated from README.Rmd. Please edit that file -->

# svrep

<!-- badges: start -->

[![R-CMD-check](https://github.com/bschneidr/svrep/workflows/R-CMD-check/badge.svg)](https://github.com/bschneidr/svrep/actions)
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

### Compare estimates from different sets of weights

In order to assess whether weighting adjustments have an impact on the
estimates we care about, we want to compare the estimates from the
different sets of weights. The function `svyby_repwts()` makes it easy
to compare estimates from different sets of weights.

``` r
# Estimate overall means (and their standard errors) from each design
svyby_repwts(
  rep_designs = list('original' = orig_rep_design,
                     'nonresponse-adjusted' = nr_adjusted_design),
  formula = ~ api00, FUN = svymean
)
#>                               Design_Name    api00       se
#> nonresponse-adjusted nonresponse-adjusted 646.4465 29.26831
#> original                         original 644.1694 26.32936

# Estimate domain means (and their standard errors) from each design
svyby_repwts(
  rep_designs = list('original' = orig_rep_design,
                     'nonresponse-adjusted' = nr_adjusted_design),
  formula = ~ api00, by = ~ stype, FUN = svymean
)
#>                                 Design_Name stype    api00       se
#> nonresponse-adjusted.E nonresponse-adjusted     E 641.9463 35.35273
#> original.E                         original     E 648.8681 25.37430
#> nonresponse-adjusted.H nonresponse-adjusted     H 699.5455 12.77214
#> original.H                         original     H 618.5714 46.34412
#> nonresponse-adjusted.M nonresponse-adjusted     M 643.3429 42.88753
#> original.M                         original     M 631.4400 33.68762
```

We can even test for differences in estimates from the two sets of
weights and calculate confidence intervals for their difference.

``` r
estimates <- svyby_repwts(
  rep_designs = list('original' = orig_rep_design,
                     'nonresponse-adjusted' = nr_adjusted_design),
  formula = ~ api00, FUN = svymean
)

diff_between_ests <- svycontrast(stat = estimates,
                                 contrasts = list(
                                   "Original vs. Adjusted" = c(-1,1)
                                 ))
print(diff_between_ests)
#>                       contrast     SE
#> Original vs. Adjusted  -2.2771 8.7637
confint(diff_between_ests)
#>                           2.5 %   97.5 %
#> Original vs. Adjusted -19.45367 14.89952
```
