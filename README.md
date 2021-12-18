
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
#>   nrows ncols degf_svy_pkg rank min_rep_wgt max_rep_wgt    max_CV SD_wgt_sums
#> 1   183    15           14   15           0    36.26464 0.5047941    403.1741
#> 
#> $`post-stratified`
#>   nrows ncols degf_svy_pkg rank min_rep_wgt max_rep_wgt    max_CV  SD_wgt_sums
#> 1   183    15           14   15           0        93.4 0.5862556 6.431099e-13
```

When carrying out nonresponse adjustments, we might want to make sure
that column sums are the same before and after the adjustments.

``` r
# Summarize each column of replicate weights,
# before and after non-response adjustments
orig_rep_col_summaries <- summarize_rep_weights(orig_rep_design, type = 'specific')
adj_rep_col_summaries <- summarize_rep_weights(nr_adjusted_design, type = 'specific')

head(adj_rep_col_summaries)
#>   Rep_Column   N N_NONZERO      SUM     MEAN       CV        L MIN      MAX
#> 1          1 183        96 6237.518 34.08480 1.086404 1.180274   0 179.4969
#> 2          2 183        97 6491.370 35.47197 1.066158 1.136692   0 158.8797
#> 3          3 183        97 6563.900 35.86830 1.086089 1.179589   0 181.0731
#> 4          4 183        94 6164.989 33.68846 1.113185 1.239180   0 177.4097
#> 5          5 183        98 6563.900 35.86830 1.074456 1.154456   0 181.0731
#> 6          6 183        97 6491.370 35.47197 1.080886 1.168314   0 180.3158

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
