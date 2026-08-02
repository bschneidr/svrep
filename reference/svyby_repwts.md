# Survey Statistics for Multiple Replicate Design Objects

A modified version of the
[`svyby()`](https://rdrr.io/pkg/survey/man/svyby.html) function from the
`survey` package. Whereas
[`svyby()`](https://rdrr.io/pkg/survey/man/svyby.html) calculates
statistics separately for each subset formed by a specified grouping
variable, `svyby_repwts()` calculates statistics separately for each
replicate design, in addition to any additional user-specified grouping
variables.

## Usage

``` r
svyby_repwts(
  rep_designs,
  formula,
  by,
  FUN,
  ...,
  deff = FALSE,
  keep.var = TRUE,
  keep.names = TRUE,
  verbose = FALSE,
  vartype = c("se", "ci", "ci", "cv", "cvpct", "var"),
  drop.empty.groups = TRUE,
  return.replicates = FALSE,
  na.rm.by = FALSE,
  na.rm.all = FALSE,
  multicore = getOption("survey.multicore")
)
```

## Arguments

- rep_designs:

  The replicate-weights survey designs to be compared. Supplied either
  as:

  - A named list of replicate-weights survey design objects, for example
    `list('nr' = nr_adjusted_design, 'ue' = ue_adjusted_design)`.

  - A 'stacked' replicate-weights survey design object created by
    [`stack_replicate_designs()`](https://bschneidr.github.io/svrep/reference/stack_replicate_designs.md).

  The designs must all have the same number of columns of replicate
  weights, of the same type (bootstrap, JKn, etc.)

- formula:

  A formula specifying the variables to pass to `FUN`

- by:

  A formula specifying factors that define subsets

- FUN:

  A function taking a formula and survey design object as its first two
  arguments. Usually a function from the `survey` package, such as
  `svytotal` or `svymean`.

- ...:

  Other arguments to `FUN`

- deff:

  A value of `TRUE` or `FALSE`, indicating whether design effects should
  be estimated if possible.

- keep.var:

  A value of `TRUE` or `FALSE`. If `FUN` returns a `svystat` object,
  indicates whether to extract standard errors from it.

- keep.names:

  Define row names based on the subsets

- verbose:

  If `TRUE`, print a label for each subset as it is processed.

- vartype:

  Report variability as one or more of standard error, confidence
  interval, coefficient of variation, percent coefficient of variation,
  or variance

- drop.empty.groups:

  If `FALSE`, report `NA` for empty groups, if `TRUE` drop them from the
  output

- return.replicates:

  If `TRUE`, return all the replicates as the "replicates" attribute of
  the result. This can be useful if you want to produce custom summaries
  of the estimates from each replicate.

- na.rm.by:

  If true, omit groups defined by `NA` values of the `by` variables

- na.rm.all:

  If true, check for groups with no non-missing observations for
  variables defined by `formula` and treat these groups as empty

- multicore:

  Use `multicore` package to distribute subsets over multiple
  processors?

## Value

An object of class `"svyby"`: a data frame showing the grouping factors
and results of `FUN` for each combination of the grouping factors. The
first grouping factor always consists of indicators for which replicate
design was used for an estimate.

## Examples

``` r
# \donttest{
data(api)

dclus1 <- svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
dclus1$variables$response_status <- sample(x = c("Respondent", "Nonrespondent",
                                                 "Ineligible", "Unknown eligibility"),
                                           size = nrow(dclus1),
                                           replace = TRUE)
orig_rep_design <- as.svrepdesign(dclus1)

# Adjust weights for cases with unknown eligibility
ue_adjusted_design <- redistribute_weights(
    design = orig_rep_design,
    reduce_if = response_status %in% c("Unknown eligibility"),
    increase_if = !response_status %in% c("Unknown eligibility"),
    by = c("stype")
)

# Adjust weights for nonresponse
nr_adjusted_design <- redistribute_weights(
    design = ue_adjusted_design,
    reduce_if = response_status %in% c("Nonrespondent"),
    increase_if = response_status == "Respondent",
    by = c("stype")
)

# Compare estimates from the three sets of replicate weights

  list_of_designs <- list('original' = orig_rep_design,
                          'unknown eligibility adjusted' = ue_adjusted_design,
                          'nonresponse adjusted' = nr_adjusted_design)

  ##_ First compare overall means for two variables
  means_by_design <- svyby_repwts(formula = ~ api00 + api99,
                                  FUN = svymean,
                                  rep_design = list_of_designs)

  print(means_by_design)
#>                                               Design_Name    api00    api99
#> nonresponse adjusted                 nonresponse adjusted 647.1131 608.7525
#> original                                         original 644.1694 606.9781
#> unknown eligibility adjusted unknown eligibility adjusted 639.1914 601.6096
#>                                   se1      se2
#> nonresponse adjusted         28.73366 28.65897
#> original                     26.32936 26.99854
#> unknown eligibility adjusted 29.30072 29.49201

  ##_ Next compare domain means for two variables
  domain_means_by_design <- svyby_repwts(formula = ~ api00 + api99,
                                         by = ~ stype,
                                         FUN = svymean,
                                         rep_design = list_of_designs)

  print(domain_means_by_design)
#>                                                 Design_Name stype    api00
#> nonresponse adjusted.E                 nonresponse adjusted     E 652.4038
#> original.E                                         original     E 648.8681
#> unknown eligibility adjusted.E unknown eligibility adjusted     E 642.8462
#> nonresponse adjusted.H                 nonresponse adjusted     H 592.7500
#> original.H                                         original     H 618.5714
#> unknown eligibility adjusted.H unknown eligibility adjusted     H 603.8889
#> nonresponse adjusted.M                 nonresponse adjusted     M 647.0818
#> original.M                                         original     M 631.4400
#> unknown eligibility adjusted.M unknown eligibility adjusted     M 637.9091
#>                                   api99      se1      se2
#> nonresponse adjusted.E         609.1058 27.77889 27.59123
#> original.E                     607.7917 25.37430 25.83542
#> unknown eligibility adjusted.E 600.5962 27.98111 27.85828
#> nonresponse adjusted.H         570.0833 57.36322 59.46811
#> original.H                     595.7143 46.34412 50.75106
#> unknown eligibility adjusted.H 583.3333 48.59150 51.31940
#> nonresponse adjusted.M         628.3727 42.78019 39.90191
#> original.M                     608.6000 33.68762 34.82521
#> unknown eligibility adjusted.M 617.6818 38.62199 39.56364

# Calculate confidence interval for difference between estimates

ests_by_design <- svyby_repwts(rep_designs = list('NR-adjusted' = nr_adjusted_design,
                                                  'Original' = orig_rep_design),
                               FUN = svymean, formula = ~ api00 + api99)

differences_in_estimates <- svycontrast(stat = ests_by_design, contrasts = list(
  'Mean of api00: NR-adjusted vs. Original' = c(1,-1,0,0),
  'Mean of api99: NR-adjusted vs. Original' = c(0,0,1,-1)
))

print(differences_in_estimates)
#>                                         contrast     SE
#> Mean of api00: NR-adjusted vs. Original   2.9437 5.5208
#> Mean of api99: NR-adjusted vs. Original   1.7744 5.5037

confint(differences_in_estimates, level = 0.95)
#>                                             2.5 %   97.5 %
#> Mean of api00: NR-adjusted vs. Original -7.876906 13.76433
#> Mean of api99: NR-adjusted vs. Original -9.012608 12.56141
# }
```
