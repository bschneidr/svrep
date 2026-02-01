# Combine Replicate Designs by Stacking

Stack replicate designs: combine rows of data, rows of replicate
weights, and the respective full-sample weights. This can be useful when
comparing estimates before and after a set of adjustments made to the
weights. Another more delicate application is when combining sets of
replicate weights from multiple years of data for a survey, although
this must be done carefully based on guidance from a data provider.

## Usage

``` r
stack_replicate_designs(..., .id = "Design_Name")
```

## Arguments

- ...:

  Replicate-weights survey design objects to combine. These can be
  supplied in one of two ways.  

  - Option 1 - A series of design objects, for example
    `'adjusted' = adjusted_design, 'orig' = orig_design`.

  - Option 2 - A list object containing design objects, for example
    `list('nr' = nr_adjusted_design, 'ue' = ue_adjusted_design)`.

  All objects must have the same specifications for `type`, `rho`,
  `mse`, `scales`, and `rscales`.

- .id:

  A single character value, which becomes the name of a new column of
  identifiers created in the output data to link each row to the design
  from which it was taken.  
  The labels used for the identifiers are taken from named arguments.

## Value

A replicate-weights survey design object, with class `svyrep.design` and
`svyrep.stacked`. The resulting survey design object always has its
value of `combined.weights` set to `TRUE`.

## Examples

``` r
# Load example data, creating a replicate design object
suppressPackageStartupMessages(library(survey))
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

# Stack the three designs, using any of the following syntax options
stacked_design <- stack_replicate_designs(orig_rep_design, ue_adjusted_design, nr_adjusted_design,
                                          .id = "which_design")
stacked_design <- stack_replicate_designs('original' = orig_rep_design,
                                          'unknown eligibility adjusted' = ue_adjusted_design,
                                          'nonresponse adjusted' = nr_adjusted_design,
                                          .id = "which_design")
list_of_designs <- list('original' = orig_rep_design,
                        'unknown eligibility adjusted' = ue_adjusted_design,
                        'nonresponse adjusted' = nr_adjusted_design)
stacked_design <- stack_replicate_designs(list_of_designs, .id = "which_design")
```
