# Convert Survey Design to Successive Differences Replicate Design

Converts a survey design object to a replicate design object with
replicate weights formed using the successive differences replication
(SDR) method. The SDR method is suitable for designs that use systematic
sampling or finely-stratified sampling designs.

## Usage

``` r
as_sdr_design(
  design,
  replicates,
  sort_variable = NULL,
  use_normal_hadamard = FALSE,
  compress = TRUE,
  mse = TRUE
)
```

## Arguments

- design:

  A survey design object created using the 'survey' (or 'srvyr')
  package, with class `'survey.design'` or `'svyimputationList'`.

- replicates:

  The target number of replicates to create. This will determine the
  order of the Hadamard matrix to use when creating replicate factors.
  If `use_normal_hadamard = TRUE`, then the actual number of replicates
  will be greater than or equal to `replicates` and determined by
  identifying the smallest available Hadamard matrix available from the
  'survey' package. If `use_normal_hadamard = FALSE`, then the actual
  number of replicates will be the smallest *power* of 4 that is greater
  or equal to the specified value of `replicates`.

- sort_variable:

  A character string specifying the name of a sorting variable. This
  variable should give the sort order used in sampling. If the design
  includes strata, then the replicate factors will be assigned after
  first sorting by the first-stage strata identifier and then sorting by
  the value of `sort_variable` within each stratum.

- use_normal_hadamard:

  Whether to use a normal Hadamard matrix: that is, a matrix whose first
  row and first column only have entries equal to 1. This means that one
  of the replicates will be an "inactive" replicate. See the "Details"
  section for more information.

- compress:

  Use a compressed representation of the replicate weights matrix. This
  reduces the computer memory required to represent the replicate
  weights and has no impact on estimates.

- mse:

  If `TRUE`, compute variances from sums of squares around the point
  estimate from the full-sample weights, If `FALSE`, compute variances
  from sums of squares around the mean estimate from the replicate
  weights.

## Value

A replicate design object, with class `svyrep.design`, which can be used
with the usual functions, such as
[`svymean()`](https://rdrr.io/pkg/survey/man/surveysummary.html) or
[`svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html).

Use `weights(..., type = 'analysis')` to extract the matrix of replicate
weights.  
Use
[`as_data_frame_with_weights()`](https://bschneidr.github.io/svrep/reference/as_data_frame_with_weights.md)
to convert the design object to a data frame with columns for the
full-sample and replicate weights.

## Statistical Overview

The successive difference replication method was proposed by Fay and
Train (1995) as a replication method appropriate for samples selected
using systematic sampling. It is designed to yield variance estimates
for totals that are equivalent to successive difference variance
estimators described in Fay and Train (1995). There are different
methods for forming the replicate factors depending on whether the
replicate variance estimator is meant to be equivalent to the SD2
variance estimator (i.e., the circular successive difference estimator)
or the SD1 variance estimator (the non-circular successive difference
estimator) described in Ash (2014). This function uses the approach
based on the SD2 variance estimator. For multistage designs, this
replication method only takes into account information about the first
stage of sampling.

The scale factor to be used for variance estimation with the replicate
weights is \\4/R\\, where \\R\\ is the number of replicates. This scale
factor will be used even when there are finite population corrections;
see the subsection below.

As an alternative to the successive difference replication estimator,
one can use a generalized replication method where the target variance
estimator is the "SD1" or "SD2" estimator. See the functions
[as_gen_boot_design](https://bschneidr.github.io/svrep/reference/as_gen_boot_design.md)
or
[as_fays_gen_rep_design](https://bschneidr.github.io/svrep/reference/as_fays_gen_rep_design.md)
for more details on generalized replication and see the help section
[variance-estimators](https://bschneidr.github.io/svrep/reference/variance-estimators.md)
for more details on the "SD1" and "SD2" variance estimators.

## Details on Stratification and Finite Population Corrections

If the design includes strata, then the replicate factors will be
assigned after first sorting by the first-stage strata identifier and
then sorting by the value of `sort_variable` within each stratum.

If there are finite population correction factors, then these finite
population correction factors will be applied to the replicate factors.
This means that variance estimates with the finite population correction
do not require any adjustment to the overall scale factor used in
variance estimation. This is the approach used by the U.S. Census Bureau
for the 5-year American Community Survey (ACS) replicate weights (U.S.
Census Bureau, 2022, p. 12-8). This approach is used regardless of
whether the design has one overall finite population correction factor
or has different finite population correction factors for different
strata.

## Details on Row Assignments for Creating Replicate Factors

The number of replicates must match the order of an available Hadamard
matrix. A Hadamard matrix can either be normal or non-normal: a normal
Hadamard matrix is one where the entries in the first row and in the
first column are all equal to one. If the user specifies
`use_normal_hadamard = TRUE`, then there are more choices of Hadamard
matrix sizes available, and so greater flexibility in choosing the
number of replicates to create. When a normal Hadamard matrix is used,
this will result in the creation of an inactive replicate (sometimes
referred to as a "dead" replicate), which is a replicate where all the
replicate factors equal one. Inactive replicates are perfectly valid for
variance estimation, though some users may find them confusing.

An important part of the process of creating replicate weights is the
assignment of rows of the Hadamard matrix to primary sampling units. The
method of Ash (2014) referred to as "RA1" is used for row assignments,
which means that the replication-based variance estimates for totals
will be equivalent to the SD2 variance estimator described by Ash
(2014). The number of cycles used with the "RA1" method is the smallest
integer greater than \\n/R\\, where \\n\\ is the number of primary
sample units and \\R\\ is the number of replicates.

## References

Ash, S. (2014). "*Using successive difference replication for estimating
variances*." **Survey Methodology**, Statistics Canada, 40(1), 47-59.

Fay, R.E. and Train, G.F. (1995). "*Aspects of Survey and Model-Based
Postcensal Estimation of Income and Poverty Characteristics for States
and Counties*." Joint Statistical Meetings, Proceedings of the Section
on Government Statistics, 154-159.

U.S. Census Bureau. (2022). "*American Community Survey and Puerto Rico
Community Survey Design and Methodology, Version 3.0.*"

## Examples

``` r
library(survey)

# Load example stratified systematic sample
data('library_stsys_sample', package = 'svrep')

## First, ensure data are sorted in same order as was used in sampling
library_stsys_sample <- library_stsys_sample |>
  sort_by(~ SAMPLING_SORT_ORDER)

## Create a survey design object
design_obj <- svydesign(
  data   = library_stsys_sample,
  strata = ~ SAMPLING_STRATUM,
  ids    = ~ 1,
  fpc    = ~ STRATUM_POP_SIZE
)

## Convert to SDR replicate design
sdr_design <- as_sdr_design(
  design              = design_obj,
  replicates          = 180,
  sort_variable       = "SAMPLING_SORT_ORDER",
  use_normal_hadamard = TRUE
)
#> Using Hadamard matrix of order 180. If `use_normal_hadamard=TRUE`, the smallest possible order is 180.
#> Finite population corrections are incorporated into the replicate factors.

## Compare to generalized bootstrap
## based on the SD2 estimator that SDR approximates
gen_boot_design <- as_gen_boot_design(
  design             = design_obj,
  variance_estimator = "SD2",
  replicates         = 180,
  exact_vcov         = TRUE
)
#> For `variance_estimator='SD2', assumes rows of data are sorted in the same order used in sampling.

## Estimate sampling variances
svytotal(x = ~ TOTSTAFF, na.rm = TRUE, design = sdr_design)
#>           total    SE
#> TOTSTAFF 180977 39077
svytotal(x = ~ TOTSTAFF, na.rm = TRUE, design = gen_boot_design)
#>           total    SE
#> TOTSTAFF 180977 38893
```
