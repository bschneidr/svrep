# Get Variance Estimator's Quadratic Form Matrix

Common variance estimators for estimated population totals can be
represented as a quadratic form. Given a choice of variance estimator
and information about the sample design, this function constructs the
matrix of the quadratic form.  
  
In notation, let \\v(\hat{Y}) =
\mathbf{\breve{y}}^{\prime}\mathbf{\Sigma}\mathbf{\breve{y}}\\, where
\\\breve{y}\\ is the vector of weighted values, \\y_i/\pi_i, \space
i=1,\dots,n\\. This function constructs the \\n \times n\\ matrix of the
quadratic form, \\\mathbf{\Sigma}\\.

## Usage

``` r
make_quad_form_matrix(
  variance_estimator = "Yates-Grundy",
  probs = NULL,
  joint_probs = NULL,
  cluster_ids = NULL,
  strata_ids = NULL,
  strata_pop_sizes = NULL,
  sort_order = NULL,
  aux_vars = NULL
)
```

## Arguments

- variance_estimator:

  The name of the variance estimator whose quadratic form matrix should
  be created. See the section "Variance Estimators" below. Options
  include:

  - **"Yates-Grundy"**: The Yates-Grundy variance estimator based on
    first-order and second-order inclusion probabilities. If this is
    used, the argument `joint_probs` must also be used.

  - **"Horvitz-Thompson"**: The Horvitz-Thompson variance estimator
    based on first-order and second-order inclusion probabilities. If
    this is used, the argument `joint_probs` must also be used.

  - **"Stratified Multistage SRS"**: The usual stratified multistage
    variance estimator based on estimating the variance of cluster
    totals within strata at each stage. If this option is used, then it
    is necessary to also use the arguments `strata_ids`, `cluster_ids`,
    `strata_pop_sizes`, and `strata_pop_sizes`.

  - **"Ultimate Cluster"**: The usual variance estimator based on
    estimating the variance of first-stage cluster totals within
    first-stage strata. If this option is used, then it is necessary to
    also use the arguments `strata_ids`, `cluster_ids`,
    `strata_pop_sizes`. Optionally, to use finite population correction
    factors, one can also use the argument `strata_pop_sizes`.

  - **"Deville-1"**: A variance estimator for unequal-probability
    sampling without replacement, described in Matei and Tillé (2005) as
    "Deville 1". If this option is used, then it is necessary to also
    use the arguments `strata_ids`, `cluster_ids`, and `probs`.

  - **"Deville-2"**: A variance estimator for unequal-probability
    sampling without replacement, described in Matei and Tillé (2005) as
    "Deville 2". If this option is used, then it is necessary to also
    use the arguments `strata_ids`, `cluster_ids`, and `probs`.

  - **"SD1"**: The non-circular successive-differences variance
    estimator described by Ash (2014), sometimes used for variance
    estimation for systematic sampling.

  - **"SD2"**: The circular successive-differences variance estimator
    described by Ash (2014). This estimator is the basis of the
    "successive-differences replication" estimator commonly used for
    variance estimation for systematic sampling.

  - **"BOSB"**:  
    The kernel-based variance estimator proposed by Breidt, Opsomer, and
    Sanchez-Borrego (2016) for use with systematic samples or other
    finely stratified designs. Uses the Epanechnikov kernel with the
    bandwidth automatically chosen to result in the smallest possible
    nonempty kernel window.

  - **"Deville-Tille"**: The estimator of Deville and Tillé (2005),
    developed for balanced sampling using the cube method.

  - **"Beaumont-Emond"**: The variance estimator of Beaumont and
    Emond (2022) for multistage unequal-probability sampling without
    replacement.

- probs:

  Required if `variance_estimator` equals `"Deville-1"`, `"Deville-2"`,
  or `"Deville-Tille"`. This should be a matrix or data frame of
  sampling probabilities. If there are multiple stages of sampling, then
  `probs` can have multiple columns, with one column for each level of
  sampling to be accounted for by the variance estimator.

- joint_probs:

  Only used if `variance_estimator = "Horvitz-Thompson"` or
  `variance_estimator = "Yates-Grundy"`. This should be a matrix of
  joint inclusion probabilities. Element `[i,i]` of the matrix is the
  first-order inclusion probability of unit `i`, while element `[i,j]`
  is the joint inclusion probability of units `i` and `j`.

- cluster_ids:

  Required unless `variance_estimator` equals `"Horvitz-Thompson"` or
  `"Yates-Grundy"`. This should be a matrix or data frame of cluster
  IDs. If there are multiple stages of sampling, then `cluster_ids` can
  have multiple columns, with one column for each level of sampling to
  be accounted for by the variance estimator.

- strata_ids:

  Required if `variance_estimator` equals `"Stratified Multistage SRS"`
  or `"Ultimate Cluster"`. This should be a matrix or data frame of
  strata IDs. If there are multiple stages of sampling, then
  `strata_ids` can have multiple columns, with one column for each level
  of sampling to be accounted for by the variance estimator.

- strata_pop_sizes:

  Required if `variance_estimator` equals `"Stratified Multistage SRS"`,
  but can optionally be used if `variance_estimator` equals
  `"Ultimate Cluster"`, `"SD1"`, or `"SD2"`. If there are multiple
  stages of sampling, then `strata_pop_sizes` can have multiple columns,
  with one column for each level of sampling to be accounted for by the
  variance estimator.

- sort_order:

  Required if `variance_estimator` equals `"SD1"` or `"SD2"`. This
  should be a vector that orders the rows of data into the order used
  for sampling.

- aux_vars:

  Required if `variance_estimator` equals `"Deville-Tille"`. A matrix of
  auxiliary variables.

## Value

The matrix of the quadratic form representing the variance estimator.

## Variance Estimators

See
[variance-estimators](https://bschneidr.github.io/svrep/reference/variance-estimators.md)
for a description of each variance estimator.

## Arguments required for each variance estimator

Below are the arguments that are required or optional for each variance
estimator.

|                           |          |             |             |            |                  |            |          |
|---------------------------|----------|-------------|-------------|------------|------------------|------------|----------|
| variance_estimator        | probs    | joint_probs | cluster_ids | strata_ids | strata_pop_sizes | sort_order | aux_vars |
| Stratified Multistage SRS |          |             | Required    | Required   | Required         |            |          |
| Ultimate Cluster          |          |             | Required    | Required   | Optional         |            |          |
| SD1                       |          |             | Required    | Optional   | Optional         | Required   |          |
| SD2                       |          |             | Required    | Optional   | Optional         | Required   |          |
| Deville-1                 | Required |             | Required    | Optional   |                  |            |          |
| Deville-2                 | Required |             | Required    | Optional   |                  |            |          |
| Beaumont-Emond            | Required |             | Required    | Optional   |                  |            |          |
| Deville-Tille             | Required |             | Required    | Optional   |                  |            | Required |
| BOSB                      |          |             | Required    | Optional   |                  |            | Required |
| Yates-Grundy              |          | Required    |             |            |                  |            |          |
| Horvitz-Thompson          |          | Required    |             |            |                  |            |          |

## See also

See
[variance-estimators](https://bschneidr.github.io/svrep/reference/variance-estimators.md)
for a description of each variance estimator.

For a two-phase design, the function
[make_twophase_quad_form](https://bschneidr.github.io/svrep/reference/make_twophase_quad_form.md)
combines the quadratic form matrix from each phase.

## Examples

``` r
# \donttest{
# Example 1: The Horvitz-Thompson Estimator
  library(survey)
  data("election", package = "survey")

  ht_quad_form_matrix <- make_quad_form_matrix(variance_estimator = "Horvitz-Thompson",
                                               joint_probs = election_jointprob)
  ##_ Produce variance estimate
  wtd_y <- as.matrix(election_pps$wt * election_pps$Bush)
  t(wtd_y) %*% ht_quad_form_matrix %*% wtd_y
#> 1 x 1 Matrix of class "dgeMatrix"
#>              [,1]
#> [1,] 6.782923e+12

  ##_ Compare against result from 'survey' package
  svytotal(x = ~ Bush,
           design = svydesign(data=election_pps,
                              variance = "HT",
                              pps = ppsmat(election_jointprob),
                              ids = ~ 1, fpc = ~ p)) |> vcov()
#>              [,1]
#> [1,] 6.782923e+12

# Example 2: Stratified multistage Sample ----

  data("mu284", package = 'survey')
  multistage_srswor_design <- svydesign(data = mu284,
                                        ids = ~ id1 + id2,
                                        fpc = ~ n1 + n2)

  multistage_srs_quad_form <- make_quad_form_matrix(
    variance_estimator = "Stratified Multistage SRS",
    cluster_ids = mu284[,c('id1', 'id2')],
    strata_ids = matrix(1, nrow = nrow(mu284), ncol = 2),
    strata_pop_sizes = mu284[,c('n1', 'n2')]
  )

  wtd_y <- as.matrix(weights(multistage_srswor_design) * mu284$y1)
  t(wtd_y) %*% multistage_srs_quad_form %*% wtd_y
#> 1 x 1 Matrix of class "dgeMatrix"
#>         [,1]
#> [1,] 5172234

  ##_ Compare against result from 'survey' package
  svytotal(x = ~ y1, design = multistage_srswor_design) |> vcov()
#>         y1
#> y1 5172234

# Example 3: Successive-differences estimator ----

  data('library_stsys_sample', package = 'svrep')

  sd1_quad_form <- make_quad_form_matrix(
    variance_estimator = 'SD1',
    cluster_ids = library_stsys_sample[,'FSCSKEY',drop=FALSE],
    strata_ids = library_stsys_sample[,'SAMPLING_STRATUM',drop=FALSE],
    strata_pop_sizes = library_stsys_sample[,'STRATUM_POP_SIZE',drop=FALSE],
    sort_order = library_stsys_sample[['SAMPLING_SORT_ORDER']]
  )

  wtd_y <- as.matrix(library_stsys_sample[['TOTCIR']] /
                      library_stsys_sample$SAMPLING_PROB)
  wtd_y[is.na(wtd_y)] <- 0

  t(wtd_y) %*% sd1_quad_form %*% wtd_y
#> 1 x 1 Matrix of class "dgeMatrix"
#>              [,1]
#> [1,] 2.217205e+17

# Example 4: Deville estimators ----

 data('library_multistage_sample', package = 'svrep')

 deville_quad_form <- make_quad_form_matrix(
     variance_estimator = 'Deville-1',
     cluster_ids = library_multistage_sample[,c("PSU_ID", "SSU_ID")],
     strata_ids = cbind(rep(1, times = nrow(library_multistage_sample)),
                        library_multistage_sample$PSU_ID),
     probs = library_multistage_sample[,c("PSU_SAMPLING_PROB",
                                          "SSU_SAMPLING_PROB")]
 )
# }
```
