# Convert Survey Design to Fay's Generalized Replication Design

Converts a survey design object to a replicate design object with
replicate weights formed using the generalized replication method of Fay
(1989). The generalized replication method forms replicate weights from
a textbook variance estimator, provided that the variance estimator can
be represented as a quadratic form whose matrix is positive semidefinite
(this covers a large class of variance estimators).

## Usage

``` r
as_fays_gen_rep_design(
  design,
  variance_estimator = NULL,
  aux_var_names = NULL,
  max_replicates = Inf,
  balanced = TRUE,
  psd_option = "warn",
  mse = TRUE,
  compress = TRUE
)
```

## Arguments

- design:

  A survey design object created using the 'survey' (or 'srvyr')
  package, with class `'survey.design'` or `'svyimputationList'`.

- variance_estimator:

  The name of the variance estimator whose quadratic form matrix should
  be created. See
  [variance-estimators](https://bschneidr.github.io/svrep/reference/variance-estimators.md)
  for a detailed description of each variance estimator. Options
  include:

  - **"Yates-Grundy"**:  
    The Yates-Grundy variance estimator based on first-order and
    second-order inclusion probabilities.

  - **"Horvitz-Thompson"**:  
    The Horvitz-Thompson variance estimator based on first-order and
    second-order inclusion probabilities.

  - **"Poisson Horvitz-Thompson"**:  
    The Horvitz-Thompson variance estimator based on assuming Poisson
    sampling, with first-order inclusion probabilities inferred from the
    sampling probabilities of the survey design object.

  - **"Stratified Multistage SRS"**:  
    The usual stratified multistage variance estimator based on
    estimating the variance of cluster totals within strata at each
    stage.

  - **"Ultimate Cluster"**:  
    The usual variance estimator based on estimating the variance of
    first-stage cluster totals within first-stage strata.

  - **"Deville-1"**:  
    A variance estimator for unequal-probability sampling without
    replacement, described in Matei and Tillé (2005) as "Deville 1".

  - **"Deville-2"**:  
    A variance estimator for unequal-probability sampling without
    replacement, described in Matei and Tillé (2005) as "Deville 2".

  - **"Deville-Tille":**  
    A variance estimator useful for balanced sampling designs, proposed
    by Deville and Tillé (2005).

  - **"SD1"**:  
    The non-circular successive-differences variance estimator described
    by Ash (2014), sometimes used for variance estimation for systematic
    sampling.

  - **"SD2"**:  
    The circular successive-differences variance estimator described by
    Ash (2014). This estimator is the basis of the
    "successive-differences replication" estimator commonly used for
    variance estimation for systematic sampling.

  - **"BOSB"**:  
    The kernel-based variance estimator proposed by Breidt, Opsomer, and
    Sanchez-Borrego (2016) for use with systematic samples or other
    finely stratified designs. Uses the Epanechnikov kernel with the
    bandwidth automatically chosen to result in the smallest possible
    nonempty kernel window.

  - **"Beaumont-Emond"**:  
    The variance estimator of Beaumont and Emond (2022) for multistage
    unequal-probability sampling without replacement.

- aux_var_names:

  (Only used if `variance_estimator = "Deville-Tille")`. A vector of the
  names of auxiliary variables used in sampling.

- max_replicates:

  The maximum number of replicates to allow (should be as large as
  possible, given computer memory/storage limitations). A
  commonly-recommended default is 500. If the number of replicates
  needed for a balanced, fully-efficient estimator is less than
  `max_replicates`, then only the number of replicates needed will be
  created. If more replicates are needed than `max_replicates`, then the
  full number of replicates needed will be created, but only a random
  subsample will be retained.

- balanced:

  If `balanced=TRUE`, the replicates will all contribute equally to
  variance estimates, but the number of replicates needed may slightly
  increase.

- psd_option:

  Either `"warn"` (the default) or `"error"`. This option specifies what
  will happen if the target variance estimator has a quadratic form
  matrix which is not positive semidefinite. This can occasionally
  happen, particularly for two-phase designs.  
  If `psd_option="error"`, then an error message will be displayed.  
  If `psd_option="warn"`, then a warning message will be displayed, and
  the quadratic form matrix will be approximated by the most similar
  positive semidefinite matrix. This approximation was suggested by
  Beaumont and Patak (2012), who note that this is conservative in the
  sense of producing overestimates of variance. Beaumont and
  Patak (2012) argue that this overestimation is expected to be small in
  magnitude. See
  [`get_nearest_psd_matrix`](https://bschneidr.github.io/svrep/reference/get_nearest_psd_matrix.md)
  for details of the approximation.

- mse:

  If `TRUE` (the default), compute variances from sums of squares around
  the point estimate from the full-sample weights. If `FALSE`, compute
  variances from sums of squares around the mean estimate from the
  replicate weights. For Fay's generalized replication method, setting
  `mse = FALSE` can potentially lead to large underestimates of
  variance.

- compress:

  This reduces the computer memory required to represent the replicate
  weights and has no impact on estimates.

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

## Statistical Details

See Fay (1989) for a full description of this replication method, or see
the documentation in
[make_fays_gen_rep_factors](https://bschneidr.github.io/svrep/reference/make_fays_gen_rep_factors.md)
for implementation details.

See
[variance-estimators](https://bschneidr.github.io/svrep/reference/variance-estimators.md)
for a description of each variance estimator available for use with this
function.

Use
[`rescale_replicates`](https://bschneidr.github.io/svrep/reference/rescale_replicates.md)
to eliminate negative adjustment factors.

## Two-Phase Designs

For a two-phase design, `variance_estimator` should be a list of
variance estimators' names, with two elements, such as
`list('Ultimate Cluster', 'Poisson Horvitz-Thompson')`. In two-phase
designs, only the following estimators may be used for the second phase:

- "Ultimate Cluster"

- "Stratified Multistage SRS"

- "Poisson Horvitz-Thompson"

For statistical details on the handling of two-phase designs, see the
documentation for
[make_twophase_quad_form](https://bschneidr.github.io/svrep/reference/make_twophase_quad_form.md).

## References

The generalized replication method was first proposed in Fay (1984). Fay
(1989) refined the generalized replication method to produce "balanced"
replicates, in the sense that each replicate contributes equally to
variance estimates. The advantage of balanced replicates is that one can
still obtain a reasonable variance estimate by using only a random
subset of the replicates.

\- Ash, S. (2014). "*Using successive difference replication for
estimating variances*." **Survey Methodology**, Statistics Canada,
40(1), 47-59.  
  
- Beaumont, J.-F.; Émond, N. (2022). "*A Bootstrap Variance Estimation
Method for Multistage Sampling and Two-Phase Sampling When Poisson
Sampling Is Used at the Second Phase*." **Stats**, *5*: 339-357.
https://doi.org/10.3390/stats5020019  
  
- Breidt, F. J., Opsomer, J. D., & Sanchez-Borrego, I. (2016).
"*Nonparametric Variance Estimation Under Fine Stratification: An
Alternative to Collapsed Strata*." **Journal of the American Statistical
Association**, 111(514), 822-833.
https://doi.org/10.1080/01621459.2015.1058264  
  
- Deville, J. C., and Tillé, Y. (2005). "*Variance approximation under
balanced sampling.*" **Journal of Statistical Planning and Inference**,
128, 569-591.  
  
- Dippo, Cathryn, Robert Fay, and David Morganstein. 1984. "Computing
Variances from Complex Samples with Replicate Weights." In, 489-94.
Alexandria, VA: American Statistical Association.
http://www.asasrms.org/Proceedings/papers/1984_094.pdf.  
  
- Fay, Robert. 1984. "Some Properties of Estimates of Variance Based on
Replication Methods." In, 495-500. Alexandria, VA: American Statistical
Association. http://www.asasrms.org/Proceedings/papers/1984_095.pdf.  
  
- Fay, Robert. 1989. "Theory And Application Of Replicate Weighting For
Variance Calculations." In, 495-500. Alexandria, VA: American
Statistical Association.
http://www.asasrms.org/Proceedings/papers/1989_033.pdf  
  
- Matei, Alina, and Yves Tillé. (2005). "*Evaluation of Variance
Approximations and Estimators in Maximum Entropy Sampling with Unequal
Probability and Fixed Sample Size.*" **Journal of Official Statistics**,
21(4):543-70.

## See also

For greater customization of the method,
[`make_quad_form_matrix`](https://bschneidr.github.io/svrep/reference/make_quad_form_matrix.md)
can be used to represent several common variance estimators as a
quadratic form's matrix, which can then be used as an input to
[`make_fays_gen_rep_factors`](https://bschneidr.github.io/svrep/reference/make_fays_gen_rep_factors.md).

## Examples

``` r
# \donttest{


  ## Load an example systematic sample ----
  data('library_stsys_sample', package = 'svrep')

  ## First, ensure data are sorted in same order as was used in sampling
  library_stsys_sample <- library_stsys_sample |>
    sort_by(~ SAMPLING_SORT_ORDER)

  ## Create a survey design object
  design_obj <- svydesign(
    data = library_stsys_sample,
    strata = ~ SAMPLING_STRATUM,
    ids = ~ 1,
    fpc = ~ STRATUM_POP_SIZE
  )

  ## Convert to generalized replicate design

  gen_rep_design_sd2 <- as_fays_gen_rep_design(
    design = design_obj,
    variance_estimator = "SD2",
    max_replicates = 250,
    mse = TRUE
  )
#> For `variance_estimator='SD2', assumes rows of data are sorted in the same order used in sampling.

  svytotal(x = ~ TOTSTAFF, na.rm = TRUE, design = gen_rep_design_sd2)
#>           total    SE
#> TOTSTAFF 180977 38893
# }
```
