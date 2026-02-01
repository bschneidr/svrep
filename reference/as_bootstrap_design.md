# Convert Survey Design to Bootstrap Replicate Design

Converts a survey design object to a replicate design object with
replicate weights formed using a bootstrap method. Supports stratified,
cluster samples with one or more stages of sampling. At each stage of
sampling, either simple random sampling (with or without replacement) or
unequal probability sampling (with or without replacement) may be used.

## Usage

``` r
as_bootstrap_design(
  design,
  type = "Rao-Wu-Yue-Beaumont",
  replicates = 500,
  compress = TRUE,
  mse = getOption("survey.replicates.mse"),
  samp_method_by_stage = NULL
)
```

## Arguments

- design:

  A survey design object created using the 'survey' (or 'srvyr')
  package, with class `'survey.design'` or `'svyimputationList'`.

- type:

  The type of bootstrap to use, which should be chosen based on its
  applicability to the sampling method used for the survey. The
  available types are the following:  

  - **"Rao-Wu-Yue-Beaumont"** (the default):  
    The bootstrap method of Beaumont and Émond (2022), which is a
    generalization of the Rao-Wu-Yue bootstrap, and is applicable to a
    wide variety of designs, including single-stage and multistage
    stratified designs. The design may have different sampling methods
    used at different stages. Each stage of sampling may potentially be
    PPS (i.e., use unequal probabilities), with or without replacement,
    and may potentially use Poisson sampling.  
      
    For a stratum with a fixed sample size of \\n\\ sampling units,
    resampling in each replicate resamples \\(n-1)\\ sampling units with
    replacement.

  - **"Rao-Wu"**:  
    The basic Rao-Wu \\(n-1)\\ bootstrap method, which is only
    applicable to single-stage designs or multistage designs where the
    first-stage sampling fractions are small (and can thus be ignored).
    Accommodates stratified designs. All sampling within a stratum must
    be simple random sampling with or without replacement, although the
    first-stage sampling is effectively treated as sampling without
    replacement.

  - **"Antal-Tille"**:  
    The doubled half bootstrap method proposed by Antal and Tillé
    (2014), which is only applicable to single-stage designs or
    multistage designs where the first-stage sampling fractions are
    small (and can thus be ignored). Accommodates stratified designs.
    Sampling within each stratum can be simple random sampling or
    unequal probability sampling with or without replacement.

  - **"Preston"**:  
    Preston's multistage rescaled bootstrap, which is applicable to
    single-stage designs or multistage designs with arbitrary sampling
    fractions. Accommodates stratified designs. All sampling within a
    stratum must be simple random sampling with or without replacement.

  - **"Canty-Davison"**:  
    The Canty-Davison bootstrap, which is only applicable to
    single-stage designs, with arbitrary sampling fractions.
    Accommodates stratified designs. All sampling with a stratum must be
    simple random sampling with or without replacement.

- replicates:

  Number of bootstrap replicates (should be as large as possible, given
  computer memory/storage limitations). A commonly-recommended default
  is 500.

- compress:

  Use a compressed representation of the replicate weights matrix. This
  reduces the computer memory required to represent the replicate
  weights and has no impact on estimates.

- mse:

  If `TRUE`, compute variances from sums of squares around the point
  estimate from the full-sample weights. If `FALSE`, compute variances
  from sums of squares around the mean estimate from the replicate
  weights.

- samp_method_by_stage:

  (Optional). By default, this function will automatically determine the
  sampling method used at each stage. However, this argument can be used
  to ensure the correct sampling method is identified for each stage.  
  Accepts a vector with length equal to the number of stages of
  sampling. Each element should be one of the following:  

  - `"SRSWOR"` - Simple random sampling, without replacement

  - `"SRSWR"` - Simple random sampling, with replacement

  - `"PPSWOR"` - Unequal probabilities of selection, without replacement

  - `"PPSWR"` - Unequal probabilities of selection, with replacement

  - `"Poisson"` - Poisson sampling: each sampling unit is selected into
    the sample at most once, with potentially different probabilities of
    inclusion for each sampling unit.

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

## References

Antal, E. and Tillé, Y. (2014). "A new resampling method for sampling
designs without replacement: The doubled half bootstrap."
**Computational Statistics**, *29*(5), 1345-1363.
https://doi.org/10.1007/s00180-014-0495-0

Beaumont, J.-F.; Émond, N. (2022). "A Bootstrap Variance Estimation
Method for Multistage Sampling and Two-Phase Sampling When Poisson
Sampling Is Used at the Second Phase." **Stats**, *5*: 339-357.
https://doi.org/10.3390/stats5020019

Canty, A.J.; Davison, A.C. (1999). "Resampling-based variance estimation
for labour force surveys." **The Statistician**, *48*: 379-391.

Preston, J. (2009). "Rescaled bootstrap for stratified multistage
sampling." **Survey Methodology**, *35*(2): 227-234.

Rao, J.N.K.; Wu, C.F.J.; Yue, K. (1992). "Some recent work on resampling
methods for complex surveys." **Survey Methodology**, *18*: 209-217.

## See also

Use
[`estimate_boot_reps_for_target_cv`](https://bschneidr.github.io/svrep/reference/estimate_boot_reps_for_target_cv.md)
to help choose the number of bootstrap replicates.

The underlying function for the Rao-Wu-Yue-Beaumont bootstrap is
[`make_rwyb_bootstrap_weights`](https://bschneidr.github.io/svrep/reference/make_rwyb_bootstrap_weights.md),
and the underlying function for the Antal-Tillé bootstrap is
[`make_doubled_half_bootstrap_weights`](https://bschneidr.github.io/svrep/reference/make_doubled_half_bootstrap_weights.md).
Other bootstrap methods are implemented using functions from the
'survey' package, including:
[`bootweights`](https://rdrr.io/pkg/survey/man/bootweights.html)
(Canty-Davison),
[`subbootweights`](https://rdrr.io/pkg/survey/man/bootweights.html)
(Rao-Wu), and
[`mrbweights`](https://rdrr.io/pkg/survey/man/bootweights.html)
(Preston).

For systematic samples, one-PSU-per-stratum designs, or other especially
complex sample designs, one can use the generalized survey bootstrap
method. See
[`as_gen_boot_design`](https://bschneidr.github.io/svrep/reference/as_gen_boot_design.md)
or
[`make_gen_boot_factors`](https://bschneidr.github.io/svrep/reference/make_gen_boot_factors.md).

## Examples

``` r
# Example 1: A multistage sample with two stages of SRSWOR

  ## Load an example dataset from a multistage sample, with two stages of SRSWOR
  data("mu284", package = 'survey')
  multistage_srswor_design <- svydesign(data = mu284,
                                        ids = ~ id1 + id2,
                                        fpc = ~ n1 + n2)

  ## Convert the survey design object to a bootstrap design
  set.seed(2022)
  bootstrap_rep_design <- as_bootstrap_design(multistage_srswor_design,
                                              replicates = 500)

  ## Compare std. error estimates from bootstrap versus linearization
  data.frame(
    'Statistic' = c('total', 'mean', 'median'),
    'SE (bootstrap)' = c(SE(svytotal(x = ~ y1, design = bootstrap_rep_design)),
                         SE(svymean(x = ~ y1, design = bootstrap_rep_design)),
                         SE(svyquantile(x = ~ y1, quantile = 0.5,
                                        design = bootstrap_rep_design))),
    'SE (linearization)' = c(SE(svytotal(x = ~ y1, design = multistage_srswor_design)),
                             SE(svymean(x = ~ y1, design = multistage_srswor_design)),
                             SE(svyquantile(x = ~ y1, quantile = 0.5,
                                            design = multistage_srswor_design))),
    check.names = FALSE
  )
#>   Statistic SE (bootstrap) SE (linearization)
#> 1     total    2311.130145        2274.254701
#> 2      mean       2.449955           2.273653
#> 3    median       2.331234           2.521210

# Example 2: A multistage-sample,
# first stage selected with unequal probabilities without replacement
# second stage selected with simple random sampling without replacement

  data("library_multistage_sample", package = "svrep")

  multistage_pps <- svydesign(data = library_multistage_sample,
                              ids = ~ PSU_ID + SSU_ID,
                              fpc = ~ PSU_SAMPLING_PROB + SSU_SAMPLING_PROB,
                              pps = "brewer")

  bootstrap_rep_design <- as_bootstrap_design(
    multistage_pps, replicates = 500,
    samp_method_by_stage = c("PPSWOR", "SRSWOR")
  )

  ## Compare std. error estimates from bootstrap versus linearization
  data.frame(
      'Statistic' = c('total', 'mean'),
      'SE (bootstrap)' = c(
          SE(svytotal(x = ~ TOTCIR, na.rm = TRUE,
                      design = bootstrap_rep_design)),
          SE(svymean(x = ~ TOTCIR, na.rm = TRUE,
                     design = bootstrap_rep_design))),
      'SE (linearization)' = c(
          SE(svytotal(x = ~ TOTCIR, na.rm = TRUE,
                      design = multistage_pps)),
          SE(svymean(x = ~ TOTCIR, na.rm = TRUE,
                     design = multistage_pps))),
      check.names = FALSE
  )
#>   Statistic SE (bootstrap) SE (linearization)
#> 1     total   266151536.55       255100437.38
#> 2      mean       45762.71           42544.16
```
