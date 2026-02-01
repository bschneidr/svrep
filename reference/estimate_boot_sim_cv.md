# Estimate Bootstrap Simulation Error

Estimates the bootstrap simulation error, expressed as a "simulation
coefficient of variation" (CV).

## Usage

``` r
estimate_boot_sim_cv(svrepstat)
```

## Arguments

- svrepstat:

  An estimate obtained from a bootstrap replicate survey design object,
  with a function such as `svymean(..., return.replicates = TRUE)` or
  `withReplicates(..., return.replicates = TRUE)`.

## Value

An object of class `sim_cv_estimate`, with methods
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html). Call
[`summary()`](https://rdrr.io/r/base/summary.html) to get a data frame
with one row for each statistic. The column `STATISTIC` gives the name
of the statistic. The column `SIMULATION_CV` gives the estimated
simulation CV of the statistic. The column `N_REPLICATES` gives the
number of bootstrap replicates.

This object has no plot method. To produce plots, use the
closely-related function
[`estimate_boot_reps_for_target_cv`](https://bschneidr.github.io/svrep/reference/estimate_boot_reps_for_target_cv.md).

## Statistical Details

Unlike other replication methods such as the jackknife or balanced
repeated replication, the bootstrap variance estimator's precision can
always be improved by using a larger number of replicates, as the use of
only a finite number of bootstrap replicates introduces simulation error
to the variance estimation process. Simulation error can be measured as
a "simulation coefficient of variation" (CV), which is the ratio of the
standard error of a bootstrap estimator to the expectation of that
bootstrap estimator, where the expectation and standard error are
evaluated with respect to the bootstrapping process given the selected
sample.

For a statistic \\\hat{\theta}\\, the simulation CV of the bootstrap
variance estimator \\v\_{B}(\hat{\theta})\\ based on \\B\\ replicate
estimates \\\hat{\theta}^{\star}\_1,\dots,\hat{\theta}^{\star}\_B\\ is
defined as follows: \$\$ CV\_{\star}(v\_{B}(\hat{\theta})) =
\frac{\sqrt{var\_{\star}(v_B(\hat{\theta}))}}{E\_{\star}(v_B(\hat{\theta}))}
= \frac{CV\_{\star}(E_2)}{\sqrt{B}} \$\$ where \$\$ E_2 =
(\hat{\theta}^{\star} - \hat{\theta})^2 \$\$ \$\$ CV\_{\star}(E_2) =
\frac{\sqrt{var\_{\star}(E_2)}}{E\_{\star}(E_2)} \$\$ and
\\var\_{\star}\\ and \\E\_{\star}\\ are evaluated with respect to the
bootstrapping process, given the selected sample.  
  
The simulation CV, denoted \\CV\_{\star}(v\_{B}(\hat{\theta}))\\, is
estimated for a given number of replicates \\B\\ by estimating
\\CV\_{\star}(E_2)\\ using observed values and dividing this by
\\\sqrt{B}\\. If the bootstrap errors are assumed to be normally
distributed, then \\CV\_{\star}(E_2)=\sqrt{2}\\ and so
\\CV\_{\star}(v\_{B}(\hat{\theta}))\\ would not need to be estimated.
Using observed replicate estimates to estimate the simulation CV instead
of assuming normality allows simulation CV to be used for a a wide array
of bootstrap methods.

## References

See Section 3.3 and Section 8 of Beaumont and Patak (2012) for details
and an example where the simulation CV is used to determine the number
of bootstrap replicates needed for various alternative bootstrap methods
in an empirical illustration.

Beaumont, J.-F. and Z. Patak. (2012), "On the Generalized Bootstrap for
Sample Surveys with Special Attention to Poisson Sampling."
**International Statistical Review**, *80*: 127-148.
[doi:10.1111/j.1751-5823.2011.00166.x](https://doi.org/10.1111/j.1751-5823.2011.00166.x)
.

## See also

Use
[`estimate_boot_reps_for_target_cv`](https://bschneidr.github.io/svrep/reference/estimate_boot_reps_for_target_cv.md)
to help choose the number of bootstrap replicates.

## Examples

``` r
# \donttest{
set.seed(2022)

# Create an example bootstrap survey design object ----
data('api', package = 'survey')

boot_design <- svydesign(id=~1,strata=~stype, weights=~pw,
                         data=apistrat, fpc=~fpc) |>
  as_bootstrap_design(replicates = 5000)

# Calculate estimates of interest and retain estimates from each replicate ----

estimated_means_and_proportions <- svymean(
  x = ~ api00 + api99 + stype, 
  design = boot_design,
  return.replicates = TRUE
)

custom_statistic <- withReplicates(
  design = boot_design,
  return.replicates = TRUE,
  theta = function(wts, data) {
     numerator   <- sum(data$api00 * wts)
     denominator <- sum(data$api99 * wts)
     statistic   <- numerator/denominator
     return(statistic)
  }
)
# Estimate simulation CV of bootstrap estimates ----

  estimate_boot_sim_cv(
    svrepstat = estimated_means_and_proportions
  )
#>   STATISTIC SIMULATION_CV N_REPLICATES
#> 1     api00    0.01153215         5000
#> 2     api99    0.01153118         5000
#> 3    stypeE    0.01735963         5000
#> 4    stypeH    0.01735958         5000
#> 5    stypeM    0.01735958         5000

  estimate_boot_sim_cv(
    svrepstat = custom_statistic
  )
#> Warning: The elements of `svrepstat` are unnamed. Placeholder names (STATISTIC_1, etc.) will be used instead.
#>     STATISTIC SIMULATION_CV N_REPLICATES
#> 1 STATISTIC_1    0.01988987         5000
# }
```
