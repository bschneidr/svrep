# Convert Survey Design to Random Group Jackknife Replicate Design

Forms a specified number of jackknife replicates based on grouping
primary sampling units (PSUs) into random, (approximately) equal-sized
groups.

## Usage

``` r
as_random_group_jackknife_design(
  design,
  replicates = 50,
  var_strat = NULL,
  var_strat_frac = NULL,
  sort_var = NULL,
  adj_method = "variance-stratum-psus",
  scale_method = "variance-stratum-psus",
  group_var_name = ".random_group",
  compress = TRUE,
  mse = getOption("survey.replicates.mse")
)
```

## Arguments

- design:

  A survey design object created using the 'survey' (or 'srvyr')
  package, with class `'survey.design'` or `'svyimputationList'`.

- replicates:

  The number of replicates to create for each variance stratum. The
  total number of replicates created is the number of variance strata
  times `replicates`. Every design stratum must have at least as many
  primary sampling units (PSUs), as `replicates`.

- var_strat:

  Specifies the name of a variable in the data that defines variance
  strata to use for the grouped jackknife. If `var_strat = NULL`, then
  there is effectively only one variance stratum.

- var_strat_frac:

  Specifies the sampling fraction to use for finite population
  corrections in each value of `var_strat`. Can use either a single
  number or a variable in the data corresponding to `var_strat`.

- sort_var:

  (Optional) Specifies the name of a variable in the data which should
  be used to sort the data before assigning random groups. If a variable
  is specified for `var_strat`, the sorting will happen within values of
  that variable.

- adj_method:

  Specifies how to calculate the replicate weight adjustment factor.
  Available options for `adj_method` include:

  - `"variance-stratum-psus"` (the default)  
    The replicate weight adjustment for a unit is based on the number of
    PSUs in its variance stratum.

  - `"variance-units"`  
    The replicate weight adjustment for a unit is based on the number of
    variance units in its variance stratum.

  See the section "Adjustment and Scale Methods" for details.

- scale_method:

  Specifies how to calculate the scale factor for each replicate.
  Available options for `scale_method` include:

  - `"variance-stratum-psus"`  
    The scale factor for a variance unit is based on its number of PSUs
    compared to the number of PSUs in its variance stratum.

  - `"variance-units"`  
    The scale factor for a variance unit is based on the number of
    variance units in its variance stratum.

  See the section "Adjustment and Scale Methods" for details.

- group_var_name:

  (Optional) The name of a new variable created to save identifiers for
  which random group each PSU was grouped into for the purpose of
  forming replicates. Specify `group_var_name = NULL` to avoid creating
  the variable in the data.

- compress:

  Use a compressed representation of the replicate weights matrix. This
  reduces the computer memory required to represent the replicate
  weights and has no impact on estimates.

- mse:

  If `TRUE`, compute variances from sums of squares around the point
  estimate from the full-sample weights. If `FALSE`, compute variances
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

## Formation of Random Groups

Within each value of `VAR_STRAT`, the data are sorted by first-stage
sampling strata, and then the PSUs in each stratum are randomly
arranged. Groups are then formed by serially placing PSUs into each
group. The first PSU in the `VAR_STRAT` is placed into the first group,
the second PSU into the second group, and so on. Once a PSU has been
assigned to the last group, the process begins again by assigning the
next PSU to the first group, the PSU after that to the second group, and
so on.

The random group that each observation is assigned to can be saved as a
variable in the data by using the function argument `group_var_name`.

## Adjustment and Scale Methods

The jackknife replication variance estimator based on \\R\\ replicates
takes the following form: \$\$ v(\hat{\theta}) = \sum\_{r=1}^{R} (1 -
f_r) \times c_r \times \left(\hat{\theta}\_r - \hat{\theta}\right)^2
\$\$ where \\r\\ indexes one of the \\R\\ sets of replicate weights,
\\c_r\\ is a corresponding scale factor for the \\r\\-th replicate, and
\\1 - f_r\\ is an optional finite population correction factor that can
potentially differ across variance strata.

To form the replicate weights, the PSUs are divided into \\\tilde{H}\\
variance strata, and the \\\tilde{h}\\-th variance stratum contains
\\G\_{\tilde{h}}\\ random groups. The number of replicates \\R\\ equals
the total number of random groups across all variance strata: \\R =
\sum\_{\tilde{h}}^{\tilde{H}} G\_{\tilde{h}}\\. In other words, each
replicate corresponds to one of the random groups from one of the
variance strata.

The weights for replicate \\r\\ corresponding to random group \\g\\
within variance stratum \\\tilde{h}\\ is defined as follows.

If case \\i\\ is not in variance stratum \\\tilde{h}\\, then
\\w\_{i}^{(r)} = w_i\\.

If case \\i\\ is in variance stratum \\\tilde{h}\\ and not in random
group \\g\\, then \\w\_{i}^{(r)} = a\_{\tilde{h}g} w_i\\.

Otherwise, if case \\i\\ is in random group \\g\\ of variance stratum
\\\tilde{h}\\, then \\w\_{i}^{(r)} = 0\\.

The R function argument `adj_method` determines how the adjustment
factor \\a\_{\tilde{h} g}\\ is calculated. When
`adj_method = "variance-units"`, then \\a\_{\tilde{h} g}\\ is calculated
based on \\G\_{\tilde{h}}\\, which is the number of random groups in
variance stratum \\\tilde{h}\\. When
`adj_method = "variance-stratum-psus"`, then \\a\_{\tilde{h} g}\\ is
calculated based on \\n\_{\tilde{h}g}\\, which is the number of PSUs in
random group \\g\\ in variance stratum \\\tilde{h}\\, as well as
\\n\_{\tilde{h}}\\, the total number of PSUs in variance stratum
\\\tilde{h}\\.

If `adj_method = "variance-units"`, then: \$\$a\_{\tilde{h}g} =
\frac{G\_{\tilde{h}}}{G\_{\tilde{h}} - 1}\$\$

If `adj_method = "variance-stratum-psus"`, then: \$\$a\_{\tilde{h}g} =
\frac{n\_{\tilde{h}}}{n\_{\tilde{h}} - n\_{\tilde{h}g}}\$\$

The scale factor \\c_r\\ for replicate \\r\\ corresponding to random
group \\g\\ within variance stratum \\\tilde{h}\\ is calculated
according to the function argument `scale_method`.

If `scale_method = "variance-units"`, then: \$\$c_r =
\frac{G\_{\tilde{h}} - 1}{G\_{\tilde{h}}}\$\$

If `scale_method = "variance-stratum-psus"`, then: \$\$c_r =
\frac{n\_{\tilde{h}} - n\_{\tilde{h}g}}{n\_{\tilde{h}}}\$\$

The sampling fraction \\f_r\\ used for finite population correction
\\1 - f_r\\ is by default assumed to equal 0. However, the user can
supply a sampling fraction for each variance stratum using the argument
`var_strat_frac`.

When variance units in a variance stratum have differing numbers of
PSUs, the combination `adj_method = "variance-stratum-psus"` and
`scale_method = "variance-units"` is recommended by Valliant, Brick, and
Dever (2008), corresponding to their method `"GJ2"`.

The random-groups jackknife method often referred to as "DAGJK"
corresponds to the options `var_strat = NULL`,
`adj_method = "variance-units"`, and `scale_method = "variance-units"`.
The DAGJK method will yield upwardly-biased variance estimates for
totals if the total number of PSUs is not a multiple of the total number
of replicates (Valliant, Brick, and Dever 2008).

## References

See Section 15.5 of Valliant, Dever, and Kreuter (2018) for an
introduction to the grouped jackknife and guidelines for creating the
random groups.

\- Valliant, R., Dever, J., Kreuter, F. (2018). "Practical Tools for
Designing and Weighting Survey Samples, 2nd edition." New York:
Springer.

See Valliant, Brick, and Dever (2008) for statistical details related to
the `adj_method` and `scale_method` arguments.

\- Valliant, Richard, Michael Brick, and Jill Dever. 2008. "Weight
Adjustments for the Grouped Jackknife Variance Estimator." *Journal of
Official Statistics*. 24: 469-88.

See Chapter 4 of Wolter (2007) for additional details of the jackknife,
including the method based on random groups.

\- Wolter, Kirk. 2007. "Introduction to Variance Estimation." New York,
NY: Springer New York. https://doi.org/10.1007/978-0-387-35099-8.

## Examples

``` r
# Load example data

 data('api', package = 'survey')

 api_strat_design <- svydesign(
   data = apistrat,
   id = ~ 1,
   strata = ~stype,
   weights = ~pw
 )

# Create a random-groups jackknife design

 jk_design <- as_random_group_jackknife_design(
   api_strat_design,
   replicates = 15
 )
 print(jk_design)
#> Call: as_random_group_jackknife_design(api_strat_design, replicates = 15)
#> with 15 replicates.
```
