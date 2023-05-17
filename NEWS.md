# svrep 0.5.1

* New Features:
  * Added argument `exact_vcov = TRUE` to `as_gen_boot_design()` and `make_gen_boot_factors()`. This argument forces the generalized bootstrap variance-covariance estimates for totals to exactly match the target variance estimator. In other words, this eliminates bootstrap simulation error for variance estimates of totals. This is similar to how, for simple survey designs, the jackknife and BRR give variance estimates for totals that exactly match the Horvitz-Thompson estimates. Using `exact_vcov` requires that the number of replicates is strictly greater than the rank of the target variance estimator.

  * Added new variance estimators ("Deville 1" and "Deville 2") available to use for the generalized bootstrap, which are particularly useful for single-stage PPSWOR designs or for multistage designs with one or more stages of PPSWOR sampling. See updated documentation for `as_gen_boot_design()` and `make_quad_form_matrix()`.

  * If the 'srvyr' package is loaded, then functions from 'svrep' that return survey design objects will always return a `tbl_svy` if their input was a `tbl_svy`. This makes it easier to use functions such as `summarize()` or `mutate()`.

* Bug Fixes:
  * Fixed bug where `as_bootstrap_design()` wouldn't create more than 50 replicates for the Rao-Wu, Preston, or Canty-Davison types.

# svrep 0.5.0

* This release adds extensive new functionality for two-phase designs. The new vignette "Replication Methods for Two-phase Sampling" describes the new functionality as well as the underlying statistical methods.

  * The function `as_gen_boot_design()` can now create generalized bootstrap weights
  for two-phase survey design objects created with the 'survey' package's `twophase()` function. 
  The user must specify a list of two variance estimators to use for each phase, e.g. `list('Stratified Multistage SRS', 'Ultimate Cluster')`.
  
  * The function `make_twophase_quad_form()` can be used to create a quadratic form
  for a two-phase variance estimator, by combining quadratic forms from each phase.
  
  * The helper function `get_nearest_psd_matrix()` can be used to approximate
  a quadratic form matrix by the nearest positive semidefinite matrix. This can be
  particularly useful for two-phase designs, since the double
  expansion estimator commonly used in practice frequently does not have a 
  variance estimator which is positive semidefinite.
  
* The function `as_gen_boot_design()` has a new argument named `psd_option`, which
controls what will happen if the target variance estimator 
has a quadratic form matrix which is not positive semi-definite. 
This can occasionally happen, particularly for two-phase designs. By default,
the function will warn the user if the quadratic form is not positive semi-definite
and then automatically approximate the matrix by the nearest positive semi-definite matrix.

* Added a new function `get_design_quad_form()`, 
which determines the quadratic form matrix of a specified variance estimator, 
by parsing the information stored in a survey design object created using the 'survey' package.
  
* Added a new function `rescale_reps()` which implements the rescaling of replicate
adjustment factors to avoid negative replicate weights. This functionality
already existed in `as_gen_boot_design()` and `make_gen_boot_factors()`, but now
it is implemented with the help of this new function.

* Added helper function `is_psd_matrix()` for checking whether a matrix is positive semi-definite,
and added a helper function `get_nearest_psd_matrix()` for approximating a square matrix
by the nearest positive semi-definite matrix.

* Minor improvements to vignettes, particularly formatting.
  
# svrep 0.4.1

* Minor Updates and Bug Fixes:
  * Fix bug in [#15](https://github.com/bschneidr/svrep/issues/15), where 
  bootstrap conversion of multistage survey design objects with `as_bootstrap_design()` 
  would throw an error when user manually specified weights in `svydesign()`.
  
  * Creation of Rao-Wu-Yue-Beaumont bootstrap replicate weights is now faster
  and takes less computer memory.

  * Typo fix in vignettes.

# svrep 0.4.0

* This release adds several functions for creating bootstrap and generalized bootstrap replicate weights. The new vignette "Bootstrap methods for surveys" provides guidance for choosing a bootstrap method and selecting the number of bootstrap replicates to use, along with statistical details and references.

  * Added function `as_bootstrap_design()` to convert a survey design
  object to a replicate design with replicate weights created
  using a bootstrap method. This is essentially a specialized version of 
  `as.svrepdesign()` that supports additional bootstrap methods
  and has detailed documentation about which bootstrap methods can be used
  for different types of sampling designs.
  
  * Added function `as_gen_boot_design()` to convert a survey design
  object to a replicate design with replicate weights created
  using the generalized survey bootstrap. The user must supply the name of
  a target variance estimator (e.g., "Horvitz-Thompson" or "Ultimate Cluster")
  used to create the generalized bootstrap factors. See the new vignette for details.

  * Added functions to help choose the number of bootstrap replicates.
  The function `estimate_boot_sim_cv()` can be used to estimate the simulation error
  in a bootstrap estimate caused by using a finite number of bootstrap replicates.
  The new function `estimate_boot_reps_for_target_cv()` estimates the number of bootstrap
  replicates needed to reduce the simulation error to a target level.

  * Added function `make_rwyb_bootstrap_weights()`, which creates
  bootstrap replicate weights for a wide range of survey designs
  using the method of Rao-Wu-Yue-Beaumont (i.e., Beaumont's 
  generalization of the Rao-Wu-Yue bootstrap method). This function
  can be used directly, or users can specify `as_bootstrap_design(type = "Rao-Wu-Yue-Beaumont")`.

  * Added function `make_gen_boot_factors()` to create replicate adjustment factors
  using the generalized survey bootstrap. The key input to `make_gen_boot_factors()`
  is the matrix of the quadratic form used to represent a variance estimator.
  The new function `make_quad_form_matrix()` can be used to represent a chosen variance
  estimator as a quadratic form, given information about the sample design. This can be
  used for stratified multistage SRS designs (with or without replacement),
  systematic samples, and PPS samples, with or without replacement.

* Minor Updates and Bug Fixes:
  * When using `as_data_frame_with_weights()`,
  ensure that the full-sample weight is named `"FULL_SAMPLE_WGT"`
  if the user does not specify something different.
  * For `calibrate_to_estimate()`, ensure that the output
  names the list of columns with perturbed control columns
  `col_selection` instead of `perturbed_control_cols`,
  so that the name matches the corresponding function argument,
  `col_selection`.
  * Improvements to documentation (formatting tweaks and typo fixes)

# svrep 0.3.0

* Added helper function `as_data_frame_with_weights()` to convert
a survey design object into a data frame with columns of 
weights (full-sample weights and, if applicable, replicate weights).
This is useful for saving data and weights to a data file.

* Added `by` argument to `summarize_rep_weights()` which allows
the specification of one or more grouping variables to use for summaries
(e.g. `by = c('stratum', 'response_status')` can be used to summarize by 
response status within each stratum).

* Added a small vignette "Nonresponse Adjustments" to illustrate how to 
conduct nonresponse adjustments using `redistribute_weights()`.

* Minor Updates and Bug Fixes:
  * Internal code update to avoid annoying but harmless warning message
  about `rho` in `calibrate_to_estimate()`.
  * Bug fix for `stack_replicate_designs()` where designs created with
  `as.svrepdesign(..., type = 'mrbbootstrap')` 
  or `as.svrepdesign(..., type = 'subbootstrap')` threw an error.

# svrep 0.2.0

* Added functions `calibrate_to_estimate()` and `calibrate_to_sample()`
for calibrating to estimated control totals with methods
that account for the sampling variance of the control totals.
For an overview of these functions, please see the new vignette
"Calibrating to Estimated Control Totals".

  * The function `calibrate_to_estimate()` requires the user
  to supply a vector of control totals and its variance-covariance matrix.
  The function applies Fuller's proposed adjustments to the replicate weights,
  in which control totals are varied across replicates by perturbing the control
  totals using a spectral decomposition of the control totals'
  variance-covariance matrix.

  * The function `calibrate_to_sample()` requires the user to supply
  a replicate design for the primary survey of interest as well as a replicate
  design for the control survey used to estimate control totals for calibration.
  The function applies Opsomer & Erciulescu's method of varying 
  the control totals across replicates of the primary survey by matching each 
  primary survey replicate to a replicate from the control survey.
  
* Added an example dataset, `lou_vax_survey`, which is a simulated survey 
measuring Covid-19 vaccination status and a handful of demographic variables,
based on a simple random sample of 1,000 residents of Louisville, Kentucky
with an approximately 50% response rate.
  * An accompanying dataset `lou_pums_microdata` provides person-level microdata
  from the American Community Survey (ACS) 2015-2019 public-use microdata sample
  (PUMS) data for Louisville, KY. The dataset `lou_pums_microdata` includes
  replicate weights to use for variance estimation and can be used to generate
  control totals for `lou_vax_survey`.

# svrep 0.1.0

* Initial release of the package.
