# Package index

## Convert a Survey Design Object to a Data Frame

- [`as_data_frame_with_weights()`](https://bschneidr.github.io/svrep/reference/as_data_frame_with_weights.md)
  : Convert Survey Design to Data Frame

## Nonresponse Adjustments

- [`redistribute_weights()`](https://bschneidr.github.io/svrep/reference/redistribute_weights.md)
  : Weighting Class Adjustment

## Creating Bootstrap Replicate Weights

- [`as_bootstrap_design()`](https://bschneidr.github.io/svrep/reference/as_bootstrap_design.md)
  : Convert Survey Design to Bootstrap Replicate Design
- [`make_doubled_half_bootstrap_weights()`](https://bschneidr.github.io/svrep/reference/make_doubled_half_bootstrap_weights.md)
  : Weights for "Doubled Half Bootstrap" of Antal and Tillé (2014)
- [`make_rwyb_bootstrap_weights()`](https://bschneidr.github.io/svrep/reference/make_rwyb_bootstrap_weights.md)
  : Weights for Rao-Wu-Yue-Beaumont Bootstrap
- [`estimate_boot_reps_for_target_cv()`](https://bschneidr.github.io/svrep/reference/estimate_boot_reps_for_target_cv.md)
  : Control Bootstrap Simulation Error
- [`estimate_boot_sim_cv()`](https://bschneidr.github.io/svrep/reference/estimate_boot_sim_cv.md)
  : Estimate Bootstrap Simulation Error

## Generalized Replication Method (including the Generalized Bootstrap)

### Fay’s Generalized Replication Method

- [`as_fays_gen_rep_design()`](https://bschneidr.github.io/svrep/reference/as_fays_gen_rep_design.md)
  : Convert Survey Design to Fay's Generalized Replication Design
- [`make_fays_gen_rep_factors()`](https://bschneidr.github.io/svrep/reference/make_fays_gen_rep_factors.md)
  : Factors for Fay's Generalized Replication Method

### Generalized Bootstrap

- [`as_gen_boot_design()`](https://bschneidr.github.io/svrep/reference/as_gen_boot_design.md)
  : Convert Survey Design to Generalized Bootstrap Replicate Design
- [`make_gen_boot_factors()`](https://bschneidr.github.io/svrep/reference/make_gen_boot_factors.md)
  : Factors for the Generalized Survey Bootstrap

### Variance Estimators Available for the Generalized Replication Methods

*The generalized replication methods (Fay’s method and the generalized
bootstrap) work by “mimicking” a target variance estimator, such as the
Horvitz-Thompson estimator. This help page describes the variance
estimators that can be used as the target for the generalized
replication methods.*

- [`variance-estimators`](https://bschneidr.github.io/svrep/reference/variance-estimators.md)
  : Variance Estimators

### Helper Functions for Working with Quadratic Forms

*These functions help the user specify the quadratic form representation
of common variance estimators and, if necessary, adjust them so that
they are positive semidefinite (a necessary prerequisite for using the
generalized replication methods).*

- [`get_design_quad_form()`](https://bschneidr.github.io/svrep/reference/get_design_quad_form.md)
  : Quadratic Form Matrix of Variance Estimator for a Survey Design
- [`make_quad_form_matrix()`](https://bschneidr.github.io/svrep/reference/make_quad_form_matrix.md)
  : Get Variance Estimator's Quadratic Form Matrix
- [`make_twophase_quad_form()`](https://bschneidr.github.io/svrep/reference/make_twophase_quad_form.md)
  : Quadratic Form Matrix for a Two-phase Design
- [`get_nearest_psd_matrix()`](https://bschneidr.github.io/svrep/reference/get_nearest_psd_matrix.md)
  : Positive Semidefinite Matrix Approximation
- [`is_psd_matrix()`](https://bschneidr.github.io/svrep/reference/is_psd_matrix.md)
  : Check if Matrix is Positive Semidefinite
- [`make_kernel_var_matrix()`](https://bschneidr.github.io/svrep/reference/make_kernel_var_matrix.md)
  : Quadratic Form Matrix of Kernel-based Variance Estimator

## Creating Jackknife Replicate Weights

- [`as_random_group_jackknife_design()`](https://bschneidr.github.io/svrep/reference/as_random_group_jackknife_design.md)
  : Convert Survey Design to Random Group Jackknife Replicate Design

## Create Replicate Weights Using Successive Differences Replication

- [`as_sdr_design()`](https://bschneidr.github.io/svrep/reference/as_sdr_design.md)
  : Convert Survey Design to Successive Differences Replicate Design
- [`make_sdr_replicate_factors()`](https://bschneidr.github.io/svrep/reference/make_sdr_replicate_factors.md)
  : Factors for the Successive Difference Replication Method

## Derive Replicate Weights for a Two-phase Sample

- [`derive_twophase_rep_design()`](https://bschneidr.github.io/svrep/reference/derive_twophase_rep_design.md)
  : Replicate Design Object for a Two-phase Sample

## Calibrating to Estimated Control Totals

- [`calibrate_to_estimate()`](https://bschneidr.github.io/svrep/reference/calibrate_to_estimate.md)
  : Sample-based Calibration to An Estimate
- [`calibrate_to_sample()`](https://bschneidr.github.io/svrep/reference/calibrate_to_sample.md)
  : Sample-based Calibration with Replicates

## General-Purpose Helper Functions

- [`shuffle_replicates()`](https://bschneidr.github.io/svrep/reference/shuffle_replicates.md)
  : Shuffle Order of Replicates in a Replicate Design Object
- [`stack_replicate_designs()`](https://bschneidr.github.io/svrep/reference/stack_replicate_designs.md)
  : Combine Replicate Designs by Stacking
- [`subsample_replicates()`](https://bschneidr.github.io/svrep/reference/subsample_replicates.md)
  : Subsample Replicates in a Replicate Design Object
- [`rescale_replicates()`](https://bschneidr.github.io/svrep/reference/rescale_replicates.md)
  : Rescale Replicate Factors
- [`add_inactive_replicates()`](https://bschneidr.github.io/svrep/reference/add_inactive_replicates.md)
  : Add Inactive Replicates to a Survey Design Object
- [`svyby_repwts()`](https://bschneidr.github.io/svrep/reference/svyby_repwts.md)
  : Survey Statistics for Multiple Replicate Design Objects
- [`get_rep_scale_coefs()`](https://bschneidr.github.io/svrep/reference/get_rep_scale_coefs.md)
  : Access Replication Scale Coefficients
- [`get_rep_type()`](https://bschneidr.github.io/svrep/reference/get_rep_type.md)
  : Access Type of Replication Method

## Diagnostic Functions to Check Replicate Weights

- [`summarize_rep_weights()`](https://bschneidr.github.io/svrep/reference/summarize_rep_weights.md)
  : Summarize Replicate Weights in a Replicate Design

## Package-level Options

- [`svrep-package-options`](https://bschneidr.github.io/svrep/reference/svrep-package-options.md)
  : Package-level Options for svrep

## Example Datasets

- [`lou_pums_microdata`](https://bschneidr.github.io/svrep/reference/lou_pums_microdata.md)
  : Auxiliary Data for Louisville Vaccination Survey
- [`lou_vax_survey`](https://bschneidr.github.io/svrep/reference/lou_vax_survey.md)
  : Data of Louisville Vaccination Survey
- [`lou_vax_survey_control_totals`](https://bschneidr.github.io/svrep/reference/lou_vax_survey_control_totals.md)
  : Control Totals for Louisville Vaccination Survey
- [`library_census`](https://bschneidr.github.io/svrep/reference/libraries.md)
  [`library_multistage_sample`](https://bschneidr.github.io/svrep/reference/libraries.md)
  [`library_stsys_sample`](https://bschneidr.github.io/svrep/reference/libraries.md)
  : U.S. Public Libraries Survey
