# svrep (development version)

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
