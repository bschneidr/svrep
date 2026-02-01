# Control Totals for Louisville Vaccination Survey

Control totals to use for raking or post-stratification for the
Louisville Vaccination Survey data. Control totals are population size
estimates from the ACS 2015-2019 5-year Public Use Microdata Sample
(PUMS) for specific demographic categories among adults in Jefferson
County, KY.  

These data were created using simulation.

## Usage

``` r
data(lou_vax_survey_control_totals)
```

## Format

A nested list object with two lists, `poststratification` and `raking`,
each of which contains two elements: `estimates` and
`variance-covariance`.

- poststratification:

  Control totals for the combination of `RACE_ETHNICITY`, `SEX`, and
  `EDUC_ATTAINMENT`.

  - estimates: A numeric vector of estimated population totals.

  - variance-covariance: A variance-covariance matrix for the estimated
    population totals.

- raking:

  Separate control totals for each of `RACE_ETHNICITY`, `SEX`, and
  `EDUC_ATTAINMENT`.

  - estimates: A numeric vector of estimated population totals.

  - variance-covariance: A variance-covariance matrix for the estimated
    population totals.
