# Data of Louisville Vaccination Survey

A survey measuring Covid-19 vaccination status and a handful of
demographic variables, based on a simple random sample of 1,000
residents of Louisville, Kentucky with an approximately 50% response
rate.  

These data were created using simulation.

## Usage

``` r
data(lou_vax_survey)
```

## Format

A data frame with 1,000 rows and 6 variables

- RESPONSE_STATUS:

  Response status to the survey ('Respondent' or 'Nonrespondent')

- RACE_ETHNICITY:

  Race and Hispanic/Latino ethnicity derived from RAC1P and HISP
  variables of ACS microdata and collapsed to a smaller number of
  categories.

- SEX:

  Male or Female

- EDUC_ATTAINMENT:

  Highest level of education attained ('Less than high school' or 'High
  school or beyond') derived from SCHL variable in ACS microdata and
  collapsed to a smaller number of categories.

- VAX_STATUS:

  Covid-19 vaccination status ('Vaccinated' or 'Unvaccinated')

- SAMPLING_WEIGHT:

  Sampling weight: equal for all cases since data come from a simple
  random sample
