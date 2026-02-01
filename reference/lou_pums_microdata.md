# Auxiliary Data for Louisville Vaccination Survey

Auxiliary data that can be used to produce control totals for the
Louisville Vaccination Survey. Consists of person-level microdata from
the American Community Survey (ACS) 2015-2019 public-use microdata
sample (PUMS) data for Louisville, KY. This microdata sample represents
all adults (persons aged 18 or over) in Louisville, KY.  

These data include replicate weights to use for variance estimation.

## Usage

``` r
data(lou_pums_microdata)
```

## Format

A data frame with 80 rows and 85 variables

- UNIQUE_ID: Unique identifier for records

- AGE: Age in years (copied from the AGEP variable in the ACS microdata)

- RACE_ETHNICITY: Race and Hispanic/Latino ethnicity derived from RAC1P
  and HISP variables of ACS microdata and collapsed to a smaller number
  of categories.

- SEX: Male or Female

- EDUC_ATTAINMENT: Highest level of education attained ('Less than high
  school' or 'High school or beyond') derived from SCHL variable in ACS
  microdata and collapsed to a smaller number of categories.

- PWGTP: Weights for the full-sample

- PWGTP1-PWGTP80: 80 columns of replicate weights created using the
  Successive Differences Replication (SDR) method.

## Examples

``` r
# \donttest{
data(lou_pums_microdata)

# Prepare the data for analysis with the survey package

  lou_pums_rep_design <- svrepdesign(
    data       = lou_pums_microdata,
    variables  = ~ UNIQUE_ID + AGE + SEX + RACE_ETHNICITY + EDUC_ATTAINMENT,
    weights    = ~ PWGTP, 
    repweights = "PWGTP\\d{1,2}",
    type       = "successive-difference",
    mse        = TRUE
  )

# Estimate population proportions
  svymean(~ SEX, design = lou_pums_rep_design)
#>              mean    SE
#> SEXMale   0.47543 7e-04
#> SEXFemale 0.52457 7e-04
# }
```
