# Convert Survey Design to Data Frame

Convert a survey design object to a data frame with weights stored as
columns

## Usage

``` r
as_data_frame_with_weights(
  design,
  full_wgt_name = "FULL_SAMPLE_WGT",
  rep_wgt_prefix = "REP_WGT_",
  vars_to_keep = NULL
)
```

## Arguments

- design:

  A survey design object, created with either the `survey` or `srvyr`
  packages.

- full_wgt_name:

  The column name to use for the full-sample weights

- rep_wgt_prefix:

  For replicate design objects, a prefix to use for the column names of
  the replicate weights. The column names will be created by appending
  the replicate number after the prefix.

- vars_to_keep:

  By default, all variables in the data will be kept. To select only a
  subset of the non-weight variables, you can supply a character vector
  of variable names to keep.

## Value

A data frame, with new columns containing the weights from the survey
design object

## Examples

``` r
data("lou_vax_survey", package = 'svrep')

# Create a survey design object
survey_design <- svydesign(data = lou_vax_survey,
                           weights = ~ SAMPLING_WEIGHT,
                           ids = ~ 1)

rep_survey_design <- as.svrepdesign(survey_design,
                                    type = "boot",
                                    replicates = 10)

# Adjust the weights for nonresponse
nr_adjusted_design <- redistribute_weights(
  design = rep_survey_design,
  reduce_if = RESPONSE_STATUS == "Nonrespondent",
  increase_if = RESPONSE_STATUS == "Respondent",
  by = c("RACE_ETHNICITY", "EDUC_ATTAINMENT")
)

# Save the survey design object as a data frame
nr_adjusted_data <- as_data_frame_with_weights(
  nr_adjusted_design,
  full_wgt_name = "NR_ADJUSTED_WGT",
  rep_wgt_prefix = "NR_ADJUSTED_REP_WGT_"
)
head(nr_adjusted_data)
#>   RESPONSE_STATUS                                          RACE_ETHNICITY
#> 1   Nonrespondent                     White alone, not Hispanic or Latino
#> 2   Nonrespondent Black or African American alone, not Hispanic or Latino
#> 3      Respondent                     White alone, not Hispanic or Latino
#> 4   Nonrespondent                     White alone, not Hispanic or Latino
#> 5   Nonrespondent                     White alone, not Hispanic or Latino
#> 6      Respondent                     White alone, not Hispanic or Latino
#>      SEX       EDUC_ATTAINMENT VAX_STATUS SAMPLING_WEIGHT NR_ADJUSTED_WGT
#> 1 Female Less than high school       <NA>         596.702           0.000
#> 2 Female High school or beyond       <NA>         596.702           0.000
#> 3 Female Less than high school Vaccinated         596.702        1223.239
#> 4 Female Less than high school       <NA>         596.702           0.000
#> 5 Female High school or beyond       <NA>         596.702           0.000
#> 6 Female High school or beyond Vaccinated         596.702        1059.068
#>   NR_ADJUSTED_REP_WGT_1 NR_ADJUSTED_REP_WGT_2 NR_ADJUSTED_REP_WGT_3
#> 1                 0.000                 0.000                 0.000
#> 2                 0.000                 0.000                 0.000
#> 3              2418.717              2411.163                 0.000
#> 4                 0.000                 0.000                 0.000
#> 5                 0.000                 0.000                 0.000
#> 6              1063.961                 0.000              1069.499
#>   NR_ADJUSTED_REP_WGT_4 NR_ADJUSTED_REP_WGT_5 NR_ADJUSTED_REP_WGT_6
#> 1                 0.000                 0.000                 0.000
#> 2                 0.000                 0.000                 0.000
#> 3              1186.883                 0.000              1113.019
#> 4                 0.000                 0.000                 0.000
#> 5                 0.000                 0.000                 0.000
#> 6                 0.000              2331.301              1075.340
#>   NR_ADJUSTED_REP_WGT_7 NR_ADJUSTED_REP_WGT_8 NR_ADJUSTED_REP_WGT_9
#> 1                 0.000                 0.000                 0.000
#> 2                 0.000                 0.000                 0.000
#> 3              2320.171              1130.892              1108.161
#> 4                 0.000                 0.000                 0.000
#> 5                 0.000                 0.000                 0.000
#> 6              1001.025              1968.179              2003.435
#>   NR_ADJUSTED_REP_WGT_10
#> 1                  0.000
#> 2                  0.000
#> 3               2412.065
#> 4                  0.000
#> 5                  0.000
#> 6               1054.884

# Check the column names of the result
colnames(nr_adjusted_data)
#>  [1] "RESPONSE_STATUS"        "RACE_ETHNICITY"         "SEX"                   
#>  [4] "EDUC_ATTAINMENT"        "VAX_STATUS"             "SAMPLING_WEIGHT"       
#>  [7] "NR_ADJUSTED_WGT"        "NR_ADJUSTED_REP_WGT_1"  "NR_ADJUSTED_REP_WGT_2" 
#> [10] "NR_ADJUSTED_REP_WGT_3"  "NR_ADJUSTED_REP_WGT_4"  "NR_ADJUSTED_REP_WGT_5" 
#> [13] "NR_ADJUSTED_REP_WGT_6"  "NR_ADJUSTED_REP_WGT_7"  "NR_ADJUSTED_REP_WGT_8" 
#> [16] "NR_ADJUSTED_REP_WGT_9"  "NR_ADJUSTED_REP_WGT_10"
```
