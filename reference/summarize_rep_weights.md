# Summarize Replicate Weights in a Replicate Design

Summarize the replicate weights of a design

## Usage

``` r
summarize_rep_weights(rep_design, type = "both", by)
```

## Arguments

- rep_design:

  A replicate design object, created with either the `survey` or `srvyr`
  packages.

- type:

  Default is `"both"`. Use `type = "overall"`, for an overall summary of
  the replicate weights. Use `type = "specific"` for a summary of each
  column of replicate weights, with each column of replicate weights
  summarized in a given row of the summary.  
    
  Use `type = "both"` for a list containing both summaries, with the
  list containing the names `"overall"` and `"both"`.

- by:

  (Optional) A character vector with the names of variables used to
  group the summaries.

## Value

If `type = "both"` (the default), the result is a list of data frames
with names `"overall"` and `"specific"`. If `type = "overall"`, the
result is a data frame providing an overall summary of the replicate
weights.  
  
The contents of the `"overall"` summary are the following:

- "nrows": Number of rows for the weights

- "ncols": Number of columns of replicate weights

- "degf_svy_pkg": The degrees of freedom according to the survey package
  in R

- "rank": The matrix rank as determined by a QR decomposition

- "avg_wgt_sum": The average column sum

- "sd_wgt_sums": The standard deviation of the column sums

- "min_rep_wgt": The minimum value of any replicate weight

- "max_rep_wgt": The maximum value of any replicate weight

If `type = "specific"`, the result is a data frame providing a summary
of each column of replicate weights, with each column of replicate
weights described in a given row of the data frame. The contents of the
`"specific"` summary are the following:

- "Rep_Column": The name of a given column of replicate weights. If
  columns are unnamed, the column number is used instead

- "N": The number of entries

- "N_NONZERO": The number of nonzero entries

- "SUM": The sum of the weights

- "MEAN": The average of the weights

- "CV": The coefficient of variation of the weights (standard deviation
  divided by mean)

- "MIN": The minimum weight

- "MAX": The maximum weight

## Examples

``` r
# Load example data
data(api)

dclus1 <- svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
dclus1$variables$response_status <- sample(x = c("Respondent", "Nonrespondent",
                                                 "Ineligible", "Unknown eligibility"),
                                           size = nrow(dclus1),
                                           replace = TRUE)
rep_design <- as.svrepdesign(dclus1)

# Adjust weights for cases with unknown eligibility
ue_adjusted_design <- redistribute_weights(
    design = rep_design,
    reduce_if = response_status %in% c("Unknown eligibility"),
    increase_if = !response_status %in% c("Unknown eligibility"),
    by = c("stype")
)

# Summarize replicate weights

summarize_rep_weights(rep_design, type = "both")
#> $overall
#>   nrows ncols degf_svy_pkg rank avg_wgt_sum sd_wgt_sums min_rep_wgt max_rep_wgt
#> 1   183    15           14   15        6194    403.1741           0    36.26464
#> 
#> $specific
#>    Rep_Column   N N_NONZERO      SUM     MEAN         CV MIN      MAX
#> 1           1 183       172 6237.518 34.08480 0.25358407   0 36.26464
#> 2           2 183       179 6491.370 35.47197 0.14989713   0 36.26464
#> 3           3 183       181 6563.900 35.86830 0.10540606   0 36.26464
#> 4           4 183       170 6164.989 33.68846 0.27729183   0 36.26464
#> 5           5 183       181 6563.900 35.86830 0.10540606   0 36.26464
#> 6           6 183       179 6491.370 35.47197 0.14989713   0 36.26464
#> 7           7 183       179 6491.370 35.47197 0.14989713   0 36.26464
#> 8           8 183       167 6056.195 33.09396 0.31037848   0 36.26464
#> 9           9 183       174 6310.047 34.48113 0.22805336   0 36.26464
#> 10         10 183       149 5403.431 29.52695 0.47900073   0 36.26464
#> 11         11 183       162 5874.872 32.10312 0.36102892   0 36.26464
#> 12         12 183       146 5294.637 28.93244 0.50479412   0 36.26464
#> 13         13 183       170 6164.989 33.68846 0.27729183   0 36.26464
#> 14         14 183       182 6600.164 36.06647 0.07432829   0 36.26464
#> 15         15 183       171 6201.253 33.88663 0.26563324   0 36.26464
#> 

# Summarize replicate weights by grouping variables

summarize_rep_weights(ue_adjusted_design, type = 'overall',
                      by = c("response_status"))
#>       response_status nrows ncols degf_svy_pkg rank avg_wgt_sum sd_wgt_sums
#> 1          Ineligible    39    15           14   15    1883.610    164.0280
#> 2       Nonrespondent    47    15           14   15    2291.801    117.9820
#> 3          Respondent    42    15           13   14    2018.589    145.4286
#> 4 Unknown eligibility    55    15           -1    0       0.000      0.0000
#>   min_rep_wgt max_rep_wgt
#> 1           0    66.48517
#> 2           0    66.48517
#> 3           0    66.48517
#> 4           0     0.00000

summarize_rep_weights(ue_adjusted_design, type = 'overall',
                      by = c("stype", "response_status"))
#>    stype     response_status nrows ncols degf_svy_pkg rank avg_wgt_sum
#> 1      E          Ineligible    29    15            9   10   1385.9323
#> 2      H          Ineligible     5    15            4    5    215.5218
#> 3      M          Ineligible     5    15            4    5    282.1557
#> 4      E       Nonrespondent    37    15           14   15   1767.8095
#> 5      H       Nonrespondent     3    15            2    3    129.0337
#> 6      M       Nonrespondent     7    15            6    7    394.9578
#> 7      E          Respondent    36    15            9   10   1720.2257
#> 8      H          Respondent     3    15            2    3    129.3024
#> 9      M          Respondent     3    15            2    3    169.0614
#> 10     E Unknown eligibility    42    15           -1    0      0.0000
#> 11     H Unknown eligibility     3    15           -1    0      0.0000
#> 12     M Unknown eligibility    10    15           -1    0      0.0000
#>    sd_wgt_sums min_rep_wgt max_rep_wgt
#> 1    157.68757           0    52.64222
#> 2     26.35779           0    48.35285
#> 3     29.62667           0    66.48517
#> 4    103.43058           0    52.64222
#> 5     19.24917           0    48.35285
#> 6     30.12000           0    66.48517
#> 7    136.71653           0    52.64222
#> 8     20.44378           0    48.35285
#> 9     23.88955           0    66.48517
#> 10     0.00000           0     0.00000
#> 11     0.00000           0     0.00000
#> 12     0.00000           0     0.00000

# Compare replicate weights

rep_wt_summaries <- lapply(list('original' = rep_design,
                                'adjusted' = ue_adjusted_design),
                           summarize_rep_weights,
                           type = "overall")
print(rep_wt_summaries)
#> $original
#>   nrows ncols degf_svy_pkg rank avg_wgt_sum sd_wgt_sums min_rep_wgt max_rep_wgt
#> 1   183    15           14   15        6194    403.1741           0    36.26464
#> 
#> $adjusted
#>   nrows ncols degf_svy_pkg rank avg_wgt_sum sd_wgt_sums min_rep_wgt max_rep_wgt
#> 1   183    15           14   15        6194    403.1741           0    66.48517
#> 
```
