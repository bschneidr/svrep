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
#> 1          Ineligible    39    15           13   14    1896.620    164.4527
#> 2       Nonrespondent    47    15           14   15    2296.912    133.9403
#> 3          Respondent    41    15           13   14    2000.468    130.3750
#> 4 Unknown eligibility    56    15           -1    0       0.000      0.0000
#>   min_rep_wgt max_rep_wgt
#> 1           0    56.98729
#> 2           0    56.98729
#> 3           0    56.98729
#> 4           0     0.00000

summarize_rep_weights(ue_adjusted_design, type = 'overall',
                      by = c("stype", "response_status"))
#>    stype     response_status nrows ncols degf_svy_pkg rank avg_wgt_sum
#> 1      E          Ineligible    29    15            7    8  1413.77685
#> 2      H          Ineligible     6    15            2    3   283.80822
#> 3      M          Ineligible     4    15            3    4   199.03463
#> 4      E       Nonrespondent    36    15           12   13  1753.97013
#> 5      H       Nonrespondent     2    15            1    2    95.02487
#> 6      M       Nonrespondent     9    15            5    6   447.91713
#> 7      E          Respondent    35    15           10   11  1706.22048
#> 8      H          Respondent     2    15            1    2    95.02487
#> 9      M          Respondent     4    15            2    3   199.22315
#> 10     E Unknown eligibility    44    15           -1    0     0.00000
#> 11     H Unknown eligibility     4    15           -1    0     0.00000
#> 12     M Unknown eligibility     8    15           -1    0     0.00000
#>    sd_wgt_sums min_rep_wgt max_rep_wgt
#> 1    149.65790           0    53.40068
#> 2     40.20252           0    56.98729
#> 3     23.57306           0    56.98729
#> 4    121.79339           0    53.40068
#> 5     18.18279           0    56.98729
#> 6     46.50443           0    56.98729
#> 7    135.21908           0    53.40068
#> 8     18.18279           0    56.98729
#> 9     32.77474           0    56.98729
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
#> 1   183    15           14   15        6194    403.1741           0    56.98729
#> 
```
