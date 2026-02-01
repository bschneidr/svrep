# Nonresponse Adjustments

In this short vignette, we’ll demonstrate how the svrep package can be
used to implement weighting adjustments for nonresponse in survey
samples. For illustration, we’ll use an example survey measuring
Covid-19 vaccination status and a handful of demographic variables,
based on a simple random sample of 1,000 residents of Louisville,
Kentucky.

``` r
library(dplyr) # For data manipulation
library(survey) # For complex survey analysis
library(srvyr) # For complex survey analysis with dplyr syntax
library(svrep)

# Load and inspect the data
data("lou_vax_survey", package = 'svrep')
head(lou_vax_survey)
#> # A tibble: 6 × 6
#>   RESPONSE_STATUS RACE_ETHNICITY                SEX   EDUC_ATTAINMENT VAX_STATUS
#>   <chr>           <chr>                         <chr> <chr>           <chr>     
#> 1 Nonrespondent   White alone, not Hispanic or… Fema… Less than high… NA        
#> 2 Nonrespondent   Black or African American al… Fema… High school or… NA        
#> 3 Respondent      White alone, not Hispanic or… Fema… Less than high… Vaccinated
#> 4 Nonrespondent   White alone, not Hispanic or… Fema… Less than high… NA        
#> 5 Nonrespondent   White alone, not Hispanic or… Fema… High school or… NA        
#> 6 Respondent      White alone, not Hispanic or… Fema… High school or… Vaccinated
#> # ℹ 1 more variable: SAMPLING_WEIGHT <dbl>
colnames(lou_vax_survey)
#> [1] "RESPONSE_STATUS" "RACE_ETHNICITY"  "SEX"             "EDUC_ATTAINMENT"
#> [5] "VAX_STATUS"      "SAMPLING_WEIGHT"
```

This vaccination survey has an overall response rate of 50.2%, which
means that estimated vaccination rates may be significantly biased by
nonresponse. We’ll use nonresponse weighting adjustments to try and
reduce potential nonresponse bias.

``` r
lou_vax_survey |> count(RESPONSE_STATUS) |> mutate(pct = n/sum(n))
#> # A tibble: 2 × 3
#>   RESPONSE_STATUS     n   pct
#>   <chr>           <int> <dbl>
#> 1 Nonrespondent     498 0.498
#> 2 Respondent        502 0.502
```

## Creating initial replicate weights

To begin with, we’ll create bootstrap replicate weights. In most cases,
we can do this simply by describing the survey design using the
[`svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html) function
and then using a function to create appropriate replicate weights. The
function
[`as.svrepdesign()`](https://rdrr.io/pkg/survey/man/as.svrepdesign.html)
from the ‘survey’ package can be used to create several types of
replicate weights, using the argument `type` (with options `'JK1'`,
`'JKn'`, `'bootstrap'`, `'BRR'`, `'Fay'`, etc.) In addition, the
function
[`as_bootstrap_design()`](https://bschneidr.github.io/svrep/reference/as_bootstrap_design.md)
can be used to create bootstrap weights using additional methods not
supported in the ‘survey’ package.

``` r
# Describe the survey design
lou_vax_survey <- svydesign(ids = ~ 1, weights = ~ SAMPLING_WEIGHT,
                            data = lou_vax_survey)

print(lou_vax_survey)
#> Independent Sampling design (with replacement)
#> svydesign(ids = ~1, weights = ~SAMPLING_WEIGHT, data = lou_vax_survey)

# Create appropriate replicate weights
lou_vax_survey <- lou_vax_survey |>
  as_bootstrap_design(replicates = 100, mse = TRUE,
                      type = "Rao-Wu-Yue-Beaumont")

print(lou_vax_survey)
#> Call: as_bootstrap_design(lou_vax_survey, replicates = 100, mse = TRUE, 
#>     type = "Rao-Wu-Yue-Beaumont")
#> Survey bootstrap with 100 replicates and MSE variances.
```

For convenience, we’ll convert this survey design object into an object
with class `tbl_svy`, which allows us to use convenient tidyverse/dplyr
syntax
([`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html),
[`summarize()`](https://dplyr.tidyverse.org/reference/summarise.html),
etc.) as well as other helpful functions from the srvyr package.

``` r
lou_vax_survey <- lou_vax_survey |> as_survey()

print(lou_vax_survey)
#> Call: Called via srvyr
#> Survey bootstrap with 100 replicates and MSE variances.
#> Data variables: 
#>   - RESPONSE_STATUS (chr), RACE_ETHNICITY (chr), SEX (chr), EDUC_ATTAINMENT
#>     (chr), VAX_STATUS (chr), SAMPLING_WEIGHT (dbl)
```

## Redistributing weight from nonrespondents to respondents

A common form of nonresponse adjustment is to simply ‘redistribute’
weight from the nonrespondents to the respondents. In other words, the
weight for each nonrespondent is set to $0$, and the weight for each
respondent is increased by a factor greater than one so that the sum of
adjusted weights in the sample of respondents equals the sum of
unadjusted weights from the full sample. For example, if the sum of
weights among respondents is $299,544.4$ and the sum of weights among
nonrespondents is $297,157.6$, then a basic nonresponse adjustment would
set the weights among nonrespondents to $0$ and multiply the weight for
each respondent by an adjustment factor equal to
$1 + (297,157.6/299,544.4)$. This type of adjustment is succinctly
described in mathematical notation below.

$$\begin{aligned}
w_{i} & {= {\textit{𝑂𝑟𝑖𝑔𝑖𝑛𝑎𝑙 𝑠𝑎𝑚𝑝𝑙𝑖𝑛𝑔 𝑤𝑒𝑖𝑔h𝑡 𝑓𝑜𝑟 𝑐𝑎𝑠𝑒}\mspace{6mu}}i} \\
 & {= 1/\pi_{i},{\mspace{6mu}\textit{𝑤h𝑒𝑟𝑒}\mspace{6mu}}\pi_{i}{\mspace{6mu}\textit{𝑖𝑠 𝑡h𝑒 𝑝𝑟𝑜𝑏𝑎𝑏𝑖𝑙𝑖𝑡𝑦 𝑐𝑎𝑠𝑒 𝑖}\mspace{6mu}}\textit{h𝑎𝑑 𝑜𝑓 𝑏𝑒𝑖𝑛𝑔 𝑠𝑎𝑚𝑝𝑙𝑒𝑑}} \\
f_{NR,i} & {= \textit{𝑁𝑜𝑛𝑟𝑒𝑠𝑝𝑜𝑛𝑠𝑒 𝑎𝑑𝑗𝑢𝑠𝑡𝑚𝑒𝑛𝑡 𝑓𝑎𝑐𝑡𝑜𝑟 𝑓𝑜𝑟 𝑐𝑎𝑠𝑒 𝑖}} \\
w_{NR,i} & {= w_{i} \times f_{NR,i} = {\textit{𝑊𝑒𝑖𝑔h𝑡 𝑓𝑜𝑟 𝑐𝑎𝑠𝑒}\mspace{6mu}}i{\mspace{6mu}\textit{𝑎𝑓𝑡𝑒𝑟 𝑛𝑜𝑛𝑟𝑒𝑠𝑝𝑜𝑛𝑠𝑒 𝑎𝑑𝑗𝑢𝑠𝑡𝑚𝑒𝑛𝑡}}} \\
 & \\
{\sum\limits_{i \in s_{resp}}w_{i}} & {= \textit{𝑆𝑢𝑚 𝑜𝑓 𝑠𝑎𝑚𝑝𝑙𝑖𝑛𝑔 𝑤𝑒𝑖𝑔h𝑡𝑠 𝑎𝑚𝑜𝑛𝑔 𝑟𝑒𝑠𝑝𝑜𝑛𝑑𝑒𝑛𝑡𝑠}} \\
{\sum\limits_{i \in s_{nonresp}}w_{i}} & {= \textit{𝑆𝑢𝑚 𝑜𝑓 𝑠𝑎𝑚𝑝𝑙𝑖𝑛𝑔 𝑤𝑒𝑖𝑔h𝑡𝑠 𝑎𝑚𝑜𝑛𝑔 𝑛𝑜𝑛𝑟𝑒𝑠𝑝𝑜𝑛𝑑𝑒𝑛𝑡𝑠}} \\
 & \\
f_{NR_{i}} & {= \begin{cases}
0 & {{\textit{𝑖𝑓 𝑐𝑎𝑠𝑒}\mspace{6mu}}i{\mspace{6mu}\textit{𝑖𝑠 𝑎 𝑛𝑜𝑛𝑟𝑒𝑠𝑝𝑜𝑛𝑑𝑒𝑛𝑡}}} \\
{1 + \frac{\sum\limits_{i \in s_{nonresp}}w_{i}}{\sum\limits_{i \in s_{resp}}w_{i}}} & {{\textit{𝑖𝑓 𝑐𝑎𝑠𝑒}\mspace{6mu}}i{\mspace{6mu}\textit{𝑖𝑠 𝑎 𝑟𝑒𝑠𝑝𝑜𝑛𝑑𝑒𝑛𝑡}}}
\end{cases}}
\end{aligned}$$

We’ll illustrate this type of adjustment with the Louisville vaccination
survey. First, we’ll inspect the sum of the sampling weights for
respondents, nonrespondents, and the overall sample.

``` r
# Weights before adjustment
lou_vax_survey |>
  group_by(RESPONSE_STATUS) |>
  cascade(
    `Sum of Weights` = sum(cur_svy_wts()),
    .fill = "TOTAL"
  )
#> # A tibble: 3 × 2
#>   RESPONSE_STATUS `Sum of Weights`
#>   <chr>                      <dbl>
#> 1 Nonrespondent            297158.
#> 2 Respondent               299544.
#> 3 TOTAL                    596702
```

Next, we’ll redistribute weight from nonrespondents to respondents using
the
[`redistribute_weights()`](https://bschneidr.github.io/svrep/reference/redistribute_weights.md)
function, which adjusts the full-sample weights as well as each set of
replicate weights. To specify which subset of data should have its
weights reduced, we supply a logical expression to the argument
`reduce_if`. To specify which subset of data should have its weights
increased, we supply a logical expression to the argument `increase_if`.

``` r
# Conduct a basic nonresponse adjustment
nr_adjusted_survey <- lou_vax_survey |>
  redistribute_weights(
    reduce_if = RESPONSE_STATUS == "Nonrespondent",
    increase_if = RESPONSE_STATUS == "Respondent"
  )
```

After making the adjustment, we can check that all of the weight from
nonrespondents has been redistributed to respondents.

``` r
# Check the sum of full-sample weights by response status
nr_adjusted_survey |>
  group_by(RESPONSE_STATUS) |>
  cascade(
    `Sum of Weights` = sum(cur_svy_wts()),
    .fill = "TOTAL"
  )
#> # A tibble: 3 × 2
#>   RESPONSE_STATUS `Sum of Weights`
#>   <chr>                      <dbl>
#> 1 Nonrespondent                  0
#> 2 Respondent                596702
#> 3 TOTAL                     596702
```

``` r
# Check sums of replicate weights by response status
nr_adjusted_survey |>
  summarize_rep_weights(
    type = "specific",
    by = "RESPONSE_STATUS"
  ) |> 
  arrange(Rep_Column, RESPONSE_STATUS) |>
  head(10)
#>    RESPONSE_STATUS Rep_Column   N N_NONZERO    SUM     MEAN        CV MIN
#> 1    Nonrespondent          1 498         0      0    0.000       NaN   0
#> 2       Respondent          1 502       309 596702 1188.649 1.0131335   0
#> 3    Nonrespondent          2 498         0      0    0.000       NaN   0
#> 4       Respondent          2 502       323 596702 1188.649 0.9866873   0
#> 5    Nonrespondent          3 498         0      0    0.000       NaN   0
#> 6       Respondent          3 502       319 596702 1188.649 0.9799219   0
#> 7    Nonrespondent          4 498         0      0    0.000       NaN   0
#> 8       Respondent          4 502       322 596702 1188.649 0.9697204   0
#> 9    Nonrespondent          5 498         0      0    0.000       NaN   0
#> 10      Respondent          5 502       323 596702 1188.649 0.9403356   0
#>         MAX
#> 1     0.000
#> 2  4962.179
#> 3     0.000
#> 4  5737.519
#> 5     0.000
#> 6  7247.393
#> 7     0.000
#> 8  4745.145
#> 9     0.000
#> 10 4792.787
```

## Conducting weighting class adjustments

Nonresponse bias is liable to occur if different subpopulations
systematically differ in terms of their response rates to the survey and
also differ in terms of what the survey is trying to measure (in this
case, vaccination status). In our example, we can see some fairly large
differences in response rates across different race/ethnicity groups.

``` r
lou_vax_survey |>
  group_by(RACE_ETHNICITY) |>
  summarize(Response_Rate = mean(RESPONSE_STATUS == "Respondent"),
            Sample_Size = n(),
            n_Respondents = sum(RESPONSE_STATUS == "Respondent"))
#> # A tibble: 4 × 4
#>   RACE_ETHNICITY                         Response_Rate Sample_Size n_Respondents
#>   <chr>                                          <dbl>       <int>         <int>
#> 1 Black or African American alone, not …         0.452         188            85
#> 2 Hispanic or Latino                             0.378          45            17
#> 3 Other Race, not Hispanic or Latino             0.492          59            29
#> 4 White alone, not Hispanic or Latino            0.524         708           371
```

Weighting adjustments may be able to help reduce nonresponse bias caused
by these differences in response rates. One standard form of adjustment
known as **weighting class adjustment** is to redistribute weights from
nonrespondents to respondents separately by different categories of
auxiliary variables (such as race/ethnicity). The survey textbook
Heeringa, West, and Berglund (2017) provides an excellent overview of
weighting class adjustments. To implement a weighting class adjustment
with the svrep package, we can simply use the `by` argument of
[`redistribute_weights()`](https://bschneidr.github.io/svrep/reference/redistribute_weights.md).

``` r
nr_adjusted_survey <- lou_vax_survey |>
  redistribute_weights(
    reduce_if = RESPONSE_STATUS == "Nonrespondent",
    increase_if = RESPONSE_STATUS == "Respondent",
    by = c("RACE_ETHNICITY")
  )
```

Multiple grouping variables may be supplied to the `by` argument. For
example, one can specify `by = c("STRATUM", "RACE_ETHNICITY")` to
redistribute weights separately by combinations of stratum and
race/ethnicity category.

### Propensity cell adjustment

The popular method of forming weighting classes based on estimated
response propensities (known as **propensity cell adjustment**) can also
be used, for example by adding a variable `PROPENSITY_CELL` to the data
and using `redistribute_weights(..., by = "PROPENSITY_CELL")`.

``` r
# Fit a response propensity model
response_propensity_model <- lou_vax_survey |>
  mutate(IS_RESPONDENT = ifelse(RESPONSE_STATUS == "Respondent", 1, 0)) |>
  svyglm(formula = IS_RESPONDENT ~ RACE_ETHNICITY + EDUC_ATTAINMENT,
         family = quasibinomial(link = 'logit'))

# Predict response propensities for individual cases
lou_vax_survey <- lou_vax_survey |>
  mutate(
    RESPONSE_PROPENSITY = predict(response_propensity_model,
                                  newdata = cur_svy(),
                                  type = "response")
  )

# Divide sample into propensity classes
lou_vax_survey <- lou_vax_survey |>
  mutate(PROPENSITY_CELL = ntile(x = RESPONSE_PROPENSITY, n = 5))

lou_vax_survey |>
    group_by(PROPENSITY_CELL) |>
    summarize(n = n(),
              min = min(RESPONSE_PROPENSITY),
              mean = mean(RESPONSE_PROPENSITY),
              max = max(RESPONSE_PROPENSITY))
#> # A tibble: 5 × 5
#>   PROPENSITY_CELL     n   min  mean   max
#>             <int> <int> <dbl> <dbl> <dbl>
#> 1               1   200 0.357 0.424 0.459
#> 2               2   200 0.459 0.484 0.488
#> 3               3   200 0.488 0.488 0.512
#> 4               4   200 0.512 0.551 0.564
#> 5               5   200 0.564 0.564 0.564

# Redistribute weights by propensity class
nr_adjusted_survey <- lou_vax_survey |>
  redistribute_weights(
    reduce_if = RESPONSE_STATUS == "Nonrespondent",
    increase_if = RESPONSE_STATUS == "Respondent",
    by = "PROPENSITY_CELL"
  )

# Inspect weights before adjustment

lou_vax_survey |>
  summarize_rep_weights(type = "specific",
                        by = c("PROPENSITY_CELL")) |>
  arrange(Rep_Column, PROPENSITY_CELL) |>
  select(PROPENSITY_CELL, Rep_Column,
         N_NONZERO, SUM) |>
  head(10)
#>    PROPENSITY_CELL Rep_Column N_NONZERO      SUM
#> 1                1          1       129 116473.4
#> 2                2          1       132 128419.3
#> 3                3          1       130 119459.9
#> 4                4          1       119 117668.0
#> 5                5          1       122 114681.5
#> 6                1          2       122 115876.1
#> 7                2          2       118 110500.4
#> 8                3          2       129 123641.0
#> 9                4          2       130 117070.7
#> 10               5          2       138 129613.9

# Inspect weights after adjustment
nr_adjusted_survey |>
  summarize_rep_weights(type = "specific",
                        by = c("PROPENSITY_CELL", "RESPONSE_STATUS")) |>
  arrange(Rep_Column, PROPENSITY_CELL, RESPONSE_STATUS) |>
  select(PROPENSITY_CELL, RESPONSE_STATUS, Rep_Column,
         N_NONZERO, SUM) |>
  head(10)
#>    PROPENSITY_CELL RESPONSE_STATUS Rep_Column N_NONZERO      SUM
#> 1                1   Nonrespondent          1         0      0.0
#> 2                1      Respondent          1        56 116473.4
#> 3                2   Nonrespondent          1         0      0.0
#> 4                2      Respondent          1        63 128419.3
#> 5                3   Nonrespondent          1         0      0.0
#> 6                3      Respondent          1        58 119459.9
#> 7                4   Nonrespondent          1         0      0.0
#> 8                4      Respondent          1        59 117668.0
#> 9                5   Nonrespondent          1         0      0.0
#> 10               5      Respondent          1        73 114681.5
```

## Saving the final weights to a data file

Once we’re satisfied with the weights, we can create a data frame with
the analysis variables and columns of replicate weights. This format is
easy to export to data files that can be loaded into R or other software
later.

``` r
data_frame_with_nr_adjusted_weights <- nr_adjusted_survey |>
  as_data_frame_with_weights(
    full_wgt_name = "NR_ADJ_WGT",
    rep_wgt_prefix = "NR_ADJ_REP_WGT_"
  )

# Preview first few column names
colnames(data_frame_with_nr_adjusted_weights) |> head(12)
#>  [1] "RESPONSE_STATUS"     "RACE_ETHNICITY"      "SEX"                
#>  [4] "EDUC_ATTAINMENT"     "VAX_STATUS"          "SAMPLING_WEIGHT"    
#>  [7] "RESPONSE_PROPENSITY" "PROPENSITY_CELL"     "NR_ADJ_WGT"         
#> [10] "NR_ADJ_REP_WGT_1"    "NR_ADJ_REP_WGT_2"    "NR_ADJ_REP_WGT_3"
```

``` r
# Write the data to a CSV file
write.csv(
  x = data_frame_with_nr_adjusted_weights,
  file = "survey-data-with-nonresponse-adjusted-weights.csv"
)
```

## Statistical background

The motivation for making this adjustment is that standard methods of
statistical inference assume that every person in the population has a
known, nonzero probability of participating in the survey (i.e. has a
nonzero chance of being sampled and a nonzero chance of responding if
they are sampled), denoted $p_{i,overall}$. Basic results in survey
sampling theory guarantee that if this assumption is true, we can
produce unbiased estimates of population means and totals by weighting
data from each respondent with the weight $1/p_{i,overall}$. Crucially,
the overall probability of participation $p_{i,overall}$ is the product
of two components: the probability that a person is sampled (denoted
$\pi_{i}$), and the probability that a person would respond to the
survey if they are sampled (denoted $p_{i}$ and referred to as the
**“response propensity”**). The sampling probability $\pi_{i}$ is known
since we can control the method of sampling, but the response propensity
$p_{i}$ is unknown and can only be estimated.

$$\begin{aligned}
w_{i}^{*} & {= 1/p_{i,overall}{\mspace{6mu}\text{(weights needed for unbiased estimation)}}} \\
p_{i,overall} & {= \pi_{i} \times p_{i}} \\
\pi_{i} & {= \textbf{𝐒𝐚𝐦𝐩𝐥𝐢𝐧𝐠 𝐩𝐫𝐨𝐛𝐚𝐛𝐢𝐥𝐢𝐭𝐲}} \\
 & {{\textit{𝑖.𝑒. 𝑡h𝑒 𝑝𝑟𝑜𝑏𝑎𝑏𝑖𝑙𝑖𝑡𝑦 𝑡h𝑎𝑡 𝑐𝑎𝑠𝑒}\mspace{6mu}}i{\mspace{6mu}\textit{𝑖𝑠 𝑟𝑎𝑛𝑑𝑜𝑚𝑙𝑦 𝑠𝑎𝑚𝑝𝑙𝑒𝑑}\mspace{6mu}}{\mspace{6mu}\text{(}}\textit{𝐾𝑛𝑜𝑤𝑛}\text{)}} \\
p_{i} & {= \textbf{𝐑𝐞𝐬𝐩𝐨𝐧𝐬𝐞 𝐩𝐫𝐨𝐩𝐞𝐧𝐬𝐢𝐭𝐲}} \\
 & {{\textit{𝑖.𝑒. 𝑡h𝑒 𝑝𝑟𝑜𝑏𝑎𝑏𝑖𝑙𝑖𝑡𝑦 𝑡h𝑎𝑡 𝑐𝑎𝑠𝑒}\mspace{6mu}}i{\mspace{6mu}\textit{𝑟𝑒𝑠𝑝𝑜𝑛𝑑𝑠, 𝑖𝑓 𝑠𝑎𝑚𝑝𝑙𝑒𝑑}\mspace{6mu}}{\mspace{6mu}\text{(}}\textit{𝑈𝑛𝑘𝑛𝑜𝑤𝑛}\text{)}} \\
 & 
\end{aligned}$$

The component $p_{i}$ must be estimated using data (with estimate
${\widehat{p}}_{i}$) and then nonresponse-adjusted weights for
respondents can be formed as
$w_{NR,i} = 1/\left( \pi_{i} \times {\widehat{p}}_{i} \right)$ and used
to obtain approximately unbiased estimates of population means and
totals. To use our earlier notation, the nonresponse adjustment factor
for respondents $f_{NR,i}$ is actually defined using
$1/{\widehat{p}}_{i}$.

$$\begin{aligned}
w_{i} & {= {\textit{𝑂𝑟𝑖𝑔𝑖𝑛𝑎𝑙 𝑠𝑎𝑚𝑝𝑙𝑖𝑛𝑔 𝑤𝑒𝑖𝑔h𝑡 𝑓𝑜𝑟 𝑐𝑎𝑠𝑒}\mspace{6mu}}i} \\
 & {= 1/\pi_{i},{\mspace{6mu}\textit{𝑤h𝑒𝑟𝑒}\mspace{6mu}}\pi_{i}{\mspace{6mu}\textit{𝑖𝑠 𝑡h𝑒 𝑝𝑟𝑜𝑏𝑎𝑏𝑖𝑙𝑖𝑡𝑦 𝑐𝑎𝑠𝑒 𝑖}\mspace{6mu}}\textit{h𝑎𝑑 𝑜𝑓 𝑏𝑒𝑖𝑛𝑔 𝑠𝑎𝑚𝑝𝑙𝑒𝑑}} \\
w_{NR,i} & {= w_{i} \times f_{NR,i} = {\textit{𝑊𝑒𝑖𝑔h𝑡 𝑓𝑜𝑟 𝑐𝑎𝑠𝑒}\mspace{6mu}}i{\mspace{6mu}\textit{𝑎𝑓𝑡𝑒𝑟 𝑛𝑜𝑛𝑟𝑒𝑠𝑝𝑜𝑛𝑠𝑒 𝑎𝑑𝑗𝑢𝑠𝑡𝑚𝑒𝑛𝑡}}} \\
 & \\
f_{NR,i} & {= \begin{cases}
0 & {{\text{if case}\mspace{6mu}}i{\mspace{6mu}\text{is a nonrespondent}}} \\
{1/{\widehat{p}}_{i}} & {{\text{if case}\mspace{6mu}}i{\mspace{6mu}\text{is a respondent}}} \\
 & 
\end{cases}} \\
{\widehat{p}}_{i} & {= \textbf{𝐄𝐬𝐭𝐢𝐦𝐚𝐭𝐞𝐝 𝐫𝐞𝐬𝐩𝐨𝐧𝐬𝐞 𝐩𝐫𝐨𝐩𝐞𝐧𝐬𝐢𝐭𝐲}}
\end{aligned}$$

In essence, different methods of nonresponse weighting adjustments vary
in terms of how they estimate ${\widehat{p}}_{i}$. The basic weight
redistribution method in effect estimates $p_{i}$ as constant across all
$i$, equal to the overall weighted response rate, and uses that to form
the weights. In other words, the basic weight redistribution essentially
is a way of forming an adjustment factor $f_{NR,i}$ based on the
estimated response propensity
${\widehat{p}}_{i} = \frac{\sum_{i \in s_{resp}}w_{i}}{\sum_{i \in s}w_{i}}$.

Weighting class adjustments and propensity cell adjustments are
essentially more refined ways of forming $f_{NR,i}$ by estimating
$p_{i}$ with a more realistic model, where $p_{i}$ is not constant
across the entire sample but instead varies among weighting classes or
propensity cells.

The reason for conducting weighting adjustments not only in the
full-sample weights but in the replicate weights is to account for the
nonresponse adjustment process when estimating sampling variances and
inferential statistics such as confidence intervals. Because of random
sampling, the adjustment factors used in the nonresponse adjustment will
vary from one sample to the next, and applying weighting adjustments
separately for each replicate reflects this variability. As we’ve seen
in this vignette, the
[`redistribute_weights()`](https://bschneidr.github.io/svrep/reference/redistribute_weights.md)
function handles this for us: after nonresponse adjustment, all of the
weight in each replicate is redistributed in the same manner that weight
was redistributed for the full-sample weights.

## Recommended Reading

- See Chapter 2, Section 2.7.3 of “Applied Survey Data Analysis” for a
  statistical explanation of the weighting adjustments described in this
  vignette.

> Heeringa, S., West, B., Berglund, P. (2017). Applied Survey Data
> Analysis, 2nd edition. Boca Raton, FL: CRC Press.

- Chapter 13 of “Practical Tools for Designing and Weighting Survey
  Samples” also provides an excellent overview of nonresponse adjustment
  methods.

> Valliant, R., Dever, J., Kreuter, F. (2018). Practical Tools for
> Designing and Weighting Survey Samples, 2nd edition. New York:
> Springer.
