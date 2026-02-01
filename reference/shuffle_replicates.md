# Shuffle Order of Replicates in a Replicate Design Object

Shuffle the order of replicates in a survey design object. In other
words, the order of the columns of replicate weights is randomly
permuted.

## Usage

``` r
shuffle_replicates(design)
```

## Arguments

- design:

  A survey design object, created with either the `survey` or `srvyr`
  packages.

## Value

An updated survey design object, where the order of the replicates has
been shuffled (i.e., the order has been randomly permuted).

## Examples

``` r
library(survey)
set.seed(2023)

# Create an example survey design object

  sample_data <- data.frame(
    STRATUM = c(1,1,1,1,2,2,2,2),
    PSU     = c(1,2,3,4,5,6,7,8)
  )

  survey_design <- svydesign(
    data = sample_data,
    strata = ~ STRATUM,
    ids = ~ PSU,
    weights = ~ 1
  )

  rep_design <- survey_design |>
    as_fays_gen_rep_design(variance_estimator = "Ultimate Cluster")

# Inspect replicates before shuffling

  rep_design |> getElement("repweights")
#>          REP_1     REP_2     REP_3     REP_4     REP_5     REP_6     REP_7
#> [1,] 1.3535534 0.6464466 1.3535534 0.6464466 0.6464466 1.3535534 0.6464466
#> [2,] 1.2154822 0.7845178 1.2154822 0.7845178 1.4511845 0.5488155 1.4511845
#> [3,] 0.4268071 0.9958427 0.4268071 0.9958427 0.6625093 0.7601404 0.6625093
#> [4,] 1.0041573 1.5731929 1.0041573 1.5731929 1.2398596 1.3374907 1.2398596
#> [5,] 1.3535534 1.3535534 1.3535534 1.3535534 0.6464466 0.6464466 0.6464466
#> [6,] 0.5488155 1.2154822 1.2154822 0.5488155 1.4511845 0.7845178 0.7845178
#> [7,] 1.3374907 1.0041573 0.4268071 0.7601404 0.6625093 0.9958427 1.5731929
#> [8,] 0.7601404 0.4268071 1.0041573 1.3374907 1.2398596 1.5731929 0.9958427
#>          REP_8
#> [1,] 1.3535534
#> [2,] 0.5488155
#> [3,] 0.7601404
#> [4,] 1.3374907
#> [5,] 0.6464466
#> [6,] 1.4511845
#> [7,] 1.2398596
#> [8,] 0.6625093
#> attr(,"scale")
#> [1] 1
#> attr(,"rscales")
#> [1] 1 1 1 1 1 1 1 1

# Inspect replicates after shuffling

  rep_design |>
    shuffle_replicates() |>
    getElement("repweights")
#>          REP_5     REP_1     REP_7     REP_8     REP_6     REP_3     REP_2
#> [1,] 0.6464466 1.3535534 0.6464466 1.3535534 1.3535534 1.3535534 0.6464466
#> [2,] 1.4511845 1.2154822 1.4511845 0.5488155 0.5488155 1.2154822 0.7845178
#> [3,] 0.6625093 0.4268071 0.6625093 0.7601404 0.7601404 0.4268071 0.9958427
#> [4,] 1.2398596 1.0041573 1.2398596 1.3374907 1.3374907 1.0041573 1.5731929
#> [5,] 0.6464466 1.3535534 0.6464466 0.6464466 0.6464466 1.3535534 1.3535534
#> [6,] 1.4511845 0.5488155 0.7845178 1.4511845 0.7845178 1.2154822 1.2154822
#> [7,] 0.6625093 1.3374907 1.5731929 1.2398596 0.9958427 0.4268071 1.0041573
#> [8,] 1.2398596 0.7601404 0.9958427 0.6625093 1.5731929 1.0041573 0.4268071
#>          REP_4
#> [1,] 0.6464466
#> [2,] 0.7845178
#> [3,] 0.9958427
#> [4,] 1.5731929
#> [5,] 1.3535534
#> [6,] 0.5488155
#> [7,] 0.7601404
#> [8,] 1.3374907
#> attr(,"scale")
#> [1] 1
#> attr(,"rscales")
#> [1] 1 1 1 1 1 1 1 1
```
