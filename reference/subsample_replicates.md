# Subsample Replicates in a Replicate Design Object

Randomly subsamples the replicates of a survey design object, to keep
only a subset. The scale factor used in estimation is increased to
account for the subsampling.

## Usage

``` r
subsample_replicates(design, n_reps)
```

## Arguments

- design:

  A survey design object, created with either the `survey` or `srvyr`
  packages.

- n_reps:

  The number of replicates to keep after subsampling

## Value

An updated survey design object, where only a random selection of the
replicates has been retained. The overall 'scale' factor for the design
(accessed with `design$scale`) is increased to account for the sampling
of replicates.

## Statistical Details

Suppose the initial replicate design has \\L\\ replicates, with
respective constants \\c_k\\ for \\k=1,\dots,L\\ used to estimate
variance with the formula \$\$v\_{R} = \sum\_{k=1}^L
c_k\left(\hat{T}\_y^{(k)}-\hat{T}\_y\right)^2\$\$

With subsampling of replicates, \\L_0\\ of the original \\L\\ replicates
are randomly selected, and then variances are estimated using the
formula: \$\$v\_{R} = \frac{L}{L_0} \sum\_{k=1}^{L_0}
c_k\left(\hat{T}\_y^{(k)}-\hat{T}\_y\right)^2\$\$

This subsampling is suggested for certain replicate designs in Fay
(1989). Kim and Wu (2013) provide a detailed theoretical justification
and also propose alternative methods of subsampling replicates.

## References

Fay, Robert. 1989. "Theory And Application Of Replicate Weighting For
Variance Calculations." In, 495-500. Alexandria, VA: American
Statistical Association.
http://www.asasrms.org/Proceedings/papers/1989_033.pdf

Kim, J.K. and Wu, C. 2013. "Sparse and Efficient Replication Variance
Estimation for Complex Surveys." **Survey Methodology**, Statistics
Canada, 39(1), 91-120.

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

# Inspect replicates before subsampling

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

# Inspect replicates after subsampling

  rep_design |>
    subsample_replicates(n_reps = 4) |>
    getElement("repweights")
#>          REP_5     REP_1     REP_7     REP_8
#> [1,] 0.6464466 1.3535534 0.6464466 1.3535534
#> [2,] 1.4511845 1.2154822 1.4511845 0.5488155
#> [3,] 0.6625093 0.4268071 0.6625093 0.7601404
#> [4,] 1.2398596 1.0041573 1.2398596 1.3374907
#> [5,] 0.6464466 1.3535534 0.6464466 0.6464466
#> [6,] 1.4511845 0.5488155 0.7845178 1.4511845
#> [7,] 0.6625093 1.3374907 1.5731929 1.2398596
#> [8,] 1.2398596 0.7601404 0.9958427 0.6625093
#> attr(,"scale")
#> [1] 4
#> attr(,"rscales")
#> [1] 1 1 1 1
```
