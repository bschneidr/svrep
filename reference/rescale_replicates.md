# Rescale Replicate Factors

Rescale replicate factors. The main application of this rescaling is to
ensure that all replicate weights are strictly positive.

Note that this rescaling has no impact on variance estimates for totals
(or other linear statistics), but variance estimates for nonlinear
statistics will be affected by the rescaling.

## Usage

``` r
rescale_replicates(x, new_scale = NULL, min_wgt = 0.01, digits = 2)
```

## Arguments

- x:

  Either a replicate survey design object, or a numeric matrix of
  replicate weights. If `x` is a matrix, then it must have an attribute
  `scale` specifying the scale factor for variance estimation.

- new_scale:

  Either a single positive number, or `NULL`. If supplied, `new_scale`
  will be the new scale factor used for estimating variances. For
  example, with `B=10` bootstrap replicates, specifying `new_scale=0.2`
  will change the scale factor for variance estimation from
  \\B^{-1}=0.1\\ to \\0.2\\.  
  If `new_scale=NULL` or is left unspecified, then the argument
  `min_wgt` should be used instead.

- min_wgt:

  Should only be used if `new_scale=NULL` or `new_scale` is left
  unspecified. Specifies the minimum acceptable value for the rescaled
  weights, which will be used to automatically determine the new scale.
  Must be at least zero and must be less than one.

- digits:

  Only used if the argument `min_wgt` is used. Specifies the number of
  decimal places to use for the scale adjustment factor, so that the
  ratio `new_scale/orig_scale` will be a number with at most `digits`
  decimal places. For example, let `orig_scale` denote the original
  scale and `new_scale` denote the new scale. Then if the user specifies
  `digits=2`, the ratio `new_scale/orig_scale` will be a number such as
  0.25 or 0.10. This results in documentation that is easier to read.

## Value

If the input is a numeric matrix, returns the rescaled matrix. If the
input is a replicate survey design object, returns an updated replicate
survey design object.

For a replicate survey design object, results depend on whether the
object has a matrix of replicate factors rather than a matrix of
replicate weights (which are the product of replicate factors and
sampling weights). If the design object has `combined.weights=FALSE`,
then the replication factors are adjusted. If the design object has
`combined.weights=TRUE`, then the replicate weights are adjusted. It is
strongly recommended to only use the rescaling method for replication
factors rather than the weights.

For a replicate survey design object, the `scale` element of the design
object will be updated appropriately.

## Details

Let \\\mathbf{A} = \left\[ \mathbf{a}^{(1)} \cdots \mathbf{a}^{(b)}
\cdots \mathbf{a}^{(B)} \right\]\\ denote the \\(n \times B)\\ matrix of
replicate adjustment factors. The overall scale factor \\C\\ is used for
estimating variances with the formula \$\$v(\hat{y}) = C \sum\_{b=1}^{B}
(\hat{y}\_b - \hat{y})^2\$\$

The scale factor is changed from \\C\\ to \\C^{\prime}\\ by rescaling
the replicate factor matrix \\\mathbf{A}\\ using the transformation \$\$
1 + \sqrt{\frac{C}{C^{\prime}}}(\mathbf{A}-1) \$\$

## References

Rescaling was suggested by Fay (1989) for the specific application of
creating replicate factors using his generalized replication method.
This kind of rescaling is commony used in balanced repeated replication
to implement Fay's method of balanced repeated replication. Beaumont and
Patak (2012) provided an extended discussion of rescaling methods in the
context of rescaling generalized bootstrap replication factors to avoid
negative replicate weights.

\- Beaumont, Jean-François, and Zdenek Patak. 2012. "On the Generalized
Bootstrap for Sample Surveys with Special Attention to Poisson Sampling:
Generalized Bootstrap for Sample Surveys." International Statistical
Review 80 (1): 127-48.
https://doi.org/10.1111/j.1751-5823.2011.00166.x.  
  
- Fay, Robert. 1989. "Theory And Application Of Replicate Weighting For
Variance Calculations." In, 495-500. Alexandria, VA: American
Statistical Association.
http://www.asasrms.org/Proceedings/papers/1989_033.pdf

## Examples

``` r
# Example 1: Rescaling a matrix of replicate weights to avoid negative weights

 rep_wgts <- matrix(
   c(1.69742746694909, -0.230761178913411, 1.53333377634192,
     0.0495043413294782, 1.81820367441039, 1.13229198793703,
     1.62482013925955, 1.0866133494029, 0.28856654131668,
     0.581930729719006, 0.91827012312825, 1.49979905894482,
     1.26281337410693, 1.99327362761477, -0.25608700039304),
   nrow = 3, ncol = 5
 )
 attr(rep_wgts, 'scale') <- 1/5

 rescaled_wgts <- rescale_replicates(rep_wgts, min_wgt = 0.01)

 print(rep_wgts)
#>            [,1]       [,2]      [,3]      [,4]      [,5]
#> [1,]  1.6974275 0.04950434 1.6248201 0.5819307  1.262813
#> [2,] -0.2307612 1.81820367 1.0866133 0.9182701  1.993274
#> [3,]  1.5333338 1.13229199 0.2885665 1.4997991 -0.256087
#> attr(,"scale")
#> [1] 0.2
 print(rescaled_wgts)
#>            [,1]      [,2]      [,3]      [,4]       [,5]
#> [1,] 1.54964984 0.2509045 1.4924273 0.6705153 1.20712596
#> [2,] 0.03002431 1.6448348 1.0682609 0.9355878 1.78280928
#> [3,] 1.42032590 1.1042607 0.4393119 1.3938968 0.01006476
#> attr(,"scale")
#> [1] 0.322
 
 # Example 2: Rescale with a new, specific scale value
 
 rescaled_wgts <- rescale_replicates(rep_wgts, new_scale = 1/10)
 print(rescaled_wgts)
#>            [,1]       [,2]         [,3]      [,4]       [,5]
#> [1,]  1.9863114 -0.3442039  1.883629115 0.4087608  1.3716742
#> [2,] -0.7405592  2.1571147  1.122489773 0.8844165  2.4047010
#> [3,]  1.7542479  1.1870891 -0.006118846 1.7068226 -0.7763753
#> attr(,"scale")
#> [1] 0.1

 # Example 3: Rescaling replicate weights of a survey design object
 set.seed(2023)
 data('mu284', package = 'survey')

 ## First create a bootstrap design object
 svy_design_object <- svydesign(
   data = mu284,
   ids = ~ id1 + id2,
   fpc = ~  n1 +  n2
 )

 boot_design <- as_gen_boot_design(
   design = svy_design_object,
   variance_estimator = "Stratified Multistage SRS",
   replicates = 5
 )

 ## Rescale the weights
 rescaled_boot_design <- boot_design |>
   rescale_replicates(min_wgt = 0.01)

 boot_wgts <- weights(boot_design, "analysis")
 rescaled_boot_wgts <- weights(rescaled_boot_design, 'analysis')

 print(boot_wgts)
#>           REP_1      REP_2     REP_3       REP_4     REP_5
#>  [1,] 34.071074  -3.352195  7.031013  35.4547243 18.681422
#>  [2,] -3.271131  12.579037 57.474328   9.3992012 25.014379
#>  [3,] 12.204302  16.611770 14.029208   6.9869037 -8.727739
#>  [4,] 40.124053  62.587721 29.834150  31.6263954 10.057763
#>  [5,]  6.857688  48.936835  5.029175  42.1974204 67.126670
#>  [6,] 38.866284  -7.883877  6.363613  35.3323661 14.104502
#>  [7,] -2.705981   5.310800 51.191780 -18.8838184 34.232137
#>  [8,] 23.948409  19.740921 21.950039   0.8683187 -2.397135
#>  [9,] 38.102201  56.396305 39.516036  39.6713935 31.130900
#> [10,]  7.987330  41.986885  8.545987  47.8769538 66.314653
#> [11,] 35.747939 -13.746937  9.901870  41.9315736  8.610798
#> [12,]  1.384506   2.579634 50.469377 -26.8411850 19.800463
#> [13,] 22.153736  11.250766 19.117806   0.9281633 -1.226728
#> [14,] 48.183146  68.452257 28.322524  31.3003309 12.972211
#> [15,]  7.066647  63.713091 11.462660  41.8092990 64.604278
#> attr(,"tau")
#> [1] 1
#> attr(,"scale")
#> [1] 0.2
#> attr(,"rscales")
#> [1] 1 1 1 1 1
 print(rescaled_boot_wgts)
#>          REP_1     REP_2    REP_3      REP_4     REP_5
#>  [1,] 25.25163  6.792090 11.91375 25.9341321 17.660472
#>  [2,] 11.89944 19.717751 41.86297 18.1492541 25.851653
#>  [3,] 14.46554 16.639588 15.36570 11.8919917  4.140530
#>  [4,] 34.99383 46.074337 29.91819 30.8022417 20.163220
#>  [5,] 15.20650 35.962609 14.30456 32.6383013 44.934993
#>  [6,] 27.61693  4.556777 11.58454 25.8737772 15.402844
#>  [7,] 12.17821 16.132594 38.76401  4.1982693 30.398438
#>  [8,] 20.25849 18.183084 19.27276  8.8739163  7.263187
#>  [9,] 33.99652 43.020337 34.69391 34.7705466 30.557832
#> [10,] 15.76371 32.534452 16.03927 35.4398085 44.534455
#> [11,] 26.07876  1.664742 13.32984 29.1289262 12.692999
#> [12,] 14.19590 14.785410 38.40768  0.2731898 23.279818
#> [13,] 19.37324 13.995200 17.87572  8.9034355  7.840506
#> [14,] 38.96908 48.967100 29.17256 30.6414059 21.600811
#> [15,] 15.30957 43.251199 17.47796 32.4468551 43.690789
#> attr(,"scale")
#> [1] 0.822
#> attr(,"rscales")
#> [1] 1 1 1 1 1
```
