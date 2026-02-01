# Quadratic Form Matrix of Variance Estimator for a Survey Design

Determines the quadratic form matrix of a specified variance estimator,
by parsing the information stored in a survey design object created
using the 'survey' package.

## Usage

``` r
get_design_quad_form(
  design,
  variance_estimator,
  ensure_psd = FALSE,
  aux_var_names = NULL
)
```

## Arguments

- design:

  A survey design object created using the 'survey' (or 'srvyr')
  package, with class `'survey.design'` or `'svyimputationList'`. Also
  accepts two-phase design objects with class `'twophase2'`; see the
  section below titled "Two-Phase Designs" for more information about
  handling of two-phase designs.

- variance_estimator:

  The name of the variance estimator whose quadratic form matrix should
  be created.  
  See the section "Variance Estimators" below. Options include:

  - **"Yates-Grundy"**:  
    The Yates-Grundy variance estimator based on first-order and
    second-order inclusion probabilities.

  - **"Horvitz-Thompson"**:  
    The Horvitz-Thompson variance estimator based on first-order and
    second-order inclusion probabilities.

  - **"Poisson Horvitz-Thompson"**:  
    The Horvitz-Thompson variance estimator based on assuming Poisson
    sampling, with first-order inclusion probabilities inferred from the
    sampling probabilities of the survey design object.

  - **"Stratified Multistage SRS"**:  
    The usual stratified multistage variance estimator based on
    estimating the variance of cluster totals within strata at each
    stage.

  - **"Ultimate Cluster"**:  
    The usual variance estimator based on estimating the variance of
    first-stage cluster totals within first-stage strata.

  - **"Deville-1"**:  
    A variance estimator for unequal-probability sampling without
    replacement, described in Matei and Tillé (2005) as "Deville 1".

  - **"Deville-2"**:  
    A variance estimator for unequal-probability sampling without
    replacement, described in Matei and Tillé (2005) as "Deville 2".

  - **"Deville-Tille":**  
    A variance estimator useful for balanced sampling designs, proposed
    by Deville and Tillé (2005).

  - **"SD1"**:  
    The non-circular successive-differences variance estimator described
    by Ash (2014), sometimes used for variance estimation for systematic
    sampling.

  - **"SD2"**:  
    The circular successive-differences variance estimator described by
    Ash (2014). This estimator is the basis of the
    "successive-differences replication" estimator commonly used for
    variance estimation for systematic sampling.

  - **"Beaumont-Emond"**:  
    The variance estimator of Beaumont and Emond (2022) for multistage
    unequal-probability sampling without replacement.

  - **"BOSB"**:  
    The kernel-based variance estimator proposed by Breidt, Opsomer, and
    Sanchez-Borrego (2016) for use with systematic samples or other
    finely stratified designs. Uses the Epanechnikov kernel with the
    bandwidth automatically chosen to result in the smallest possible
    nonempty kernel window.

- ensure_psd:

  If `TRUE` (the default), ensures that the result is a positive
  semidefinite matrix. This is necessary if the quadratic form is used
  as an input for replication methods such as the generalized bootstrap.
  For mathematical details, please see the documentation for the
  function
  [`get_nearest_psd_matrix()`](https://bschneidr.github.io/svrep/reference/get_nearest_psd_matrix.md).
  The approximation method is discussed by Beaumont and Patak (2012) in
  the context of forming replicate weights for two-phase samples. The
  authors argue that this approximation should lead to only a small
  overestimation of variance.

- aux_var_names:

  Only required if `variance_estimator = "Deville-Tille"` or if
  `variance_estimator = "BOSB"`. For the Deville-Tillé estimator, this
  should be a character vector of variable names for auxiliary variables
  to be used in the variance estimator. For the BOSB estimator, this
  should be a string giving a single variable name to use as an
  auxiliary variable in the kernel-based variance estimator of Breidt,
  Opsomer, and Sanchez-Borrego (2016).

## Value

A matrix representing the quadratic form of a specified variance
estimator, based on extracting information about clustering,
stratification, and selection probabilities from the survey design
object.

## Variance Estimators

See
[variance-estimators](https://bschneidr.github.io/svrep/reference/variance-estimators.md)
for a description of each variance estimator.

## Two-Phase Designs

For a two-phase design, `variance_estimator` should be a list of
variance estimators' names, with two elements, such as
`list('Ultimate Cluster', 'Poisson Horvitz-Thompson')`. In two-phase
designs, only the following estimators may be used for the second phase:

- "Ultimate Cluster"

- "Stratified Multistage SRS"

- "Poisson Horvitz-Thompson"

For statistical details on the handling of two-phase designs, see the
documentation for
[make_twophase_quad_form](https://bschneidr.github.io/svrep/reference/make_twophase_quad_form.md).

## References

\- Ash, S. (2014). "*Using successive difference replication for
estimating variances*." **Survey Methodology**, Statistics Canada,
40(1), 47-59.  
  
- Beaumont, Jean-François, and Zdenek Patak. (2012). "*On the
Generalized Bootstrap for Sample Surveys with Special Attention to
Poisson Sampling: Generalized Bootstrap for Sample Surveys.*"
**International Statistical Review** 80 (1): 127-48.  
  
- Bellhouse, D.R. (1985). "*Computing Methods for Variance Estimation in
Complex Surveys*." **Journal of Official Statistics**, Vol.1, No.3.  
  
- Breidt, F. J., Opsomer, J. D., & Sanchez-Borrego, I. (2016).
"*Nonparametric Variance Estimation Under Fine Stratification: An
Alternative to Collapsed Strata*." **Journal of the American Statistical
Association**, 111(514), 822-833.
https://doi.org/10.1080/01621459.2015.1058264  
  
- Deville, J.C., and Tillé, Y. (2005). "*Variance approximation under
balanced sampling.*" **Journal of Statistical Planning and Inference**,
128, 569-591.  
  
- Särndal, C.E., Swensson, B., & Wretman, J. (1992). "*Model Assisted
Survey Sampling*." Springer New York.

## Examples

``` r
# \donttest{

# Example 1: Quadratic form for successive-difference variance estimator ----

   data('library_stsys_sample', package = 'svrep')

   ## First, ensure data are sorted in same order as was used in sampling
   library_stsys_sample <- library_stsys_sample |>
     sort_by(~ SAMPLING_SORT_ORDER)

   ## Create a survey design object
   design_obj <- svydesign(
     data = library_stsys_sample,
     strata = ~ SAMPLING_STRATUM,
     ids = ~ 1,
     fpc = ~ STRATUM_POP_SIZE
   )

   ## Obtain quadratic form
   quad_form_matrix <- get_design_quad_form(
     design = design_obj,
     variance_estimator = "SD2"
   )
#> For `variance_estimator='SD2', assumes rows of data are sorted in the same order used in sampling.

   ## Estimate variance of estimated population total
   y <- design_obj$variables$LIBRARIA
   wts <- weights(design_obj, type = 'sampling')
   y_wtd <- as.matrix(y) * wts
   y_wtd[is.na(y_wtd)] <- 0

   pop_total <- sum(y_wtd)

   var_est <- t(y_wtd) %*% quad_form_matrix %*% y_wtd
   std_error <- sqrt(var_est)

   print(pop_total); print(std_error)
#> [1] 65641.55
#> 1 x 1 Matrix of class "dgeMatrix"
#>          [,1]
#> [1,] 13972.04

   # Compare to estimate from assuming SRS
   svytotal(x = ~ LIBRARIA, na.rm = TRUE,
            design = design_obj)
#>          total    SE
#> LIBRARIA 65642 14008

# Example 2: Two-phase design (second phase is nonresponse) ----

  ## Estimate response propensities, separately by stratum
  library_stsys_sample[['RESPONSE_PROB']] <- svyglm(
    design = design_obj,
    formula = I(RESPONSE_STATUS == "Survey Respondent") ~ SAMPLING_STRATUM,
    family = quasibinomial(link = 'logit')
  ) |> predict(type = 'response')

  ## Create a survey design object,
  ## where nonresponse is treated as a second phase of sampling
  twophase_design <- twophase(
    data = library_stsys_sample,
    strata = list(~ SAMPLING_STRATUM, NULL),
    id = list(~ 1, ~ 1),
    fpc = list(~ STRATUM_POP_SIZE, NULL),
    probs = list(NULL, ~ RESPONSE_PROB),
    subset = ~ I(RESPONSE_STATUS == "Survey Respondent")
  )

  ## Obtain quadratic form for the two-phase variance estimator,
  ## where first phase variance contribution estimated
  ## using the successive differences estimator
  ## and second phase variance contribution estimated
  ## using the Horvitz-Thompson estimator
  ## (with joint probabilities based on assumption of Poisson sampling)
  get_design_quad_form(
    design = twophase_design,
    variance_estimator = list(
      "SD2",
      "Poisson Horvitz-Thompson"
    )
  )
#> For `variance_estimator='SD2', assumes rows of data are sorted in the same order used in sampling.
#> 211 x 211 sparse Matrix of class "dgCMatrix"
#>                                                                         
#>   [1,]  0.9753086 -0.9753086  .          .          .          .        
#>   [2,] -0.9753086  0.9753086  .          .          .          .        
#>   [3,]  .          .          0.9777778 -0.4888889  .          .        
#>   [4,]  .          .         -0.4888889  0.9777778 -0.4888889  .        
#>   [5,]  .          .          .         -0.4888889  0.9777778 -0.4888889
#>   [6,]  .          .          .          .         -0.4888889  0.9777778
#>   [7,]  .          .         -0.4888889  .          .         -0.4888889
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                          
#>   [1,]  .          .          .         .            .          .        
#>   [2,]  .          .          .         .            .          .        
#>   [3,] -0.4888889  .          .         .            .          .        
#>   [4,]  .          .          .         .            .          .        
#>   [5,]  .          .          .         .            .          .        
#>   [6,] -0.4888889  .          .         .            .          .        
#>   [7,]  0.9777778  .          .         .            .          .        
#>   [8,]  .          0.9661017 -0.9661017 .            .          .        
#>   [9,]  .         -0.9661017  0.9661017 .            .          .        
#>  [10,]  .          .          .         6.18486e-10  .          .        
#>  [11,]  .          .          .         .            0.9775281 -0.9775281
#>  [12,]  .          .          .         .           -0.9775281  0.9775281
#>  [13,]  .          .          .         .            .          .        
#>  [14,]  .          .          .         .            .          .        
#>  [15,]  .          .          .         .            .          .        
#>  [16,]  .          .          .         .            .          .        
#>  [17,]  .          .          .         .            .          .        
#>  [18,]  .          .          .         .            .          .        
#>  [19,]  .          .          .         .            .          .        
#>  [20,]  .          .          .         .            .          .        
#>  [21,]  .          .          .         .            .          .        
#>  [22,]  .          .          .         .            .          .        
#>  [23,]  .          .          .         .            .          .        
#>  [24,]  .          .          .         .            .          .        
#>  [25,]  .          .          .         .            .          .        
#>  [26,]  .          .          .         .            .          .        
#>  [27,]  .          .          .         .            .          .        
#>  [28,]  .          .          .         .            .          .        
#>  [29,]  .          .          .         .            .          .        
#>  [30,]  .          .          .         .            .          .        
#>  [31,]  .          .          .         .            .          .        
#>  [32,]  .          .          .         .            .          .        
#>  [33,]  .          .          .         .            .          .        
#>  [34,]  .          .          .         .            .          .        
#>  [35,]  .          .          .         .            .          .        
#>  [36,]  .          .          .         .            .          .        
#>  [37,]  .          .          .         .            .          .        
#>  [38,]  .          .          .         .            .          .        
#>  [39,]  .          .          .         .            .          .        
#>  [40,]  .          .          .         .            .          .        
#>  [41,]  .          .          .         .            .          .        
#>  [42,]  .          .          .         .            .          .        
#>  [43,]  .          .          .         .            .          .        
#>  [44,]  .          .          .         .            .          .        
#>  [45,]  .          .          .         .            .          .        
#>  [46,]  .          .          .         .            .          .        
#>  [47,]  .          .          .         .            .          .        
#>  [48,]  .          .          .         .            .          .        
#>  [49,]  .          .          .         .            .          .        
#>  [50,]  .          .          .         .            .          .        
#>  [51,]  .          .          .         .            .          .        
#>  [52,]  .          .          .         .            .          .        
#>  [53,]  .          .          .         .            .          .        
#>  [54,]  .          .          .         .            .          .        
#>  [55,]  .          .          .         .            .          .        
#>  [56,]  .          .          .         .            .          .        
#>  [57,]  .          .          .         .            .          .        
#>  [58,]  .          .          .         .            .          .        
#>  [59,]  .          .          .         .            .          .        
#>  [60,]  .          .          .         .            .          .        
#>  [61,]  .          .          .         .            .          .        
#>  [62,]  .          .          .         .            .          .        
#>  [63,]  .          .          .         .            .          .        
#>  [64,]  .          .          .         .            .          .        
#>  [65,]  .          .          .         .            .          .        
#>  [66,]  .          .          .         .            .          .        
#>  [67,]  .          .          .         .            .          .        
#>  [68,]  .          .          .         .            .          .        
#>  [69,]  .          .          .         .            .          .        
#>  [70,]  .          .          .         .            .          .        
#>  [71,]  .          .          .         .            .          .        
#>  [72,]  .          .          .         .            .          .        
#>  [73,]  .          .          .         .            .          .        
#>  [74,]  .          .          .         .            .          .        
#>  [75,]  .          .          .         .            .          .        
#>  [76,]  .          .          .         .            .          .        
#>  [77,]  .          .          .         .            .          .        
#>  [78,]  .          .          .         .            .          .        
#>  [79,]  .          .          .         .            .          .        
#>  [80,]  .          .          .         .            .          .        
#>  [81,]  .          .          .         .            .          .        
#>  [82,]  .          .          .         .            .          .        
#>  [83,]  .          .          .         .            .          .        
#>  [84,]  .          .          .         .            .          .        
#>  [85,]  .          .          .         .            .          .        
#>  [86,]  .          .          .         .            .          .        
#>  [87,]  .          .          .         .            .          .        
#>  [88,]  .          .          .         .            .          .        
#>  [89,]  .          .          .         .            .          .        
#>  [90,]  .          .          .         .            .          .        
#>  [91,]  .          .          .         .            .          .        
#>  [92,]  .          .          .         .            .          .        
#>  [93,]  .          .          .         .            .          .        
#>  [94,]  .          .          .         .            .          .        
#>  [95,]  .          .          .         .            .          .        
#>  [96,]  .          .          .         .            .          .        
#>  [97,]  .          .          .         .            .          .        
#>  [98,]  .          .          .         .            .          .        
#>  [99,]  .          .          .         .            .          .        
#> [100,]  .          .          .         .            .          .        
#> [101,]  .          .          .         .            .          .        
#> [102,]  .          .          .         .            .          .        
#> [103,]  .          .          .         .            .          .        
#> [104,]  .          .          .         .            .          .        
#> [105,]  .          .          .         .            .          .        
#> [106,]  .          .          .         .            .          .        
#> [107,]  .          .          .         .            .          .        
#> [108,]  .          .          .         .            .          .        
#> [109,]  .          .          .         .            .          .        
#> [110,]  .          .          .         .            .          .        
#> [111,]  .          .          .         .            .          .        
#> [112,]  .          .          .         .            .          .        
#> [113,]  .          .          .         .            .          .        
#> [114,]  .          .          .         .            .          .        
#> [115,]  .          .          .         .            .          .        
#> [116,]  .          .          .         .            .          .        
#> [117,]  .          .          .         .            .          .        
#> [118,]  .          .          .         .            .          .        
#> [119,]  .          .          .         .            .          .        
#> [120,]  .          .          .         .            .          .        
#> [121,]  .          .          .         .            .          .        
#> [122,]  .          .          .         .            .          .        
#> [123,]  .          .          .         .            .          .        
#> [124,]  .          .          .         .            .          .        
#> [125,]  .          .          .         .            .          .        
#> [126,]  .          .          .         .            .          .        
#> [127,]  .          .          .         .            .          .        
#> [128,]  .          .          .         .            .          .        
#> [129,]  .          .          .         .            .          .        
#> [130,]  .          .          .         .            .          .        
#> [131,]  .          .          .         .            .          .        
#> [132,]  .          .          .         .            .          .        
#> [133,]  .          .          .         .            .          .        
#> [134,]  .          .          .         .            .          .        
#> [135,]  .          .          .         .            .          .        
#> [136,]  .          .          .         .            .          .        
#> [137,]  .          .          .         .            .          .        
#> [138,]  .          .          .         .            .          .        
#> [139,]  .          .          .         .            .          .        
#> [140,]  .          .          .         .            .          .        
#> [141,]  .          .          .         .            .          .        
#> [142,]  .          .          .         .            .          .        
#> [143,]  .          .          .         .            .          .        
#> [144,]  .          .          .         .            .          .        
#> [145,]  .          .          .         .            .          .        
#> [146,]  .          .          .         .            .          .        
#> [147,]  .          .          .         .            .          .        
#> [148,]  .          .          .         .            .          .        
#> [149,]  .          .          .         .            .          .        
#> [150,]  .          .          .         .            .          .        
#> [151,]  .          .          .         .            .          .        
#> [152,]  .          .          .         .            .          .        
#> [153,]  .          .          .         .            .          .        
#> [154,]  .          .          .         .            .          .        
#> [155,]  .          .          .         .            .          .        
#> [156,]  .          .          .         .            .          .        
#> [157,]  .          .          .         .            .          .        
#> [158,]  .          .          .         .            .          .        
#> [159,]  .          .          .         .            .          .        
#> [160,]  .          .          .         .            .          .        
#> [161,]  .          .          .         .            .          .        
#> [162,]  .          .          .         .            .          .        
#> [163,]  .          .          .         .            .          .        
#> [164,]  .          .          .         .            .          .        
#> [165,]  .          .          .         .            .          .        
#> [166,]  .          .          .         .            .          .        
#> [167,]  .          .          .         .            .          .        
#> [168,]  .          .          .         .            .          .        
#> [169,]  .          .          .         .            .          .        
#> [170,]  .          .          .         .            .          .        
#> [171,]  .          .          .         .            .          .        
#> [172,]  .          .          .         .            .          .        
#> [173,]  .          .          .         .            .          .        
#> [174,]  .          .          .         .            .          .        
#> [175,]  .          .          .         .            .          .        
#> [176,]  .          .          .         .            .          .        
#> [177,]  .          .          .         .            .          .        
#> [178,]  .          .          .         .            .          .        
#> [179,]  .          .          .         .            .          .        
#> [180,]  .          .          .         .            .          .        
#> [181,]  .          .          .         .            .          .        
#> [182,]  .          .          .         .            .          .        
#> [183,]  .          .          .         .            .          .        
#> [184,]  .          .          .         .            .          .        
#> [185,]  .          .          .         .            .          .        
#> [186,]  .          .          .         .            .          .        
#> [187,]  .          .          .         .            .          .        
#> [188,]  .          .          .         .            .          .        
#> [189,]  .          .          .         .            .          .        
#> [190,]  .          .          .         .            .          .        
#> [191,]  .          .          .         .            .          .        
#> [192,]  .          .          .         .            .          .        
#> [193,]  .          .          .         .            .          .        
#> [194,]  .          .          .         .            .          .        
#> [195,]  .          .          .         .            .          .        
#> [196,]  .          .          .         .            .          .        
#> [197,]  .          .          .         .            .          .        
#> [198,]  .          .          .         .            .          .        
#> [199,]  .          .          .         .            .          .        
#> [200,]  .          .          .         .            .          .        
#> [201,]  .          .          .         .            .          .        
#> [202,]  .          .          .         .            .          .        
#> [203,]  .          .          .         .            .          .        
#> [204,]  .          .          .         .            .          .        
#> [205,]  .          .          .         .            .          .        
#> [206,]  .          .          .         .            .          .        
#> [207,]  .          .          .         .            .          .        
#> [208,]  .          .          .         .            .          .        
#> [209,]  .          .          .         .            .          .        
#> [210,]  .          .          .         .            .          .        
#> [211,]  .          .          .         .            .          .        
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  0.9817352 -0.4885845  .          .          .          .        
#>  [14,] -0.4885845  0.9817352 -0.4885845  .          .          .        
#>  [15,]  .         -0.4885845  0.9817352 -0.4885845  .          .        
#>  [16,]  .          .         -0.4885845  0.9817352  .          .        
#>  [17,]  .          .          .          .          0.9734513 -0.4867257
#>  [18,]  .          .          .          .         -0.4867257  0.9734513
#>  [19,]  .          .          .          .         -0.4867257 -0.4867257
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                          
#>   [1,]  .          .          .          .          .         .          
#>   [2,]  .          .          .          .          .         .          
#>   [3,]  .          .          .          .          .         .          
#>   [4,]  .          .          .          .          .         .          
#>   [5,]  .          .          .          .          .         .          
#>   [6,]  .          .          .          .          .         .          
#>   [7,]  .          .          .          .          .         .          
#>   [8,]  .          .          .          .          .         .          
#>   [9,]  .          .          .          .          .         .          
#>  [10,]  .          .          .          .          .         .          
#>  [11,]  .          .          .          .          .         .          
#>  [12,]  .          .          .          .          .         .          
#>  [13,]  .          .          .          .          .         .          
#>  [14,]  .          .          .          .          .         .          
#>  [15,]  .          .          .          .          .         .          
#>  [16,]  .          .          .          .          .         .          
#>  [17,] -0.4867257  .          .          .          .         .          
#>  [18,] -0.4867257  .          .          .          .         .          
#>  [19,]  0.9734513  .          .          .          .         .          
#>  [20,]  .          0.9791667 -0.4895833  .         -0.4895833 .          
#>  [21,]  .         -0.4895833  0.9791667 -0.4895833  .         .          
#>  [22,]  .          .         -0.4895833  0.9791667 -0.4895833 .          
#>  [23,]  .         -0.4895833  .         -0.4895833  0.9791667 .          
#>  [24,]  .          .          .          .          .         6.18486e-10
#>  [25,]  .          .          .          .          .         .          
#>  [26,]  .          .          .          .          .         .          
#>  [27,]  .          .          .          .          .         .          
#>  [28,]  .          .          .          .          .         .          
#>  [29,]  .          .          .          .          .         .          
#>  [30,]  .          .          .          .          .         .          
#>  [31,]  .          .          .          .          .         .          
#>  [32,]  .          .          .          .          .         .          
#>  [33,]  .          .          .          .          .         .          
#>  [34,]  .          .          .          .          .         .          
#>  [35,]  .          .          .          .          .         .          
#>  [36,]  .          .          .          .          .         .          
#>  [37,]  .          .          .          .          .         .          
#>  [38,]  .          .          .          .          .         .          
#>  [39,]  .          .          .          .          .         .          
#>  [40,]  .          .          .          .          .         .          
#>  [41,]  .          .          .          .          .         .          
#>  [42,]  .          .          .          .          .         .          
#>  [43,]  .          .          .          .          .         .          
#>  [44,]  .          .          .          .          .         .          
#>  [45,]  .          .          .          .          .         .          
#>  [46,]  .          .          .          .          .         .          
#>  [47,]  .          .          .          .          .         .          
#>  [48,]  .          .          .          .          .         .          
#>  [49,]  .          .          .          .          .         .          
#>  [50,]  .          .          .          .          .         .          
#>  [51,]  .          .          .          .          .         .          
#>  [52,]  .          .          .          .          .         .          
#>  [53,]  .          .          .          .          .         .          
#>  [54,]  .          .          .          .          .         .          
#>  [55,]  .          .          .          .          .         .          
#>  [56,]  .          .          .          .          .         .          
#>  [57,]  .          .          .          .          .         .          
#>  [58,]  .          .          .          .          .         .          
#>  [59,]  .          .          .          .          .         .          
#>  [60,]  .          .          .          .          .         .          
#>  [61,]  .          .          .          .          .         .          
#>  [62,]  .          .          .          .          .         .          
#>  [63,]  .          .          .          .          .         .          
#>  [64,]  .          .          .          .          .         .          
#>  [65,]  .          .          .          .          .         .          
#>  [66,]  .          .          .          .          .         .          
#>  [67,]  .          .          .          .          .         .          
#>  [68,]  .          .          .          .          .         .          
#>  [69,]  .          .          .          .          .         .          
#>  [70,]  .          .          .          .          .         .          
#>  [71,]  .          .          .          .          .         .          
#>  [72,]  .          .          .          .          .         .          
#>  [73,]  .          .          .          .          .         .          
#>  [74,]  .          .          .          .          .         .          
#>  [75,]  .          .          .          .          .         .          
#>  [76,]  .          .          .          .          .         .          
#>  [77,]  .          .          .          .          .         .          
#>  [78,]  .          .          .          .          .         .          
#>  [79,]  .          .          .          .          .         .          
#>  [80,]  .          .          .          .          .         .          
#>  [81,]  .          .          .          .          .         .          
#>  [82,]  .          .          .          .          .         .          
#>  [83,]  .          .          .          .          .         .          
#>  [84,]  .          .          .          .          .         .          
#>  [85,]  .          .          .          .          .         .          
#>  [86,]  .          .          .          .          .         .          
#>  [87,]  .          .          .          .          .         .          
#>  [88,]  .          .          .          .          .         .          
#>  [89,]  .          .          .          .          .         .          
#>  [90,]  .          .          .          .          .         .          
#>  [91,]  .          .          .          .          .         .          
#>  [92,]  .          .          .          .          .         .          
#>  [93,]  .          .          .          .          .         .          
#>  [94,]  .          .          .          .          .         .          
#>  [95,]  .          .          .          .          .         .          
#>  [96,]  .          .          .          .          .         .          
#>  [97,]  .          .          .          .          .         .          
#>  [98,]  .          .          .          .          .         .          
#>  [99,]  .          .          .          .          .         .          
#> [100,]  .          .          .          .          .         .          
#> [101,]  .          .          .          .          .         .          
#> [102,]  .          .          .          .          .         .          
#> [103,]  .          .          .          .          .         .          
#> [104,]  .          .          .          .          .         .          
#> [105,]  .          .          .          .          .         .          
#> [106,]  .          .          .          .          .         .          
#> [107,]  .          .          .          .          .         .          
#> [108,]  .          .          .          .          .         .          
#> [109,]  .          .          .          .          .         .          
#> [110,]  .          .          .          .          .         .          
#> [111,]  .          .          .          .          .         .          
#> [112,]  .          .          .          .          .         .          
#> [113,]  .          .          .          .          .         .          
#> [114,]  .          .          .          .          .         .          
#> [115,]  .          .          .          .          .         .          
#> [116,]  .          .          .          .          .         .          
#> [117,]  .          .          .          .          .         .          
#> [118,]  .          .          .          .          .         .          
#> [119,]  .          .          .          .          .         .          
#> [120,]  .          .          .          .          .         .          
#> [121,]  .          .          .          .          .         .          
#> [122,]  .          .          .          .          .         .          
#> [123,]  .          .          .          .          .         .          
#> [124,]  .          .          .          .          .         .          
#> [125,]  .          .          .          .          .         .          
#> [126,]  .          .          .          .          .         .          
#> [127,]  .          .          .          .          .         .          
#> [128,]  .          .          .          .          .         .          
#> [129,]  .          .          .          .          .         .          
#> [130,]  .          .          .          .          .         .          
#> [131,]  .          .          .          .          .         .          
#> [132,]  .          .          .          .          .         .          
#> [133,]  .          .          .          .          .         .          
#> [134,]  .          .          .          .          .         .          
#> [135,]  .          .          .          .          .         .          
#> [136,]  .          .          .          .          .         .          
#> [137,]  .          .          .          .          .         .          
#> [138,]  .          .          .          .          .         .          
#> [139,]  .          .          .          .          .         .          
#> [140,]  .          .          .          .          .         .          
#> [141,]  .          .          .          .          .         .          
#> [142,]  .          .          .          .          .         .          
#> [143,]  .          .          .          .          .         .          
#> [144,]  .          .          .          .          .         .          
#> [145,]  .          .          .          .          .         .          
#> [146,]  .          .          .          .          .         .          
#> [147,]  .          .          .          .          .         .          
#> [148,]  .          .          .          .          .         .          
#> [149,]  .          .          .          .          .         .          
#> [150,]  .          .          .          .          .         .          
#> [151,]  .          .          .          .          .         .          
#> [152,]  .          .          .          .          .         .          
#> [153,]  .          .          .          .          .         .          
#> [154,]  .          .          .          .          .         .          
#> [155,]  .          .          .          .          .         .          
#> [156,]  .          .          .          .          .         .          
#> [157,]  .          .          .          .          .         .          
#> [158,]  .          .          .          .          .         .          
#> [159,]  .          .          .          .          .         .          
#> [160,]  .          .          .          .          .         .          
#> [161,]  .          .          .          .          .         .          
#> [162,]  .          .          .          .          .         .          
#> [163,]  .          .          .          .          .         .          
#> [164,]  .          .          .          .          .         .          
#> [165,]  .          .          .          .          .         .          
#> [166,]  .          .          .          .          .         .          
#> [167,]  .          .          .          .          .         .          
#> [168,]  .          .          .          .          .         .          
#> [169,]  .          .          .          .          .         .          
#> [170,]  .          .          .          .          .         .          
#> [171,]  .          .          .          .          .         .          
#> [172,]  .          .          .          .          .         .          
#> [173,]  .          .          .          .          .         .          
#> [174,]  .          .          .          .          .         .          
#> [175,]  .          .          .          .          .         .          
#> [176,]  .          .          .          .          .         .          
#> [177,]  .          .          .          .          .         .          
#> [178,]  .          .          .          .          .         .          
#> [179,]  .          .          .          .          .         .          
#> [180,]  .          .          .          .          .         .          
#> [181,]  .          .          .          .          .         .          
#> [182,]  .          .          .          .          .         .          
#> [183,]  .          .          .          .          .         .          
#> [184,]  .          .          .          .          .         .          
#> [185,]  .          .          .          .          .         .          
#> [186,]  .          .          .          .          .         .          
#> [187,]  .          .          .          .          .         .          
#> [188,]  .          .          .          .          .         .          
#> [189,]  .          .          .          .          .         .          
#> [190,]  .          .          .          .          .         .          
#> [191,]  .          .          .          .          .         .          
#> [192,]  .          .          .          .          .         .          
#> [193,]  .          .          .          .          .         .          
#> [194,]  .          .          .          .          .         .          
#> [195,]  .          .          .          .          .         .          
#> [196,]  .          .          .          .          .         .          
#> [197,]  .          .          .          .          .         .          
#> [198,]  .          .          .          .          .         .          
#> [199,]  .          .          .          .          .         .          
#> [200,]  .          .          .          .          .         .          
#> [201,]  .          .          .          .          .         .          
#> [202,]  .          .          .          .          .         .          
#> [203,]  .          .          .          .          .         .          
#> [204,]  .          .          .          .          .         .          
#> [205,]  .          .          .          .          .         .          
#> [206,]  .          .          .          .          .         .          
#> [207,]  .          .          .          .          .         .          
#> [208,]  .          .          .          .          .         .          
#> [209,]  .          .          .          .          .         .          
#> [210,]  .          .          .          .          .         .          
#> [211,]  .          .          .          .          .         .          
#>                                                                             
#>   [1,]  .          .          .      .      .          .         .          
#>   [2,]  .          .          .      .      .          .         .          
#>   [3,]  .          .          .      .      .          .         .          
#>   [4,]  .          .          .      .      .          .         .          
#>   [5,]  .          .          .      .      .          .         .          
#>   [6,]  .          .          .      .      .          .         .          
#>   [7,]  .          .          .      .      .          .         .          
#>   [8,]  .          .          .      .      .          .         .          
#>   [9,]  .          .          .      .      .          .         .          
#>  [10,]  .          .          .      .      .          .         .          
#>  [11,]  .          .          .      .      .          .         .          
#>  [12,]  .          .          .      .      .          .         .          
#>  [13,]  .          .          .      .      .          .         .          
#>  [14,]  .          .          .      .      .          .         .          
#>  [15,]  .          .          .      .      .          .         .          
#>  [16,]  .          .          .      .      .          .         .          
#>  [17,]  .          .          .      .      .          .         .          
#>  [18,]  .          .          .      .      .          .         .          
#>  [19,]  .          .          .      .      .          .         .          
#>  [20,]  .          .          .      .      .          .         .          
#>  [21,]  .          .          .      .      .          .         .          
#>  [22,]  .          .          .      .      .          .         .          
#>  [23,]  .          .          .      .      .          .         .          
#>  [24,]  .          .          .      .      .          .         .          
#>  [25,]  0.9047619 -0.9047619  .      .      .          .         .          
#>  [26,] -0.9047619  0.9047619  .      .      .          .         .          
#>  [27,]  .          .          0.975 -0.975  .          .         .          
#>  [28,]  .          .         -0.975  0.975  .          .         .          
#>  [29,]  .          .          .      .      0.9677419 -0.9677419 .          
#>  [30,]  .          .          .      .     -0.9677419  0.9677419 .          
#>  [31,]  .          .          .      .      .          .         6.18486e-10
#>  [32,]  .          .          .      .      .          .         .          
#>  [33,]  .          .          .      .      .          .         .          
#>  [34,]  .          .          .      .      .          .         .          
#>  [35,]  .          .          .      .      .          .         .          
#>  [36,]  .          .          .      .      .          .         .          
#>  [37,]  .          .          .      .      .          .         .          
#>  [38,]  .          .          .      .      .          .         .          
#>  [39,]  .          .          .      .      .          .         .          
#>  [40,]  .          .          .      .      .          .         .          
#>  [41,]  .          .          .      .      .          .         .          
#>  [42,]  .          .          .      .      .          .         .          
#>  [43,]  .          .          .      .      .          .         .          
#>  [44,]  .          .          .      .      .          .         .          
#>  [45,]  .          .          .      .      .          .         .          
#>  [46,]  .          .          .      .      .          .         .          
#>  [47,]  .          .          .      .      .          .         .          
#>  [48,]  .          .          .      .      .          .         .          
#>  [49,]  .          .          .      .      .          .         .          
#>  [50,]  .          .          .      .      .          .         .          
#>  [51,]  .          .          .      .      .          .         .          
#>  [52,]  .          .          .      .      .          .         .          
#>  [53,]  .          .          .      .      .          .         .          
#>  [54,]  .          .          .      .      .          .         .          
#>  [55,]  .          .          .      .      .          .         .          
#>  [56,]  .          .          .      .      .          .         .          
#>  [57,]  .          .          .      .      .          .         .          
#>  [58,]  .          .          .      .      .          .         .          
#>  [59,]  .          .          .      .      .          .         .          
#>  [60,]  .          .          .      .      .          .         .          
#>  [61,]  .          .          .      .      .          .         .          
#>  [62,]  .          .          .      .      .          .         .          
#>  [63,]  .          .          .      .      .          .         .          
#>  [64,]  .          .          .      .      .          .         .          
#>  [65,]  .          .          .      .      .          .         .          
#>  [66,]  .          .          .      .      .          .         .          
#>  [67,]  .          .          .      .      .          .         .          
#>  [68,]  .          .          .      .      .          .         .          
#>  [69,]  .          .          .      .      .          .         .          
#>  [70,]  .          .          .      .      .          .         .          
#>  [71,]  .          .          .      .      .          .         .          
#>  [72,]  .          .          .      .      .          .         .          
#>  [73,]  .          .          .      .      .          .         .          
#>  [74,]  .          .          .      .      .          .         .          
#>  [75,]  .          .          .      .      .          .         .          
#>  [76,]  .          .          .      .      .          .         .          
#>  [77,]  .          .          .      .      .          .         .          
#>  [78,]  .          .          .      .      .          .         .          
#>  [79,]  .          .          .      .      .          .         .          
#>  [80,]  .          .          .      .      .          .         .          
#>  [81,]  .          .          .      .      .          .         .          
#>  [82,]  .          .          .      .      .          .         .          
#>  [83,]  .          .          .      .      .          .         .          
#>  [84,]  .          .          .      .      .          .         .          
#>  [85,]  .          .          .      .      .          .         .          
#>  [86,]  .          .          .      .      .          .         .          
#>  [87,]  .          .          .      .      .          .         .          
#>  [88,]  .          .          .      .      .          .         .          
#>  [89,]  .          .          .      .      .          .         .          
#>  [90,]  .          .          .      .      .          .         .          
#>  [91,]  .          .          .      .      .          .         .          
#>  [92,]  .          .          .      .      .          .         .          
#>  [93,]  .          .          .      .      .          .         .          
#>  [94,]  .          .          .      .      .          .         .          
#>  [95,]  .          .          .      .      .          .         .          
#>  [96,]  .          .          .      .      .          .         .          
#>  [97,]  .          .          .      .      .          .         .          
#>  [98,]  .          .          .      .      .          .         .          
#>  [99,]  .          .          .      .      .          .         .          
#> [100,]  .          .          .      .      .          .         .          
#> [101,]  .          .          .      .      .          .         .          
#> [102,]  .          .          .      .      .          .         .          
#> [103,]  .          .          .      .      .          .         .          
#> [104,]  .          .          .      .      .          .         .          
#> [105,]  .          .          .      .      .          .         .          
#> [106,]  .          .          .      .      .          .         .          
#> [107,]  .          .          .      .      .          .         .          
#> [108,]  .          .          .      .      .          .         .          
#> [109,]  .          .          .      .      .          .         .          
#> [110,]  .          .          .      .      .          .         .          
#> [111,]  .          .          .      .      .          .         .          
#> [112,]  .          .          .      .      .          .         .          
#> [113,]  .          .          .      .      .          .         .          
#> [114,]  .          .          .      .      .          .         .          
#> [115,]  .          .          .      .      .          .         .          
#> [116,]  .          .          .      .      .          .         .          
#> [117,]  .          .          .      .      .          .         .          
#> [118,]  .          .          .      .      .          .         .          
#> [119,]  .          .          .      .      .          .         .          
#> [120,]  .          .          .      .      .          .         .          
#> [121,]  .          .          .      .      .          .         .          
#> [122,]  .          .          .      .      .          .         .          
#> [123,]  .          .          .      .      .          .         .          
#> [124,]  .          .          .      .      .          .         .          
#> [125,]  .          .          .      .      .          .         .          
#> [126,]  .          .          .      .      .          .         .          
#> [127,]  .          .          .      .      .          .         .          
#> [128,]  .          .          .      .      .          .         .          
#> [129,]  .          .          .      .      .          .         .          
#> [130,]  .          .          .      .      .          .         .          
#> [131,]  .          .          .      .      .          .         .          
#> [132,]  .          .          .      .      .          .         .          
#> [133,]  .          .          .      .      .          .         .          
#> [134,]  .          .          .      .      .          .         .          
#> [135,]  .          .          .      .      .          .         .          
#> [136,]  .          .          .      .      .          .         .          
#> [137,]  .          .          .      .      .          .         .          
#> [138,]  .          .          .      .      .          .         .          
#> [139,]  .          .          .      .      .          .         .          
#> [140,]  .          .          .      .      .          .         .          
#> [141,]  .          .          .      .      .          .         .          
#> [142,]  .          .          .      .      .          .         .          
#> [143,]  .          .          .      .      .          .         .          
#> [144,]  .          .          .      .      .          .         .          
#> [145,]  .          .          .      .      .          .         .          
#> [146,]  .          .          .      .      .          .         .          
#> [147,]  .          .          .      .      .          .         .          
#> [148,]  .          .          .      .      .          .         .          
#> [149,]  .          .          .      .      .          .         .          
#> [150,]  .          .          .      .      .          .         .          
#> [151,]  .          .          .      .      .          .         .          
#> [152,]  .          .          .      .      .          .         .          
#> [153,]  .          .          .      .      .          .         .          
#> [154,]  .          .          .      .      .          .         .          
#> [155,]  .          .          .      .      .          .         .          
#> [156,]  .          .          .      .      .          .         .          
#> [157,]  .          .          .      .      .          .         .          
#> [158,]  .          .          .      .      .          .         .          
#> [159,]  .          .          .      .      .          .         .          
#> [160,]  .          .          .      .      .          .         .          
#> [161,]  .          .          .      .      .          .         .          
#> [162,]  .          .          .      .      .          .         .          
#> [163,]  .          .          .      .      .          .         .          
#> [164,]  .          .          .      .      .          .         .          
#> [165,]  .          .          .      .      .          .         .          
#> [166,]  .          .          .      .      .          .         .          
#> [167,]  .          .          .      .      .          .         .          
#> [168,]  .          .          .      .      .          .         .          
#> [169,]  .          .          .      .      .          .         .          
#> [170,]  .          .          .      .      .          .         .          
#> [171,]  .          .          .      .      .          .         .          
#> [172,]  .          .          .      .      .          .         .          
#> [173,]  .          .          .      .      .          .         .          
#> [174,]  .          .          .      .      .          .         .          
#> [175,]  .          .          .      .      .          .         .          
#> [176,]  .          .          .      .      .          .         .          
#> [177,]  .          .          .      .      .          .         .          
#> [178,]  .          .          .      .      .          .         .          
#> [179,]  .          .          .      .      .          .         .          
#> [180,]  .          .          .      .      .          .         .          
#> [181,]  .          .          .      .      .          .         .          
#> [182,]  .          .          .      .      .          .         .          
#> [183,]  .          .          .      .      .          .         .          
#> [184,]  .          .          .      .      .          .         .          
#> [185,]  .          .          .      .      .          .         .          
#> [186,]  .          .          .      .      .          .         .          
#> [187,]  .          .          .      .      .          .         .          
#> [188,]  .          .          .      .      .          .         .          
#> [189,]  .          .          .      .      .          .         .          
#> [190,]  .          .          .      .      .          .         .          
#> [191,]  .          .          .      .      .          .         .          
#> [192,]  .          .          .      .      .          .         .          
#> [193,]  .          .          .      .      .          .         .          
#> [194,]  .          .          .      .      .          .         .          
#> [195,]  .          .          .      .      .          .         .          
#> [196,]  .          .          .      .      .          .         .          
#> [197,]  .          .          .      .      .          .         .          
#> [198,]  .          .          .      .      .          .         .          
#> [199,]  .          .          .      .      .          .         .          
#> [200,]  .          .          .      .      .          .         .          
#> [201,]  .          .          .      .      .          .         .          
#> [202,]  .          .          .      .      .          .         .          
#> [203,]  .          .          .      .      .          .         .          
#> [204,]  .          .          .      .      .          .         .          
#> [205,]  .          .          .      .      .          .         .          
#> [206,]  .          .          .      .      .          .         .          
#> [207,]  .          .          .      .      .          .         .          
#> [208,]  .          .          .      .      .          .         .          
#> [209,]  .          .          .      .      .          .         .          
#> [210,]  .          .          .      .      .          .         .          
#> [211,]  .          .          .      .      .          .         .          
#>                                                                          
#>   [1,] .            .          .          .          .          .        
#>   [2,] .            .          .          .          .          .        
#>   [3,] .            .          .          .          .          .        
#>   [4,] .            .          .          .          .          .        
#>   [5,] .            .          .          .          .          .        
#>   [6,] .            .          .          .          .          .        
#>   [7,] .            .          .          .          .          .        
#>   [8,] .            .          .          .          .          .        
#>   [9,] .            .          .          .          .          .        
#>  [10,] .            .          .          .          .          .        
#>  [11,] .            .          .          .          .          .        
#>  [12,] .            .          .          .          .          .        
#>  [13,] .            .          .          .          .          .        
#>  [14,] .            .          .          .          .          .        
#>  [15,] .            .          .          .          .          .        
#>  [16,] .            .          .          .          .          .        
#>  [17,] .            .          .          .          .          .        
#>  [18,] .            .          .          .          .          .        
#>  [19,] .            .          .          .          .          .        
#>  [20,] .            .          .          .          .          .        
#>  [21,] .            .          .          .          .          .        
#>  [22,] .            .          .          .          .          .        
#>  [23,] .            .          .          .          .          .        
#>  [24,] .            .          .          .          .          .        
#>  [25,] .            .          .          .          .          .        
#>  [26,] .            .          .          .          .          .        
#>  [27,] .            .          .          .          .          .        
#>  [28,] .            .          .          .          .          .        
#>  [29,] .            .          .          .          .          .        
#>  [30,] .            .          .          .          .          .        
#>  [31,] .            .          .          .          .          .        
#>  [32,] 6.18486e-10  .          .          .          .          .        
#>  [33,] .            0.9797422 -0.4898711  .          .          .        
#>  [34,] .           -0.4898711  0.9797422 -0.4898711  .          .        
#>  [35,] .            .         -0.4898711  0.9797422 -0.4898711  .        
#>  [36,] .            .          .         -0.4898711  0.9797422 -0.4898711
#>  [37,] .            .          .          .         -0.4898711  0.9797422
#>  [38,] .            .          .          .          .         -0.4898711
#>  [39,] .            .          .          .          .          .        
#>  [40,] .            .          .          .          .          .        
#>  [41,] .            .          .          .          .          .        
#>  [42,] .            .          .          .          .          .        
#>  [43,] .           -0.4898711  .          .          .          .        
#>  [44,] .            .          .          .          .          .        
#>  [45,] .            .          .          .          .          .        
#>  [46,] .            .          .          .          .          .        
#>  [47,] .            .          .          .          .          .        
#>  [48,] .            .          .          .          .          .        
#>  [49,] .            .          .          .          .          .        
#>  [50,] .            .          .          .          .          .        
#>  [51,] .            .          .          .          .          .        
#>  [52,] .            .          .          .          .          .        
#>  [53,] .            .          .          .          .          .        
#>  [54,] .            .          .          .          .          .        
#>  [55,] .            .          .          .          .          .        
#>  [56,] .            .          .          .          .          .        
#>  [57,] .            .          .          .          .          .        
#>  [58,] .            .          .          .          .          .        
#>  [59,] .            .          .          .          .          .        
#>  [60,] .            .          .          .          .          .        
#>  [61,] .            .          .          .          .          .        
#>  [62,] .            .          .          .          .          .        
#>  [63,] .            .          .          .          .          .        
#>  [64,] .            .          .          .          .          .        
#>  [65,] .            .          .          .          .          .        
#>  [66,] .            .          .          .          .          .        
#>  [67,] .            .          .          .          .          .        
#>  [68,] .            .          .          .          .          .        
#>  [69,] .            .          .          .          .          .        
#>  [70,] .            .          .          .          .          .        
#>  [71,] .            .          .          .          .          .        
#>  [72,] .            .          .          .          .          .        
#>  [73,] .            .          .          .          .          .        
#>  [74,] .            .          .          .          .          .        
#>  [75,] .            .          .          .          .          .        
#>  [76,] .            .          .          .          .          .        
#>  [77,] .            .          .          .          .          .        
#>  [78,] .            .          .          .          .          .        
#>  [79,] .            .          .          .          .          .        
#>  [80,] .            .          .          .          .          .        
#>  [81,] .            .          .          .          .          .        
#>  [82,] .            .          .          .          .          .        
#>  [83,] .            .          .          .          .          .        
#>  [84,] .            .          .          .          .          .        
#>  [85,] .            .          .          .          .          .        
#>  [86,] .            .          .          .          .          .        
#>  [87,] .            .          .          .          .          .        
#>  [88,] .            .          .          .          .          .        
#>  [89,] .            .          .          .          .          .        
#>  [90,] .            .          .          .          .          .        
#>  [91,] .            .          .          .          .          .        
#>  [92,] .            .          .          .          .          .        
#>  [93,] .            .          .          .          .          .        
#>  [94,] .            .          .          .          .          .        
#>  [95,] .            .          .          .          .          .        
#>  [96,] .            .          .          .          .          .        
#>  [97,] .            .          .          .          .          .        
#>  [98,] .            .          .          .          .          .        
#>  [99,] .            .          .          .          .          .        
#> [100,] .            .          .          .          .          .        
#> [101,] .            .          .          .          .          .        
#> [102,] .            .          .          .          .          .        
#> [103,] .            .          .          .          .          .        
#> [104,] .            .          .          .          .          .        
#> [105,] .            .          .          .          .          .        
#> [106,] .            .          .          .          .          .        
#> [107,] .            .          .          .          .          .        
#> [108,] .            .          .          .          .          .        
#> [109,] .            .          .          .          .          .        
#> [110,] .            .          .          .          .          .        
#> [111,] .            .          .          .          .          .        
#> [112,] .            .          .          .          .          .        
#> [113,] .            .          .          .          .          .        
#> [114,] .            .          .          .          .          .        
#> [115,] .            .          .          .          .          .        
#> [116,] .            .          .          .          .          .        
#> [117,] .            .          .          .          .          .        
#> [118,] .            .          .          .          .          .        
#> [119,] .            .          .          .          .          .        
#> [120,] .            .          .          .          .          .        
#> [121,] .            .          .          .          .          .        
#> [122,] .            .          .          .          .          .        
#> [123,] .            .          .          .          .          .        
#> [124,] .            .          .          .          .          .        
#> [125,] .            .          .          .          .          .        
#> [126,] .            .          .          .          .          .        
#> [127,] .            .          .          .          .          .        
#> [128,] .            .          .          .          .          .        
#> [129,] .            .          .          .          .          .        
#> [130,] .            .          .          .          .          .        
#> [131,] .            .          .          .          .          .        
#> [132,] .            .          .          .          .          .        
#> [133,] .            .          .          .          .          .        
#> [134,] .            .          .          .          .          .        
#> [135,] .            .          .          .          .          .        
#> [136,] .            .          .          .          .          .        
#> [137,] .            .          .          .          .          .        
#> [138,] .            .          .          .          .          .        
#> [139,] .            .          .          .          .          .        
#> [140,] .            .          .          .          .          .        
#> [141,] .            .          .          .          .          .        
#> [142,] .            .          .          .          .          .        
#> [143,] .            .          .          .          .          .        
#> [144,] .            .          .          .          .          .        
#> [145,] .            .          .          .          .          .        
#> [146,] .            .          .          .          .          .        
#> [147,] .            .          .          .          .          .        
#> [148,] .            .          .          .          .          .        
#> [149,] .            .          .          .          .          .        
#> [150,] .            .          .          .          .          .        
#> [151,] .            .          .          .          .          .        
#> [152,] .            .          .          .          .          .        
#> [153,] .            .          .          .          .          .        
#> [154,] .            .          .          .          .          .        
#> [155,] .            .          .          .          .          .        
#> [156,] .            .          .          .          .          .        
#> [157,] .            .          .          .          .          .        
#> [158,] .            .          .          .          .          .        
#> [159,] .            .          .          .          .          .        
#> [160,] .            .          .          .          .          .        
#> [161,] .            .          .          .          .          .        
#> [162,] .            .          .          .          .          .        
#> [163,] .            .          .          .          .          .        
#> [164,] .            .          .          .          .          .        
#> [165,] .            .          .          .          .          .        
#> [166,] .            .          .          .          .          .        
#> [167,] .            .          .          .          .          .        
#> [168,] .            .          .          .          .          .        
#> [169,] .            .          .          .          .          .        
#> [170,] .            .          .          .          .          .        
#> [171,] .            .          .          .          .          .        
#> [172,] .            .          .          .          .          .        
#> [173,] .            .          .          .          .          .        
#> [174,] .            .          .          .          .          .        
#> [175,] .            .          .          .          .          .        
#> [176,] .            .          .          .          .          .        
#> [177,] .            .          .          .          .          .        
#> [178,] .            .          .          .          .          .        
#> [179,] .            .          .          .          .          .        
#> [180,] .            .          .          .          .          .        
#> [181,] .            .          .          .          .          .        
#> [182,] .            .          .          .          .          .        
#> [183,] .            .          .          .          .          .        
#> [184,] .            .          .          .          .          .        
#> [185,] .            .          .          .          .          .        
#> [186,] .            .          .          .          .          .        
#> [187,] .            .          .          .          .          .        
#> [188,] .            .          .          .          .          .        
#> [189,] .            .          .          .          .          .        
#> [190,] .            .          .          .          .          .        
#> [191,] .            .          .          .          .          .        
#> [192,] .            .          .          .          .          .        
#> [193,] .            .          .          .          .          .        
#> [194,] .            .          .          .          .          .        
#> [195,] .            .          .          .          .          .        
#> [196,] .            .          .          .          .          .        
#> [197,] .            .          .          .          .          .        
#> [198,] .            .          .          .          .          .        
#> [199,] .            .          .          .          .          .        
#> [200,] .            .          .          .          .          .        
#> [201,] .            .          .          .          .          .        
#> [202,] .            .          .          .          .          .        
#> [203,] .            .          .          .          .          .        
#> [204,] .            .          .          .          .          .        
#> [205,] .            .          .          .          .          .        
#> [206,] .            .          .          .          .          .        
#> [207,] .            .          .          .          .          .        
#> [208,] .            .          .          .          .          .        
#> [209,] .            .          .          .          .          .        
#> [210,] .            .          .          .          .          .        
#> [211,] .            .          .          .          .          .        
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .         -0.4898711
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,] -0.4898711  .          .          .          .          .        
#>  [38,]  0.9797422 -0.4898711  .          .          .          .        
#>  [39,] -0.4898711  0.9797422 -0.4898711  .          .          .        
#>  [40,]  .         -0.4898711  0.9797422 -0.4898711  .          .        
#>  [41,]  .          .         -0.4898711  0.9797422 -0.4898711  .        
#>  [42,]  .          .          .         -0.4898711  0.9797422 -0.4898711
#>  [43,]  .          .          .          .         -0.4898711  0.9797422
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  0.9711538 -0.4855769 -0.4855769  .          .          .        
#>  [45,] -0.4855769  0.9711538 -0.4855769  .          .          .        
#>  [46,] -0.4855769 -0.4855769  0.9711538  .          .          .        
#>  [47,]  .          .          .          0.9791332 -0.4895666  .        
#>  [48,]  .          .          .         -0.4895666  0.9791332 -0.4895666
#>  [49,]  .          .          .          .         -0.4895666  0.9791332
#>  [50,]  .          .          .          .          .         -0.4895666
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .         -0.4895666  .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,] -0.4895666  .          .          .          .          .        
#>  [50,]  0.9791332 -0.4895666  .          .          .          .        
#>  [51,] -0.4895666  0.9791332 -0.4895666  .          .          .        
#>  [52,]  .         -0.4895666  0.9791332 -0.4895666  .          .        
#>  [53,]  .          .         -0.4895666  0.9791332 -0.4895666  .        
#>  [54,]  .          .          .         -0.4895666  0.9791332 -0.4895666
#>  [55,]  .          .          .          .         -0.4895666  0.9791332
#>  [56,]  .          .          .          .          .         -0.4895666
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .         -0.4895666  .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,] -0.4895666  .          .          .          .          .        
#>  [56,]  0.9791332 -0.4895666  .          .          .          .        
#>  [57,] -0.4895666  0.9791332 -0.4895666  .          .          .        
#>  [58,]  .         -0.4895666  0.9791332 -0.4895666  .          .        
#>  [59,]  .          .         -0.4895666  0.9791332  .          .        
#>  [60,]  .          .          .          .          0.9788136 -0.4894068
#>  [61,]  .          .          .          .         -0.4894068  0.9788136
#>  [62,]  .          .          .          .          .         -0.4894068
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .         -0.4894068  .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .         -0.4894068  .          .          .        
#>  [61,] -0.4894068  .          .          .          .          .        
#>  [62,]  0.9788136 -0.4894068  .          .          .          .        
#>  [63,] -0.4894068  0.9788136 -0.4894068  .          .          .        
#>  [64,]  .         -0.4894068  0.9788136  .          .          .        
#>  [65,]  .          .          .          0.9784615 -0.4892308  .        
#>  [66,]  .          .          .         -0.4892308  0.9784615 -0.4892308
#>  [67,]  .          .          .          .         -0.4892308  0.9784615
#>  [68,]  .          .          .          .          .         -0.4892308
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .         -0.4892308  .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                           
#>   [1,]  .          .          .          .          .       .       .     
#>   [2,]  .          .          .          .          .       .       .     
#>   [3,]  .          .          .          .          .       .       .     
#>   [4,]  .          .          .          .          .       .       .     
#>   [5,]  .          .          .          .          .       .       .     
#>   [6,]  .          .          .          .          .       .       .     
#>   [7,]  .          .          .          .          .       .       .     
#>   [8,]  .          .          .          .          .       .       .     
#>   [9,]  .          .          .          .          .       .       .     
#>  [10,]  .          .          .          .          .       .       .     
#>  [11,]  .          .          .          .          .       .       .     
#>  [12,]  .          .          .          .          .       .       .     
#>  [13,]  .          .          .          .          .       .       .     
#>  [14,]  .          .          .          .          .       .       .     
#>  [15,]  .          .          .          .          .       .       .     
#>  [16,]  .          .          .          .          .       .       .     
#>  [17,]  .          .          .          .          .       .       .     
#>  [18,]  .          .          .          .          .       .       .     
#>  [19,]  .          .          .          .          .       .       .     
#>  [20,]  .          .          .          .          .       .       .     
#>  [21,]  .          .          .          .          .       .       .     
#>  [22,]  .          .          .          .          .       .       .     
#>  [23,]  .          .          .          .          .       .       .     
#>  [24,]  .          .          .          .          .       .       .     
#>  [25,]  .          .          .          .          .       .       .     
#>  [26,]  .          .          .          .          .       .       .     
#>  [27,]  .          .          .          .          .       .       .     
#>  [28,]  .          .          .          .          .       .       .     
#>  [29,]  .          .          .          .          .       .       .     
#>  [30,]  .          .          .          .          .       .       .     
#>  [31,]  .          .          .          .          .       .       .     
#>  [32,]  .          .          .          .          .       .       .     
#>  [33,]  .          .          .          .          .       .       .     
#>  [34,]  .          .          .          .          .       .       .     
#>  [35,]  .          .          .          .          .       .       .     
#>  [36,]  .          .          .          .          .       .       .     
#>  [37,]  .          .          .          .          .       .       .     
#>  [38,]  .          .          .          .          .       .       .     
#>  [39,]  .          .          .          .          .       .       .     
#>  [40,]  .          .          .          .          .       .       .     
#>  [41,]  .          .          .          .          .       .       .     
#>  [42,]  .          .          .          .          .       .       .     
#>  [43,]  .          .          .          .          .       .       .     
#>  [44,]  .          .          .          .          .       .       .     
#>  [45,]  .          .          .          .          .       .       .     
#>  [46,]  .          .          .          .          .       .       .     
#>  [47,]  .          .          .          .          .       .       .     
#>  [48,]  .          .          .          .          .       .       .     
#>  [49,]  .          .          .          .          .       .       .     
#>  [50,]  .          .          .          .          .       .       .     
#>  [51,]  .          .          .          .          .       .       .     
#>  [52,]  .          .          .          .          .       .       .     
#>  [53,]  .          .          .          .          .       .       .     
#>  [54,]  .          .          .          .          .       .       .     
#>  [55,]  .          .          .          .          .       .       .     
#>  [56,]  .          .          .          .          .       .       .     
#>  [57,]  .          .          .          .          .       .       .     
#>  [58,]  .          .          .          .          .       .       .     
#>  [59,]  .          .          .          .          .       .       .     
#>  [60,]  .          .          .          .          .       .       .     
#>  [61,]  .          .          .          .          .       .       .     
#>  [62,]  .          .          .          .          .       .       .     
#>  [63,]  .          .          .          .          .       .       .     
#>  [64,]  .          .          .          .          .       .       .     
#>  [65,]  .          .          .         -0.4892308  .       .       .     
#>  [66,]  .          .          .          .          .       .       .     
#>  [67,] -0.4892308  .          .          .          .       .       .     
#>  [68,]  0.9784615 -0.4892308  .          .          .       .       .     
#>  [69,] -0.4892308  0.9784615 -0.4892308  .          .       .       .     
#>  [70,]  .         -0.4892308  0.9784615 -0.4892308  .       .       .     
#>  [71,]  .          .         -0.4892308  0.9784615  .       .       .     
#>  [72,]  .          .          .          .          0.9750 -0.4875 -0.4875
#>  [73,]  .          .          .          .         -0.4875  0.9750 -0.4875
#>  [74,]  .          .          .          .         -0.4875 -0.4875  0.9750
#>  [75,]  .          .          .          .          .       .       .     
#>  [76,]  .          .          .          .          .       .       .     
#>  [77,]  .          .          .          .          .       .       .     
#>  [78,]  .          .          .          .          .       .       .     
#>  [79,]  .          .          .          .          .       .       .     
#>  [80,]  .          .          .          .          .       .       .     
#>  [81,]  .          .          .          .          .       .       .     
#>  [82,]  .          .          .          .          .       .       .     
#>  [83,]  .          .          .          .          .       .       .     
#>  [84,]  .          .          .          .          .       .       .     
#>  [85,]  .          .          .          .          .       .       .     
#>  [86,]  .          .          .          .          .       .       .     
#>  [87,]  .          .          .          .          .       .       .     
#>  [88,]  .          .          .          .          .       .       .     
#>  [89,]  .          .          .          .          .       .       .     
#>  [90,]  .          .          .          .          .       .       .     
#>  [91,]  .          .          .          .          .       .       .     
#>  [92,]  .          .          .          .          .       .       .     
#>  [93,]  .          .          .          .          .       .       .     
#>  [94,]  .          .          .          .          .       .       .     
#>  [95,]  .          .          .          .          .       .       .     
#>  [96,]  .          .          .          .          .       .       .     
#>  [97,]  .          .          .          .          .       .       .     
#>  [98,]  .          .          .          .          .       .       .     
#>  [99,]  .          .          .          .          .       .       .     
#> [100,]  .          .          .          .          .       .       .     
#> [101,]  .          .          .          .          .       .       .     
#> [102,]  .          .          .          .          .       .       .     
#> [103,]  .          .          .          .          .       .       .     
#> [104,]  .          .          .          .          .       .       .     
#> [105,]  .          .          .          .          .       .       .     
#> [106,]  .          .          .          .          .       .       .     
#> [107,]  .          .          .          .          .       .       .     
#> [108,]  .          .          .          .          .       .       .     
#> [109,]  .          .          .          .          .       .       .     
#> [110,]  .          .          .          .          .       .       .     
#> [111,]  .          .          .          .          .       .       .     
#> [112,]  .          .          .          .          .       .       .     
#> [113,]  .          .          .          .          .       .       .     
#> [114,]  .          .          .          .          .       .       .     
#> [115,]  .          .          .          .          .       .       .     
#> [116,]  .          .          .          .          .       .       .     
#> [117,]  .          .          .          .          .       .       .     
#> [118,]  .          .          .          .          .       .       .     
#> [119,]  .          .          .          .          .       .       .     
#> [120,]  .          .          .          .          .       .       .     
#> [121,]  .          .          .          .          .       .       .     
#> [122,]  .          .          .          .          .       .       .     
#> [123,]  .          .          .          .          .       .       .     
#> [124,]  .          .          .          .          .       .       .     
#> [125,]  .          .          .          .          .       .       .     
#> [126,]  .          .          .          .          .       .       .     
#> [127,]  .          .          .          .          .       .       .     
#> [128,]  .          .          .          .          .       .       .     
#> [129,]  .          .          .          .          .       .       .     
#> [130,]  .          .          .          .          .       .       .     
#> [131,]  .          .          .          .          .       .       .     
#> [132,]  .          .          .          .          .       .       .     
#> [133,]  .          .          .          .          .       .       .     
#> [134,]  .          .          .          .          .       .       .     
#> [135,]  .          .          .          .          .       .       .     
#> [136,]  .          .          .          .          .       .       .     
#> [137,]  .          .          .          .          .       .       .     
#> [138,]  .          .          .          .          .       .       .     
#> [139,]  .          .          .          .          .       .       .     
#> [140,]  .          .          .          .          .       .       .     
#> [141,]  .          .          .          .          .       .       .     
#> [142,]  .          .          .          .          .       .       .     
#> [143,]  .          .          .          .          .       .       .     
#> [144,]  .          .          .          .          .       .       .     
#> [145,]  .          .          .          .          .       .       .     
#> [146,]  .          .          .          .          .       .       .     
#> [147,]  .          .          .          .          .       .       .     
#> [148,]  .          .          .          .          .       .       .     
#> [149,]  .          .          .          .          .       .       .     
#> [150,]  .          .          .          .          .       .       .     
#> [151,]  .          .          .          .          .       .       .     
#> [152,]  .          .          .          .          .       .       .     
#> [153,]  .          .          .          .          .       .       .     
#> [154,]  .          .          .          .          .       .       .     
#> [155,]  .          .          .          .          .       .       .     
#> [156,]  .          .          .          .          .       .       .     
#> [157,]  .          .          .          .          .       .       .     
#> [158,]  .          .          .          .          .       .       .     
#> [159,]  .          .          .          .          .       .       .     
#> [160,]  .          .          .          .          .       .       .     
#> [161,]  .          .          .          .          .       .       .     
#> [162,]  .          .          .          .          .       .       .     
#> [163,]  .          .          .          .          .       .       .     
#> [164,]  .          .          .          .          .       .       .     
#> [165,]  .          .          .          .          .       .       .     
#> [166,]  .          .          .          .          .       .       .     
#> [167,]  .          .          .          .          .       .       .     
#> [168,]  .          .          .          .          .       .       .     
#> [169,]  .          .          .          .          .       .       .     
#> [170,]  .          .          .          .          .       .       .     
#> [171,]  .          .          .          .          .       .       .     
#> [172,]  .          .          .          .          .       .       .     
#> [173,]  .          .          .          .          .       .       .     
#> [174,]  .          .          .          .          .       .       .     
#> [175,]  .          .          .          .          .       .       .     
#> [176,]  .          .          .          .          .       .       .     
#> [177,]  .          .          .          .          .       .       .     
#> [178,]  .          .          .          .          .       .       .     
#> [179,]  .          .          .          .          .       .       .     
#> [180,]  .          .          .          .          .       .       .     
#> [181,]  .          .          .          .          .       .       .     
#> [182,]  .          .          .          .          .       .       .     
#> [183,]  .          .          .          .          .       .       .     
#> [184,]  .          .          .          .          .       .       .     
#> [185,]  .          .          .          .          .       .       .     
#> [186,]  .          .          .          .          .       .       .     
#> [187,]  .          .          .          .          .       .       .     
#> [188,]  .          .          .          .          .       .       .     
#> [189,]  .          .          .          .          .       .       .     
#> [190,]  .          .          .          .          .       .       .     
#> [191,]  .          .          .          .          .       .       .     
#> [192,]  .          .          .          .          .       .       .     
#> [193,]  .          .          .          .          .       .       .     
#> [194,]  .          .          .          .          .       .       .     
#> [195,]  .          .          .          .          .       .       .     
#> [196,]  .          .          .          .          .       .       .     
#> [197,]  .          .          .          .          .       .       .     
#> [198,]  .          .          .          .          .       .       .     
#> [199,]  .          .          .          .          .       .       .     
#> [200,]  .          .          .          .          .       .       .     
#> [201,]  .          .          .          .          .       .       .     
#> [202,]  .          .          .          .          .       .       .     
#> [203,]  .          .          .          .          .       .       .     
#> [204,]  .          .          .          .          .       .       .     
#> [205,]  .          .          .          .          .       .       .     
#> [206,]  .          .          .          .          .       .       .     
#> [207,]  .          .          .          .          .       .       .     
#> [208,]  .          .          .          .          .       .       .     
#> [209,]  .          .          .          .          .       .       .     
#> [210,]  .          .          .          .          .       .       .     
#> [211,]  .          .          .          .          .       .       .     
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  0.9701493 -0.9701493  .          .          .          .        
#>  [76,] -0.9701493  0.9701493  .          .          .          .        
#>  [77,]  .          .          0.9783198 -0.4891599  .          .        
#>  [78,]  .          .         -0.4891599  0.9783198 -0.4891599  .        
#>  [79,]  .          .          .         -0.4891599  0.9783198 -0.4891599
#>  [80,]  .          .          .          .         -0.4891599  0.9783198
#>  [81,]  .          .          .          .          .         -0.4891599
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .         -0.4891599  .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .         -0.4891599  .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,] -0.4891599  .          .          .          .          .        
#>  [81,]  0.9783198 -0.4891599  .          .          .          .        
#>  [82,] -0.4891599  0.9783198 -0.4891599  .          .          .        
#>  [83,]  .         -0.4891599  0.9783198 -0.4891599  .          .        
#>  [84,]  .          .         -0.4891599  0.9783198  .          .        
#>  [85,]  .          .          .          .          0.9166667 -0.9166667
#>  [86,]  .          .          .          .         -0.9166667  0.9166667
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  0.9765625 -0.4882813  .          .          .         -0.4882813
#>  [88,] -0.4882813  0.9765625 -0.4882813  .          .          .        
#>  [89,]  .         -0.4882813  0.9765625 -0.4882813  .          .        
#>  [90,]  .          .         -0.4882813  0.9765625 -0.4882813  .        
#>  [91,]  .          .          .         -0.4882813  0.9765625 -0.4882813
#>  [92,] -0.4882813  .          .          .         -0.4882813  0.9765625
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                             
#>   [1,]  .         .         .         .         .         .         .       
#>   [2,]  .         .         .         .         .         .         .       
#>   [3,]  .         .         .         .         .         .         .       
#>   [4,]  .         .         .         .         .         .         .       
#>   [5,]  .         .         .         .         .         .         .       
#>   [6,]  .         .         .         .         .         .         .       
#>   [7,]  .         .         .         .         .         .         .       
#>   [8,]  .         .         .         .         .         .         .       
#>   [9,]  .         .         .         .         .         .         .       
#>  [10,]  .         .         .         .         .         .         .       
#>  [11,]  .         .         .         .         .         .         .       
#>  [12,]  .         .         .         .         .         .         .       
#>  [13,]  .         .         .         .         .         .         .       
#>  [14,]  .         .         .         .         .         .         .       
#>  [15,]  .         .         .         .         .         .         .       
#>  [16,]  .         .         .         .         .         .         .       
#>  [17,]  .         .         .         .         .         .         .       
#>  [18,]  .         .         .         .         .         .         .       
#>  [19,]  .         .         .         .         .         .         .       
#>  [20,]  .         .         .         .         .         .         .       
#>  [21,]  .         .         .         .         .         .         .       
#>  [22,]  .         .         .         .         .         .         .       
#>  [23,]  .         .         .         .         .         .         .       
#>  [24,]  .         .         .         .         .         .         .       
#>  [25,]  .         .         .         .         .         .         .       
#>  [26,]  .         .         .         .         .         .         .       
#>  [27,]  .         .         .         .         .         .         .       
#>  [28,]  .         .         .         .         .         .         .       
#>  [29,]  .         .         .         .         .         .         .       
#>  [30,]  .         .         .         .         .         .         .       
#>  [31,]  .         .         .         .         .         .         .       
#>  [32,]  .         .         .         .         .         .         .       
#>  [33,]  .         .         .         .         .         .         .       
#>  [34,]  .         .         .         .         .         .         .       
#>  [35,]  .         .         .         .         .         .         .       
#>  [36,]  .         .         .         .         .         .         .       
#>  [37,]  .         .         .         .         .         .         .       
#>  [38,]  .         .         .         .         .         .         .       
#>  [39,]  .         .         .         .         .         .         .       
#>  [40,]  .         .         .         .         .         .         .       
#>  [41,]  .         .         .         .         .         .         .       
#>  [42,]  .         .         .         .         .         .         .       
#>  [43,]  .         .         .         .         .         .         .       
#>  [44,]  .         .         .         .         .         .         .       
#>  [45,]  .         .         .         .         .         .         .       
#>  [46,]  .         .         .         .         .         .         .       
#>  [47,]  .         .         .         .         .         .         .       
#>  [48,]  .         .         .         .         .         .         .       
#>  [49,]  .         .         .         .         .         .         .       
#>  [50,]  .         .         .         .         .         .         .       
#>  [51,]  .         .         .         .         .         .         .       
#>  [52,]  .         .         .         .         .         .         .       
#>  [53,]  .         .         .         .         .         .         .       
#>  [54,]  .         .         .         .         .         .         .       
#>  [55,]  .         .         .         .         .         .         .       
#>  [56,]  .         .         .         .         .         .         .       
#>  [57,]  .         .         .         .         .         .         .       
#>  [58,]  .         .         .         .         .         .         .       
#>  [59,]  .         .         .         .         .         .         .       
#>  [60,]  .         .         .         .         .         .         .       
#>  [61,]  .         .         .         .         .         .         .       
#>  [62,]  .         .         .         .         .         .         .       
#>  [63,]  .         .         .         .         .         .         .       
#>  [64,]  .         .         .         .         .         .         .       
#>  [65,]  .         .         .         .         .         .         .       
#>  [66,]  .         .         .         .         .         .         .       
#>  [67,]  .         .         .         .         .         .         .       
#>  [68,]  .         .         .         .         .         .         .       
#>  [69,]  .         .         .         .         .         .         .       
#>  [70,]  .         .         .         .         .         .         .       
#>  [71,]  .         .         .         .         .         .         .       
#>  [72,]  .         .         .         .         .         .         .       
#>  [73,]  .         .         .         .         .         .         .       
#>  [74,]  .         .         .         .         .         .         .       
#>  [75,]  .         .         .         .         .         .         .       
#>  [76,]  .         .         .         .         .         .         .       
#>  [77,]  .         .         .         .         .         .         .       
#>  [78,]  .         .         .         .         .         .         .       
#>  [79,]  .         .         .         .         .         .         .       
#>  [80,]  .         .         .         .         .         .         .       
#>  [81,]  .         .         .         .         .         .         .       
#>  [82,]  .         .         .         .         .         .         .       
#>  [83,]  .         .         .         .         .         .         .       
#>  [84,]  .         .         .         .         .         .         .       
#>  [85,]  .         .         .         .         .         .         .       
#>  [86,]  .         .         .         .         .         .         .       
#>  [87,]  .         .         .         .         .         .         .       
#>  [88,]  .         .         .         .         .         .         .       
#>  [89,]  .         .         .         .         .         .         .       
#>  [90,]  .         .         .         .         .         .         .       
#>  [91,]  .         .         .         .         .         .         .       
#>  [92,]  .         .         .         .         .         .         .       
#>  [93,]  0.979798 -0.489899  .         .         .         .         .       
#>  [94,] -0.489899  0.979798 -0.489899  .         .         .         .       
#>  [95,]  .        -0.489899  0.979798 -0.489899  .         .         .       
#>  [96,]  .         .        -0.489899  0.979798 -0.489899  .         .       
#>  [97,]  .         .         .        -0.489899  0.979798 -0.489899  .       
#>  [98,]  .         .         .         .        -0.489899  0.979798 -0.489899
#>  [99,]  .         .         .         .         .        -0.489899  0.979798
#> [100,] -0.489899  .         .         .         .         .        -0.489899
#> [101,]  .         .         .         .         .         .         .       
#> [102,]  .         .         .         .         .         .         .       
#> [103,]  .         .         .         .         .         .         .       
#> [104,]  .         .         .         .         .         .         .       
#> [105,]  .         .         .         .         .         .         .       
#> [106,]  .         .         .         .         .         .         .       
#> [107,]  .         .         .         .         .         .         .       
#> [108,]  .         .         .         .         .         .         .       
#> [109,]  .         .         .         .         .         .         .       
#> [110,]  .         .         .         .         .         .         .       
#> [111,]  .         .         .         .         .         .         .       
#> [112,]  .         .         .         .         .         .         .       
#> [113,]  .         .         .         .         .         .         .       
#> [114,]  .         .         .         .         .         .         .       
#> [115,]  .         .         .         .         .         .         .       
#> [116,]  .         .         .         .         .         .         .       
#> [117,]  .         .         .         .         .         .         .       
#> [118,]  .         .         .         .         .         .         .       
#> [119,]  .         .         .         .         .         .         .       
#> [120,]  .         .         .         .         .         .         .       
#> [121,]  .         .         .         .         .         .         .       
#> [122,]  .         .         .         .         .         .         .       
#> [123,]  .         .         .         .         .         .         .       
#> [124,]  .         .         .         .         .         .         .       
#> [125,]  .         .         .         .         .         .         .       
#> [126,]  .         .         .         .         .         .         .       
#> [127,]  .         .         .         .         .         .         .       
#> [128,]  .         .         .         .         .         .         .       
#> [129,]  .         .         .         .         .         .         .       
#> [130,]  .         .         .         .         .         .         .       
#> [131,]  .         .         .         .         .         .         .       
#> [132,]  .         .         .         .         .         .         .       
#> [133,]  .         .         .         .         .         .         .       
#> [134,]  .         .         .         .         .         .         .       
#> [135,]  .         .         .         .         .         .         .       
#> [136,]  .         .         .         .         .         .         .       
#> [137,]  .         .         .         .         .         .         .       
#> [138,]  .         .         .         .         .         .         .       
#> [139,]  .         .         .         .         .         .         .       
#> [140,]  .         .         .         .         .         .         .       
#> [141,]  .         .         .         .         .         .         .       
#> [142,]  .         .         .         .         .         .         .       
#> [143,]  .         .         .         .         .         .         .       
#> [144,]  .         .         .         .         .         .         .       
#> [145,]  .         .         .         .         .         .         .       
#> [146,]  .         .         .         .         .         .         .       
#> [147,]  .         .         .         .         .         .         .       
#> [148,]  .         .         .         .         .         .         .       
#> [149,]  .         .         .         .         .         .         .       
#> [150,]  .         .         .         .         .         .         .       
#> [151,]  .         .         .         .         .         .         .       
#> [152,]  .         .         .         .         .         .         .       
#> [153,]  .         .         .         .         .         .         .       
#> [154,]  .         .         .         .         .         .         .       
#> [155,]  .         .         .         .         .         .         .       
#> [156,]  .         .         .         .         .         .         .       
#> [157,]  .         .         .         .         .         .         .       
#> [158,]  .         .         .         .         .         .         .       
#> [159,]  .         .         .         .         .         .         .       
#> [160,]  .         .         .         .         .         .         .       
#> [161,]  .         .         .         .         .         .         .       
#> [162,]  .         .         .         .         .         .         .       
#> [163,]  .         .         .         .         .         .         .       
#> [164,]  .         .         .         .         .         .         .       
#> [165,]  .         .         .         .         .         .         .       
#> [166,]  .         .         .         .         .         .         .       
#> [167,]  .         .         .         .         .         .         .       
#> [168,]  .         .         .         .         .         .         .       
#> [169,]  .         .         .         .         .         .         .       
#> [170,]  .         .         .         .         .         .         .       
#> [171,]  .         .         .         .         .         .         .       
#> [172,]  .         .         .         .         .         .         .       
#> [173,]  .         .         .         .         .         .         .       
#> [174,]  .         .         .         .         .         .         .       
#> [175,]  .         .         .         .         .         .         .       
#> [176,]  .         .         .         .         .         .         .       
#> [177,]  .         .         .         .         .         .         .       
#> [178,]  .         .         .         .         .         .         .       
#> [179,]  .         .         .         .         .         .         .       
#> [180,]  .         .         .         .         .         .         .       
#> [181,]  .         .         .         .         .         .         .       
#> [182,]  .         .         .         .         .         .         .       
#> [183,]  .         .         .         .         .         .         .       
#> [184,]  .         .         .         .         .         .         .       
#> [185,]  .         .         .         .         .         .         .       
#> [186,]  .         .         .         .         .         .         .       
#> [187,]  .         .         .         .         .         .         .       
#> [188,]  .         .         .         .         .         .         .       
#> [189,]  .         .         .         .         .         .         .       
#> [190,]  .         .         .         .         .         .         .       
#> [191,]  .         .         .         .         .         .         .       
#> [192,]  .         .         .         .         .         .         .       
#> [193,]  .         .         .         .         .         .         .       
#> [194,]  .         .         .         .         .         .         .       
#> [195,]  .         .         .         .         .         .         .       
#> [196,]  .         .         .         .         .         .         .       
#> [197,]  .         .         .         .         .         .         .       
#> [198,]  .         .         .         .         .         .         .       
#> [199,]  .         .         .         .         .         .         .       
#> [200,]  .         .         .         .         .         .         .       
#> [201,]  .         .         .         .         .         .         .       
#> [202,]  .         .         .         .         .         .         .       
#> [203,]  .         .         .         .         .         .         .       
#> [204,]  .         .         .         .         .         .         .       
#> [205,]  .         .         .         .         .         .         .       
#> [206,]  .         .         .         .         .         .         .       
#> [207,]  .         .         .         .         .         .         .       
#> [208,]  .         .         .         .         .         .         .       
#> [209,]  .         .         .         .         .         .         .       
#> [210,]  .         .         .         .         .         .         .       
#> [211,]  .         .         .         .         .         .         .       
#>                                                                                
#>   [1,]  .         .          .          .          .     .     .    .          
#>   [2,]  .         .          .          .          .     .     .    .          
#>   [3,]  .         .          .          .          .     .     .    .          
#>   [4,]  .         .          .          .          .     .     .    .          
#>   [5,]  .         .          .          .          .     .     .    .          
#>   [6,]  .         .          .          .          .     .     .    .          
#>   [7,]  .         .          .          .          .     .     .    .          
#>   [8,]  .         .          .          .          .     .     .    .          
#>   [9,]  .         .          .          .          .     .     .    .          
#>  [10,]  .         .          .          .          .     .     .    .          
#>  [11,]  .         .          .          .          .     .     .    .          
#>  [12,]  .         .          .          .          .     .     .    .          
#>  [13,]  .         .          .          .          .     .     .    .          
#>  [14,]  .         .          .          .          .     .     .    .          
#>  [15,]  .         .          .          .          .     .     .    .          
#>  [16,]  .         .          .          .          .     .     .    .          
#>  [17,]  .         .          .          .          .     .     .    .          
#>  [18,]  .         .          .          .          .     .     .    .          
#>  [19,]  .         .          .          .          .     .     .    .          
#>  [20,]  .         .          .          .          .     .     .    .          
#>  [21,]  .         .          .          .          .     .     .    .          
#>  [22,]  .         .          .          .          .     .     .    .          
#>  [23,]  .         .          .          .          .     .     .    .          
#>  [24,]  .         .          .          .          .     .     .    .          
#>  [25,]  .         .          .          .          .     .     .    .          
#>  [26,]  .         .          .          .          .     .     .    .          
#>  [27,]  .         .          .          .          .     .     .    .          
#>  [28,]  .         .          .          .          .     .     .    .          
#>  [29,]  .         .          .          .          .     .     .    .          
#>  [30,]  .         .          .          .          .     .     .    .          
#>  [31,]  .         .          .          .          .     .     .    .          
#>  [32,]  .         .          .          .          .     .     .    .          
#>  [33,]  .         .          .          .          .     .     .    .          
#>  [34,]  .         .          .          .          .     .     .    .          
#>  [35,]  .         .          .          .          .     .     .    .          
#>  [36,]  .         .          .          .          .     .     .    .          
#>  [37,]  .         .          .          .          .     .     .    .          
#>  [38,]  .         .          .          .          .     .     .    .          
#>  [39,]  .         .          .          .          .     .     .    .          
#>  [40,]  .         .          .          .          .     .     .    .          
#>  [41,]  .         .          .          .          .     .     .    .          
#>  [42,]  .         .          .          .          .     .     .    .          
#>  [43,]  .         .          .          .          .     .     .    .          
#>  [44,]  .         .          .          .          .     .     .    .          
#>  [45,]  .         .          .          .          .     .     .    .          
#>  [46,]  .         .          .          .          .     .     .    .          
#>  [47,]  .         .          .          .          .     .     .    .          
#>  [48,]  .         .          .          .          .     .     .    .          
#>  [49,]  .         .          .          .          .     .     .    .          
#>  [50,]  .         .          .          .          .     .     .    .          
#>  [51,]  .         .          .          .          .     .     .    .          
#>  [52,]  .         .          .          .          .     .     .    .          
#>  [53,]  .         .          .          .          .     .     .    .          
#>  [54,]  .         .          .          .          .     .     .    .          
#>  [55,]  .         .          .          .          .     .     .    .          
#>  [56,]  .         .          .          .          .     .     .    .          
#>  [57,]  .         .          .          .          .     .     .    .          
#>  [58,]  .         .          .          .          .     .     .    .          
#>  [59,]  .         .          .          .          .     .     .    .          
#>  [60,]  .         .          .          .          .     .     .    .          
#>  [61,]  .         .          .          .          .     .     .    .          
#>  [62,]  .         .          .          .          .     .     .    .          
#>  [63,]  .         .          .          .          .     .     .    .          
#>  [64,]  .         .          .          .          .     .     .    .          
#>  [65,]  .         .          .          .          .     .     .    .          
#>  [66,]  .         .          .          .          .     .     .    .          
#>  [67,]  .         .          .          .          .     .     .    .          
#>  [68,]  .         .          .          .          .     .     .    .          
#>  [69,]  .         .          .          .          .     .     .    .          
#>  [70,]  .         .          .          .          .     .     .    .          
#>  [71,]  .         .          .          .          .     .     .    .          
#>  [72,]  .         .          .          .          .     .     .    .          
#>  [73,]  .         .          .          .          .     .     .    .          
#>  [74,]  .         .          .          .          .     .     .    .          
#>  [75,]  .         .          .          .          .     .     .    .          
#>  [76,]  .         .          .          .          .     .     .    .          
#>  [77,]  .         .          .          .          .     .     .    .          
#>  [78,]  .         .          .          .          .     .     .    .          
#>  [79,]  .         .          .          .          .     .     .    .          
#>  [80,]  .         .          .          .          .     .     .    .          
#>  [81,]  .         .          .          .          .     .     .    .          
#>  [82,]  .         .          .          .          .     .     .    .          
#>  [83,]  .         .          .          .          .     .     .    .          
#>  [84,]  .         .          .          .          .     .     .    .          
#>  [85,]  .         .          .          .          .     .     .    .          
#>  [86,]  .         .          .          .          .     .     .    .          
#>  [87,]  .         .          .          .          .     .     .    .          
#>  [88,]  .         .          .          .          .     .     .    .          
#>  [89,]  .         .          .          .          .     .     .    .          
#>  [90,]  .         .          .          .          .     .     .    .          
#>  [91,]  .         .          .          .          .     .     .    .          
#>  [92,]  .         .          .          .          .     .     .    .          
#>  [93,] -0.489899  .          .          .          .     .     .    .          
#>  [94,]  .         .          .          .          .     .     .    .          
#>  [95,]  .         .          .          .          .     .     .    .          
#>  [96,]  .         .          .          .          .     .     .    .          
#>  [97,]  .         .          .          .          .     .     .    .          
#>  [98,]  .         .          .          .          .     .     .    .          
#>  [99,] -0.489899  .          .          .          .     .     .    .          
#> [100,]  0.979798  .          .          .          .     .     .    .          
#> [101,]  .         0.9779412 -0.4889706 -0.4889706  .     .     .    .          
#> [102,]  .        -0.4889706  0.9779412 -0.4889706  .     .     .    .          
#> [103,]  .        -0.4889706 -0.4889706  0.9779412  .     .     .    .          
#> [104,]  .         .          .          .          0.98 -0.49 -0.49 .          
#> [105,]  .         .          .          .         -0.49  0.98 -0.49 .          
#> [106,]  .         .          .          .         -0.49 -0.49  0.98 .          
#> [107,]  .         .          .          .          .     .     .    6.18486e-10
#> [108,]  .         .          .          .          .     .     .    .          
#> [109,]  .         .          .          .          .     .     .    .          
#> [110,]  .         .          .          .          .     .     .    .          
#> [111,]  .         .          .          .          .     .     .    .          
#> [112,]  .         .          .          .          .     .     .    .          
#> [113,]  .         .          .          .          .     .     .    .          
#> [114,]  .         .          .          .          .     .     .    .          
#> [115,]  .         .          .          .          .     .     .    .          
#> [116,]  .         .          .          .          .     .     .    .          
#> [117,]  .         .          .          .          .     .     .    .          
#> [118,]  .         .          .          .          .     .     .    .          
#> [119,]  .         .          .          .          .     .     .    .          
#> [120,]  .         .          .          .          .     .     .    .          
#> [121,]  .         .          .          .          .     .     .    .          
#> [122,]  .         .          .          .          .     .     .    .          
#> [123,]  .         .          .          .          .     .     .    .          
#> [124,]  .         .          .          .          .     .     .    .          
#> [125,]  .         .          .          .          .     .     .    .          
#> [126,]  .         .          .          .          .     .     .    .          
#> [127,]  .         .          .          .          .     .     .    .          
#> [128,]  .         .          .          .          .     .     .    .          
#> [129,]  .         .          .          .          .     .     .    .          
#> [130,]  .         .          .          .          .     .     .    .          
#> [131,]  .         .          .          .          .     .     .    .          
#> [132,]  .         .          .          .          .     .     .    .          
#> [133,]  .         .          .          .          .     .     .    .          
#> [134,]  .         .          .          .          .     .     .    .          
#> [135,]  .         .          .          .          .     .     .    .          
#> [136,]  .         .          .          .          .     .     .    .          
#> [137,]  .         .          .          .          .     .     .    .          
#> [138,]  .         .          .          .          .     .     .    .          
#> [139,]  .         .          .          .          .     .     .    .          
#> [140,]  .         .          .          .          .     .     .    .          
#> [141,]  .         .          .          .          .     .     .    .          
#> [142,]  .         .          .          .          .     .     .    .          
#> [143,]  .         .          .          .          .     .     .    .          
#> [144,]  .         .          .          .          .     .     .    .          
#> [145,]  .         .          .          .          .     .     .    .          
#> [146,]  .         .          .          .          .     .     .    .          
#> [147,]  .         .          .          .          .     .     .    .          
#> [148,]  .         .          .          .          .     .     .    .          
#> [149,]  .         .          .          .          .     .     .    .          
#> [150,]  .         .          .          .          .     .     .    .          
#> [151,]  .         .          .          .          .     .     .    .          
#> [152,]  .         .          .          .          .     .     .    .          
#> [153,]  .         .          .          .          .     .     .    .          
#> [154,]  .         .          .          .          .     .     .    .          
#> [155,]  .         .          .          .          .     .     .    .          
#> [156,]  .         .          .          .          .     .     .    .          
#> [157,]  .         .          .          .          .     .     .    .          
#> [158,]  .         .          .          .          .     .     .    .          
#> [159,]  .         .          .          .          .     .     .    .          
#> [160,]  .         .          .          .          .     .     .    .          
#> [161,]  .         .          .          .          .     .     .    .          
#> [162,]  .         .          .          .          .     .     .    .          
#> [163,]  .         .          .          .          .     .     .    .          
#> [164,]  .         .          .          .          .     .     .    .          
#> [165,]  .         .          .          .          .     .     .    .          
#> [166,]  .         .          .          .          .     .     .    .          
#> [167,]  .         .          .          .          .     .     .    .          
#> [168,]  .         .          .          .          .     .     .    .          
#> [169,]  .         .          .          .          .     .     .    .          
#> [170,]  .         .          .          .          .     .     .    .          
#> [171,]  .         .          .          .          .     .     .    .          
#> [172,]  .         .          .          .          .     .     .    .          
#> [173,]  .         .          .          .          .     .     .    .          
#> [174,]  .         .          .          .          .     .     .    .          
#> [175,]  .         .          .          .          .     .     .    .          
#> [176,]  .         .          .          .          .     .     .    .          
#> [177,]  .         .          .          .          .     .     .    .          
#> [178,]  .         .          .          .          .     .     .    .          
#> [179,]  .         .          .          .          .     .     .    .          
#> [180,]  .         .          .          .          .     .     .    .          
#> [181,]  .         .          .          .          .     .     .    .          
#> [182,]  .         .          .          .          .     .     .    .          
#> [183,]  .         .          .          .          .     .     .    .          
#> [184,]  .         .          .          .          .     .     .    .          
#> [185,]  .         .          .          .          .     .     .    .          
#> [186,]  .         .          .          .          .     .     .    .          
#> [187,]  .         .          .          .          .     .     .    .          
#> [188,]  .         .          .          .          .     .     .    .          
#> [189,]  .         .          .          .          .     .     .    .          
#> [190,]  .         .          .          .          .     .     .    .          
#> [191,]  .         .          .          .          .     .     .    .          
#> [192,]  .         .          .          .          .     .     .    .          
#> [193,]  .         .          .          .          .     .     .    .          
#> [194,]  .         .          .          .          .     .     .    .          
#> [195,]  .         .          .          .          .     .     .    .          
#> [196,]  .         .          .          .          .     .     .    .          
#> [197,]  .         .          .          .          .     .     .    .          
#> [198,]  .         .          .          .          .     .     .    .          
#> [199,]  .         .          .          .          .     .     .    .          
#> [200,]  .         .          .          .          .     .     .    .          
#> [201,]  .         .          .          .          .     .     .    .          
#> [202,]  .         .          .          .          .     .     .    .          
#> [203,]  .         .          .          .          .     .     .    .          
#> [204,]  .         .          .          .          .     .     .    .          
#> [205,]  .         .          .          .          .     .     .    .          
#> [206,]  .         .          .          .          .     .     .    .          
#> [207,]  .         .          .          .          .     .     .    .          
#> [208,]  .         .          .          .          .     .     .    .          
#> [209,]  .         .          .          .          .     .     .    .          
#> [210,]  .         .          .          .          .     .     .    .          
#> [211,]  .         .          .          .          .     .     .    .          
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  0.9622642 -0.9622642  .          .          .          .        
#> [109,] -0.9622642  0.9622642  .          .          .          .        
#> [110,]  .          .          0.9756098 -0.9756098  .          .        
#> [111,]  .          .         -0.9756098  0.9756098  .          .        
#> [112,]  .          .          .          .          0.9756098 -0.9756098
#> [113,]  .          .          .          .         -0.9756098  0.9756098
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  0.9746835 -0.9746835  .          .          .          .        
#> [115,] -0.9746835  0.9746835  .          .          .          .        
#> [116,]  .          .          0.9836066 -0.4897541  .          .        
#> [117,]  .          .         -0.4897541  0.9836066 -0.4897541  .        
#> [118,]  .          .          .         -0.4897541  0.9836066 -0.4897541
#> [119,]  .          .          .          .         -0.4897541  0.9836066
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  0.9821429 -0.4888393  .          .          .          .        
#> [121,] -0.4888393  0.9821429 -0.4888393  .          .          .        
#> [122,]  .         -0.4888393  0.9821429 -0.4888393  .          .        
#> [123,]  .          .         -0.4888393  0.9821429  .          .        
#> [124,]  .          .          .          .          0.9797297 -0.4898649
#> [125,]  .          .          .          .         -0.4898649  0.9797297
#> [126,]  .          .          .          .          .         -0.4898649
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .         -0.4898649  .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .         -0.4898649  .          .        
#> [125,] -0.4898649  .          .          .          .          .        
#> [126,]  0.9797297 -0.4898649  .          .          .          .        
#> [127,] -0.4898649  0.9797297 -0.4898649  .          .          .        
#> [128,]  .         -0.4898649  0.9797297 -0.4898649  .          .        
#> [129,]  .          .         -0.4898649  0.9797297  .          .        
#> [130,]  .          .          .          .          0.9795918 -0.9795918
#> [131,]  .          .          .          .         -0.9795918  0.9795918
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                               
#>   [1,]  .          .          .         .         .         .         .       
#>   [2,]  .          .          .         .         .         .         .       
#>   [3,]  .          .          .         .         .         .         .       
#>   [4,]  .          .          .         .         .         .         .       
#>   [5,]  .          .          .         .         .         .         .       
#>   [6,]  .          .          .         .         .         .         .       
#>   [7,]  .          .          .         .         .         .         .       
#>   [8,]  .          .          .         .         .         .         .       
#>   [9,]  .          .          .         .         .         .         .       
#>  [10,]  .          .          .         .         .         .         .       
#>  [11,]  .          .          .         .         .         .         .       
#>  [12,]  .          .          .         .         .         .         .       
#>  [13,]  .          .          .         .         .         .         .       
#>  [14,]  .          .          .         .         .         .         .       
#>  [15,]  .          .          .         .         .         .         .       
#>  [16,]  .          .          .         .         .         .         .       
#>  [17,]  .          .          .         .         .         .         .       
#>  [18,]  .          .          .         .         .         .         .       
#>  [19,]  .          .          .         .         .         .         .       
#>  [20,]  .          .          .         .         .         .         .       
#>  [21,]  .          .          .         .         .         .         .       
#>  [22,]  .          .          .         .         .         .         .       
#>  [23,]  .          .          .         .         .         .         .       
#>  [24,]  .          .          .         .         .         .         .       
#>  [25,]  .          .          .         .         .         .         .       
#>  [26,]  .          .          .         .         .         .         .       
#>  [27,]  .          .          .         .         .         .         .       
#>  [28,]  .          .          .         .         .         .         .       
#>  [29,]  .          .          .         .         .         .         .       
#>  [30,]  .          .          .         .         .         .         .       
#>  [31,]  .          .          .         .         .         .         .       
#>  [32,]  .          .          .         .         .         .         .       
#>  [33,]  .          .          .         .         .         .         .       
#>  [34,]  .          .          .         .         .         .         .       
#>  [35,]  .          .          .         .         .         .         .       
#>  [36,]  .          .          .         .         .         .         .       
#>  [37,]  .          .          .         .         .         .         .       
#>  [38,]  .          .          .         .         .         .         .       
#>  [39,]  .          .          .         .         .         .         .       
#>  [40,]  .          .          .         .         .         .         .       
#>  [41,]  .          .          .         .         .         .         .       
#>  [42,]  .          .          .         .         .         .         .       
#>  [43,]  .          .          .         .         .         .         .       
#>  [44,]  .          .          .         .         .         .         .       
#>  [45,]  .          .          .         .         .         .         .       
#>  [46,]  .          .          .         .         .         .         .       
#>  [47,]  .          .          .         .         .         .         .       
#>  [48,]  .          .          .         .         .         .         .       
#>  [49,]  .          .          .         .         .         .         .       
#>  [50,]  .          .          .         .         .         .         .       
#>  [51,]  .          .          .         .         .         .         .       
#>  [52,]  .          .          .         .         .         .         .       
#>  [53,]  .          .          .         .         .         .         .       
#>  [54,]  .          .          .         .         .         .         .       
#>  [55,]  .          .          .         .         .         .         .       
#>  [56,]  .          .          .         .         .         .         .       
#>  [57,]  .          .          .         .         .         .         .       
#>  [58,]  .          .          .         .         .         .         .       
#>  [59,]  .          .          .         .         .         .         .       
#>  [60,]  .          .          .         .         .         .         .       
#>  [61,]  .          .          .         .         .         .         .       
#>  [62,]  .          .          .         .         .         .         .       
#>  [63,]  .          .          .         .         .         .         .       
#>  [64,]  .          .          .         .         .         .         .       
#>  [65,]  .          .          .         .         .         .         .       
#>  [66,]  .          .          .         .         .         .         .       
#>  [67,]  .          .          .         .         .         .         .       
#>  [68,]  .          .          .         .         .         .         .       
#>  [69,]  .          .          .         .         .         .         .       
#>  [70,]  .          .          .         .         .         .         .       
#>  [71,]  .          .          .         .         .         .         .       
#>  [72,]  .          .          .         .         .         .         .       
#>  [73,]  .          .          .         .         .         .         .       
#>  [74,]  .          .          .         .         .         .         .       
#>  [75,]  .          .          .         .         .         .         .       
#>  [76,]  .          .          .         .         .         .         .       
#>  [77,]  .          .          .         .         .         .         .       
#>  [78,]  .          .          .         .         .         .         .       
#>  [79,]  .          .          .         .         .         .         .       
#>  [80,]  .          .          .         .         .         .         .       
#>  [81,]  .          .          .         .         .         .         .       
#>  [82,]  .          .          .         .         .         .         .       
#>  [83,]  .          .          .         .         .         .         .       
#>  [84,]  .          .          .         .         .         .         .       
#>  [85,]  .          .          .         .         .         .         .       
#>  [86,]  .          .          .         .         .         .         .       
#>  [87,]  .          .          .         .         .         .         .       
#>  [88,]  .          .          .         .         .         .         .       
#>  [89,]  .          .          .         .         .         .         .       
#>  [90,]  .          .          .         .         .         .         .       
#>  [91,]  .          .          .         .         .         .         .       
#>  [92,]  .          .          .         .         .         .         .       
#>  [93,]  .          .          .         .         .         .         .       
#>  [94,]  .          .          .         .         .         .         .       
#>  [95,]  .          .          .         .         .         .         .       
#>  [96,]  .          .          .         .         .         .         .       
#>  [97,]  .          .          .         .         .         .         .       
#>  [98,]  .          .          .         .         .         .         .       
#>  [99,]  .          .          .         .         .         .         .       
#> [100,]  .          .          .         .         .         .         .       
#> [101,]  .          .          .         .         .         .         .       
#> [102,]  .          .          .         .         .         .         .       
#> [103,]  .          .          .         .         .         .         .       
#> [104,]  .          .          .         .         .         .         .       
#> [105,]  .          .          .         .         .         .         .       
#> [106,]  .          .          .         .         .         .         .       
#> [107,]  .          .          .         .         .         .         .       
#> [108,]  .          .          .         .         .         .         .       
#> [109,]  .          .          .         .         .         .         .       
#> [110,]  .          .          .         .         .         .         .       
#> [111,]  .          .          .         .         .         .         .       
#> [112,]  .          .          .         .         .         .         .       
#> [113,]  .          .          .         .         .         .         .       
#> [114,]  .          .          .         .         .         .         .       
#> [115,]  .          .          .         .         .         .         .       
#> [116,]  .          .          .         .         .         .         .       
#> [117,]  .          .          .         .         .         .         .       
#> [118,]  .          .          .         .         .         .         .       
#> [119,]  .          .          .         .         .         .         .       
#> [120,]  .          .          .         .         .         .         .       
#> [121,]  .          .          .         .         .         .         .       
#> [122,]  .          .          .         .         .         .         .       
#> [123,]  .          .          .         .         .         .         .       
#> [124,]  .          .          .         .         .         .         .       
#> [125,]  .          .          .         .         .         .         .       
#> [126,]  .          .          .         .         .         .         .       
#> [127,]  .          .          .         .         .         .         .       
#> [128,]  .          .          .         .         .         .         .       
#> [129,]  .          .          .         .         .         .         .       
#> [130,]  .          .          .         .         .         .         .       
#> [131,]  .          .          .         .         .         .         .       
#> [132,]  0.9090909 -0.9090909  .         .         .         .         .       
#> [133,] -0.9090909  0.9090909  .         .         .         .         .       
#> [134,]  .          .          0.978836 -0.489418  .         .         .       
#> [135,]  .          .         -0.489418  0.978836 -0.489418  .         .       
#> [136,]  .          .          .        -0.489418  0.978836 -0.489418  .       
#> [137,]  .          .          .         .        -0.489418  0.978836 -0.489418
#> [138,]  .          .          .         .         .        -0.489418  0.978836
#> [139,]  .          .          .         .         .         .        -0.489418
#> [140,]  .          .          .         .         .         .         .       
#> [141,]  .          .          .         .         .         .         .       
#> [142,]  .          .          .         .         .         .         .       
#> [143,]  .          .          .         .         .         .         .       
#> [144,]  .          .          .         .         .         .         .       
#> [145,]  .          .          .         .         .         .         .       
#> [146,]  .          .          .         .         .         .         .       
#> [147,]  .          .          .         .         .         .         .       
#> [148,]  .          .          .         .         .         .         .       
#> [149,]  .          .         -0.489418  .         .         .         .       
#> [150,]  .          .          .         .         .         .         .       
#> [151,]  .          .          .         .         .         .         .       
#> [152,]  .          .          .         .         .         .         .       
#> [153,]  .          .          .         .         .         .         .       
#> [154,]  .          .          .         .         .         .         .       
#> [155,]  .          .          .         .         .         .         .       
#> [156,]  .          .          .         .         .         .         .       
#> [157,]  .          .          .         .         .         .         .       
#> [158,]  .          .          .         .         .         .         .       
#> [159,]  .          .          .         .         .         .         .       
#> [160,]  .          .          .         .         .         .         .       
#> [161,]  .          .          .         .         .         .         .       
#> [162,]  .          .          .         .         .         .         .       
#> [163,]  .          .          .         .         .         .         .       
#> [164,]  .          .          .         .         .         .         .       
#> [165,]  .          .          .         .         .         .         .       
#> [166,]  .          .          .         .         .         .         .       
#> [167,]  .          .          .         .         .         .         .       
#> [168,]  .          .          .         .         .         .         .       
#> [169,]  .          .          .         .         .         .         .       
#> [170,]  .          .          .         .         .         .         .       
#> [171,]  .          .          .         .         .         .         .       
#> [172,]  .          .          .         .         .         .         .       
#> [173,]  .          .          .         .         .         .         .       
#> [174,]  .          .          .         .         .         .         .       
#> [175,]  .          .          .         .         .         .         .       
#> [176,]  .          .          .         .         .         .         .       
#> [177,]  .          .          .         .         .         .         .       
#> [178,]  .          .          .         .         .         .         .       
#> [179,]  .          .          .         .         .         .         .       
#> [180,]  .          .          .         .         .         .         .       
#> [181,]  .          .          .         .         .         .         .       
#> [182,]  .          .          .         .         .         .         .       
#> [183,]  .          .          .         .         .         .         .       
#> [184,]  .          .          .         .         .         .         .       
#> [185,]  .          .          .         .         .         .         .       
#> [186,]  .          .          .         .         .         .         .       
#> [187,]  .          .          .         .         .         .         .       
#> [188,]  .          .          .         .         .         .         .       
#> [189,]  .          .          .         .         .         .         .       
#> [190,]  .          .          .         .         .         .         .       
#> [191,]  .          .          .         .         .         .         .       
#> [192,]  .          .          .         .         .         .         .       
#> [193,]  .          .          .         .         .         .         .       
#> [194,]  .          .          .         .         .         .         .       
#> [195,]  .          .          .         .         .         .         .       
#> [196,]  .          .          .         .         .         .         .       
#> [197,]  .          .          .         .         .         .         .       
#> [198,]  .          .          .         .         .         .         .       
#> [199,]  .          .          .         .         .         .         .       
#> [200,]  .          .          .         .         .         .         .       
#> [201,]  .          .          .         .         .         .         .       
#> [202,]  .          .          .         .         .         .         .       
#> [203,]  .          .          .         .         .         .         .       
#> [204,]  .          .          .         .         .         .         .       
#> [205,]  .          .          .         .         .         .         .       
#> [206,]  .          .          .         .         .         .         .       
#> [207,]  .          .          .         .         .         .         .       
#> [208,]  .          .          .         .         .         .         .       
#> [209,]  .          .          .         .         .         .         .       
#> [210,]  .          .          .         .         .         .         .       
#> [211,]  .          .          .         .         .         .         .       
#>                                                                             
#>   [1,]  .         .         .         .         .         .         .       
#>   [2,]  .         .         .         .         .         .         .       
#>   [3,]  .         .         .         .         .         .         .       
#>   [4,]  .         .         .         .         .         .         .       
#>   [5,]  .         .         .         .         .         .         .       
#>   [6,]  .         .         .         .         .         .         .       
#>   [7,]  .         .         .         .         .         .         .       
#>   [8,]  .         .         .         .         .         .         .       
#>   [9,]  .         .         .         .         .         .         .       
#>  [10,]  .         .         .         .         .         .         .       
#>  [11,]  .         .         .         .         .         .         .       
#>  [12,]  .         .         .         .         .         .         .       
#>  [13,]  .         .         .         .         .         .         .       
#>  [14,]  .         .         .         .         .         .         .       
#>  [15,]  .         .         .         .         .         .         .       
#>  [16,]  .         .         .         .         .         .         .       
#>  [17,]  .         .         .         .         .         .         .       
#>  [18,]  .         .         .         .         .         .         .       
#>  [19,]  .         .         .         .         .         .         .       
#>  [20,]  .         .         .         .         .         .         .       
#>  [21,]  .         .         .         .         .         .         .       
#>  [22,]  .         .         .         .         .         .         .       
#>  [23,]  .         .         .         .         .         .         .       
#>  [24,]  .         .         .         .         .         .         .       
#>  [25,]  .         .         .         .         .         .         .       
#>  [26,]  .         .         .         .         .         .         .       
#>  [27,]  .         .         .         .         .         .         .       
#>  [28,]  .         .         .         .         .         .         .       
#>  [29,]  .         .         .         .         .         .         .       
#>  [30,]  .         .         .         .         .         .         .       
#>  [31,]  .         .         .         .         .         .         .       
#>  [32,]  .         .         .         .         .         .         .       
#>  [33,]  .         .         .         .         .         .         .       
#>  [34,]  .         .         .         .         .         .         .       
#>  [35,]  .         .         .         .         .         .         .       
#>  [36,]  .         .         .         .         .         .         .       
#>  [37,]  .         .         .         .         .         .         .       
#>  [38,]  .         .         .         .         .         .         .       
#>  [39,]  .         .         .         .         .         .         .       
#>  [40,]  .         .         .         .         .         .         .       
#>  [41,]  .         .         .         .         .         .         .       
#>  [42,]  .         .         .         .         .         .         .       
#>  [43,]  .         .         .         .         .         .         .       
#>  [44,]  .         .         .         .         .         .         .       
#>  [45,]  .         .         .         .         .         .         .       
#>  [46,]  .         .         .         .         .         .         .       
#>  [47,]  .         .         .         .         .         .         .       
#>  [48,]  .         .         .         .         .         .         .       
#>  [49,]  .         .         .         .         .         .         .       
#>  [50,]  .         .         .         .         .         .         .       
#>  [51,]  .         .         .         .         .         .         .       
#>  [52,]  .         .         .         .         .         .         .       
#>  [53,]  .         .         .         .         .         .         .       
#>  [54,]  .         .         .         .         .         .         .       
#>  [55,]  .         .         .         .         .         .         .       
#>  [56,]  .         .         .         .         .         .         .       
#>  [57,]  .         .         .         .         .         .         .       
#>  [58,]  .         .         .         .         .         .         .       
#>  [59,]  .         .         .         .         .         .         .       
#>  [60,]  .         .         .         .         .         .         .       
#>  [61,]  .         .         .         .         .         .         .       
#>  [62,]  .         .         .         .         .         .         .       
#>  [63,]  .         .         .         .         .         .         .       
#>  [64,]  .         .         .         .         .         .         .       
#>  [65,]  .         .         .         .         .         .         .       
#>  [66,]  .         .         .         .         .         .         .       
#>  [67,]  .         .         .         .         .         .         .       
#>  [68,]  .         .         .         .         .         .         .       
#>  [69,]  .         .         .         .         .         .         .       
#>  [70,]  .         .         .         .         .         .         .       
#>  [71,]  .         .         .         .         .         .         .       
#>  [72,]  .         .         .         .         .         .         .       
#>  [73,]  .         .         .         .         .         .         .       
#>  [74,]  .         .         .         .         .         .         .       
#>  [75,]  .         .         .         .         .         .         .       
#>  [76,]  .         .         .         .         .         .         .       
#>  [77,]  .         .         .         .         .         .         .       
#>  [78,]  .         .         .         .         .         .         .       
#>  [79,]  .         .         .         .         .         .         .       
#>  [80,]  .         .         .         .         .         .         .       
#>  [81,]  .         .         .         .         .         .         .       
#>  [82,]  .         .         .         .         .         .         .       
#>  [83,]  .         .         .         .         .         .         .       
#>  [84,]  .         .         .         .         .         .         .       
#>  [85,]  .         .         .         .         .         .         .       
#>  [86,]  .         .         .         .         .         .         .       
#>  [87,]  .         .         .         .         .         .         .       
#>  [88,]  .         .         .         .         .         .         .       
#>  [89,]  .         .         .         .         .         .         .       
#>  [90,]  .         .         .         .         .         .         .       
#>  [91,]  .         .         .         .         .         .         .       
#>  [92,]  .         .         .         .         .         .         .       
#>  [93,]  .         .         .         .         .         .         .       
#>  [94,]  .         .         .         .         .         .         .       
#>  [95,]  .         .         .         .         .         .         .       
#>  [96,]  .         .         .         .         .         .         .       
#>  [97,]  .         .         .         .         .         .         .       
#>  [98,]  .         .         .         .         .         .         .       
#>  [99,]  .         .         .         .         .         .         .       
#> [100,]  .         .         .         .         .         .         .       
#> [101,]  .         .         .         .         .         .         .       
#> [102,]  .         .         .         .         .         .         .       
#> [103,]  .         .         .         .         .         .         .       
#> [104,]  .         .         .         .         .         .         .       
#> [105,]  .         .         .         .         .         .         .       
#> [106,]  .         .         .         .         .         .         .       
#> [107,]  .         .         .         .         .         .         .       
#> [108,]  .         .         .         .         .         .         .       
#> [109,]  .         .         .         .         .         .         .       
#> [110,]  .         .         .         .         .         .         .       
#> [111,]  .         .         .         .         .         .         .       
#> [112,]  .         .         .         .         .         .         .       
#> [113,]  .         .         .         .         .         .         .       
#> [114,]  .         .         .         .         .         .         .       
#> [115,]  .         .         .         .         .         .         .       
#> [116,]  .         .         .         .         .         .         .       
#> [117,]  .         .         .         .         .         .         .       
#> [118,]  .         .         .         .         .         .         .       
#> [119,]  .         .         .         .         .         .         .       
#> [120,]  .         .         .         .         .         .         .       
#> [121,]  .         .         .         .         .         .         .       
#> [122,]  .         .         .         .         .         .         .       
#> [123,]  .         .         .         .         .         .         .       
#> [124,]  .         .         .         .         .         .         .       
#> [125,]  .         .         .         .         .         .         .       
#> [126,]  .         .         .         .         .         .         .       
#> [127,]  .         .         .         .         .         .         .       
#> [128,]  .         .         .         .         .         .         .       
#> [129,]  .         .         .         .         .         .         .       
#> [130,]  .         .         .         .         .         .         .       
#> [131,]  .         .         .         .         .         .         .       
#> [132,]  .         .         .         .         .         .         .       
#> [133,]  .         .         .         .         .         .         .       
#> [134,]  .         .         .         .         .         .         .       
#> [135,]  .         .         .         .         .         .         .       
#> [136,]  .         .         .         .         .         .         .       
#> [137,]  .         .         .         .         .         .         .       
#> [138,] -0.489418  .         .         .         .         .         .       
#> [139,]  0.978836 -0.489418  .         .         .         .         .       
#> [140,] -0.489418  0.978836 -0.489418  .         .         .         .       
#> [141,]  .        -0.489418  0.978836 -0.489418  .         .         .       
#> [142,]  .         .        -0.489418  0.978836 -0.489418  .         .       
#> [143,]  .         .         .        -0.489418  0.978836 -0.489418  .       
#> [144,]  .         .         .         .        -0.489418  0.978836 -0.489418
#> [145,]  .         .         .         .         .        -0.489418  0.978836
#> [146,]  .         .         .         .         .         .        -0.489418
#> [147,]  .         .         .         .         .         .         .       
#> [148,]  .         .         .         .         .         .         .       
#> [149,]  .         .         .         .         .         .         .       
#> [150,]  .         .         .         .         .         .         .       
#> [151,]  .         .         .         .         .         .         .       
#> [152,]  .         .         .         .         .         .         .       
#> [153,]  .         .         .         .         .         .         .       
#> [154,]  .         .         .         .         .         .         .       
#> [155,]  .         .         .         .         .         .         .       
#> [156,]  .         .         .         .         .         .         .       
#> [157,]  .         .         .         .         .         .         .       
#> [158,]  .         .         .         .         .         .         .       
#> [159,]  .         .         .         .         .         .         .       
#> [160,]  .         .         .         .         .         .         .       
#> [161,]  .         .         .         .         .         .         .       
#> [162,]  .         .         .         .         .         .         .       
#> [163,]  .         .         .         .         .         .         .       
#> [164,]  .         .         .         .         .         .         .       
#> [165,]  .         .         .         .         .         .         .       
#> [166,]  .         .         .         .         .         .         .       
#> [167,]  .         .         .         .         .         .         .       
#> [168,]  .         .         .         .         .         .         .       
#> [169,]  .         .         .         .         .         .         .       
#> [170,]  .         .         .         .         .         .         .       
#> [171,]  .         .         .         .         .         .         .       
#> [172,]  .         .         .         .         .         .         .       
#> [173,]  .         .         .         .         .         .         .       
#> [174,]  .         .         .         .         .         .         .       
#> [175,]  .         .         .         .         .         .         .       
#> [176,]  .         .         .         .         .         .         .       
#> [177,]  .         .         .         .         .         .         .       
#> [178,]  .         .         .         .         .         .         .       
#> [179,]  .         .         .         .         .         .         .       
#> [180,]  .         .         .         .         .         .         .       
#> [181,]  .         .         .         .         .         .         .       
#> [182,]  .         .         .         .         .         .         .       
#> [183,]  .         .         .         .         .         .         .       
#> [184,]  .         .         .         .         .         .         .       
#> [185,]  .         .         .         .         .         .         .       
#> [186,]  .         .         .         .         .         .         .       
#> [187,]  .         .         .         .         .         .         .       
#> [188,]  .         .         .         .         .         .         .       
#> [189,]  .         .         .         .         .         .         .       
#> [190,]  .         .         .         .         .         .         .       
#> [191,]  .         .         .         .         .         .         .       
#> [192,]  .         .         .         .         .         .         .       
#> [193,]  .         .         .         .         .         .         .       
#> [194,]  .         .         .         .         .         .         .       
#> [195,]  .         .         .         .         .         .         .       
#> [196,]  .         .         .         .         .         .         .       
#> [197,]  .         .         .         .         .         .         .       
#> [198,]  .         .         .         .         .         .         .       
#> [199,]  .         .         .         .         .         .         .       
#> [200,]  .         .         .         .         .         .         .       
#> [201,]  .         .         .         .         .         .         .       
#> [202,]  .         .         .         .         .         .         .       
#> [203,]  .         .         .         .         .         .         .       
#> [204,]  .         .         .         .         .         .         .       
#> [205,]  .         .         .         .         .         .         .       
#> [206,]  .         .         .         .         .         .         .       
#> [207,]  .         .         .         .         .         .         .       
#> [208,]  .         .         .         .         .         .         .       
#> [209,]  .         .         .         .         .         .         .       
#> [210,]  .         .         .         .         .         .         .       
#> [211,]  .         .         .         .         .         .         .       
#>                                                                                
#>   [1,]  .         .         .         .         .          .          .        
#>   [2,]  .         .         .         .         .          .          .        
#>   [3,]  .         .         .         .         .          .          .        
#>   [4,]  .         .         .         .         .          .          .        
#>   [5,]  .         .         .         .         .          .          .        
#>   [6,]  .         .         .         .         .          .          .        
#>   [7,]  .         .         .         .         .          .          .        
#>   [8,]  .         .         .         .         .          .          .        
#>   [9,]  .         .         .         .         .          .          .        
#>  [10,]  .         .         .         .         .          .          .        
#>  [11,]  .         .         .         .         .          .          .        
#>  [12,]  .         .         .         .         .          .          .        
#>  [13,]  .         .         .         .         .          .          .        
#>  [14,]  .         .         .         .         .          .          .        
#>  [15,]  .         .         .         .         .          .          .        
#>  [16,]  .         .         .         .         .          .          .        
#>  [17,]  .         .         .         .         .          .          .        
#>  [18,]  .         .         .         .         .          .          .        
#>  [19,]  .         .         .         .         .          .          .        
#>  [20,]  .         .         .         .         .          .          .        
#>  [21,]  .         .         .         .         .          .          .        
#>  [22,]  .         .         .         .         .          .          .        
#>  [23,]  .         .         .         .         .          .          .        
#>  [24,]  .         .         .         .         .          .          .        
#>  [25,]  .         .         .         .         .          .          .        
#>  [26,]  .         .         .         .         .          .          .        
#>  [27,]  .         .         .         .         .          .          .        
#>  [28,]  .         .         .         .         .          .          .        
#>  [29,]  .         .         .         .         .          .          .        
#>  [30,]  .         .         .         .         .          .          .        
#>  [31,]  .         .         .         .         .          .          .        
#>  [32,]  .         .         .         .         .          .          .        
#>  [33,]  .         .         .         .         .          .          .        
#>  [34,]  .         .         .         .         .          .          .        
#>  [35,]  .         .         .         .         .          .          .        
#>  [36,]  .         .         .         .         .          .          .        
#>  [37,]  .         .         .         .         .          .          .        
#>  [38,]  .         .         .         .         .          .          .        
#>  [39,]  .         .         .         .         .          .          .        
#>  [40,]  .         .         .         .         .          .          .        
#>  [41,]  .         .         .         .         .          .          .        
#>  [42,]  .         .         .         .         .          .          .        
#>  [43,]  .         .         .         .         .          .          .        
#>  [44,]  .         .         .         .         .          .          .        
#>  [45,]  .         .         .         .         .          .          .        
#>  [46,]  .         .         .         .         .          .          .        
#>  [47,]  .         .         .         .         .          .          .        
#>  [48,]  .         .         .         .         .          .          .        
#>  [49,]  .         .         .         .         .          .          .        
#>  [50,]  .         .         .         .         .          .          .        
#>  [51,]  .         .         .         .         .          .          .        
#>  [52,]  .         .         .         .         .          .          .        
#>  [53,]  .         .         .         .         .          .          .        
#>  [54,]  .         .         .         .         .          .          .        
#>  [55,]  .         .         .         .         .          .          .        
#>  [56,]  .         .         .         .         .          .          .        
#>  [57,]  .         .         .         .         .          .          .        
#>  [58,]  .         .         .         .         .          .          .        
#>  [59,]  .         .         .         .         .          .          .        
#>  [60,]  .         .         .         .         .          .          .        
#>  [61,]  .         .         .         .         .          .          .        
#>  [62,]  .         .         .         .         .          .          .        
#>  [63,]  .         .         .         .         .          .          .        
#>  [64,]  .         .         .         .         .          .          .        
#>  [65,]  .         .         .         .         .          .          .        
#>  [66,]  .         .         .         .         .          .          .        
#>  [67,]  .         .         .         .         .          .          .        
#>  [68,]  .         .         .         .         .          .          .        
#>  [69,]  .         .         .         .         .          .          .        
#>  [70,]  .         .         .         .         .          .          .        
#>  [71,]  .         .         .         .         .          .          .        
#>  [72,]  .         .         .         .         .          .          .        
#>  [73,]  .         .         .         .         .          .          .        
#>  [74,]  .         .         .         .         .          .          .        
#>  [75,]  .         .         .         .         .          .          .        
#>  [76,]  .         .         .         .         .          .          .        
#>  [77,]  .         .         .         .         .          .          .        
#>  [78,]  .         .         .         .         .          .          .        
#>  [79,]  .         .         .         .         .          .          .        
#>  [80,]  .         .         .         .         .          .          .        
#>  [81,]  .         .         .         .         .          .          .        
#>  [82,]  .         .         .         .         .          .          .        
#>  [83,]  .         .         .         .         .          .          .        
#>  [84,]  .         .         .         .         .          .          .        
#>  [85,]  .         .         .         .         .          .          .        
#>  [86,]  .         .         .         .         .          .          .        
#>  [87,]  .         .         .         .         .          .          .        
#>  [88,]  .         .         .         .         .          .          .        
#>  [89,]  .         .         .         .         .          .          .        
#>  [90,]  .         .         .         .         .          .          .        
#>  [91,]  .         .         .         .         .          .          .        
#>  [92,]  .         .         .         .         .          .          .        
#>  [93,]  .         .         .         .         .          .          .        
#>  [94,]  .         .         .         .         .          .          .        
#>  [95,]  .         .         .         .         .          .          .        
#>  [96,]  .         .         .         .         .          .          .        
#>  [97,]  .         .         .         .         .          .          .        
#>  [98,]  .         .         .         .         .          .          .        
#>  [99,]  .         .         .         .         .          .          .        
#> [100,]  .         .         .         .         .          .          .        
#> [101,]  .         .         .         .         .          .          .        
#> [102,]  .         .         .         .         .          .          .        
#> [103,]  .         .         .         .         .          .          .        
#> [104,]  .         .         .         .         .          .          .        
#> [105,]  .         .         .         .         .          .          .        
#> [106,]  .         .         .         .         .          .          .        
#> [107,]  .         .         .         .         .          .          .        
#> [108,]  .         .         .         .         .          .          .        
#> [109,]  .         .         .         .         .          .          .        
#> [110,]  .         .         .         .         .          .          .        
#> [111,]  .         .         .         .         .          .          .        
#> [112,]  .         .         .         .         .          .          .        
#> [113,]  .         .         .         .         .          .          .        
#> [114,]  .         .         .         .         .          .          .        
#> [115,]  .         .         .         .         .          .          .        
#> [116,]  .         .         .         .         .          .          .        
#> [117,]  .         .         .         .         .          .          .        
#> [118,]  .         .         .         .         .          .          .        
#> [119,]  .         .         .         .         .          .          .        
#> [120,]  .         .         .         .         .          .          .        
#> [121,]  .         .         .         .         .          .          .        
#> [122,]  .         .         .         .         .          .          .        
#> [123,]  .         .         .         .         .          .          .        
#> [124,]  .         .         .         .         .          .          .        
#> [125,]  .         .         .         .         .          .          .        
#> [126,]  .         .         .         .         .          .          .        
#> [127,]  .         .         .         .         .          .          .        
#> [128,]  .         .         .         .         .          .          .        
#> [129,]  .         .         .         .         .          .          .        
#> [130,]  .         .         .         .         .          .          .        
#> [131,]  .         .         .         .         .          .          .        
#> [132,]  .         .         .         .         .          .          .        
#> [133,]  .         .         .         .         .          .          .        
#> [134,]  .         .         .        -0.489418  .          .          .        
#> [135,]  .         .         .         .         .          .          .        
#> [136,]  .         .         .         .         .          .          .        
#> [137,]  .         .         .         .         .          .          .        
#> [138,]  .         .         .         .         .          .          .        
#> [139,]  .         .         .         .         .          .          .        
#> [140,]  .         .         .         .         .          .          .        
#> [141,]  .         .         .         .         .          .          .        
#> [142,]  .         .         .         .         .          .          .        
#> [143,]  .         .         .         .         .          .          .        
#> [144,]  .         .         .         .         .          .          .        
#> [145,] -0.489418  .         .         .         .          .          .        
#> [146,]  0.978836 -0.489418  .         .         .          .          .        
#> [147,] -0.489418  0.978836 -0.489418  .         .          .          .        
#> [148,]  .        -0.489418  0.978836 -0.489418  .          .          .        
#> [149,]  .         .        -0.489418  0.978836  .          .          .        
#> [150,]  .         .         .         .         0.9760956 -0.4880478  .        
#> [151,]  .         .         .         .        -0.4880478  0.9760956 -0.4880478
#> [152,]  .         .         .         .         .         -0.4880478  0.9760956
#> [153,]  .         .         .         .         .          .         -0.4880478
#> [154,]  .         .         .         .         .          .          .        
#> [155,]  .         .         .         .        -0.4880478  .          .        
#> [156,]  .         .         .         .         .          .          .        
#> [157,]  .         .         .         .         .          .          .        
#> [158,]  .         .         .         .         .          .          .        
#> [159,]  .         .         .         .         .          .          .        
#> [160,]  .         .         .         .         .          .          .        
#> [161,]  .         .         .         .         .          .          .        
#> [162,]  .         .         .         .         .          .          .        
#> [163,]  .         .         .         .         .          .          .        
#> [164,]  .         .         .         .         .          .          .        
#> [165,]  .         .         .         .         .          .          .        
#> [166,]  .         .         .         .         .          .          .        
#> [167,]  .         .         .         .         .          .          .        
#> [168,]  .         .         .         .         .          .          .        
#> [169,]  .         .         .         .         .          .          .        
#> [170,]  .         .         .         .         .          .          .        
#> [171,]  .         .         .         .         .          .          .        
#> [172,]  .         .         .         .         .          .          .        
#> [173,]  .         .         .         .         .          .          .        
#> [174,]  .         .         .         .         .          .          .        
#> [175,]  .         .         .         .         .          .          .        
#> [176,]  .         .         .         .         .          .          .        
#> [177,]  .         .         .         .         .          .          .        
#> [178,]  .         .         .         .         .          .          .        
#> [179,]  .         .         .         .         .          .          .        
#> [180,]  .         .         .         .         .          .          .        
#> [181,]  .         .         .         .         .          .          .        
#> [182,]  .         .         .         .         .          .          .        
#> [183,]  .         .         .         .         .          .          .        
#> [184,]  .         .         .         .         .          .          .        
#> [185,]  .         .         .         .         .          .          .        
#> [186,]  .         .         .         .         .          .          .        
#> [187,]  .         .         .         .         .          .          .        
#> [188,]  .         .         .         .         .          .          .        
#> [189,]  .         .         .         .         .          .          .        
#> [190,]  .         .         .         .         .          .          .        
#> [191,]  .         .         .         .         .          .          .        
#> [192,]  .         .         .         .         .          .          .        
#> [193,]  .         .         .         .         .          .          .        
#> [194,]  .         .         .         .         .          .          .        
#> [195,]  .         .         .         .         .          .          .        
#> [196,]  .         .         .         .         .          .          .        
#> [197,]  .         .         .         .         .          .          .        
#> [198,]  .         .         .         .         .          .          .        
#> [199,]  .         .         .         .         .          .          .        
#> [200,]  .         .         .         .         .          .          .        
#> [201,]  .         .         .         .         .          .          .        
#> [202,]  .         .         .         .         .          .          .        
#> [203,]  .         .         .         .         .          .          .        
#> [204,]  .         .         .         .         .          .          .        
#> [205,]  .         .         .         .         .          .          .        
#> [206,]  .         .         .         .         .          .          .        
#> [207,]  .         .         .         .         .          .          .        
#> [208,]  .         .         .         .         .          .          .        
#> [209,]  .         .         .         .         .          .          .        
#> [210,]  .         .         .         .         .          .          .        
#> [211,]  .         .         .         .         .          .          .        
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .         -0.4880478  .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,] -0.4880478  .          .          .          .          .        
#> [153,]  0.9760956 -0.4880478  .          .          .          .        
#> [154,] -0.4880478  0.9760956 -0.4880478  .          .          .        
#> [155,]  .         -0.4880478  0.9760956  .          .          .        
#> [156,]  .          .          .          0.9747899 -0.4873950 -0.4873950
#> [157,]  .          .          .         -0.4873950  0.9747899 -0.4873950
#> [158,]  .          .          .         -0.4873950 -0.4873950  0.9747899
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                            
#>   [1,]  .          .          .          .     .     .     .     .     .   
#>   [2,]  .          .          .          .     .     .     .     .     .   
#>   [3,]  .          .          .          .     .     .     .     .     .   
#>   [4,]  .          .          .          .     .     .     .     .     .   
#>   [5,]  .          .          .          .     .     .     .     .     .   
#>   [6,]  .          .          .          .     .     .     .     .     .   
#>   [7,]  .          .          .          .     .     .     .     .     .   
#>   [8,]  .          .          .          .     .     .     .     .     .   
#>   [9,]  .          .          .          .     .     .     .     .     .   
#>  [10,]  .          .          .          .     .     .     .     .     .   
#>  [11,]  .          .          .          .     .     .     .     .     .   
#>  [12,]  .          .          .          .     .     .     .     .     .   
#>  [13,]  .          .          .          .     .     .     .     .     .   
#>  [14,]  .          .          .          .     .     .     .     .     .   
#>  [15,]  .          .          .          .     .     .     .     .     .   
#>  [16,]  .          .          .          .     .     .     .     .     .   
#>  [17,]  .          .          .          .     .     .     .     .     .   
#>  [18,]  .          .          .          .     .     .     .     .     .   
#>  [19,]  .          .          .          .     .     .     .     .     .   
#>  [20,]  .          .          .          .     .     .     .     .     .   
#>  [21,]  .          .          .          .     .     .     .     .     .   
#>  [22,]  .          .          .          .     .     .     .     .     .   
#>  [23,]  .          .          .          .     .     .     .     .     .   
#>  [24,]  .          .          .          .     .     .     .     .     .   
#>  [25,]  .          .          .          .     .     .     .     .     .   
#>  [26,]  .          .          .          .     .     .     .     .     .   
#>  [27,]  .          .          .          .     .     .     .     .     .   
#>  [28,]  .          .          .          .     .     .     .     .     .   
#>  [29,]  .          .          .          .     .     .     .     .     .   
#>  [30,]  .          .          .          .     .     .     .     .     .   
#>  [31,]  .          .          .          .     .     .     .     .     .   
#>  [32,]  .          .          .          .     .     .     .     .     .   
#>  [33,]  .          .          .          .     .     .     .     .     .   
#>  [34,]  .          .          .          .     .     .     .     .     .   
#>  [35,]  .          .          .          .     .     .     .     .     .   
#>  [36,]  .          .          .          .     .     .     .     .     .   
#>  [37,]  .          .          .          .     .     .     .     .     .   
#>  [38,]  .          .          .          .     .     .     .     .     .   
#>  [39,]  .          .          .          .     .     .     .     .     .   
#>  [40,]  .          .          .          .     .     .     .     .     .   
#>  [41,]  .          .          .          .     .     .     .     .     .   
#>  [42,]  .          .          .          .     .     .     .     .     .   
#>  [43,]  .          .          .          .     .     .     .     .     .   
#>  [44,]  .          .          .          .     .     .     .     .     .   
#>  [45,]  .          .          .          .     .     .     .     .     .   
#>  [46,]  .          .          .          .     .     .     .     .     .   
#>  [47,]  .          .          .          .     .     .     .     .     .   
#>  [48,]  .          .          .          .     .     .     .     .     .   
#>  [49,]  .          .          .          .     .     .     .     .     .   
#>  [50,]  .          .          .          .     .     .     .     .     .   
#>  [51,]  .          .          .          .     .     .     .     .     .   
#>  [52,]  .          .          .          .     .     .     .     .     .   
#>  [53,]  .          .          .          .     .     .     .     .     .   
#>  [54,]  .          .          .          .     .     .     .     .     .   
#>  [55,]  .          .          .          .     .     .     .     .     .   
#>  [56,]  .          .          .          .     .     .     .     .     .   
#>  [57,]  .          .          .          .     .     .     .     .     .   
#>  [58,]  .          .          .          .     .     .     .     .     .   
#>  [59,]  .          .          .          .     .     .     .     .     .   
#>  [60,]  .          .          .          .     .     .     .     .     .   
#>  [61,]  .          .          .          .     .     .     .     .     .   
#>  [62,]  .          .          .          .     .     .     .     .     .   
#>  [63,]  .          .          .          .     .     .     .     .     .   
#>  [64,]  .          .          .          .     .     .     .     .     .   
#>  [65,]  .          .          .          .     .     .     .     .     .   
#>  [66,]  .          .          .          .     .     .     .     .     .   
#>  [67,]  .          .          .          .     .     .     .     .     .   
#>  [68,]  .          .          .          .     .     .     .     .     .   
#>  [69,]  .          .          .          .     .     .     .     .     .   
#>  [70,]  .          .          .          .     .     .     .     .     .   
#>  [71,]  .          .          .          .     .     .     .     .     .   
#>  [72,]  .          .          .          .     .     .     .     .     .   
#>  [73,]  .          .          .          .     .     .     .     .     .   
#>  [74,]  .          .          .          .     .     .     .     .     .   
#>  [75,]  .          .          .          .     .     .     .     .     .   
#>  [76,]  .          .          .          .     .     .     .     .     .   
#>  [77,]  .          .          .          .     .     .     .     .     .   
#>  [78,]  .          .          .          .     .     .     .     .     .   
#>  [79,]  .          .          .          .     .     .     .     .     .   
#>  [80,]  .          .          .          .     .     .     .     .     .   
#>  [81,]  .          .          .          .     .     .     .     .     .   
#>  [82,]  .          .          .          .     .     .     .     .     .   
#>  [83,]  .          .          .          .     .     .     .     .     .   
#>  [84,]  .          .          .          .     .     .     .     .     .   
#>  [85,]  .          .          .          .     .     .     .     .     .   
#>  [86,]  .          .          .          .     .     .     .     .     .   
#>  [87,]  .          .          .          .     .     .     .     .     .   
#>  [88,]  .          .          .          .     .     .     .     .     .   
#>  [89,]  .          .          .          .     .     .     .     .     .   
#>  [90,]  .          .          .          .     .     .     .     .     .   
#>  [91,]  .          .          .          .     .     .     .     .     .   
#>  [92,]  .          .          .          .     .     .     .     .     .   
#>  [93,]  .          .          .          .     .     .     .     .     .   
#>  [94,]  .          .          .          .     .     .     .     .     .   
#>  [95,]  .          .          .          .     .     .     .     .     .   
#>  [96,]  .          .          .          .     .     .     .     .     .   
#>  [97,]  .          .          .          .     .     .     .     .     .   
#>  [98,]  .          .          .          .     .     .     .     .     .   
#>  [99,]  .          .          .          .     .     .     .     .     .   
#> [100,]  .          .          .          .     .     .     .     .     .   
#> [101,]  .          .          .          .     .     .     .     .     .   
#> [102,]  .          .          .          .     .     .     .     .     .   
#> [103,]  .          .          .          .     .     .     .     .     .   
#> [104,]  .          .          .          .     .     .     .     .     .   
#> [105,]  .          .          .          .     .     .     .     .     .   
#> [106,]  .          .          .          .     .     .     .     .     .   
#> [107,]  .          .          .          .     .     .     .     .     .   
#> [108,]  .          .          .          .     .     .     .     .     .   
#> [109,]  .          .          .          .     .     .     .     .     .   
#> [110,]  .          .          .          .     .     .     .     .     .   
#> [111,]  .          .          .          .     .     .     .     .     .   
#> [112,]  .          .          .          .     .     .     .     .     .   
#> [113,]  .          .          .          .     .     .     .     .     .   
#> [114,]  .          .          .          .     .     .     .     .     .   
#> [115,]  .          .          .          .     .     .     .     .     .   
#> [116,]  .          .          .          .     .     .     .     .     .   
#> [117,]  .          .          .          .     .     .     .     .     .   
#> [118,]  .          .          .          .     .     .     .     .     .   
#> [119,]  .          .          .          .     .     .     .     .     .   
#> [120,]  .          .          .          .     .     .     .     .     .   
#> [121,]  .          .          .          .     .     .     .     .     .   
#> [122,]  .          .          .          .     .     .     .     .     .   
#> [123,]  .          .          .          .     .     .     .     .     .   
#> [124,]  .          .          .          .     .     .     .     .     .   
#> [125,]  .          .          .          .     .     .     .     .     .   
#> [126,]  .          .          .          .     .     .     .     .     .   
#> [127,]  .          .          .          .     .     .     .     .     .   
#> [128,]  .          .          .          .     .     .     .     .     .   
#> [129,]  .          .          .          .     .     .     .     .     .   
#> [130,]  .          .          .          .     .     .     .     .     .   
#> [131,]  .          .          .          .     .     .     .     .     .   
#> [132,]  .          .          .          .     .     .     .     .     .   
#> [133,]  .          .          .          .     .     .     .     .     .   
#> [134,]  .          .          .          .     .     .     .     .     .   
#> [135,]  .          .          .          .     .     .     .     .     .   
#> [136,]  .          .          .          .     .     .     .     .     .   
#> [137,]  .          .          .          .     .     .     .     .     .   
#> [138,]  .          .          .          .     .     .     .     .     .   
#> [139,]  .          .          .          .     .     .     .     .     .   
#> [140,]  .          .          .          .     .     .     .     .     .   
#> [141,]  .          .          .          .     .     .     .     .     .   
#> [142,]  .          .          .          .     .     .     .     .     .   
#> [143,]  .          .          .          .     .     .     .     .     .   
#> [144,]  .          .          .          .     .     .     .     .     .   
#> [145,]  .          .          .          .     .     .     .     .     .   
#> [146,]  .          .          .          .     .     .     .     .     .   
#> [147,]  .          .          .          .     .     .     .     .     .   
#> [148,]  .          .          .          .     .     .     .     .     .   
#> [149,]  .          .          .          .     .     .     .     .     .   
#> [150,]  .          .          .          .     .     .     .     .     .   
#> [151,]  .          .          .          .     .     .     .     .     .   
#> [152,]  .          .          .          .     .     .     .     .     .   
#> [153,]  .          .          .          .     .     .     .     .     .   
#> [154,]  .          .          .          .     .     .     .     .     .   
#> [155,]  .          .          .          .     .     .     .     .     .   
#> [156,]  .          .          .          .     .     .     .     .     .   
#> [157,]  .          .          .          .     .     .     .     .     .   
#> [158,]  .          .          .          .     .     .     .     .     .   
#> [159,]  0.9787234 -0.4893617 -0.4893617  .     .     .     .     .     .   
#> [160,] -0.4893617  0.9787234 -0.4893617  .     .     .     .     .     .   
#> [161,] -0.4893617 -0.4893617  0.9787234  .     .     .     .     .     .   
#> [162,]  .          .          .          0.98 -0.49  .     .     .     .   
#> [163,]  .          .          .         -0.49  0.98 -0.49  .     .     .   
#> [164,]  .          .          .          .    -0.49  0.98 -0.49  .     .   
#> [165,]  .          .          .          .     .    -0.49  0.98 -0.49  .   
#> [166,]  .          .          .          .     .     .    -0.49  0.98 -0.49
#> [167,]  .          .          .          .     .     .     .    -0.49  0.98
#> [168,]  .          .          .          .     .     .     .     .    -0.49
#> [169,]  .          .          .          .     .     .     .     .     .   
#> [170,]  .          .          .         -0.49  .     .     .     .     .   
#> [171,]  .          .          .          .     .     .     .     .     .   
#> [172,]  .          .          .          .     .     .     .     .     .   
#> [173,]  .          .          .          .     .     .     .     .     .   
#> [174,]  .          .          .          .     .     .     .     .     .   
#> [175,]  .          .          .          .     .     .     .     .     .   
#> [176,]  .          .          .          .     .     .     .     .     .   
#> [177,]  .          .          .          .     .     .     .     .     .   
#> [178,]  .          .          .          .     .     .     .     .     .   
#> [179,]  .          .          .          .     .     .     .     .     .   
#> [180,]  .          .          .          .     .     .     .     .     .   
#> [181,]  .          .          .          .     .     .     .     .     .   
#> [182,]  .          .          .          .     .     .     .     .     .   
#> [183,]  .          .          .          .     .     .     .     .     .   
#> [184,]  .          .          .          .     .     .     .     .     .   
#> [185,]  .          .          .          .     .     .     .     .     .   
#> [186,]  .          .          .          .     .     .     .     .     .   
#> [187,]  .          .          .          .     .     .     .     .     .   
#> [188,]  .          .          .          .     .     .     .     .     .   
#> [189,]  .          .          .          .     .     .     .     .     .   
#> [190,]  .          .          .          .     .     .     .     .     .   
#> [191,]  .          .          .          .     .     .     .     .     .   
#> [192,]  .          .          .          .     .     .     .     .     .   
#> [193,]  .          .          .          .     .     .     .     .     .   
#> [194,]  .          .          .          .     .     .     .     .     .   
#> [195,]  .          .          .          .     .     .     .     .     .   
#> [196,]  .          .          .          .     .     .     .     .     .   
#> [197,]  .          .          .          .     .     .     .     .     .   
#> [198,]  .          .          .          .     .     .     .     .     .   
#> [199,]  .          .          .          .     .     .     .     .     .   
#> [200,]  .          .          .          .     .     .     .     .     .   
#> [201,]  .          .          .          .     .     .     .     .     .   
#> [202,]  .          .          .          .     .     .     .     .     .   
#> [203,]  .          .          .          .     .     .     .     .     .   
#> [204,]  .          .          .          .     .     .     .     .     .   
#> [205,]  .          .          .          .     .     .     .     .     .   
#> [206,]  .          .          .          .     .     .     .     .     .   
#> [207,]  .          .          .          .     .     .     .     .     .   
#> [208,]  .          .          .          .     .     .     .     .     .   
#> [209,]  .          .          .          .     .     .     .     .     .   
#> [210,]  .          .          .          .     .     .     .     .     .   
#> [211,]  .          .          .          .     .     .     .     .     .   
#>                                                                              
#>   [1,]  .     .     .     .          .          .         .         .        
#>   [2,]  .     .     .     .          .          .         .         .        
#>   [3,]  .     .     .     .          .          .         .         .        
#>   [4,]  .     .     .     .          .          .         .         .        
#>   [5,]  .     .     .     .          .          .         .         .        
#>   [6,]  .     .     .     .          .          .         .         .        
#>   [7,]  .     .     .     .          .          .         .         .        
#>   [8,]  .     .     .     .          .          .         .         .        
#>   [9,]  .     .     .     .          .          .         .         .        
#>  [10,]  .     .     .     .          .          .         .         .        
#>  [11,]  .     .     .     .          .          .         .         .        
#>  [12,]  .     .     .     .          .          .         .         .        
#>  [13,]  .     .     .     .          .          .         .         .        
#>  [14,]  .     .     .     .          .          .         .         .        
#>  [15,]  .     .     .     .          .          .         .         .        
#>  [16,]  .     .     .     .          .          .         .         .        
#>  [17,]  .     .     .     .          .          .         .         .        
#>  [18,]  .     .     .     .          .          .         .         .        
#>  [19,]  .     .     .     .          .          .         .         .        
#>  [20,]  .     .     .     .          .          .         .         .        
#>  [21,]  .     .     .     .          .          .         .         .        
#>  [22,]  .     .     .     .          .          .         .         .        
#>  [23,]  .     .     .     .          .          .         .         .        
#>  [24,]  .     .     .     .          .          .         .         .        
#>  [25,]  .     .     .     .          .          .         .         .        
#>  [26,]  .     .     .     .          .          .         .         .        
#>  [27,]  .     .     .     .          .          .         .         .        
#>  [28,]  .     .     .     .          .          .         .         .        
#>  [29,]  .     .     .     .          .          .         .         .        
#>  [30,]  .     .     .     .          .          .         .         .        
#>  [31,]  .     .     .     .          .          .         .         .        
#>  [32,]  .     .     .     .          .          .         .         .        
#>  [33,]  .     .     .     .          .          .         .         .        
#>  [34,]  .     .     .     .          .          .         .         .        
#>  [35,]  .     .     .     .          .          .         .         .        
#>  [36,]  .     .     .     .          .          .         .         .        
#>  [37,]  .     .     .     .          .          .         .         .        
#>  [38,]  .     .     .     .          .          .         .         .        
#>  [39,]  .     .     .     .          .          .         .         .        
#>  [40,]  .     .     .     .          .          .         .         .        
#>  [41,]  .     .     .     .          .          .         .         .        
#>  [42,]  .     .     .     .          .          .         .         .        
#>  [43,]  .     .     .     .          .          .         .         .        
#>  [44,]  .     .     .     .          .          .         .         .        
#>  [45,]  .     .     .     .          .          .         .         .        
#>  [46,]  .     .     .     .          .          .         .         .        
#>  [47,]  .     .     .     .          .          .         .         .        
#>  [48,]  .     .     .     .          .          .         .         .        
#>  [49,]  .     .     .     .          .          .         .         .        
#>  [50,]  .     .     .     .          .          .         .         .        
#>  [51,]  .     .     .     .          .          .         .         .        
#>  [52,]  .     .     .     .          .          .         .         .        
#>  [53,]  .     .     .     .          .          .         .         .        
#>  [54,]  .     .     .     .          .          .         .         .        
#>  [55,]  .     .     .     .          .          .         .         .        
#>  [56,]  .     .     .     .          .          .         .         .        
#>  [57,]  .     .     .     .          .          .         .         .        
#>  [58,]  .     .     .     .          .          .         .         .        
#>  [59,]  .     .     .     .          .          .         .         .        
#>  [60,]  .     .     .     .          .          .         .         .        
#>  [61,]  .     .     .     .          .          .         .         .        
#>  [62,]  .     .     .     .          .          .         .         .        
#>  [63,]  .     .     .     .          .          .         .         .        
#>  [64,]  .     .     .     .          .          .         .         .        
#>  [65,]  .     .     .     .          .          .         .         .        
#>  [66,]  .     .     .     .          .          .         .         .        
#>  [67,]  .     .     .     .          .          .         .         .        
#>  [68,]  .     .     .     .          .          .         .         .        
#>  [69,]  .     .     .     .          .          .         .         .        
#>  [70,]  .     .     .     .          .          .         .         .        
#>  [71,]  .     .     .     .          .          .         .         .        
#>  [72,]  .     .     .     .          .          .         .         .        
#>  [73,]  .     .     .     .          .          .         .         .        
#>  [74,]  .     .     .     .          .          .         .         .        
#>  [75,]  .     .     .     .          .          .         .         .        
#>  [76,]  .     .     .     .          .          .         .         .        
#>  [77,]  .     .     .     .          .          .         .         .        
#>  [78,]  .     .     .     .          .          .         .         .        
#>  [79,]  .     .     .     .          .          .         .         .        
#>  [80,]  .     .     .     .          .          .         .         .        
#>  [81,]  .     .     .     .          .          .         .         .        
#>  [82,]  .     .     .     .          .          .         .         .        
#>  [83,]  .     .     .     .          .          .         .         .        
#>  [84,]  .     .     .     .          .          .         .         .        
#>  [85,]  .     .     .     .          .          .         .         .        
#>  [86,]  .     .     .     .          .          .         .         .        
#>  [87,]  .     .     .     .          .          .         .         .        
#>  [88,]  .     .     .     .          .          .         .         .        
#>  [89,]  .     .     .     .          .          .         .         .        
#>  [90,]  .     .     .     .          .          .         .         .        
#>  [91,]  .     .     .     .          .          .         .         .        
#>  [92,]  .     .     .     .          .          .         .         .        
#>  [93,]  .     .     .     .          .          .         .         .        
#>  [94,]  .     .     .     .          .          .         .         .        
#>  [95,]  .     .     .     .          .          .         .         .        
#>  [96,]  .     .     .     .          .          .         .         .        
#>  [97,]  .     .     .     .          .          .         .         .        
#>  [98,]  .     .     .     .          .          .         .         .        
#>  [99,]  .     .     .     .          .          .         .         .        
#> [100,]  .     .     .     .          .          .         .         .        
#> [101,]  .     .     .     .          .          .         .         .        
#> [102,]  .     .     .     .          .          .         .         .        
#> [103,]  .     .     .     .          .          .         .         .        
#> [104,]  .     .     .     .          .          .         .         .        
#> [105,]  .     .     .     .          .          .         .         .        
#> [106,]  .     .     .     .          .          .         .         .        
#> [107,]  .     .     .     .          .          .         .         .        
#> [108,]  .     .     .     .          .          .         .         .        
#> [109,]  .     .     .     .          .          .         .         .        
#> [110,]  .     .     .     .          .          .         .         .        
#> [111,]  .     .     .     .          .          .         .         .        
#> [112,]  .     .     .     .          .          .         .         .        
#> [113,]  .     .     .     .          .          .         .         .        
#> [114,]  .     .     .     .          .          .         .         .        
#> [115,]  .     .     .     .          .          .         .         .        
#> [116,]  .     .     .     .          .          .         .         .        
#> [117,]  .     .     .     .          .          .         .         .        
#> [118,]  .     .     .     .          .          .         .         .        
#> [119,]  .     .     .     .          .          .         .         .        
#> [120,]  .     .     .     .          .          .         .         .        
#> [121,]  .     .     .     .          .          .         .         .        
#> [122,]  .     .     .     .          .          .         .         .        
#> [123,]  .     .     .     .          .          .         .         .        
#> [124,]  .     .     .     .          .          .         .         .        
#> [125,]  .     .     .     .          .          .         .         .        
#> [126,]  .     .     .     .          .          .         .         .        
#> [127,]  .     .     .     .          .          .         .         .        
#> [128,]  .     .     .     .          .          .         .         .        
#> [129,]  .     .     .     .          .          .         .         .        
#> [130,]  .     .     .     .          .          .         .         .        
#> [131,]  .     .     .     .          .          .         .         .        
#> [132,]  .     .     .     .          .          .         .         .        
#> [133,]  .     .     .     .          .          .         .         .        
#> [134,]  .     .     .     .          .          .         .         .        
#> [135,]  .     .     .     .          .          .         .         .        
#> [136,]  .     .     .     .          .          .         .         .        
#> [137,]  .     .     .     .          .          .         .         .        
#> [138,]  .     .     .     .          .          .         .         .        
#> [139,]  .     .     .     .          .          .         .         .        
#> [140,]  .     .     .     .          .          .         .         .        
#> [141,]  .     .     .     .          .          .         .         .        
#> [142,]  .     .     .     .          .          .         .         .        
#> [143,]  .     .     .     .          .          .         .         .        
#> [144,]  .     .     .     .          .          .         .         .        
#> [145,]  .     .     .     .          .          .         .         .        
#> [146,]  .     .     .     .          .          .         .         .        
#> [147,]  .     .     .     .          .          .         .         .        
#> [148,]  .     .     .     .          .          .         .         .        
#> [149,]  .     .     .     .          .          .         .         .        
#> [150,]  .     .     .     .          .          .         .         .        
#> [151,]  .     .     .     .          .          .         .         .        
#> [152,]  .     .     .     .          .          .         .         .        
#> [153,]  .     .     .     .          .          .         .         .        
#> [154,]  .     .     .     .          .          .         .         .        
#> [155,]  .     .     .     .          .          .         .         .        
#> [156,]  .     .     .     .          .          .         .         .        
#> [157,]  .     .     .     .          .          .         .         .        
#> [158,]  .     .     .     .          .          .         .         .        
#> [159,]  .     .     .     .          .          .         .         .        
#> [160,]  .     .     .     .          .          .         .         .        
#> [161,]  .     .     .     .          .          .         .         .        
#> [162,]  .     .    -0.49  .          .          .         .         .        
#> [163,]  .     .     .     .          .          .         .         .        
#> [164,]  .     .     .     .          .          .         .         .        
#> [165,]  .     .     .     .          .          .         .         .        
#> [166,]  .     .     .     .          .          .         .         .        
#> [167,] -0.49  .     .     .          .          .         .         .        
#> [168,]  0.98 -0.49  .     .          .          .         .         .        
#> [169,] -0.49  0.98 -0.49  .          .          .         .         .        
#> [170,]  .    -0.49  0.98  .          .          .         .         .        
#> [171,]  .     .     .     0.9583333 -0.9583333  .         .         .        
#> [172,]  .     .     .    -0.9583333  0.9583333  .         .         .        
#> [173,]  .     .     .     .          .          0.952381 -0.952381  .        
#> [174,]  .     .     .     .          .         -0.952381  0.952381  .        
#> [175,]  .     .     .     .          .          .         .         0.9819820
#> [176,]  .     .     .     .          .          .         .        -0.4864865
#> [177,]  .     .     .     .          .          .         .         .        
#> [178,]  .     .     .     .          .          .         .         .        
#> [179,]  .     .     .     .          .          .         .         .        
#> [180,]  .     .     .     .          .          .         .         .        
#> [181,]  .     .     .     .          .          .         .         .        
#> [182,]  .     .     .     .          .          .         .         .        
#> [183,]  .     .     .     .          .          .         .         .        
#> [184,]  .     .     .     .          .          .         .         .        
#> [185,]  .     .     .     .          .          .         .         .        
#> [186,]  .     .     .     .          .          .         .         .        
#> [187,]  .     .     .     .          .          .         .         .        
#> [188,]  .     .     .     .          .          .         .         .        
#> [189,]  .     .     .     .          .          .         .         .        
#> [190,]  .     .     .     .          .          .         .         .        
#> [191,]  .     .     .     .          .          .         .         .        
#> [192,]  .     .     .     .          .          .         .         .        
#> [193,]  .     .     .     .          .          .         .         .        
#> [194,]  .     .     .     .          .          .         .         .        
#> [195,]  .     .     .     .          .          .         .         .        
#> [196,]  .     .     .     .          .          .         .         .        
#> [197,]  .     .     .     .          .          .         .         .        
#> [198,]  .     .     .     .          .          .         .         .        
#> [199,]  .     .     .     .          .          .         .         .        
#> [200,]  .     .     .     .          .          .         .         .        
#> [201,]  .     .     .     .          .          .         .         .        
#> [202,]  .     .     .     .          .          .         .         .        
#> [203,]  .     .     .     .          .          .         .         .        
#> [204,]  .     .     .     .          .          .         .         .        
#> [205,]  .     .     .     .          .          .         .         .        
#> [206,]  .     .     .     .          .          .         .         .        
#> [207,]  .     .     .     .          .          .         .         .        
#> [208,]  .     .     .     .          .          .         .         .        
#> [209,]  .     .     .     .          .          .         .         .        
#> [210,]  .     .     .     .          .          .         .         .        
#> [211,]  .     .     .     .          .          .         .         .        
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,] -0.4864865  .          .          .          .          .        
#> [176,]  0.9819820  .          .          .          .          .        
#> [177,]  .          0.9784946 -0.4892473  .         -0.4892473  .        
#> [178,]  .         -0.4892473  0.9784946 -0.4892473  .          .        
#> [179,]  .          .         -0.4892473  0.9784946 -0.4892473  .        
#> [180,]  .         -0.4892473  .         -0.4892473  0.9784946  .        
#> [181,]  .          .          .          .          .          0.9820467
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .         -0.4892280
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  0.9820467 -0.4892280  .          .          .          .        
#> [183,] -0.4892280  0.9820467 -0.4892280  .          .          .        
#> [184,]  .         -0.4892280  0.9820467  .          .          .        
#> [185,]  .          .          .          0.9820467 -0.4892280  .        
#> [186,]  .          .          .         -0.4892280  0.9820467 -0.4892280
#> [187,]  .          .          .          .         -0.4892280  0.9820467
#> [188,]  .          .          .          .          .         -0.4892280
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,]  .          .          .          .          .          .        
#> [206,]  .          .          .          .          .          .        
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                       
#>   [1,]  .          .          .          .         .         .        
#>   [2,]  .          .          .          .         .         .        
#>   [3,]  .          .          .          .         .         .        
#>   [4,]  .          .          .          .         .         .        
#>   [5,]  .          .          .          .         .         .        
#>   [6,]  .          .          .          .         .         .        
#>   [7,]  .          .          .          .         .         .        
#>   [8,]  .          .          .          .         .         .        
#>   [9,]  .          .          .          .         .         .        
#>  [10,]  .          .          .          .         .         .        
#>  [11,]  .          .          .          .         .         .        
#>  [12,]  .          .          .          .         .         .        
#>  [13,]  .          .          .          .         .         .        
#>  [14,]  .          .          .          .         .         .        
#>  [15,]  .          .          .          .         .         .        
#>  [16,]  .          .          .          .         .         .        
#>  [17,]  .          .          .          .         .         .        
#>  [18,]  .          .          .          .         .         .        
#>  [19,]  .          .          .          .         .         .        
#>  [20,]  .          .          .          .         .         .        
#>  [21,]  .          .          .          .         .         .        
#>  [22,]  .          .          .          .         .         .        
#>  [23,]  .          .          .          .         .         .        
#>  [24,]  .          .          .          .         .         .        
#>  [25,]  .          .          .          .         .         .        
#>  [26,]  .          .          .          .         .         .        
#>  [27,]  .          .          .          .         .         .        
#>  [28,]  .          .          .          .         .         .        
#>  [29,]  .          .          .          .         .         .        
#>  [30,]  .          .          .          .         .         .        
#>  [31,]  .          .          .          .         .         .        
#>  [32,]  .          .          .          .         .         .        
#>  [33,]  .          .          .          .         .         .        
#>  [34,]  .          .          .          .         .         .        
#>  [35,]  .          .          .          .         .         .        
#>  [36,]  .          .          .          .         .         .        
#>  [37,]  .          .          .          .         .         .        
#>  [38,]  .          .          .          .         .         .        
#>  [39,]  .          .          .          .         .         .        
#>  [40,]  .          .          .          .         .         .        
#>  [41,]  .          .          .          .         .         .        
#>  [42,]  .          .          .          .         .         .        
#>  [43,]  .          .          .          .         .         .        
#>  [44,]  .          .          .          .         .         .        
#>  [45,]  .          .          .          .         .         .        
#>  [46,]  .          .          .          .         .         .        
#>  [47,]  .          .          .          .         .         .        
#>  [48,]  .          .          .          .         .         .        
#>  [49,]  .          .          .          .         .         .        
#>  [50,]  .          .          .          .         .         .        
#>  [51,]  .          .          .          .         .         .        
#>  [52,]  .          .          .          .         .         .        
#>  [53,]  .          .          .          .         .         .        
#>  [54,]  .          .          .          .         .         .        
#>  [55,]  .          .          .          .         .         .        
#>  [56,]  .          .          .          .         .         .        
#>  [57,]  .          .          .          .         .         .        
#>  [58,]  .          .          .          .         .         .        
#>  [59,]  .          .          .          .         .         .        
#>  [60,]  .          .          .          .         .         .        
#>  [61,]  .          .          .          .         .         .        
#>  [62,]  .          .          .          .         .         .        
#>  [63,]  .          .          .          .         .         .        
#>  [64,]  .          .          .          .         .         .        
#>  [65,]  .          .          .          .         .         .        
#>  [66,]  .          .          .          .         .         .        
#>  [67,]  .          .          .          .         .         .        
#>  [68,]  .          .          .          .         .         .        
#>  [69,]  .          .          .          .         .         .        
#>  [70,]  .          .          .          .         .         .        
#>  [71,]  .          .          .          .         .         .        
#>  [72,]  .          .          .          .         .         .        
#>  [73,]  .          .          .          .         .         .        
#>  [74,]  .          .          .          .         .         .        
#>  [75,]  .          .          .          .         .         .        
#>  [76,]  .          .          .          .         .         .        
#>  [77,]  .          .          .          .         .         .        
#>  [78,]  .          .          .          .         .         .        
#>  [79,]  .          .          .          .         .         .        
#>  [80,]  .          .          .          .         .         .        
#>  [81,]  .          .          .          .         .         .        
#>  [82,]  .          .          .          .         .         .        
#>  [83,]  .          .          .          .         .         .        
#>  [84,]  .          .          .          .         .         .        
#>  [85,]  .          .          .          .         .         .        
#>  [86,]  .          .          .          .         .         .        
#>  [87,]  .          .          .          .         .         .        
#>  [88,]  .          .          .          .         .         .        
#>  [89,]  .          .          .          .         .         .        
#>  [90,]  .          .          .          .         .         .        
#>  [91,]  .          .          .          .         .         .        
#>  [92,]  .          .          .          .         .         .        
#>  [93,]  .          .          .          .         .         .        
#>  [94,]  .          .          .          .         .         .        
#>  [95,]  .          .          .          .         .         .        
#>  [96,]  .          .          .          .         .         .        
#>  [97,]  .          .          .          .         .         .        
#>  [98,]  .          .          .          .         .         .        
#>  [99,]  .          .          .          .         .         .        
#> [100,]  .          .          .          .         .         .        
#> [101,]  .          .          .          .         .         .        
#> [102,]  .          .          .          .         .         .        
#> [103,]  .          .          .          .         .         .        
#> [104,]  .          .          .          .         .         .        
#> [105,]  .          .          .          .         .         .        
#> [106,]  .          .          .          .         .         .        
#> [107,]  .          .          .          .         .         .        
#> [108,]  .          .          .          .         .         .        
#> [109,]  .          .          .          .         .         .        
#> [110,]  .          .          .          .         .         .        
#> [111,]  .          .          .          .         .         .        
#> [112,]  .          .          .          .         .         .        
#> [113,]  .          .          .          .         .         .        
#> [114,]  .          .          .          .         .         .        
#> [115,]  .          .          .          .         .         .        
#> [116,]  .          .          .          .         .         .        
#> [117,]  .          .          .          .         .         .        
#> [118,]  .          .          .          .         .         .        
#> [119,]  .          .          .          .         .         .        
#> [120,]  .          .          .          .         .         .        
#> [121,]  .          .          .          .         .         .        
#> [122,]  .          .          .          .         .         .        
#> [123,]  .          .          .          .         .         .        
#> [124,]  .          .          .          .         .         .        
#> [125,]  .          .          .          .         .         .        
#> [126,]  .          .          .          .         .         .        
#> [127,]  .          .          .          .         .         .        
#> [128,]  .          .          .          .         .         .        
#> [129,]  .          .          .          .         .         .        
#> [130,]  .          .          .          .         .         .        
#> [131,]  .          .          .          .         .         .        
#> [132,]  .          .          .          .         .         .        
#> [133,]  .          .          .          .         .         .        
#> [134,]  .          .          .          .         .         .        
#> [135,]  .          .          .          .         .         .        
#> [136,]  .          .          .          .         .         .        
#> [137,]  .          .          .          .         .         .        
#> [138,]  .          .          .          .         .         .        
#> [139,]  .          .          .          .         .         .        
#> [140,]  .          .          .          .         .         .        
#> [141,]  .          .          .          .         .         .        
#> [142,]  .          .          .          .         .         .        
#> [143,]  .          .          .          .         .         .        
#> [144,]  .          .          .          .         .         .        
#> [145,]  .          .          .          .         .         .        
#> [146,]  .          .          .          .         .         .        
#> [147,]  .          .          .          .         .         .        
#> [148,]  .          .          .          .         .         .        
#> [149,]  .          .          .          .         .         .        
#> [150,]  .          .          .          .         .         .        
#> [151,]  .          .          .          .         .         .        
#> [152,]  .          .          .          .         .         .        
#> [153,]  .          .          .          .         .         .        
#> [154,]  .          .          .          .         .         .        
#> [155,]  .          .          .          .         .         .        
#> [156,]  .          .          .          .         .         .        
#> [157,]  .          .          .          .         .         .        
#> [158,]  .          .          .          .         .         .        
#> [159,]  .          .          .          .         .         .        
#> [160,]  .          .          .          .         .         .        
#> [161,]  .          .          .          .         .         .        
#> [162,]  .          .          .          .         .         .        
#> [163,]  .          .          .          .         .         .        
#> [164,]  .          .          .          .         .         .        
#> [165,]  .          .          .          .         .         .        
#> [166,]  .          .          .          .         .         .        
#> [167,]  .          .          .          .         .         .        
#> [168,]  .          .          .          .         .         .        
#> [169,]  .          .          .          .         .         .        
#> [170,]  .          .          .          .         .         .        
#> [171,]  .          .          .          .         .         .        
#> [172,]  .          .          .          .         .         .        
#> [173,]  .          .          .          .         .         .        
#> [174,]  .          .          .          .         .         .        
#> [175,]  .          .          .          .         .         .        
#> [176,]  .          .          .          .         .         .        
#> [177,]  .          .          .          .         .         .        
#> [178,]  .          .          .          .         .         .        
#> [179,]  .          .          .          .         .         .        
#> [180,]  .          .          .          .         .         .        
#> [181,]  .          .         -0.4892280  .         .         .        
#> [182,]  .          .          .          .         .         .        
#> [183,]  .          .          .          .         .         .        
#> [184,]  .          .          .          .         .         .        
#> [185,]  .          .          .          .         .         .        
#> [186,]  .          .          .          .         .         .        
#> [187,] -0.4892280  .          .          .         .         .        
#> [188,]  0.9820467 -0.4892280  .          .         .         .        
#> [189,] -0.4892280  0.9820467 -0.4892280  .         .         .        
#> [190,]  .         -0.4892280  0.9820467  .         .         .        
#> [191,]  .          .          .          0.971831 -0.971831  .        
#> [192,]  .          .          .         -0.971831  0.971831  .        
#> [193,]  .          .          .          .         .         0.9784946
#> [194,]  .          .          .          .         .        -0.9784946
#> [195,]  .          .          .          .         .         .        
#> [196,]  .          .          .          .         .         .        
#> [197,]  .          .          .          .         .         .        
#> [198,]  .          .          .          .         .         .        
#> [199,]  .          .          .          .         .         .        
#> [200,]  .          .          .          .         .         .        
#> [201,]  .          .          .          .         .         .        
#> [202,]  .          .          .          .         .         .        
#> [203,]  .          .          .          .         .         .        
#> [204,]  .          .          .          .         .         .        
#> [205,]  .          .          .          .         .         .        
#> [206,]  .          .          .          .         .         .        
#> [207,]  .          .          .          .         .         .        
#> [208,]  .          .          .          .         .         .        
#> [209,]  .          .          .          .         .         .        
#> [210,]  .          .          .          .         .         .        
#> [211,]  .          .          .          .         .         .        
#>                                                                          
#>   [1,]  .         .            .          .          .          .        
#>   [2,]  .         .            .          .          .          .        
#>   [3,]  .         .            .          .          .          .        
#>   [4,]  .         .            .          .          .          .        
#>   [5,]  .         .            .          .          .          .        
#>   [6,]  .         .            .          .          .          .        
#>   [7,]  .         .            .          .          .          .        
#>   [8,]  .         .            .          .          .          .        
#>   [9,]  .         .            .          .          .          .        
#>  [10,]  .         .            .          .          .          .        
#>  [11,]  .         .            .          .          .          .        
#>  [12,]  .         .            .          .          .          .        
#>  [13,]  .         .            .          .          .          .        
#>  [14,]  .         .            .          .          .          .        
#>  [15,]  .         .            .          .          .          .        
#>  [16,]  .         .            .          .          .          .        
#>  [17,]  .         .            .          .          .          .        
#>  [18,]  .         .            .          .          .          .        
#>  [19,]  .         .            .          .          .          .        
#>  [20,]  .         .            .          .          .          .        
#>  [21,]  .         .            .          .          .          .        
#>  [22,]  .         .            .          .          .          .        
#>  [23,]  .         .            .          .          .          .        
#>  [24,]  .         .            .          .          .          .        
#>  [25,]  .         .            .          .          .          .        
#>  [26,]  .         .            .          .          .          .        
#>  [27,]  .         .            .          .          .          .        
#>  [28,]  .         .            .          .          .          .        
#>  [29,]  .         .            .          .          .          .        
#>  [30,]  .         .            .          .          .          .        
#>  [31,]  .         .            .          .          .          .        
#>  [32,]  .         .            .          .          .          .        
#>  [33,]  .         .            .          .          .          .        
#>  [34,]  .         .            .          .          .          .        
#>  [35,]  .         .            .          .          .          .        
#>  [36,]  .         .            .          .          .          .        
#>  [37,]  .         .            .          .          .          .        
#>  [38,]  .         .            .          .          .          .        
#>  [39,]  .         .            .          .          .          .        
#>  [40,]  .         .            .          .          .          .        
#>  [41,]  .         .            .          .          .          .        
#>  [42,]  .         .            .          .          .          .        
#>  [43,]  .         .            .          .          .          .        
#>  [44,]  .         .            .          .          .          .        
#>  [45,]  .         .            .          .          .          .        
#>  [46,]  .         .            .          .          .          .        
#>  [47,]  .         .            .          .          .          .        
#>  [48,]  .         .            .          .          .          .        
#>  [49,]  .         .            .          .          .          .        
#>  [50,]  .         .            .          .          .          .        
#>  [51,]  .         .            .          .          .          .        
#>  [52,]  .         .            .          .          .          .        
#>  [53,]  .         .            .          .          .          .        
#>  [54,]  .         .            .          .          .          .        
#>  [55,]  .         .            .          .          .          .        
#>  [56,]  .         .            .          .          .          .        
#>  [57,]  .         .            .          .          .          .        
#>  [58,]  .         .            .          .          .          .        
#>  [59,]  .         .            .          .          .          .        
#>  [60,]  .         .            .          .          .          .        
#>  [61,]  .         .            .          .          .          .        
#>  [62,]  .         .            .          .          .          .        
#>  [63,]  .         .            .          .          .          .        
#>  [64,]  .         .            .          .          .          .        
#>  [65,]  .         .            .          .          .          .        
#>  [66,]  .         .            .          .          .          .        
#>  [67,]  .         .            .          .          .          .        
#>  [68,]  .         .            .          .          .          .        
#>  [69,]  .         .            .          .          .          .        
#>  [70,]  .         .            .          .          .          .        
#>  [71,]  .         .            .          .          .          .        
#>  [72,]  .         .            .          .          .          .        
#>  [73,]  .         .            .          .          .          .        
#>  [74,]  .         .            .          .          .          .        
#>  [75,]  .         .            .          .          .          .        
#>  [76,]  .         .            .          .          .          .        
#>  [77,]  .         .            .          .          .          .        
#>  [78,]  .         .            .          .          .          .        
#>  [79,]  .         .            .          .          .          .        
#>  [80,]  .         .            .          .          .          .        
#>  [81,]  .         .            .          .          .          .        
#>  [82,]  .         .            .          .          .          .        
#>  [83,]  .         .            .          .          .          .        
#>  [84,]  .         .            .          .          .          .        
#>  [85,]  .         .            .          .          .          .        
#>  [86,]  .         .            .          .          .          .        
#>  [87,]  .         .            .          .          .          .        
#>  [88,]  .         .            .          .          .          .        
#>  [89,]  .         .            .          .          .          .        
#>  [90,]  .         .            .          .          .          .        
#>  [91,]  .         .            .          .          .          .        
#>  [92,]  .         .            .          .          .          .        
#>  [93,]  .         .            .          .          .          .        
#>  [94,]  .         .            .          .          .          .        
#>  [95,]  .         .            .          .          .          .        
#>  [96,]  .         .            .          .          .          .        
#>  [97,]  .         .            .          .          .          .        
#>  [98,]  .         .            .          .          .          .        
#>  [99,]  .         .            .          .          .          .        
#> [100,]  .         .            .          .          .          .        
#> [101,]  .         .            .          .          .          .        
#> [102,]  .         .            .          .          .          .        
#> [103,]  .         .            .          .          .          .        
#> [104,]  .         .            .          .          .          .        
#> [105,]  .         .            .          .          .          .        
#> [106,]  .         .            .          .          .          .        
#> [107,]  .         .            .          .          .          .        
#> [108,]  .         .            .          .          .          .        
#> [109,]  .         .            .          .          .          .        
#> [110,]  .         .            .          .          .          .        
#> [111,]  .         .            .          .          .          .        
#> [112,]  .         .            .          .          .          .        
#> [113,]  .         .            .          .          .          .        
#> [114,]  .         .            .          .          .          .        
#> [115,]  .         .            .          .          .          .        
#> [116,]  .         .            .          .          .          .        
#> [117,]  .         .            .          .          .          .        
#> [118,]  .         .            .          .          .          .        
#> [119,]  .         .            .          .          .          .        
#> [120,]  .         .            .          .          .          .        
#> [121,]  .         .            .          .          .          .        
#> [122,]  .         .            .          .          .          .        
#> [123,]  .         .            .          .          .          .        
#> [124,]  .         .            .          .          .          .        
#> [125,]  .         .            .          .          .          .        
#> [126,]  .         .            .          .          .          .        
#> [127,]  .         .            .          .          .          .        
#> [128,]  .         .            .          .          .          .        
#> [129,]  .         .            .          .          .          .        
#> [130,]  .         .            .          .          .          .        
#> [131,]  .         .            .          .          .          .        
#> [132,]  .         .            .          .          .          .        
#> [133,]  .         .            .          .          .          .        
#> [134,]  .         .            .          .          .          .        
#> [135,]  .         .            .          .          .          .        
#> [136,]  .         .            .          .          .          .        
#> [137,]  .         .            .          .          .          .        
#> [138,]  .         .            .          .          .          .        
#> [139,]  .         .            .          .          .          .        
#> [140,]  .         .            .          .          .          .        
#> [141,]  .         .            .          .          .          .        
#> [142,]  .         .            .          .          .          .        
#> [143,]  .         .            .          .          .          .        
#> [144,]  .         .            .          .          .          .        
#> [145,]  .         .            .          .          .          .        
#> [146,]  .         .            .          .          .          .        
#> [147,]  .         .            .          .          .          .        
#> [148,]  .         .            .          .          .          .        
#> [149,]  .         .            .          .          .          .        
#> [150,]  .         .            .          .          .          .        
#> [151,]  .         .            .          .          .          .        
#> [152,]  .         .            .          .          .          .        
#> [153,]  .         .            .          .          .          .        
#> [154,]  .         .            .          .          .          .        
#> [155,]  .         .            .          .          .          .        
#> [156,]  .         .            .          .          .          .        
#> [157,]  .         .            .          .          .          .        
#> [158,]  .         .            .          .          .          .        
#> [159,]  .         .            .          .          .          .        
#> [160,]  .         .            .          .          .          .        
#> [161,]  .         .            .          .          .          .        
#> [162,]  .         .            .          .          .          .        
#> [163,]  .         .            .          .          .          .        
#> [164,]  .         .            .          .          .          .        
#> [165,]  .         .            .          .          .          .        
#> [166,]  .         .            .          .          .          .        
#> [167,]  .         .            .          .          .          .        
#> [168,]  .         .            .          .          .          .        
#> [169,]  .         .            .          .          .          .        
#> [170,]  .         .            .          .          .          .        
#> [171,]  .         .            .          .          .          .        
#> [172,]  .         .            .          .          .          .        
#> [173,]  .         .            .          .          .          .        
#> [174,]  .         .            .          .          .          .        
#> [175,]  .         .            .          .          .          .        
#> [176,]  .         .            .          .          .          .        
#> [177,]  .         .            .          .          .          .        
#> [178,]  .         .            .          .          .          .        
#> [179,]  .         .            .          .          .          .        
#> [180,]  .         .            .          .          .          .        
#> [181,]  .         .            .          .          .          .        
#> [182,]  .         .            .          .          .          .        
#> [183,]  .         .            .          .          .          .        
#> [184,]  .         .            .          .          .          .        
#> [185,]  .         .            .          .          .          .        
#> [186,]  .         .            .          .          .          .        
#> [187,]  .         .            .          .          .          .        
#> [188,]  .         .            .          .          .          .        
#> [189,]  .         .            .          .          .          .        
#> [190,]  .         .            .          .          .          .        
#> [191,]  .         .            .          .          .          .        
#> [192,]  .         .            .          .          .          .        
#> [193,] -0.9784946 .            .          .          .          .        
#> [194,]  0.9784946 .            .          .          .          .        
#> [195,]  .         6.18486e-10  .          .          .          .        
#> [196,]  .         .            0.9836957 -0.4891304  .          .        
#> [197,]  .         .           -0.4891304  0.9836957 -0.4891304  .        
#> [198,]  .         .            .         -0.4891304  0.9836957  .        
#> [199,]  .         .            .          .          .          0.9666667
#> [200,]  .         .            .          .          .         -0.9666667
#> [201,]  .         .            .          .          .          .        
#> [202,]  .         .            .          .          .          .        
#> [203,]  .         .            .          .          .          .        
#> [204,]  .         .            .          .          .          .        
#> [205,]  .         .            .          .          .          .        
#> [206,]  .         .            .          .          .          .        
#> [207,]  .         .            .          .          .          .        
#> [208,]  .         .            .          .          .          .        
#> [209,]  .         .            .          .          .          .        
#> [210,]  .         .            .          .          .          .        
#> [211,]  .         .            .          .          .          .        
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,] -0.9666667  .          .          .          .          .        
#> [200,]  0.9666667  .          .          .          .          .        
#> [201,]  .          0.9816273 -0.4895013  .          .          .        
#> [202,]  .         -0.4895013  0.9816273 -0.4895013  .          .        
#> [203,]  .          .         -0.4895013  0.9816273 -0.4895013  .        
#> [204,]  .          .          .         -0.4895013  0.9816273 -0.4895013
#> [205,]  .          .          .          .         -0.4895013  0.9816273
#> [206,]  .          .          .          .          .         -0.4895013
#> [207,]  .          .          .          .          .          .        
#> [208,]  .          .          .          .          .          .        
#> [209,]  .          .          .          .          .          .        
#> [210,]  .          .          .          .          .          .        
#> [211,]  .          .          .          .          .          .        
#>                                                                         
#>   [1,]  .          .          .          .          .          .        
#>   [2,]  .          .          .          .          .          .        
#>   [3,]  .          .          .          .          .          .        
#>   [4,]  .          .          .          .          .          .        
#>   [5,]  .          .          .          .          .          .        
#>   [6,]  .          .          .          .          .          .        
#>   [7,]  .          .          .          .          .          .        
#>   [8,]  .          .          .          .          .          .        
#>   [9,]  .          .          .          .          .          .        
#>  [10,]  .          .          .          .          .          .        
#>  [11,]  .          .          .          .          .          .        
#>  [12,]  .          .          .          .          .          .        
#>  [13,]  .          .          .          .          .          .        
#>  [14,]  .          .          .          .          .          .        
#>  [15,]  .          .          .          .          .          .        
#>  [16,]  .          .          .          .          .          .        
#>  [17,]  .          .          .          .          .          .        
#>  [18,]  .          .          .          .          .          .        
#>  [19,]  .          .          .          .          .          .        
#>  [20,]  .          .          .          .          .          .        
#>  [21,]  .          .          .          .          .          .        
#>  [22,]  .          .          .          .          .          .        
#>  [23,]  .          .          .          .          .          .        
#>  [24,]  .          .          .          .          .          .        
#>  [25,]  .          .          .          .          .          .        
#>  [26,]  .          .          .          .          .          .        
#>  [27,]  .          .          .          .          .          .        
#>  [28,]  .          .          .          .          .          .        
#>  [29,]  .          .          .          .          .          .        
#>  [30,]  .          .          .          .          .          .        
#>  [31,]  .          .          .          .          .          .        
#>  [32,]  .          .          .          .          .          .        
#>  [33,]  .          .          .          .          .          .        
#>  [34,]  .          .          .          .          .          .        
#>  [35,]  .          .          .          .          .          .        
#>  [36,]  .          .          .          .          .          .        
#>  [37,]  .          .          .          .          .          .        
#>  [38,]  .          .          .          .          .          .        
#>  [39,]  .          .          .          .          .          .        
#>  [40,]  .          .          .          .          .          .        
#>  [41,]  .          .          .          .          .          .        
#>  [42,]  .          .          .          .          .          .        
#>  [43,]  .          .          .          .          .          .        
#>  [44,]  .          .          .          .          .          .        
#>  [45,]  .          .          .          .          .          .        
#>  [46,]  .          .          .          .          .          .        
#>  [47,]  .          .          .          .          .          .        
#>  [48,]  .          .          .          .          .          .        
#>  [49,]  .          .          .          .          .          .        
#>  [50,]  .          .          .          .          .          .        
#>  [51,]  .          .          .          .          .          .        
#>  [52,]  .          .          .          .          .          .        
#>  [53,]  .          .          .          .          .          .        
#>  [54,]  .          .          .          .          .          .        
#>  [55,]  .          .          .          .          .          .        
#>  [56,]  .          .          .          .          .          .        
#>  [57,]  .          .          .          .          .          .        
#>  [58,]  .          .          .          .          .          .        
#>  [59,]  .          .          .          .          .          .        
#>  [60,]  .          .          .          .          .          .        
#>  [61,]  .          .          .          .          .          .        
#>  [62,]  .          .          .          .          .          .        
#>  [63,]  .          .          .          .          .          .        
#>  [64,]  .          .          .          .          .          .        
#>  [65,]  .          .          .          .          .          .        
#>  [66,]  .          .          .          .          .          .        
#>  [67,]  .          .          .          .          .          .        
#>  [68,]  .          .          .          .          .          .        
#>  [69,]  .          .          .          .          .          .        
#>  [70,]  .          .          .          .          .          .        
#>  [71,]  .          .          .          .          .          .        
#>  [72,]  .          .          .          .          .          .        
#>  [73,]  .          .          .          .          .          .        
#>  [74,]  .          .          .          .          .          .        
#>  [75,]  .          .          .          .          .          .        
#>  [76,]  .          .          .          .          .          .        
#>  [77,]  .          .          .          .          .          .        
#>  [78,]  .          .          .          .          .          .        
#>  [79,]  .          .          .          .          .          .        
#>  [80,]  .          .          .          .          .          .        
#>  [81,]  .          .          .          .          .          .        
#>  [82,]  .          .          .          .          .          .        
#>  [83,]  .          .          .          .          .          .        
#>  [84,]  .          .          .          .          .          .        
#>  [85,]  .          .          .          .          .          .        
#>  [86,]  .          .          .          .          .          .        
#>  [87,]  .          .          .          .          .          .        
#>  [88,]  .          .          .          .          .          .        
#>  [89,]  .          .          .          .          .          .        
#>  [90,]  .          .          .          .          .          .        
#>  [91,]  .          .          .          .          .          .        
#>  [92,]  .          .          .          .          .          .        
#>  [93,]  .          .          .          .          .          .        
#>  [94,]  .          .          .          .          .          .        
#>  [95,]  .          .          .          .          .          .        
#>  [96,]  .          .          .          .          .          .        
#>  [97,]  .          .          .          .          .          .        
#>  [98,]  .          .          .          .          .          .        
#>  [99,]  .          .          .          .          .          .        
#> [100,]  .          .          .          .          .          .        
#> [101,]  .          .          .          .          .          .        
#> [102,]  .          .          .          .          .          .        
#> [103,]  .          .          .          .          .          .        
#> [104,]  .          .          .          .          .          .        
#> [105,]  .          .          .          .          .          .        
#> [106,]  .          .          .          .          .          .        
#> [107,]  .          .          .          .          .          .        
#> [108,]  .          .          .          .          .          .        
#> [109,]  .          .          .          .          .          .        
#> [110,]  .          .          .          .          .          .        
#> [111,]  .          .          .          .          .          .        
#> [112,]  .          .          .          .          .          .        
#> [113,]  .          .          .          .          .          .        
#> [114,]  .          .          .          .          .          .        
#> [115,]  .          .          .          .          .          .        
#> [116,]  .          .          .          .          .          .        
#> [117,]  .          .          .          .          .          .        
#> [118,]  .          .          .          .          .          .        
#> [119,]  .          .          .          .          .          .        
#> [120,]  .          .          .          .          .          .        
#> [121,]  .          .          .          .          .          .        
#> [122,]  .          .          .          .          .          .        
#> [123,]  .          .          .          .          .          .        
#> [124,]  .          .          .          .          .          .        
#> [125,]  .          .          .          .          .          .        
#> [126,]  .          .          .          .          .          .        
#> [127,]  .          .          .          .          .          .        
#> [128,]  .          .          .          .          .          .        
#> [129,]  .          .          .          .          .          .        
#> [130,]  .          .          .          .          .          .        
#> [131,]  .          .          .          .          .          .        
#> [132,]  .          .          .          .          .          .        
#> [133,]  .          .          .          .          .          .        
#> [134,]  .          .          .          .          .          .        
#> [135,]  .          .          .          .          .          .        
#> [136,]  .          .          .          .          .          .        
#> [137,]  .          .          .          .          .          .        
#> [138,]  .          .          .          .          .          .        
#> [139,]  .          .          .          .          .          .        
#> [140,]  .          .          .          .          .          .        
#> [141,]  .          .          .          .          .          .        
#> [142,]  .          .          .          .          .          .        
#> [143,]  .          .          .          .          .          .        
#> [144,]  .          .          .          .          .          .        
#> [145,]  .          .          .          .          .          .        
#> [146,]  .          .          .          .          .          .        
#> [147,]  .          .          .          .          .          .        
#> [148,]  .          .          .          .          .          .        
#> [149,]  .          .          .          .          .          .        
#> [150,]  .          .          .          .          .          .        
#> [151,]  .          .          .          .          .          .        
#> [152,]  .          .          .          .          .          .        
#> [153,]  .          .          .          .          .          .        
#> [154,]  .          .          .          .          .          .        
#> [155,]  .          .          .          .          .          .        
#> [156,]  .          .          .          .          .          .        
#> [157,]  .          .          .          .          .          .        
#> [158,]  .          .          .          .          .          .        
#> [159,]  .          .          .          .          .          .        
#> [160,]  .          .          .          .          .          .        
#> [161,]  .          .          .          .          .          .        
#> [162,]  .          .          .          .          .          .        
#> [163,]  .          .          .          .          .          .        
#> [164,]  .          .          .          .          .          .        
#> [165,]  .          .          .          .          .          .        
#> [166,]  .          .          .          .          .          .        
#> [167,]  .          .          .          .          .          .        
#> [168,]  .          .          .          .          .          .        
#> [169,]  .          .          .          .          .          .        
#> [170,]  .          .          .          .          .          .        
#> [171,]  .          .          .          .          .          .        
#> [172,]  .          .          .          .          .          .        
#> [173,]  .          .          .          .          .          .        
#> [174,]  .          .          .          .          .          .        
#> [175,]  .          .          .          .          .          .        
#> [176,]  .          .          .          .          .          .        
#> [177,]  .          .          .          .          .          .        
#> [178,]  .          .          .          .          .          .        
#> [179,]  .          .          .          .          .          .        
#> [180,]  .          .          .          .          .          .        
#> [181,]  .          .          .          .          .          .        
#> [182,]  .          .          .          .          .          .        
#> [183,]  .          .          .          .          .          .        
#> [184,]  .          .          .          .          .          .        
#> [185,]  .          .          .          .          .          .        
#> [186,]  .          .          .          .          .          .        
#> [187,]  .          .          .          .          .          .        
#> [188,]  .          .          .          .          .          .        
#> [189,]  .          .          .          .          .          .        
#> [190,]  .          .          .          .          .          .        
#> [191,]  .          .          .          .          .          .        
#> [192,]  .          .          .          .          .          .        
#> [193,]  .          .          .          .          .          .        
#> [194,]  .          .          .          .          .          .        
#> [195,]  .          .          .          .          .          .        
#> [196,]  .          .          .          .          .          .        
#> [197,]  .          .          .          .          .          .        
#> [198,]  .          .          .          .          .          .        
#> [199,]  .          .          .          .          .          .        
#> [200,]  .          .          .          .          .          .        
#> [201,]  .          .          .          .          .          .        
#> [202,]  .          .          .          .          .          .        
#> [203,]  .          .          .          .          .          .        
#> [204,]  .          .          .          .          .          .        
#> [205,] -0.4895013  .          .          .          .          .        
#> [206,]  0.9816273 -0.4895013  .          .          .          .        
#> [207,] -0.4895013  0.9816273  .          .          .          .        
#> [208,]  .          .          0.9791667 -0.9791667  .          .        
#> [209,]  .          .         -0.9791667  0.9791667  .          .        
#> [210,]  .          .          .          .          0.9130435 -0.9130435
#> [211,]  .          .          .          .         -0.9130435  0.9130435
# }
```
