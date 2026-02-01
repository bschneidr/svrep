# Quadratic Form Matrix for a Two-phase Design

This function combines quadratic forms from each phase of a two phase
design, so that the combined variance of the entire two-phase sampling
design can be estimated.

## Usage

``` r
make_twophase_quad_form(
  sigma_1,
  sigma_2,
  phase_2_joint_probs,
  ensure_psd = TRUE
)
```

## Arguments

- sigma_1:

  The quadratic form for the first phase variance estimator, subsetted
  to only include cases selected in the phase two sample.

- sigma_2:

  The quadratic form for the second phase variance estimator,
  conditional on the selection of the first phase sample.

- phase_2_joint_probs:

  The matrix of conditional joint inclusion probabilities for the second
  phase, given the selected first phase sample.

- ensure_psd:

  If `TRUE` (the default), ensures that the result is a positive
  semidefinite matrix. This is necessary if the quadratic form is used
  as an input for replication methods such as the generalized bootstrap.
  For details, see the help section entitled "Ensuring the Result is
  Positive Semidefinite".

## Value

A quadratic form matrix that can be used to estimate the sampling
variance from a two-phase sample design.

## Statistical Details

The two-phase variance estimator has a quadratic form matrix
\\\boldsymbol{\Sigma}\_{ab}\\ given by: \$\$ \boldsymbol{\Sigma}\_{ab} =
{W}^{-1}\_b(\boldsymbol{\Sigma}\_{a^\prime} \circ D_b ){W}^{-1}\_b +
\boldsymbol{\Sigma}\_b \$\$ The first term estimates the variance
contribution from the first phase of sampling, while the second term
estimates the variance contribution from the second phase of sampling.  

The full quadratic form of the variance estimator is: \$\$ v(\hat{t_y})
= \breve{\breve{y^{'}}} \boldsymbol{\Sigma}\_{ab} \breve{\breve{y}} \$\$
where the weighted variable \\\breve{\breve{y}}\_k =
\frac{y_k}{\pi\_{ak}\pi\_{bk}}\\, is formed using the first phase
inclusion probability, denoted \\\pi\_{ak}\\, and the conditional second
phase inclusion probability (given the selected first phase sample),
denoted \\\pi\_{bk}\\.  

The notation for this estimator is as follows:  

- \\n_a\\ denotes the first phase sample size.

- \\n_b\\ denotes the second phase sample size.

- \\\boldsymbol{\Sigma}\_a\\ denotes the matrix of dimension \\n_a
  \times n_a\\ representing the quadratic form for the variance
  estimator used for the full first-phase design.

- \\\boldsymbol{\Sigma}\_{a^\prime}\\ denotes the matrix of dimension
  \\n_b \times n_b\\ formed by subsetting the rows and columns of
  \\\boldsymbol{\Sigma}\_a\\ to only include cases selected in the
  second-phase sample.

- \\\boldsymbol{\Sigma}\_{b}\\ denotes the matrix of dimension \\n_b
  \times n_b\\ representing the Horvitz-Thompson estimator of variance
  for the second-phase sample, conditional on the selected first-phase
  sample.

- \\\boldsymbol{D}\_b\\ denotes the \\n_b \times n_b\\ matrix of weights
  formed by the inverses of the second-phase joint inclusion
  probabilities, with element \\kl\\ equal to \\\pi\_{bkl}^{-1}\\, where
  \\\pi\_{bkl}\\ is the conditional probability that units \\k\\ and
  \\l\\ are included in the second-phase sample, given the selected
  first-phase sample. Note that this matrix will often not be positive
  semidefinite, and so the two-phase variance estimator has a quadratic
  form which is not necessarily positive semidefinite.

- \\\boldsymbol{W}\_b\\ denotes the diagonal \\n_b \times n_b\\ matrix
  whose \\k\\-th diagonal entry is the second-phase weight
  \\\pi\_{bk}^{-1}\\, where \\\pi\_{bk}\\ is the conditional probability
  that unit \\k\\ is included in the second-phase sample, given the
  selected first-phase sample.

## Ensuring the Result is Positive Semidefinite

Note that the matrix \\(\boldsymbol{\Sigma}\_{a^\prime} \circ D_b )\\
may not be positive semidefinite, since the matrix \\D_b\\ is not
guaranteed to be positive semidefinite. If
\\(\boldsymbol{\Sigma}\_{a^\prime} \circ D_b )\\ is found not to be
positive semidefinite, then it is approximated by the nearest positive
semidefinite matrix in the Frobenius norm, using the method of Higham
(1988).  
  
This approximation is discussed by Beaumont and Patak (2012) in the
context of forming replicate weights for two-phase samples. The authors
argue that this approximation should lead to only a small overestimation
of variance.  
  
Since \\(\boldsymbol{\Sigma}\_{a^\prime} \circ D_b )\\ is a real,
symmetric matrix, this is equivalent to "zeroing out" negative
eigenvalues. To be more precise, denote
\\A=(\boldsymbol{\Sigma}\_{a^\prime} \circ D_b )\\. Then we can form the
spectral decomposition \\A=\Gamma \Lambda \Gamma^{\prime}\\, where
\\\Lambda\\ is the diagonal matrix whose entries are eigenvalues of
\\A\\. The method of Higham (1988) is to approximate \\A\\ with
\\\tilde{A} = \Gamma \Lambda\_{+} \Gamma^{\prime}\\, where the \\ii\\-th
entry of \\\Lambda\_{+}\\ is \\\max(\Lambda\_{ii}, 0)\\.

## References

See Section 7.5 of Tillé (2020) or Section 9.3 of Särndal, Swensson, and
Wretman (1992) for an overview of variance estimation for two-phase
sampling. In the case where the Horvitz-Thompson variance estimator is
used for both phases, the method used in this function is equivalent to
equation (9.3.8) of Särndal, Swensson, and Wretman (1992) and equation
(7.7) of Tillé (2020). However, this function can be used for any
combination of first-phase and second-phase variance estimators,
provided that the joint inclusion probabilities from the second-phase
design are available and are all nonzero.  
  

- Beaumont, Jean-François, and Zdenek Patak. (2012). "On the Generalized
  Bootstrap for Sample Surveys with Special Attention to Poisson
  Sampling: Generalized Bootstrap for Sample Surveys." International
  Statistical Review 80 (1): 127-48.  
    

- Higham, N. J. (1988). "*Computing a nearest symmetric positive
  semidefinite matrix.*" Linear Algebra and Its Applications, 103,
  103-118.  
    

- Särndal, C.-E., Swensson, B., & Wretman, J. (1992). "*Model Assisted
  Survey Sampling*." Springer New York.  
    

- Tillé, Y. (2020). "*Sampling and estimation from finite populations*."
  (I. Hekimi, Trans.). Wiley.

## See also

For each phase of sampling, the function
[make_quad_form_matrix](https://bschneidr.github.io/svrep/reference/make_quad_form_matrix.md)
can be used to create the appropriate quadratic form matrix.

## Examples

``` r
# \donttest{
library(survey)

## ---------------------- Example 1 ------------------------##
## First phase is a stratified multistage sample            ##
## Second phase is a simple random sample                   ##
##----------------------------------------------------------##
data('library_multistage_sample', package = 'svrep')

# Load first-phase sample
  twophase_sample <- library_multistage_sample

# Select second-phase sample
  set.seed(2022)

  twophase_sample[['SECOND_PHASE_SELECTION']] <- sampling::srswor(
    n = 100,
    N = nrow(twophase_sample)
  ) |> as.logical()

# Declare survey design
  twophase_design <- twophase(
    method = "full",
    data = twophase_sample,
    # Identify the subset of first-phase elements
    # which were selected into the second-phase sample
    subset = ~ SECOND_PHASE_SELECTION,
    # Describe clusters, probabilities, and population sizes
    # at each phase of sampling
    id = list(~ PSU_ID + SSU_ID,
              ~ 1),
    probs = list(~ PSU_SAMPLING_PROB + SSU_SAMPLING_PROB,
                 NULL),
    fpc = list(~ PSU_POP_SIZE + SSU_POP_SIZE,
               NULL)
  )

# Get quadratic form matrix for the first phase design
  first_phase_sigma <- get_design_quad_form(
    design = twophase_design$phase1$full,
    variance_estimator = "Stratified Multistage SRS"
  )

# Subset to only include cases sampled in second phase

  first_phase_sigma <- first_phase_sigma[twophase_design$subset,
                                         twophase_design$subset]

# Get quadratic form matrix for the second-phase design
  second_phase_sigma <- get_design_quad_form(
    design = twophase_design$phase2,
    variance_estimator = "Ultimate Cluster"
  )

# Get second-phase joint probabilities
  n <- twophase_design$phase2$fpc$sampsize[1,1]
  N <- twophase_design$phase2$fpc$popsize[1,1]

  second_phase_joint_probs <- Matrix::Matrix((n/N)*((n-1)/(N-1)),
                                     nrow = n, ncol = n)
  diag(second_phase_joint_probs) <- rep(n/N, times = n)

# Get quadratic form for entire two-phase variance estimator
  twophase_quad_form <- make_twophase_quad_form(
   sigma_1 = first_phase_sigma,
   sigma_2 = second_phase_sigma,
   phase_2_joint_probs = second_phase_joint_probs
 )
#> Warning: Approximating (sigma_1/phase_2_joint_probs) with the nearest positive semidefinite matrix, since the matrix (1/phase_2_joint_probs) is not positive semidefinite. This is expected to result in a small overestimation of variance. See `help('make_twophase_quad_form', package = 'svrep')` for details.

 # Use for variance estimation

   rep_factors <- make_gen_boot_factors(
     Sigma = twophase_quad_form,
     num_replicates = 500
   )

   library(survey)

   combined_weights <- 1/twophase_design$prob

   twophase_rep_design <- svrepdesign(
     data = twophase_sample |>
       subset(SECOND_PHASE_SELECTION),
     type = 'other',
     repweights = rep_factors,
     weights = combined_weights,
     combined.weights = FALSE,
     scale = attr(rep_factors, 'scale'),
     rscales = attr(rep_factors, 'rscales')
   )

   svymean(x = ~ LIBRARIA, design = twophase_rep_design)
#>            mean     SE
#> LIBRARIA 5.4293 2.0079


## ---------------------- Example 2 ------------------------##
## First phase is a stratified systematic sample            ##
## Second phase is nonresponse, modeled as Poisson sampling ##
##----------------------------------------------------------##

data('library_stsys_sample', package = 'svrep')

# Determine quadratic form for full first-phase sample variance estimator

  full_phase1_quad_form <- make_quad_form_matrix(
    variance_estimator = "SD2",
    cluster_ids = library_stsys_sample[,'FSCSKEY',drop=FALSE],
    strata_ids = library_stsys_sample[,'SAMPLING_STRATUM',drop=FALSE],
    strata_pop_sizes = library_stsys_sample[,'STRATUM_POP_SIZE',drop=FALSE],
    sort_order = library_stsys_sample$SAMPLING_SORT_ORDER
  )

# Identify cases included in phase two sample
# (in this example, respondents)
  phase2_inclusion <- (
    library_stsys_sample$RESPONSE_STATUS == "Survey Respondent"
  )
  phase2_sample <- library_stsys_sample[phase2_inclusion,]

# Estimate response propensities

  response_propensities <- glm(
    data = library_stsys_sample,
    family = quasibinomial('logit'),
    formula = phase2_inclusion ~ 1,
    weights = 1/library_stsys_sample$SAMPLING_PROB
  ) |>
    predict(type = "response",
            newdata = phase2_sample)

# Estimate conditional joint inclusion probabilities for second phase

  phase2_joint_probs <- outer(response_propensities, response_propensities)
  diag(phase2_joint_probs) <- response_propensities

# Determine quadratic form for variance estimator of second phase
# (Horvitz-Thompson estimator for nonresponse modeled as Poisson sampling)

  phase2_quad_form <- make_quad_form_matrix(
    variance_estimator = "Horvitz-Thompson",
    joint_probs = phase2_joint_probs
  )

# Create combined quadratic form for entire design

 twophase_quad_form <- make_twophase_quad_form(
   sigma_1 = full_phase1_quad_form[phase2_inclusion, phase2_inclusion],
   sigma_2 = phase2_quad_form,
   phase_2_joint_probs = phase2_joint_probs
 )
#> Warning: Approximating (sigma_1/phase_2_joint_probs) with the nearest positive semidefinite matrix, since the matrix (1/phase_2_joint_probs) is not positive semidefinite. This is expected to result in a small overestimation of variance. See `help('make_twophase_quad_form', package = 'svrep')` for details.

 combined_weights <- 1/(phase2_sample$SAMPLING_PROB * response_propensities)

# Use for variance estimation

  rep_factors <- make_gen_boot_factors(
    Sigma = twophase_quad_form,
    num_replicates = 500
  )

  twophase_rep_design <- svrepdesign(
    data = phase2_sample,
    type = 'other',
    repweights = rep_factors,
    weights = combined_weights,
    combined.weights = FALSE,
    scale = attr(rep_factors, 'scale'),
    rscales = attr(rep_factors, 'rscales')
  )

  svymean(x = ~ LIBRARIA, design = twophase_rep_design)
#>            mean     SE
#> LIBRARIA 7.3534 1.5718
# }
```
