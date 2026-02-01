# Replication Methods for Two-phase Sampling

This vignette provides an overview of two-phase sampling and how the
‘svrep’ package can be used to estimate sampling variances for
estimators that are commonly-used in two-phase sampling. Readers who are
already familiar with two-phase sampling and its applications are
encouraged to skip to Section 3 of this vignette for a description of
variance estimators and their implementation in the ‘svrep’ package.

## Two-phase Sampling vs. Multistage Sampling

Two-phase sampling (also known as “double sampling”) is a common feature
in surveys. In a two-phase sample, a large first-phase sample is
selected, and a smaller second-phase sample is selected from the
first-phase sample. Multistage cluster sampling is a special case of
two-phase sampling, where a second-phase sample of secondary sampling
units (SSUs) are selected from a first-phase sample of primary sampling
units (PSUs). In the specific case of multistage sampling, the
second-phase sampling of SSUs must sample at least one SSU from within
each PSU and must sample independently across PSUs (in other words, each
PSU is treated as a stratum in second-phase sampling). Two-phase
sampling in general does not have these restrictions: the second-phase
sample design can be arbitrary, and some primary sampling units might
not appear at all in the second-phase sample.

## Applications of Two-Phase Sampling

The flexibility of two-phase sampling can be quite valuable, and for
this reason two-phase samples are commonly-used in practice. We
highlight two common applications of two-phase sampling below:

- Any given survey conducted using an online panel is necessarily a
  two-phase sample, where the panel recruitment represents the first
  phase of sampling and the process of requesting panelists to
  participate in a specific survey represents the second phase of
  sampling. Often, the recruitment sampling is quite complex (e.g.,
  three-stage stratified cluster sampling), but the sampling of
  panelists for a given survey is conducted using simple random sampling
  or stratified simple random sampling from the list of panelists.

- Statistical agencies often reduce the cost of a small survey by
  drawing its sample from respondents to a larger survey that’s already
  being conducted. For example, the U.S. Census Bureau conducts the
  [National Survey of College Graduates
  (NSCG)](https://ncses.nsf.gov/surveys/national-survey-college-graduates/2023)
  by sampling from households that responded to the [American Community
  Survey (ACS)](https://www.census.gov/programs-surveys/acs/).
  Similarly, the [National Study of Caregiving
  (NSOC)](https://nhats.org/researcher/nsoc) is conducted by sampling
  respondents to the [National Health and Aging Trends Study
  (NHATS)](https://nhats.org/researcher/nhats).

The information from the first-phase sample is useful for both design
and analysis of the second-phase sample. From a design standpoint,
information collected in the first-phase sample can be used to stratify
units or assign unequal sampling probabilities for the second-phase
sampling, which can result in more precise estimates relative to using
simple random sampling. From an analysis standpoint, the information
collected in the first-phase sample can also be used to improve
estimators, by using raking, post-stratification, or generalized
regression (GREG) to calibrate the small second-phase sample to the
large first-phase sample.

## Replicate Variance Estimation with the ‘svrep’ Package

In this vignette, we’ll show how to use the generalized bootstrap to
estimate sampling variances for estimates based on two-phase sample
designs. Other types of replication such as the jackknife or balanced
repeated replication (BRR) can theoretically be used, but the ‘svrep’
package only implements two-phase replication methods for the
generalized bootstrap and for Fay’s generalized replication method. In
theory, other replication methods can be used for two-phase samples, but
their applicability is much more limited.

### Overview of the Generalized Bootstrap

The basic idea of the generalized bootstrap is to “mimic” a target
variance estimator for population totals, where the target variance
estimator is appropriate to the particular sampling design and can be
written down as a quadratic form. For example, the generalized bootstrap
can mimic the Horvitz-Thompson estimator or the usual variance estimator
used for simple random sampling. To be more precise, by “mimic”, we mean
that the generalized bootstrap variance estimate for a population total
on average exactly matches the variance estimate produced by the target
variance estimator.

In order to mimic a target variance estimator, we have to specify the
target variance estimator for a population total
$\widehat{Y} = \sum_{i = 1}^{n}\left( y_{i}/\pi_{i} \right)$ as a
quadratic form. That is, we have to specify a variance estimator
$v\left( \widehat{Y} \right)$ as
$v\left( \widehat{Y} \right) = \sum_{i = 1}^{n}\sum_{i = 1}^{n}\sigma_{ij}\left( w_{i}y_{i} \right)\left( w_{j}y_{j} \right)$,
for some set of values $\sigma_{ij},i,j \in \{ 1,\ldots,n\}$. In matrix
notation, we write
$v\left( \widehat{Y} \right) = {\breve{y}}^{\prime}\Sigma\breve{y}$,
where $\Sigma$ is the symmetric, positive semi-definite matrix of
dimension $n \times n$, with element $ij$ equal to $\sigma_{ij}$, and
$\breve{y}$ is a vector whose $i$-th element is $w_{i}y_{i}$.

When using the generalized bootstrap, the difficult part of the variance
estimation process is simply identifying the quadratic form. Once the
quadratic form has been written down, it is easy to create replicate
weights using the generalized bootstrap. Fortunately, the ‘svrep’
package can automatically identify the appropriate quadratic form to use
for variance estimators for many single-phase and two-phase sample
designs. The user simply needs to supply the necessary data, describe
the survey design, and select a target variance estimator to use for
each phase of sampling.

For a broad overview of the generalized survey bootstrap and its use in
the ‘svrep’ package, the reader is encouraged to read the ‘svrep’
package vignette titled “Bootstrap Methods for Surveys”. For a thorough
overview of the generalized survey bootstrap and its theory, Beaumont
and Patak (2012) provide a clear introduction and several useful
suggestions for its implementation in practice. The present vignette
simply describes the application of the generalized bootstrap to
two-phase samples and how it can be implemented with the ‘svrep’
package.

### Creating Example Data

In the example below, we create a two-phase survey design:

- The first phase is a stratified multistage sample, where the first
  stage sample of PSUs was selected using unequal probability sampling
  without replacement (PPSWOR) and the second stage sample was selected
  using simple random sampling without replacement (SRSWOR).

- The second phase sample is a simple random sample without replacement
  from the first phase sample.

This type of design would be fairly typical for a survey conducted on an
online panel, where the panel recruitment uses a complex design but the
sampling of panelists for a given survey uses simple random sampling of
panelists.

The particular dataset we’ll use comes from the Public Libraries Survey
(PLS), an annual survey of public libraries in the U.S, with data from
FY2020.

``` r
data('library_multistage_sample', package = 'svrep')

# Load first-phase sample
  twophase_sample <- library_multistage_sample

# Select second-phase sample
  set.seed(2020)
  
  twophase_sample[['SECOND_PHASE_SELECTION']] <- sampling::srswor(
    n = 100,
    N = nrow(twophase_sample)
  ) |> as.logical()
```

### Describing the Two-phase Survey Design

Next, we use the ‘survey’ package’s function
[`twophase()`](https://rdrr.io/pkg/survey/man/twophase.html) to describe
the sample design at each phase, in terms of stratification, clustering,
probabilities, and population sizes. Note that we use a
[`list()`](https://rdrr.io/r/base/list.html) for most arguments, where
the first element of the list describes the first phase of sampling, and
the second element of the list describes the second phase of sampling.

``` r
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
```

### Creating Generalized Bootstrap Replicates

Once the two-phase design has been described, we can use the
[`as_gen_boot_design()`](https://bschneidr.github.io/svrep/reference/as_gen_boot_design.md)
function to create generalized bootstrap replicate weights. This
requires us to specify the desired number of replicates and the target
variance estimator for each phase of sampling. Note that different
target variance estimators may be used for each phase, since each phase
might have a very different design.

``` r
# Obtain a generalized bootstrap replicates
# based on 
#   - The phase 1 estimator is the usual variance estimator
#     for stratified multistage simple random sampling
#   - The phase 2 estimator is the usual variance estimator
#     for single-stage simple random sampling

twophase_boot_design <- as_gen_boot_design(
  design = twophase_design,
  variance_estimator = list(
    "Phase 1" = "Stratified Multistage SRS",
    "Phase 2" = "Ultimate Cluster"
  ),
  replicates = 1000
)
```

The result is a replicate survey design object which can be used for
estimation with the usual functions from the ‘survey’ and ‘srvyr’
packages.

``` r
twophase_boot_design |> svymean(x = ~ LIBRARIA, na.rm = TRUE)
#>            mean     SE
#> LIBRARIA 7.6044 1.9818
```

When using
[`as_gen_boot_design()`](https://bschneidr.github.io/svrep/reference/as_gen_boot_design.md)
for two-phase designs, it’s useful to know that you will often see a
warning message about needing to approximate the first-phase variance
estimator’s quadratic form.

``` r
twophase_boot_design <- as_gen_boot_design(
  design = twophase_design,
  variance_estimator = list(
    "Phase 1" = "Stratified Multistage SRS",
    "Phase 2" = "Ultimate Cluster"
  )
)
#> Warning in as_gen_boot_design.twophase2(design = twophase_design,
#> variance_estimator = list(`Phase 1` = "Stratified Multistage SRS", : The sample
#> quadratic form matrix for this design and variance estimator is not positive
#> semidefinite. It will be approximated by the nearest positive semidefinite
#> matrix.
```

As you can see from the output above, the function emitted a warning
message. The generalized bootstrap works by mimicking a variance
estimator but requires that variance estimator to be represented as a
positive semidefinite quadratic form. In two-phase designs, however, it
is often the case that the usual variance estimator cannot be
represented exactly as a positive semidefinite quadratic form. In such
cases, Beaumont and Patak (2012) suggest using an approximation of the
actual quadratic form matrix by the most similar positive semidefinite
matrix. This approximation will in general never lead to an
underestimation of variance, and Beaumont and Patak (2012) argue that
this should only produce a small overestimate of variance in practice.
Section 5 of this vignette provides more details of this approximation.

### Create Replicates Using Fay’s Generalized Replication Method

Instead of the generalized bootstrap, we can instead use Fay’s
generalized replication method. The R code looks almost exactly the same
as for the generalized bootstrap.

``` r
twophase_genrep_design <- as_fays_gen_rep_design(
  design = twophase_design,
  variance_estimator = list(
    "Phase 1" = "Stratified Multistage SRS",
    "Phase 2" = "Ultimate Cluster"
  ),
  max_replicates = 500
)
#> Warning in as_fays_gen_rep_design.twophase2(design = twophase_design,
#> variance_estimator = list(`Phase 1` = "Stratified Multistage SRS", : The sample
#> quadratic form matrix for this design and variance estimator is not positive
#> semidefinite. It will be approximated by the nearest positive semidefinite
#> matrix.
```

The key difference from a programming standpoint is that we use the
argument `max_replicates` to specify the maximum number of replicates
that can be created. If the function determines that fewer than
`max_replicates` are needed to obtain a fully-efficient variance
estimator, then the actual number of replicates created will be less
than `max_replicates`.

### Calibrating Second-phase Weights to First-phase Estimates

In two-phase sampling, it can be helpful to calibrate the weights from
the small second-phase sample using estimates produced from the larger,
more reliable first-phase sample. The main reason for doing this is to
produce more precise estimates for variables only measured in the
second-phase sample, and calibration is effective at this if the
calibration variables are associated with the second-phase variables of
interest. But calibration is also nice because it forces second-phase
estimates for calibration variables to match the first-phase estimates,
thus improving the consistency of the two sets of estimates.

Calibrating the weights for the second-phase sample is straightforward
and can be done using usual software and methods. However, care is
needed to ensure that resulting variance estimates appropriately reflect
the fact that we are calibrating to *estimates* rather than to known
population values. This is fairly easy when replication methods are used
for variance estimation, but requires the use of the appropriate
functions from the ‘svrep’ package. Section 4.3.1 of this memo discusses
the theory of replicate variance estimation for two-phase calibration,
based on more detailed treatments of this topic by Fuller (1998) and
Lohr (2022).

Below is a general process for using the ‘svrep’ package to calibrate a
second-phase sample to first-phase estimates while ensuring that
replicate weights are adjusted appropriately for the purpose of variance
estimation. There are two useful functions from the ‘svrep’ package for
this purpose, which we present as “Option 1” and “Option 2” in the
following overview.

#### Preliminaries

**Ensure that the calibration variables do not have any missing values**

First, we need to ensure that the variables we want to use for
calibration do not have any missing values in either the first-phase or
second-phase sample. Some imputation might be necessary.

``` r
# Impute missing values (if necessary)
twophase_sample <- twophase_sample |>
  mutate(
    TOTCIR = ifelse(
      is.na(TOTCIR),
      stats::weighted.mean(TOTCIR, na.rm = TRUE,
                           w = 1/SAMPLING_PROB),
      TOTCIR
    ),
    TOTSTAFF = ifelse(
      is.na(TOTSTAFF),
      stats::weighted.mean(TOTSTAFF, na.rm = TRUE,
                           w = 1/SAMPLING_PROB),
      TOTSTAFF
    )
  )
```

**(If you haven’t already) Create replicate weights for the second-phase
sample**

Before calibration, we need to create replicate weights for the
second-phase sample that appropriately reflect the sampling variance of
the entire two-phase design. We did that already in this document, but
we’ll repeat that code again below.

``` r
# Describe the two-phase survey design
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

# Create replicate weights for the second-phase sample
# (meant to reflect variance of the entire two-phase design)
  twophase_boot_design <- as_gen_boot_design(
    design = twophase_design,
    variance_estimator = list(
      "Phase 1" = "Stratified Multistage SRS",
      "Phase 2" = "Ultimate Cluster"
    ),
    replicates = 1000,
    mse = TRUE
  )
```

#### Option 1: Calibrate to a set of estimates and their variance-covariance matrix

In this approach, we use data from the first-phase sample to produce
estimated totals to use for calibration of the second-phase sample. To
ensure that the calibration of the second-phase sample appropriately
reflects the variance of the first-phase estimated totals, we also need
to estimate the variance of the first-phase totals.

There are many ways to estimate the first-phase variance, but for
convenience we’ll use the generalized bootstrap.

``` r
# Extract a survey design object representing the first phase sample
  first_phase_design <- twophase_design$phase1$full

# Create replicate weights for the first-phase sample
  first_phase_gen_boot <- as_gen_boot_design(
    design = first_phase_design,
    variance_estimator = "Stratified Multistage SRS",
    replicates = 1000
  )
  
# Estimate first-phase totals and their sampling-covariance
  first_phase_estimates <- svytotal(
    x = ~ TOTCIR + TOTSTAFF,
    design = first_phase_gen_boot
  )
  
  first_phase_totals <- coef(first_phase_estimates)
  first_phase_vcov <- vcov(first_phase_estimates)
  
  print(first_phase_totals)
#>       TOTCIR     TOTSTAFF 
#> 1648795905.4     152846.6
  print(first_phase_vcov)
#>                TOTCIR     TOTSTAFF
#> TOTCIR   6.606150e+16 5.853993e+12
#> TOTSTAFF 5.853993e+12 5.747174e+08
#> attr(,"means")
#> [1] 1645727222.2     152190.5
```

After we’ve estimated the first-phase totals, we can use the function
[`calibrate_to_estimate()`](https://bschneidr.github.io/svrep/reference/calibrate_to_estimate.md)
to calibrate the two-phase survey design object to the first-phase
totals. This function is discussed in more detail in the vignette titled
“Sample-based Calibration”, and the underlying method is described in
Fuller (1998).

``` r
calibrated_twophase_design <- calibrate_to_estimate(
  rep_design = twophase_boot_design,
  # Specify the variables in the data to use for calibration
  cal_formula = ~ TOTCIR + TOTSTAFF,
  # Supply the first-phase estimates and their variance
  estimate = first_phase_totals,
  vcov_estimate = first_phase_vcov,
)
#> Selection of replicate columns whose control totals will be perturbed will be done at random.
#> For tips on reproducible selection, see `help('calibrate_to_estimate')`
```

Let’s examine the results from calibration. First, we’ll check that the
calibrated second-phase estimates match the first-phase estimates.

``` r
# Display second-phase estimates for calibration variables
svytotal(
  x = ~ TOTCIR + TOTSTAFF,
  design = calibrated_twophase_design
)
#>               total        SE
#> TOTCIR   1648795905 257024311
#> TOTSTAFF     152847     23973

# Display the original first-phase estimates (which are identical!)
print(first_phase_estimates)
#>               total        SE
#> TOTCIR   1648795905 257024311
#> TOTSTAFF     152847     23973
```

Next, we’ll inspect an estimate for a variable that wasn’t used in
calibration.

``` r
# Inspect calibrated second-phase estimate
svytotal(
  x = ~ LIBRARIA, na.rm = TRUE,
  design = calibrated_twophase_design
)
#>          total    SE
#> LIBRARIA 57355 37647

# Compare to uncalibrated second-phase estimate
svytotal(
  x = ~ LIBRARIA, na.rm = TRUE,
  design = twophase_boot_design
)
#>          total    SE
#> LIBRARIA 54368 12039

# Compare to first-phase estimate
svytotal(
  x = ~ LIBRARIA, na.rm = TRUE,
  design = first_phase_gen_boot
)
#>          total     SE
#> LIBRARIA 55696 9171.3
```

#### Option 2: Calibrate to independently-generated first-phase replicates

If we have the data for the first-phase sample available and there are
replicate weights created for the first-phase sample, then we have an
arguably better method available to handle calibration. We can simply
produce replicate estimates of the first-phase totals using each
first-phase replicate, and then we can calibrate each second-phase
replicate to one of the first-phase replicate totals.

To do this, we first create replicate weights for the first-phase design
using the generalized bootstrap (or any other replication method).

``` r
# Extract a survey design object representing the first phase sample
  first_phase_design <- twophase_design$phase1$full

# Create replicate weights for the first-phase sample
  first_phase_gen_boot <- as_gen_boot_design(
    design = first_phase_design,
    variance_estimator = "Stratified Multistage SRS",
    replicates = 1000
  )
```

After we’ve created the first-phase replicates, we can use the function
[`calibrate_to_sample()`](https://bschneidr.github.io/svrep/reference/calibrate_to_sample.md)
to calibrate the two-phase survey design object to replicate estimates
created using the first-phase replicate design. This function is
discussed in more detail in the vignette titled “Sample-based
Calibration”. See Section 4.3.1 of this vignette for the underlying
theory, which is based on Fuller (1998) and Opsomer and Erciulescu
(2021).[¹](#fn1)

``` r
calibrated_twophase_design <- calibrate_to_sample(
  primary_rep_design = twophase_boot_design,
  # Supply the first-phase replicate design
  control_rep_design = first_phase_gen_boot,
  # Specify the variables in the data to use for calibration
  cal_formula = ~ TOTCIR + TOTSTAFF
)
#> Matching between primary and control replicates will be done at random.
#> For tips on reproducible matching, see `help('calibrate_to_sample')`
```

Let’s examine the results from calibration. First, we’ll check that the
calibrated second-phase estimates match the first-phase estimates.

``` r
# Display second-phase estimates for calibration variables
calibrated_ests <- svytotal(
  x = ~ TOTCIR + TOTSTAFF,
  design = calibrated_twophase_design
)

print(calibrated_ests)
#>               total        SE
#> TOTCIR   1648795905 242527993
#> TOTSTAFF     152847     22856

# Display the original first-phase estimates (which are identical!)
first_phase_ests <- svytotal(
  x = ~ TOTCIR + TOTSTAFF,
  design = first_phase_gen_boot
)

print(first_phase_ests)
#>               total        SE
#> TOTCIR   1648795905 242515035
#> TOTSTAFF     152847     22854
```

As expected, the variance estimate for the calibrated second-phase
estimate is the same as the variance estimate for the first-phase
estimate, allowing a small tolerance for numeric differences.

``` r
ratio_of_variances <- vcov(calibrated_ests)/vcov(first_phase_ests)
ratio_of_variances
#>             TOTCIR  TOTSTAFF
#> TOTCIR   1.0001069 0.9998445
#> TOTSTAFF 0.9998445 1.0002008
#> attr(,"means")
#>       TOTCIR     TOTSTAFF 
#> 1648795905.4     152846.6
```

Next, we’ll inspect an estimate for a variable that wasn’t used in
calibration.

``` r
# Inspect calibrated second-phase estimate
svytotal(
  x = ~ LIBRARIA, na.rm = TRUE,
  design = calibrated_twophase_design
)
#>          total    SE
#> LIBRARIA 57355 42417

# Compare to uncalibrated second-phase estimate
svytotal(
  x = ~ LIBRARIA, na.rm = TRUE,
  design = twophase_boot_design
)
#>          total    SE
#> LIBRARIA 54368 12039

# Compare to first-phase estimate
svytotal(
  x = ~ LIBRARIA, na.rm = TRUE,
  design = first_phase_gen_boot
)
#>          total     SE
#> LIBRARIA 55696 8876.4
```

### Ratio Estimation

A special case of calibration that is commonly used in two-phase samples
is ratio estimation. Whether you use the function
[`calibrate_to_sample()`](https://bschneidr.github.io/svrep/reference/calibrate_to_sample.md)
or
[`calibrate_to_estimate()`](https://bschneidr.github.io/svrep/reference/calibrate_to_estimate.md),
the syntax is very similar.

``` r
ratio_calib_design <- calibrate_to_sample(
  primary_rep_design = twophase_boot_design,
  # Supply the first-phase replicate design
  control_rep_design = first_phase_gen_boot,
  # Specify the GREG formula.
  # For ratio estimation, we add `-1` to the formula 
  # (i.e., we remove the intercept from the working model)
  # and specify only a single variable
  cal_formula = ~ -1 + TOTSTAFF,
  variance = 1
)
#> Matching between primary and control replicates will be done at random.
#> For tips on reproducible matching, see `help('calibrate_to_sample')`
```

Note that for ratio estimation, the calibration formula includes `-1` to
ensure that ratio estimation is used instead of regression estimation.
This is similar to how, when fitting a regression model in R, we use
`lm(y ~ -1 + x)` to fit a linear model without an intercept. Specifying
the parameter `variance = 1` indicates that the working model used in
calibration is homoscedastic, and so the same adjustment factor will be
used for every case’s weights. This can be seen when we compare the
adjusted weights to the unadjusted weights.

``` r
ratio_adjusted_weights <- weights(ratio_calib_design, type = "sampling")
unadjusted_weights <- weights(twophase_boot_design, type = "sampling")

adjustment_factors <- ratio_adjusted_weights/unadjusted_weights
head(adjustment_factors)
#>        1        3        5        7       10       13 
#> 1.090189 1.090189 1.090189 1.090189 1.090189 1.090189
```

Note that the adjustment factor for the weights is simply the ratio of
the first-phase estimated total to the second-phase estimated total.

``` r
phase1_total <- svytotal(
  x = ~ TOTSTAFF,
  first_phase_design
) |> coef()
phase2_total <- svytotal(
  x = ~ TOTSTAFF,
  twophase_boot_design
) |> coef()

phase1_total/phase2_total
#> TOTSTAFF 
#> 1.090189
```

## Design-based Estimators for Two-phase Sampling

In the section below, we first describe the double expansion estimator
(DEE) which produces unbiased estimates for two-phase samples, using
only information about the sampling design in both phases. Next, we
describe calibration estimators which adjust the weights from the
double-expansion estimator so that sampling variances can be reduced
using information from the first-phase sample. We’ll examine both the
theoretical sampling variance of each estimator as well as approaches
for estimating variance using replication methods.

The interested reader is encouraged to consult chapter 9.3 of Särndal,
Swensson, and Wretman (1992) or chapter 12 of Lohr (2022) for a more
detailed discussion of two-phase sampling.

### Notation

We use the following notation to denote each sample and its size.

#### Notation for Samples and Sample Size

\$\$ \begin{aligned} s_a &: \text{The set of units in the first-phase
sample} \\ s_b &: \text{The set of units in the second-phase sample} \\
& \space \space \space \text{Note that }s_b \text{ is a subset of } s_a
\\ n_a &: \text{The number of units in }s_1 \\ n_b &: \text{The number
of units in }s_2 \\ \end{aligned} \$\$

#### Notation for Probabilities and Weights

We use the following notation to denote the **inclusion probability** of
each unit, for each phase:

$$\begin{aligned}
\pi_{i}^{(a)} & {:{\text{The probability unit}\mspace{6mu}}i{\mspace{6mu}\text{is included in}\mspace{6mu}}s_{a}} \\
\pi_{i}^{(b|s_{a})} & {:{\text{The conditional probability unit}\mspace{6mu}}i{\mspace{6mu}\text{is included in}\mspace{6mu}}s_{b},} \\
 & {{\mspace{6mu}\text{given the realized first-phase sample}\mspace{6mu}}s_{a}} \\
\pi_{i} & {:{\text{The}\mspace{6mu}}\textbf{𝐮𝐧𝐜𝐨𝐧𝐝𝐢𝐭𝐢𝐨𝐧𝐚𝐥}{\mspace{6mu}\text{probability unit}\mspace{6mu}}i{\mspace{6mu}\text{is included in}\mspace{6mu}}s_{b}} \\
 & 
\end{aligned}$$

In practice, the probability $\pi_{i}$ is prohibitively difficult to
calculate, because it requires us to figure out $\pi_{i}^{(b|s_{a})}$
for every possible first-phase sample $s_{a}$, not just the particular
$s_{a}$ that we actually selected. So instead, we define the useful
quantity $\pi^{*}$, which depends only on the particular first-phase
sample $s_{a}$ that we actually selected.

$$\pi_{i}^{*}:=\pi_{i}^{(b|s_{a})} \times \pi_{i}^{(a)}$$

For variance estimation, it’s also necessary to consider the **joint
inclusion probability** (sometimes referred to as “second order
probability”), which is simply the probability that a pair of units $i$
and $j$ are both included in a sample.

$$\begin{aligned}
\pi_{ij}^{(a)} & {:{\text{The probability units}\mspace{6mu}}i{\mspace{6mu}\text{and}\mspace{6mu}}j{\mspace{6mu}\text{are both included in}\mspace{6mu}}s_{a}} \\
\pi_{ij}^{(b|s_{a})} & {:{\text{The conditional probability units}\mspace{6mu}}i{\mspace{6mu}\text{and}\mspace{6mu}}j{\mspace{6mu}\text{are both included in}\mspace{6mu}}s_{b},} \\
 & {{\mspace{6mu}\text{given the realized first-phase sample}\mspace{6mu}}s_{a}} \\
 & 
\end{aligned}$$

We also define the quantity $\pi_{ij}^{*}$ similar to $\pi_{i}^{*}$.

$$\pi_{ij}^{*}:=\pi_{ij}^{(b|s_{a})} \times \pi_{ij}^{(a)}$$

The probabilities and the $\pi_{i}^{*}$ values are used to define
sampling weights for the survey.

$$\begin{aligned}
w_{i}^{(a)} & {:=1/\pi_{i}^{(a)}} \\
w_{i}^{(b|s_{a})} & {:=1/\pi_{i}^{(b|s_{a})}} \\
w_{i}^{*} & {:=1/\pi_{i}^{*} = w_{i}^{(b|s_{a})} \times w_{i}^{(a)}}
\end{aligned}$$

### The Double Expansion Estimator

Suppose we wish to estimate a population total $Y$, using observed
values $y_{i}$ in our second-phase sample, $s_{b}$. Särndal, Swensson,
and Wretman (1992) show that we can produce an unbiased estimate of $Y$
using the second-phase sample $s_{b}$, as follows:

$$\begin{aligned}
{\widehat{Y}}^{(b)} & {= \sum\limits_{i = 1}^{n_{(b)}}w_{i}^{*} \times y_{i}} \\
 & {= \sum\limits_{i = 1}^{n_{(b)}}w_{i}^{(b|s_{a})} \times w_{i}^{(a)} \times y_{i}}
\end{aligned}$$

This estimator has been dubbed the “double expansion estimator”, using
the sampling jargon that refers to weighting a sample value $y_{i}$ as
“expanding” $y_{i}$ from the sample to the population. The name “double
expansion” is used because the weight $w_{i}^{*}$ can be thought of as
first using the weight $w_{i}^{(b|s_{a})}$ to “expand” the quantity
$y_{i}$ and then using the weight $w_{i}^{(a)}$ to expand the quantity
$w_{i}^{(b|s_{a})} \times y_{i}$.

#### Variance of the Double Expansion Estimator

The sampling variance of the double expansion estimator is the sum of
two different components.

$$\begin{aligned}
{V\left( {\widehat{Y}}^{(b)} \right)} & {= V\left( {\widehat{Y}}^{(a)} \right) + E\left( V\left\lbrack {\widehat{Y}}^{(b)} \mid s_{a} \right\rbrack \right)} \\
 & \\
{\text{where:}\mspace{6mu}} & {{\widehat{Y}}^{(a)} = \sum\limits_{i = 1}^{n_{(a)}}w_{i}^{(a)} \times y_{i}} \\
{\text{and}\mspace{6mu}} & {V\left\lbrack {\widehat{Y}}^{(b)} \mid s_{a} \right\rbrack{\mspace{6mu}\text{is the variance of}\mspace{6mu}}{\widehat{Y}}^{(b)}} \\
 & {{\mspace{6mu}\text{across all samples}\mspace{6mu}}s_{b}} \\
 & {{\mspace{6mu}\text{drawn from a given}\mspace{6mu}}s_{a}}
\end{aligned}$$

The first component is the variance of the estimate
${\widehat{Y}}^{(a)}$ that we would obtain if we used the entire
first-phase sample $s_{a}$ for our estimate, rather than using the
subset $s_{b}$.

The second component is the additional variance caused by using the
subset $s_{b}$ instead of $s_{a}$. It is equal to the expected value
(across all samples $s_{a}$) of the conditional variance of
${\widehat{Y}}^{(b)}$ across all samples $s_{b}$ (conditioning on a
given first-phase sample $s_{a}$).

##### Estimating the Variance of the Double Expansion Estimator

Both variance components can be estimated using only the values $y_{i}$
observed in $s_{b}$. For the second component, we simply estimate
$V\left\lbrack {\widehat{Y}}^{(b)} \mid s_{a} \right\rbrack$, which is
an unbiased estimate for its expectation,
$E\left( V\left\lbrack {\widehat{Y}}^{(b)} \mid s_{a} \right\rbrack \right)$.

Thus, our variance estimate for the double expansion estimator takes the
following form:

$$\widehat{V}\left( {\widehat{Y}}^{(b)} \right) = \widehat{V}\left\lbrack {\widehat{Y}}^{(a)} \right\rbrack + \widehat{V}\left\lbrack {\widehat{Y}}^{(b)} \mid s_{a} \right\rbrack$$

###### Estimating the second-phase variance component

For estimating
$\widehat{V}\left\lbrack {\widehat{Y}}^{(b)} \mid s_{a} \right\rbrack$,
we simply choose a variance estimator for the second-phase design,
taking the first-phase sample as a given. We assume that this variance
estimator can be written as a quadratic form.

$$\begin{aligned}
{\widehat{V}\left\lbrack {\widehat{Y}}^{(b)} \mid s_{a} \right\rbrack} & {= \sum\limits_{i = 1}^{n_{b}}\sum\limits_{i = 1}^{n_{b}}\sigma_{ij}^{(b)}\left( w_{i}^{*}y_{i} \right)\left( w_{j}^{*}y_{j} \right)} \\
 & 
\end{aligned}$$

For the Horvitz-Thompson estimator, for instance, we would use
$\sigma_{ij}^{(b)} = \left( 1 - \frac{\pi_{i}^{b|s_{a}}\pi_{j}^{b|s_{a}}}{\pi_{ij}^{b|s_{a}}} \right)$.

This quadratic form can also be written in matrix notation:

$$\begin{aligned}
{\widehat{V}\left\lbrack {\widehat{Y}}^{(b)} \mid s_{a} \right\rbrack} & {= \left( W^{*}y \right)^{\prime}\Sigma_{b}\left( W^{*}y \right)} \\
{\text{where}\mspace{6mu}} & {\Sigma_{b}{\mspace{6mu}\text{is an}\mspace{6mu}}n_{b} \times n_{b}{\mspace{6mu}\text{symmetric matrix}}} \\
 & {{\mspace{6mu}\text{with entry}\mspace{6mu}}ij{\mspace{6mu}\text{equal to}\mspace{6mu}}\sigma_{ij}^{(b)}} \\
{\text{and}\mspace{6mu}} & {W^{*}{\mspace{6mu}\text{is the}\mspace{6mu}}n_{b} \times n_{b}{\mspace{6mu}\text{diagonal matrix}}} \\
 & {{\mspace{6mu}\text{with entry}\mspace{6mu}}ii{\mspace{6mu}\text{equal to}\mspace{6mu}}w_{i}^{*}} \\
 & {y{\mspace{6mu}\text{is the}\mspace{6mu}}n_{b} \times 1{\mspace{6mu}\text{vector of values}}} \\
 & \text{for the variable of interest}
\end{aligned}$$

###### Estimating the first-phase variance component

Estimating the first variance component,
$V\left( {\widehat{Y}}^{(a)} \right)$, is only slightly trickier. First,
we need to choose a variance estimator appropriate to the first-phase
design, which we would use if we had $y_{i}$ observed for the entire
sample $s_{a}$. We’ll denote that variance estimator
$\widetilde{V}\left\lbrack {\widehat{Y}}^{(a)} \right\rbrack$.

$$\begin{aligned}
{\widetilde{V}\left\lbrack {\widehat{Y}}^{(a)} \right\rbrack} & {= \sum\limits_{i = 1}^{n_{a}}\sum\limits_{i = 1}^{n_{a}}\sigma_{ij}^{(a)}\left( w_{i}^{(a)}y_{i} \right)\left( w_{i}^{(a)}y_{j} \right)} \\
 & 
\end{aligned}$$

In matrix notation, we can write:

$$\begin{aligned}
{\widetilde{V}\left\lbrack {\widehat{Y}}^{(a)} \right\rbrack} & {= \left( W^{(a)}y \right)^{\prime}\left( \Sigma_{a} \right)\left( W^{(a)}y \right)} \\
{\text{where}\mspace{6mu}} & {\Sigma_{a}{\mspace{6mu}\text{is an}\mspace{6mu}}n_{a} \times n_{a}{\mspace{6mu}\text{symmetric matrix}}} \\
 & {{\mspace{6mu}\text{with entry}\mspace{6mu}}ij{\mspace{6mu}\text{equal to}\mspace{6mu}}\sigma_{ij}} \\
{\text{and}\mspace{6mu}} & {W^{(a)}{\mspace{6mu}\text{is the}\mspace{6mu}}n_{a} \times n_{a}{\mspace{6mu}\text{diagonal matrix}}} \\
 & {{\mspace{6mu}\text{with entry}\mspace{6mu}}ii{\mspace{6mu}\text{equal to}\mspace{6mu}}w_{i}^{(a)}}
\end{aligned}$$

However, since we’re working with the subsample $s_{b}$ instead of
$s_{a}$, we need to estimate
$\widetilde{V}\left\lbrack {\widehat{Y}}^{(a)} \right\rbrack$ using only
the data from $s_{b}$. We can use the second-phase joint inclusion
probabilities $\pi_{ij}^{(b \mid s_{a})}$ to produce an unbiased
estimate of
$\widetilde{V}\left\lbrack {\widehat{Y}}^{(a)} \right\rbrack$ using only
the data from $s_{b}$.

$$\begin{aligned}
{\widehat{V}\left\lbrack {\widehat{Y}}^{(a)} \right\rbrack} & {= \sum\limits_{i = 1}^{n_{b}}\sum\limits_{i = 1}^{n_{b}}\frac{1}{\pi_{ij}^{(b \mid s_{a})}}\sigma_{ij}^{(a)}\left( w_{i}^{(a)}y_{i} \right)\left( w_{i}^{(a)}y_{j} \right)} \\
 & 
\end{aligned}$$

We can also write this in matrix notation:

$$\begin{aligned}
{\widehat{V}\left\lbrack {\widehat{Y}}^{(a)} \right\rbrack} & {= \left( W^{(a)}y \right)^{\prime}\left( \Sigma_{a^{\prime}} \circ D_{b} \right)\left( W^{(a)}y \right)} \\
{\text{where}\mspace{6mu}} & {\Sigma_{a^{\prime}}{\mspace{6mu}\text{is an}\mspace{6mu}}n_{b} \times n_{b}{\mspace{6mu}\text{symmetric matrix}}} \\
 & {{\mspace{6mu}\text{with entry}\mspace{6mu}}ij{\mspace{6mu}\text{equal to}\mspace{6mu}}\sigma_{ij}} \\
{\text{and}\mspace{6mu}} & {W^{(a)}{\mspace{6mu}\text{is the}\mspace{6mu}}n_{b} \times n_{b}{\mspace{6mu}\text{diagonal matrix}}} \\
 & {{\mspace{6mu}\text{with entry}\mspace{6mu}}ii{\mspace{6mu}\text{equal to}\mspace{6mu}}w_{i}^{(a)}} \\
{\mspace{6mu}\text{and}\mspace{6mu}} & {D_{b}{\mspace{6mu}\text{is an}\mspace{6mu}}n_{b} \times n_{b}{\mspace{6mu}\text{symmetric matrix}}} \\
 & {{\mspace{6mu}\text{with entry}\mspace{6mu}}ij{\mspace{6mu}\text{equal to}\mspace{6mu}}\frac{1}{\pi_{ij}^{(b \mid s_{a})}}} \\
 & 
\end{aligned}$$

As a sidenote, that matrix $D_{b}$ is very likely the source of any
warning messages you’ll see about a two-phase variance estimator not
being positive semidefinite. [²](#fn2)

###### Combining the two estimated variance components

Putting the two estimated variance components together, we thus obtain
the following unbiased variance estimator for the double expansion
estimator.

$$\begin{aligned}
{\widehat{V}\left( {\widehat{Y}}^{(b)} \right)} & {= \widehat{V}\left( {\widehat{Y}}^{(a)} \right) + \widehat{V}\left\lbrack {\widehat{Y}}^{(b)} \mid s_{a} \right\rbrack} \\
 & {= \sum\limits_{i = 1}^{n_{b}}\sum\limits_{i = 1}^{n_{b}}\frac{1}{\pi_{ij}^{(b \mid s_{a})}}\sigma_{ij}^{(a)}\left( w_{i}^{(a)}y_{i} \right)\left( w_{i}^{(a)}y_{j} \right)} \\
 & {+ \sum\limits_{i = 1}^{n_{b}}\sum\limits_{i = 1}^{n_{b}}\sigma_{ij}^{(b)}\left( w_{i}^{*}y_{i} \right)\left( w_{j}^{*}y_{j} \right)} \\
 & 
\end{aligned}$$

In matrix notation, we can write this as follows:

$$\begin{aligned}
{\widehat{V}\left( {\widehat{Y}}^{(b)} \right)} & {= \widehat{V}\left( {\widehat{Y}}^{(a)} \right) + \widehat{V}\left\lbrack {\widehat{Y}}^{(b)} \mid s_{a} \right\rbrack} \\
 & {= \left( W^{(a)}y \right)^{\prime}\left( \Sigma_{a^{\prime}} \circ D_{b} \right)\left( W^{(a)}y \right)} \\
 & {+ \left( W^{*}y \right)^{\prime}\Sigma_{b}\left( W^{*}y \right)} \\
 & 
\end{aligned}$$

Because quadratic forms are additive and because
$W^{*} = W^{(a)}W^{(b \mid s_{a})}$, we can more compactly write the
estimator as follows:

$$\begin{aligned}
{\widehat{V}\left( {\widehat{Y}}^{(b)} \right)} & {= \left( W^{*}y \right)^{\prime}\Sigma_{ab}\left( W^{*}y \right)} \\
{\text{where}\mspace{6mu}} & \\
\Sigma_{ab} & {= {W^{(b)}}^{- 1}\left( \Sigma_{a^{\prime}} \circ D_{b} \right){W^{(b)}}^{- 1} + \Sigma_{b}} \\
{\text{where}\mspace{6mu}} & {W^{(b)}{\mspace{6mu}\text{is the}\mspace{6mu}}n_{b} \times n_{b}{\mspace{6mu}\text{diagonal matrix}}} \\
 & {{\mspace{6mu}\text{with entry}\mspace{6mu}}ii{\mspace{6mu}\text{equal to}\mspace{6mu}}w_{i}^{(b \mid s_{a})}}
\end{aligned}$$

In the ‘svrep’ package, $\Sigma_{ab}$ can be constructed with the inputs
$\Sigma_{a^{\prime}}$, $\Sigma_{b}$, and $\left( 1/D_{b} \right)$, using
the function
[`make_twophase_quad_form()`](https://bschneidr.github.io/svrep/reference/make_twophase_quad_form.md).

Click to show/hide example of using `make_twophase_quad_form()`

``` r
set.seed(2022)
y <- rnorm(n = 100)

# Select first phase sample, SRS without replacement
  phase_1_sample_indicators <- sampling::srswor(n = 50, N = 100) |>
    as.logical()
  
  phase_1_sample <- y[phase_1_sample_indicators]
  
# Make variance estimator for first-phase variance component
  Sigma_a <- make_quad_form_matrix(
    variance_estimator = "Ultimate Cluster",
    cluster_ids = as.matrix(1:50),
    strata_ids = rep(1, times = 50) |> as.matrix(),
    strata_pop_sizes = rep(100, times = 50) |> as.matrix()
  )

# Select second stage sample, SRS without replacment
  phase_2_sample_indicators <- sampling::srswor(n = 5, N = 50) |>
    as.logical()
  
  phase_2_sample <- phase_1_sample[phase_2_sample_indicators]
  
# Estimate two-phase variance
  Sigma_a_prime <- Sigma_a[phase_2_sample_indicators,
                           phase_2_sample_indicators]
  
  phase_2_joint_probs <- outer(rep(5/50, times = 5),
                               rep(4/49, times = 5))
  diag(phase_2_joint_probs) <- rep(5/50, times = 5)

  Sigma_b <- make_quad_form_matrix(
    variance_estimator = "Ultimate Cluster",
    cluster_ids = as.matrix(1:5),
    strata_ids = rep(1, times = 5) |> as.matrix(),
    strata_pop_sizes = rep(50, times = 5) |> as.matrix()
  )
  
  sigma_ab <- make_twophase_quad_form(
    sigma_1 = Sigma_a_prime,
    sigma_2 = Sigma_b,
    phase_2_joint_probs = phase_2_joint_probs
  )
  
  wts <- rep(
    (50/100)^(-1) * (5/50)^(-1),
    times = 5
  )
  W_star <- diag(wts)
  
  W_star_y <- W_star %*% phase_2_sample
  t(W_star_y) %*% sigma_ab %*% (W_star_y)
#> 1 x 1 Matrix of class "dgeMatrix"
#>          [,1]
#> [1,] 2182.221
  
# Since both phases are SRS without replacement,
# variance estimate for a total should be similar to the following
  5 * var(W_star_y)
#>          [,1]
#> [1,] 2297.075
```

  

The matrix notation is useful for understanding replication methods of
variance estimation for two-phase samples. Any unbiased replication
variance estimator for two-phase samples should generate each set of
adjustment factors so that the sets of replicate weights have
expectation $\mathbf{1}_{n_{b}}$ and variance-covariance matrix
$\mathbf{\Sigma}_{ab}$.

The generalized bootstrap does this by generating draws from a
multivariate normal distribution with those parameters. For some
specific combinations of simple first-phase and second-phase designs,
there are jackknife and BRR methods which have been developed to
accomplish this same goal (see Lohr (2022) for some examples). The
generalized bootstrap however is much easier to use for the complex
designs actually encountered in most settings and also enjoys other
advantages [³](#fn3).

### Calibration Estimators

This section describes calibration estimators (such as raking,
post-stratification, or ratio estimators) commonly used for two-phase
designs. For a more detailed treatment of such estimators, see Chapter
11 of Lohr (2022) or Chapter 6 of Särndal, Swensson, and Wretman (1992).

In two-phase sampling, it can be helpful to calibrate the weights from
the small second-phase sample $s_{b}$ so that estimates for variables
$x_{1},\ldots,x_{p}$ measured in both phases match the estimates
produced using the larger, more reliable sample $s_{a}$. For a variable
$y$ measured only in the second-phase sample, this can lead to more
precise estimates if the calibration variables $x_{1},\ldots,x_{p}$ are
associated with $y$.

If generalized regression (GREG) is used, the two-phase GREG estimator
can be written as follows:

$${\widehat{Y}}_{\text{GREG}}^{(b)} = {\widehat{Y}}^{(a)} + \left( {\widehat{\mathbf{X}}}^{(a)} - {\widehat{\mathbf{X}}}^{(b)} \right){\widehat{\mathbf{B}}}^{(b)}$$

where ${\widehat{\mathbf{X}}}^{(a)}$ is the $p$-length vector of
estimated population totals for variables $x_{1},\ldots,x_{p}$ estimates
using the first-phase data, ${\widehat{\mathbf{X}}}^{(b)}$ is the vector
of estimated population totals using the second-phase data, and
${\widehat{\mathbf{B}}}^{(b)}$ is estimated using the following:

$${\widehat{\mathbf{B}}}^{(b)} = \left( \sum\limits_{i = 1}^{n_{(b)}}w_{i}^{*}\frac{1}{\sigma_{i}^{2}}\mathbf{x}_{i}\mathbf{x}_{i}^{T} \right)^{- 1}\sum\limits_{i = 1}^{n_{(b)}}w_{i}^{*}\frac{1}{\sigma_{i}^{2}}\mathbf{x}_{i}y_{i}$$

where the constants $\sigma_{i}$ are chosen based on the specific type
of calibration desired.[⁴](#fn4)

The GREG estimator can also be expressed as a weighted estimator based
on modified weights ${\widetilde{w}}_{i}^{*}:=g_{i}w_{i}^{*}$ where the
modification factor $g$ is suitably chosen for the specific method of
calibration used (post-stratification, raking, etc.)

$$\begin{aligned}
{\widehat{Y}}_{\text{GREG}}^{(b)} & {= \sum\limits_{i = 1}^{n_{(b)}}{\widetilde{w}}_{i}^{*}y_{i} = \sum\limits_{i = 1}^{n_{(b)}}\left( g_{i}w_{i}^{*} \right)y_{i}}
\end{aligned}$$

The modification factors $g_{i}$ (commonly referred to as “g-weights”)
can be expressed as:

$$g_{i} = 1 + \left( {\widehat{\mathbf{X}}}^{(a)} - {\widehat{\mathbf{X}}}^{(b)} \right)^{\prime}\left( \sum\limits_{i = 1}^{n_{(b)}}w_{i}^{*}\frac{1}{\sigma_{i}^{2}}\mathbf{x}_{i}\mathbf{x}_{i}^{T} \right)^{- 1}\sum\limits_{i = 1}^{n_{(b)}}w_{i}^{*}\frac{1}{\sigma_{i}^{2}}\mathbf{x}_{i}$$

The calibrated second-phase weights
${\widetilde{w}}_{i}^{*} = g_{i}w_{i}^{*}$ from the GREG estimator
ensure that the second-phase estimates for the variables
$x_{1},\ldots,x_{p}$ match the first-phase estimates.

$$\sum\limits_{i = 1}^{n_{(b)}}{\widetilde{w}}_{i}^{*}x_{i} = \sum\limits_{i = 1}^{n_{(a)}}w^{(a)}x_{i}$$

#### Variance of the Calibration Estimator

If we assume that the second-phase calibration estimator
${\widehat{Y}}_{GREG}^{(b)}$ is unbiased for the first-phase estimate
${\widehat{Y}}^{(a)}$ (which should at least approximately the case),
then we can decompose the calibration estimator’s variance into a
first-phase component and a second-phase component as follows:

$$\begin{aligned}
{V\left( {\widehat{Y}}_{GREG}^{(b)} \right)} & {= V\left\lbrack E\left( {\widehat{Y}}_{GREG}^{(b)} \mid \mathbf{Z} \right) \right\rbrack + E\left\lbrack V\left( {\widehat{Y}}_{GREG}^{(b)} \mid \mathbf{Z} \right) \right\rbrack} \\
 & {= V\left\lbrack {\widehat{Y}}^{(a)} \right\rbrack + E\left\lbrack V\left( {\widehat{Y}}_{GREG}^{(b)} \mid \mathbf{Z} \right) \right\rbrack}
\end{aligned}$$

where the first term is the first-phase variance component and the
second term is the second-phase variance component.

Using only the second-phase sample, the variance of the calibration
estimator can thus be estimated unbiasedly by the following estimator:

$$V\left( {\widehat{Y}}_{GREG}^{(b)} \right) = \widehat{V}\left\lbrack {\widehat{Y}}^{(a)} \right\rbrack + \widehat{V}\left\lbrack {\widehat{E}}^{(b)} \mid s_{a} \right\rbrack$$

where ${\widehat{E}}^{(b)} = \sum_{i = 1}^{n_{(b)}}w^{*}e_{i}$ and
$e_{i} = y_{i} - \mathbf{x}_{i}^{\prime}{\widehat{\mathbf{B}}}^{(b)}$ is
the “residual” of the GREG model.

This is the same as the variance estimator we saw earlier for the
uncalibrated estimator, ${\widehat{Y}}^{(b)}$, except that the
second-phase component for the GREG estimator uses ${\widehat{E}}^{(b)}$
in place of ${\widehat{Y}}^{(b)}$

$$\widehat{V}\left( {\widehat{Y}}^{(b)} \right) = \widehat{V}\left\lbrack {\widehat{Y}}^{(a)} \right\rbrack + \widehat{V}\left\lbrack {\widehat{Y}}^{(b)} \mid s_{a} \right\rbrack$$

This decomposition is useful for understanding the theoretical variance
of the calibration estimator and how it can be estimated in general.

#### Replication Variance Estimation

For variance estimation using replication methods, another (approximate)
decomposition proves to be more useful. Fuller (1998) decomposes the
two-phase calibration estimator’s variance as follows.

$$V\left( {\widehat{Y}}_{GREG}^{(b)} \right) \approx E\left\lbrack V\left( {\widetilde{E}}^{(b)} \mid s_{a} \right) \right\rbrack + \mathbf{B}^{\prime}\mathbf{V}\left( {\widehat{\mathbf{X}}}^{(a)} \right)\mathbf{B}$$

where $\mathbf{B}$ is the finite-population version of
${\widehat{\mathbf{B}}}^{(b)}$ that we could calculate if we had data
from the entire population rather than just the second-phase sample
$s_{b}$, and
${\widetilde{E}}^{(b)} = \sum_{i = 1}^{n_{(b)}}w_{i}^{*}\left( y_{i} - \mathbf{x}_{i}^{\prime}\mathbf{B} \right)$
is the weighted sum of second-phase residuals based on using
$\mathbf{B}$.

This decomposition of the variance suggests the following estimator:

$$\widehat{V}\left( {\widehat{Y}}_{GREG}^{(b)} \right):=\widehat{V}\left( {\widehat{E}}^{(b)} \mid s_{a} \right) + \left( {\widehat{\mathbf{B}}}^{(b)} \right)^{\prime}\widehat{\mathbf{V}}\left( {\widehat{\mathbf{X}}}^{(a)} \right)\left( {\widehat{\mathbf{B}}}^{(b)} \right)$$

The first component is estimated using only the second-phase data and a
conditional variance estimator for the second-phase design (taking the
selected first-phase sample as given). The second component depends on
the first-phase estimates ${\widehat{\mathbf{X}}}^{(a)}$ as well as the
first-phase variance estimate
$\widehat{V}\left( {\widehat{\mathbf{X}}}^{(a)} \right)$ and the values
$\mathbf{B}^{(b)}$ used in the calibration.

Fuller (1998) proposed a replication-based version of this estimator. To
describe this estimator, first we suppose that we have developed
two-phase replicate weights appropriate for the double-expansion
estimator.

$$\begin{aligned}
{\widehat{V}\left( {\widehat{Y}}^{(b)} \right)} & {= K_{(b)}\sum\limits_{r = 1}^{R_{(b)}}\left( {\widehat{Y}}_{(r)}^{(b)} - {\widehat{Y}}^{(b)} \right)^{2}} \\
{\text{where}\mspace{6mu}} & {{\widehat{Y}}_{(r)}^{(b)} = \sum\limits_{i = 1}^{n_{(b)}}w_{r,i}y_{i}} \\
 & {{\text{is the}\mspace{6mu}}r\text{-th}{\mspace{6mu}\text{replicate estimate}}} \\
 & {\text{for the second-phase sample}\mspace{6mu}} \\
{\text{and}\mspace{6mu}} & {K_{(b)}{\mspace{6mu}\text{is a constant specific}}} \\
 & \text{to the replication method}
\end{aligned}$$

Now suppose that we have a $k$-length vector of estimated first-phase
totals, ${\widehat{\mathbf{X}}}^{(a)}$, which will be used in
calibration of the second phase weights. And we suppose that these
estimated totals also have an estimated variance-covariance matrix,
denoted
$\widehat{\mathbf{V}}\left( {\widehat{\mathbf{X}}}^{(a)} \right)$, which
is a $k \times k$ matrix.

Then we can decompose the variance-covariance matrix as follows:

$$\widehat{\mathbf{V}}\left( {\widehat{\mathbf{X}}}^{(a)} \right) = K_{(b)}\sum\limits_{i = 1}^{R_{(b)}}{\mathbf{δ}}_{i}^{\prime}{\mathbf{δ}}_{i}$$

where ${\mathbf{δ}}_{i}$ is a vector of dimension $k$, and $K_{(b)}$ is
the constant mentioned earlier. There are multiple ways to do this
decomposition. Two particularly useful methods are to either use an
eigendecomposition, as suggested by Fuller (1998), or instead use
replicate estimates from the first-phase survey, as suggested by Opsomer
and Erciulescu (2021).

Fuller demonstrates that we can obtain a reasonable variance estimator
for the two-phase calibration estimator by using these $R_{(b)}$ vectors
${\mathbf{δ}}_{r}$ to form $R_{(b)}$ different control totals to use as
the calibration targets for the $R_{(b)}$ second-phase replicates. In
other words, we simply calibrate the $r$-th set of replicate weights to
the $r$-th control total
${\widehat{\mathbf{X}}}^{(a)} + {\mathbf{δ}}_{r}$. Crucially, the order
of the vectors ${\mathbf{δ}}_{r}$ should be totally random, so that the
vectors ${\mathbf{δ}}_{r}$ are independent of the sets of replicate
weights $\mathbf{w}_{r}$.

Fuller (1998) shows that calibrating the second-phase replicates to the
random calibration targets described above results in a variance
estimator that is consistent for the variance of the two-phase
calibration estimator.

This is the underlying estimator described in R code earlier in this
vignette through the use of the functions
[`calibrate_to_estimate()`](https://bschneidr.github.io/svrep/reference/calibrate_to_estimate.md)
or
[`calibrate_to_sample()`](https://bschneidr.github.io/svrep/reference/calibrate_to_sample.md).
The essential difference between the two functions is in how they form
the vectors ${\mathbf{δ}}_{r}$.

The function
[`calibrate_to_estimate()`](https://bschneidr.github.io/svrep/reference/calibrate_to_estimate.md)
forms the vectors ${\mathbf{δ}}_{r}$ using an eigen-decomposition of a
specified variance-covariance matrix.

``` r
# Print first phase estimates and their variance-covariance
print(first_phase_totals)
#>       TOTCIR     TOTSTAFF 
#> 1648795905.4     152846.6
print(first_phase_vcov)
#>                TOTCIR     TOTSTAFF
#> TOTCIR   6.606150e+16 5.853993e+12
#> TOTSTAFF 5.853993e+12 5.747174e+08
#> attr(,"means")
#> [1] 1645727222.2     152190.5

# Calibrate the two-phase replicate design
# to the totals estimated from the first-phase sample
calibrated_twophase_design <- calibrate_to_estimate(
  rep_design = twophase_boot_design,
  # Specify the variables in the data to use for calibration
  cal_formula = ~ TOTCIR + TOTSTAFF,
  # Supply the first-phase estimates and their variance
  estimate = first_phase_totals,
  vcov_estimate = first_phase_vcov,
)
#> Selection of replicate columns whose control totals will be perturbed will be done at random.
#> For tips on reproducible selection, see `help('calibrate_to_estimate')`
```

In contrast, the function
[`calibrate_to_sample()`](https://bschneidr.github.io/svrep/reference/calibrate_to_sample.md)
forms the vectors ${\mathbf{δ}}_{r}$ by using replicate estimates from
the first-phase sample.

``` r
calibrated_twophase_design <- calibrate_to_sample(
  primary_rep_design = twophase_boot_design,
  # Supply the first-phase replicate design
  control_rep_design = first_phase_gen_boot,
  # Specify the variables in the data to use for calibration
  cal_formula = ~ TOTCIR + TOTSTAFF
)
#> Matching between primary and control replicates will be done at random.
#> For tips on reproducible matching, see `help('calibrate_to_sample')`
```

## Ensuring the Variance Estimator is Positive Semidefinite

If you’ve made it this far in the vignette, then you’re probably now
well-aware that variance estimators for two-phase designs are often not
the positive semidefinite quadratic form we’d like them to be. Instead,
they’re usually close to but not quite a positive semidefinite quadratic
form, owing to the difficulty of estimating the first-phase variance
component.[⁵](#fn5)

One solution for handling a quadratic form matrix $\Sigma_{ab}$ which is
not positive semidefinite is to approximate it by
${\widetilde{\Sigma}}_{ab} = \Gamma\Lambda^{*}\Gamma^{\prime}$, where
$\Gamma$ is a matrix of eigenvalues of $\Sigma_{ab}$, $\Lambda$ is the
diagonal matrix of eigenvalues of $\Sigma_{ab}$, and $\Lambda^{*}$ is an
updated version of $\Lambda$ where negative eigenvalues have been
replaced by $0$. This solution is suggested by Beaumont and Patak (2012)
as a general-purpose solution for implementing the generalized bootstrap
when the target variance estimator that it’s mimicking isn’t positive
semidefinite. Beaumont and Patak (2012) argue that using
${\widetilde{\Sigma}}_{ab}$ instead of $\Sigma_{ab}$ should only result
in a small overestimation.

### Usage with the Generalized Bootstrap

When the function
[`as_gen_boot_design()`](https://bschneidr.github.io/svrep/reference/as_gen_boot_design.md)
is used to create generalized bootstrap replicate weights, it will warn
you if the target variance estimator is not positive semidefinite and
will let you know that it will therefore approximate the target variance
estimator using the method described above.

``` r
gen_boot_design <- as_gen_boot_design(
  design = twophase_design,
  variance_estimator = list(
    'Phase 1' = "Ultimate Cluster",
    'Phase 2' = "Ultimate Cluster"
  )
)
#> Warning in as_gen_boot_design.twophase2(design = twophase_design,
#> variance_estimator = list(`Phase 1` = "Ultimate Cluster", : The sample
#> quadratic form matrix for this design and variance estimator is not positive
#> semidefinite. It will be approximated by the nearest positive semidefinite
#> matrix.
```

### Helper Functions for Ensuring an Estimator is Positive Semidefinite

The ‘svrep’ package has two functions which can be helpful for dealing
with matrices which we hope are positive semidefinite but which might
not be.

The function
[`is_psd_matrix()`](https://bschneidr.github.io/svrep/reference/is_psd_matrix.md)
simply checks whether a matrix is positive semidefinite. It works by
estimating the matrix’s eigenvalues and determining whether any of them
are negative.

``` r
twophase_quad_form_matrix <- get_design_quad_form(
  design = twophase_design,
  variance_estimator = list(
    'Phase 1' = "Ultimate Cluster",
    'Phase 2' = "Ultimate Cluster"
  )
)

twophase_quad_form_matrix |> is_psd_matrix()
#> [1] FALSE
```

If the matrix isn’t positive semidefinite (but is at least symmetric),
then the function
[`get_nearest_psd_matrix()`](https://bschneidr.github.io/svrep/reference/get_nearest_psd_matrix.md)
will implement the approximation method described earlier.

Approximating this quadratic form by one which *is* positive
semidefinite leads to a very similar (but slightly larger) estimated
standard error.

``` r
approx_quad_form <- get_nearest_psd_matrix(twophase_quad_form_matrix)
```

In the example two-phase design based on the library survey from
earlier, we can see that this approximation results in a standard error
estimate which is only slightly larger than the standard error estimate
based on the quadratic form which wasn’t quite positive semidefinite.

``` r
# Extract weights and a single variable from the second-phase sample
## NOTE: To get second-phase data,
##       we use `my_design$phase1$sample$variables`.
##       To get first-phase data,
##       we use `my_design$phase1$full$variables

wts <- weights(twophase_design, type = "sampling")
y <- twophase_design$phase1$sample$variables$TOTSTAFF
wtd_y <- as.matrix(wts * y)

# Estimate standard errors
std_error <- as.numeric(
  t(wtd_y) %*% twophase_quad_form_matrix %*% wtd_y
) |> sqrt()

approx_std_error <- as.numeric(
  t(wtd_y) %*% approx_quad_form %*% wtd_y
) |> sqrt()

print(approx_std_error)
#> [1] 20498.68
print(std_error)
#> [1] 19765.59

approx_std_error / std_error
#> [1] 1.037089
```

## References

Beaumont, Jean-François, and Zdenek Patak. 2012. “On the Generalized
Bootstrap for Sample Surveys with Special Attention to Poisson Sampling:
*Generalized Bootstrap for Sample Surveys*.” *International Statistical
Review* 80 (1): 127–48.
<https://doi.org/10.1111/j.1751-5823.2011.00166.x>.

Fuller, Wayne A. 1998. “Replication Variance Estimation for Two-Phase
Samples.” *Statistica Sinica* 8 (4): 1153–64.

Lohr, Sharon L. 2022. *Sampling: Design and Analysis*. Third edition.
Chapman & Hall CRC Texts in Statistical Science. Boca Raton: CRC Press.

Opsomer, J. D., and A. L. Erciulescu. 2021. “Replication Variance
Estimation After Sample-Based Calibration.” *Survey Methodology,
Statistics Canada* Vol. 47 (No. 2).

Särndal, Carl-Erik, Bengt Swensson, and Jan Wretman. 1992. *Model
Assisted Survey Sampling*. Springer Series in Statistics. New York, NY:
Springer New York.

------------------------------------------------------------------------

1.  This method is justified by the arguments made in Fuller (1998),
    although the idea of using one sample’s replicate estimates as
    control totals for another survey’s replicate estimates appears to
    have first been proposed by Opsomer and Erciulescu (2021), in the
    context of two independent samples.

2.  When a two-phase variance estimator is not positive semidefinite,
    the culprit is usually the matrix $D_{b}$. The matrix $D_{b}$ is
    frequently not positive semidefinite and as a result causes the
    whole quadratic form not to be positive semidefinite. The matrices
    $\Sigma_{a}$ and $\Sigma_{b}$ are positive semidefinite for the
    usual single-phase variance estimators for most single-phase designs
    in practice, and so $\Sigma_{a^{\prime}}$ is too (since it is a
    principal submatrix of $\Sigma_{a}$). But when $D_{b}$ isn’t also
    positive semidefinite, then
    $\left( \Sigma_{a^{\prime}} \circ D_{b} \right)$ generally won’t be
    either. For one simple example of a $D_{b}$ which is not positive
    semidefinite, consider a simple random sample without replacement of
    size $n = 2$ from a population of $N = 4$. The matrix
    $D_{b} = (\begin{matrix}
    (1/2)^{- 1} & (1/6)^{- 1} \\
    (1/6)^{- 1} & (1/2)^{- 1}
    \end{matrix}) = (\begin{matrix}
    2 & 6 \\
    6 & 2
    \end{matrix})$ is not positive semidefinite, since it has
    eigenvalues $8$ and $- 4$.

3.  See the vignette “Bootstrap Methods for Surveys”

4.  The specific type of calibration used depends both on the number and
    nature of the variables available and on the “working model” chosen
    to characterize the relationship between the variables $y$ and
    $x_{1},\ldots,x_{p}$. The choice of calibration method affects the
    choice of predictors $x_{1},...,x_{p}$ as well as the choice of
    constants $\sigma_{i}$ used in the GREG model.

5.  Recall that the first-phase variance component is estimated by an
    estimator of an estimator. To be precise,
    $\widehat{V}\left\lbrack {\widehat{Y}}^{(a)} \right\rbrack$ is an
    estimate using the second-phase sample $s_{b}$ of the variance
    estimator
    $\widetilde{V}\left\lbrack {\widehat{Y}}^{(a)} \right\rbrack$ based
    on $s_{a}$ that we would have used if we had observed values of $y$
    for the entire first-phase sample.
