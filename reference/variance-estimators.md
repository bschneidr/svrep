# Variance Estimators

This help page describes variance estimators which are commonly used for
survey samples. These variance estimators can be used as the basis of
the generalized replication methods, implemented with the functions
[`as_fays_gen_rep_design()`](https://bschneidr.github.io/svrep/reference/as_fays_gen_rep_design.md),
[`as_gen_boot_design()`](https://bschneidr.github.io/svrep/reference/as_gen_boot_design.md),
[`make_fays_gen_rep_factors()`](https://bschneidr.github.io/svrep/reference/make_fays_gen_rep_factors.md),
or
[`make_gen_boot_factors()`](https://bschneidr.github.io/svrep/reference/make_gen_boot_factors.md)

## Shared Notation

Let \\s\\ denote the selected sample of size \\n\\, with elements
\\i=1,\dots,n\\. Element \\i\\ in the sample had probability \\\pi_i\\
of being included in the sample. The *pair* of elements \\ij\\ was
sampled with probability \\\pi\_{ij}\\.

The population total for a variable is denoted \\Y = \sum\_{i \in
U}y_i\\, and the Horvitz-Thompson estimator for \\\hat{Y}\\ is denoted
\\\hat{Y} = \sum\_{i \in s} y_i/\pi_i\\. For convenience, we denote
\\\breve{y}\_i = y_i/\pi_i\\.

The true sampling variance of \\\hat{Y}\\ is denoted \\V(\hat{Y})\\,
while an estimator of this sampling variance is denoted \\v(\hat{Y})\\.

## Horvitz-Thompson

The **Horvitz-Thompson** variance estimator: \$\$ v(\hat{Y}) = \sum\_{i
\in s}\sum\_{j \in s} (1 - \frac{\pi_i \pi_j}{\pi\_{ij}})
\frac{y_i}{\pi_i} \frac{y_j}{\pi_j} \$\$

## Yates-Grundy

The **Yates-Grundy** variance estimator: \$\$ v(\hat{Y}) =
-\frac{1}{2}\sum\_{i \in s}\sum\_{j \in s} (1 - \frac{\pi_i
\pi_j}{\pi\_{ij}}) (\frac{y_i}{\pi_i} - \frac{y_j}{\pi_j})^2 \$\$

## Poisson Horvitz-Thompson

The **Poisson Horvitz-Thompson** variance estimator is simply the
Horvitz-Thompson variance estimator, but where \\\pi\_{ij}=\pi_i \times
\pi_j\\, which is the case for Poisson sampling.  
  

## Stratified Multistage SRS

The **Stratified Multistage SRS** variance estimator is the recursive
variance estimator proposed by Bellhouse (1985) and used in the 'survey'
package's function
[svyrecvar](https://rdrr.io/pkg/survey/man/svyrecvar.html). In the case
of simple random sampling without replacement (with one or more stages),
this estimator exactly matches the Horvitz-Thompson estimator.

The estimator can be used for any number of sampling stages. For
illustration, we describe its use for two sampling stages. \$\$
v(\hat{Y}) = \hat{V}\_1 + \hat{V}\_2 \$\$ where \$\$ \hat{V}\_1 =
\sum\_{h=1}^{H} (1 - \frac{n_h}{N_h})\frac{n_h}{n_h - 1}
\sum\_{i=1}^{n_h} (y\_{hi.} - \bar{y}\_{hi.})^2 \$\$ and \$\$ \hat{V}\_2
= \sum\_{h=1}^{H} \frac{n_h}{N_h} \sum\_{i=1}^{n_h}v\_{hi}(y\_{hi.})
\$\$ where \\n_h\\ is the number of sampled clusters in stratum \\h\\,
\\N_h\\ is the number of population clusters in stratum \\h\\,
\\y\_{hi.}\\ is the weighted cluster total in cluster \\i\\ of stratum
\\h\\, \\\bar{y}\_{hi.}\\ is the mean weighted cluster total of stratum
\\h\\, (\\\bar{y}\_{hi.} = \frac{1}{n_h}\sum\_{i=1}^{n_h}y\_{hi.}\\),
and \\v\_{hi}(y\_{hi.})\\ is the estimated sampling variance of
\\y\_{hi.}\\.  
  

## Ultimate Cluster

The **Ultimate Cluster** variance estimator is simply the stratified
multistage SRS variance estimator, but ignoring variances from later
stages of sampling. \$\$ v(\hat{Y}) = \hat{V}\_1 \$\$ This is the
variance estimator used in the 'survey' package when the user specifies
`option(survey.ultimate.cluster = TRUE)` or uses
`svyrecvar(..., one.stage = TRUE)`. When the first-stage sampling
fractions are small, analysts often omit the finite population
corrections \\(1-\frac{n_h}{N_h})\\ when using the ultimate cluster
estimator.

## SD1 and SD2 (Successive Difference Estimators)

The **SD1** and **SD2** variance estimators are "successive difference"
estimators sometimes used for systematic sampling designs. Ash (2014)
describes each estimator as follows: \$\$ \hat{v}\_{S D 1}(\hat{Y}) =
\left(1-\frac{n}{N}\right) \frac{n}{2(n-1)}
\sum\_{k=2}^n\left(\breve{y}\_k-\breve{y}\_{k-1}\right)^2 \$\$ \$\$
\hat{v}\_{S D 2}(\hat{Y}) = \left(1-\frac{n}{N}\right)
\frac{1}{2}\left\[\sum\_{k=2}^n\left(\breve{y}\_k-\breve{y}\_{k-1}\right)^2+\left(\breve{y}\_n-\breve{y}\_1\right)^2\right\]
\$\$ where \\\breve{y}\_k = y_k/\pi_k\\ is the weighted value of unit
\\k\\ with selection probability \\\pi_k\\. The SD1 estimator is
recommended by Wolter (1984). The SD2 estimator is the basis of the
successive difference replication estimator commonly used for systematic
sampling designs and is more conservative. See Ash (2014) for details.  
  
For multistage samples, SD1 and SD2 are applied to the clusters at each
stage, separately by stratum. For later stages of sampling, the variance
estimate from a stratum is multiplied by the product of sampling
fractions from earlier stages of sampling. For example, at a third stage
of sampling, the variance estimate from a third-stage stratum is
multiplied by \\\frac{n_1}{N_1}\frac{n_2}{N_2}\\, which is the product
of sampling fractions from the first-stage stratum and second-stage
stratum.

## Beaumont-Emond

The **"Beaumont-Emond"** variance estimator was proposed by Beaumont and
Emond (2022), intended for designs that use fixed-size,
unequal-probability random sampling without replacement. The variance
estimator is simply the Horvitz-Thompson variance estimator with the
following approximation for the joint inclusion probabilities. \$\$
\pi\_{kl} \approx \pi_k \pi_l \frac{n - 1}{(n-1) +
\sqrt{(1-\pi_k)(1-\pi_l)}} \$\$ In the case of cluster sampling, this
approximation is applied to the clusters rather than the units within
clusters, with \\n\\ denoting the number of sampled clusters. and the
probabilities \\\pi\\ referring to the cluster's sampling probability.
For stratified samples, the joint probability for units \\k\\ and \\l\\
in different strata is simply the product of \\\pi_k\\ and \\\pi_l\\.

For multistage samples, this approximation is applied to the clusters at
each stage, separately by stratum. For later stages of sampling, the
variance estimate from a stratum is multiplied by the product of
sampling probabilities from earlier stages of sampling. For example, at
a third stage of sampling, the variance estimate from a third-stage
stratum is multiplied by \\\pi_1 \times \pi\_{(2 \| 1)}\\, where
\\\pi_1\\ is the sampling probability of the first-stage unit and
\\\pi\_{(2\|1)}\\ is the sampling probability of the second-stage unit
within the first-stage unit.

## Deville 1 and Deville 2

The **"Deville-1"** and **"Deville-2"** variance estimators are clearly
described in Matei and Tillé (2005), and are intended for designs that
use fixed-size, unequal-probability random sampling without replacement.
These variance estimators have been shown to be effective for designs
that use a fixed sample size with a high-entropy sampling method. This
includes most PPSWOR sampling methods, but unequal-probability
systematic sampling is an important exception.

These variance estimators take the following form: \$\$ \hat{v}(\hat{Y})
= \sum\_{i=1}^{n} c_i (\breve{y}\_i -
\frac{1}{\sum\_{i=k}^{n}c_k}\sum\_{k=1}^{n}c_k \breve{y}\_k)^2 \$\$
where \\\breve{y}\_i = y_i/\pi_i\\ is the weighted value of the the
variable of interest, and \\c_i\\ depend on the method used:

- **"Deville-1"**: \$\$c_i=\left(1-\pi_i\right) \frac{n}{n-1}\$\$

- **"Deville-2"**: \$\$c_i = (1-\pi_i) \left\[1 - \sum\_{k=1}^{n}
  \left(\frac{1-\pi_k}{\sum\_{k=1}^{n}(1-\pi_k)}\right)^2
  \right\]^{-1}\$\$

In the case of simple random sampling without replacement (SRSWOR),
these estimators are both identical to the usual stratified multistage
SRS estimator (which is itself a special case of the Horvitz-Thompson
estimator).

For multistage samples, "Deville-1" and "Deville-2" are applied to the
clusters at each stage, separately by stratum. For later stages of
sampling, the variance estimate from a stratum is multiplied by the
product of sampling probabilities from earlier stages of sampling. For
example, at a third stage of sampling, the variance estimate from a
third-stage stratum is multiplied by \\\pi_1 \times \pi\_{(2 \| 1)}\\,
where \\\pi_1\\ is the sampling probability of the first-stage unit and
\\\pi\_{(2\|1)}\\ is the sampling probability of the second-stage unit
within the first-stage unit.

## BOSB

This kernel-based variance estimator was proposed by Breidt, Opsomer,
and Sanchez-Borrego (2016), for use with samples selected using
systematic sampling or where only a single sampling unit is selected
from each stratum (sometimes referred to as "fine stratification").

Suppose there are \\n\\ sampled units, and for each unit \\i\\ there is
a numeric population characteristic \\x_i\\ and there is a weighted
total \\\hat{Y}\_i\\, where \\\hat{Y}\_i\\ is only observed in the
selected sample but \\x_i\\ is known prior to sampling.

The variance estimator has the following form:

\$\$ \hat{V}\_{ker}=\frac{1}{C_d} \sum\_{i=1}^n
(\hat{Y}\_i-\sum\_{j=1}^n d_j(i) \hat{Y}\_j)^2 \$\$

The terms \\d_j(i)\\ are kernel weights given by

\$\$ d_j(i)=\frac{K(\frac{x_i-x_j}{h})}{\sum\_{j=1}^n
K(\frac{x_i-x_j}{h})} \$\$

where \\K(\cdot)\\ is a symmetric, bounded kernel function and \\h\\ is
a bandwidth parameter. The normalizing constant \\C_d\\ is computed as:

\$\$ C_d=\frac{1}{n} \sum\_{i=1}^n(1-2 d_i(i)+\sum\_{j=1}^H d_j^2(i))
\$\$

For most functions in the 'svrep' package, the kernel function is the
Epanechnikov kernel and the bandwidth is automatically selected to yield
the smallest possible nonempty kernel window, as was recommended by
Breidt, Opsomer, and Sanchez-Borrego (2016). That's the case for the
functions
[`as_fays_gen_rep_design()`](https://bschneidr.github.io/svrep/reference/as_fays_gen_rep_design.md),
[`as_gen_boot_design()`](https://bschneidr.github.io/svrep/reference/as_gen_boot_design.md),
[`make_quad_form_matrix()`](https://bschneidr.github.io/svrep/reference/make_quad_form_matrix.md),
etc. However, users can construct the quadratic form matrix of this
variance estimator using a different kernel and a different bandwidth by
directly working with the function
[`make_kernel_var_matrix()`](https://bschneidr.github.io/svrep/reference/make_kernel_var_matrix.md).

## Deville-Tillé

See Section 6.8 of Tillé (2020) for more detail on this estimator,
including an explanation of its quadratic form. See Deville and Tillé
(2005) for the results of a simulation study comparing this and other
alternative estimators for balanced sampling.

The estimator can be written as follows: \$\$ v(\hat{Y})=\sum\_{k \in S}
\frac{c_k}{\pi_k^2}\left(y_k-\hat{y}\_k^\*\right)^2, \$\$ where \$\$
\hat{y}\_k^\*=\mathbf{z}\_k^{\top}\left(\sum\_{\ell \in S} c\_{\ell}
\frac{\mathbf{z}\_{\ell}
\mathbf{z}\_{\ell}^{\prime}}{\pi\_{\ell}^2}\right)^{-1} \sum\_{\ell \in
S} c\_{\ell} \frac{\mathbf{z}\_{\ell} y\_{\ell}}{\pi\_{\ell}^2} \$\$ and
\\\mathbf{z}\_k\\ denotes the vector of auxiliary variables for
observation \\k\\ included in sample \\S\\, with inclusion probability
\\\pi_k\\. The value \\c_k\\ is set to \\\frac{n}{n-q}(1-\pi_k)\\, where
\\n\\ is the number of observations and \\q\\ is the number of auxiliary
variables.

## References

Ash, S. (2014). "*Using successive difference replication for estimating
variances*." **Survey Methodology**, Statistics Canada, 40(1), 47-59.

Beaumont, J.F.; Émond, N. (2022). "*A Bootstrap Variance Estimation
Method for Multistage Sampling and Two-Phase Sampling When Poisson
Sampling Is Used at the Second Phase*." **Stats**, *5*: 339-357.
https://doi.org/10.3390/stats5020019

Bellhouse, D.R. (1985). "*Computing Methods for Variance Estimation in
Complex Surveys*." **Journal of Official Statistics**, Vol.1, No.3.

Breidt, F. J., Opsomer, J. D., & Sanchez-Borrego, I. (2016).
"*Nonparametric Variance Estimation Under Fine Stratification: An
Alternative to Collapsed Strata*." **Journal of the American Statistical
Association**, 111(514), 822-833.
https://doi.org/10.1080/01621459.2015.1058264

Deville, J.C., and Tillé, Y. (2005). "*Variance approximation under
balanced sampling.*" **Journal of Statistical Planning and Inference**,
128, 569-591.

Tillé, Y. (2020). "*Sampling and estimation from finite populations*."
(I. Hekimi, Trans.). Wiley.

Matei, Alina, and Yves Tillé. (2005). "*Evaluation of Variance
Approximations and Estimators in Maximum Entropy Sampling with Unequal
Probability and Fixed Sample Size.*" **Journal of Official Statistics**,
21(4):543-70.
