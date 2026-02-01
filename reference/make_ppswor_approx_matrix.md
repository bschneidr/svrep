# Create a quadratic form's matrix to represent a variance estimator for PPSWOR designs, based on commonly-used approximations

Several variance estimators for designs that use unequal probability
sampling without replacement (i.e., PPSWOR), variance estimation tends
to be more accurate when using an approximation estimator that uses the
first-order inclusion probabilities (i.e., the basic sampling weights)
and ignores the joint inclusion probabilities. This function returns the
matrix of the quadratic form used to represent such variance estimators.

## Usage

``` r
make_ppswor_approx_matrix(probs, method = "Deville-1")
```

## Arguments

- probs:

  A vector of first-order inclusion probabilities

- method:

  A string specifying the approximation method to use. See the "Details"
  section below. Options include:

  - "Deville-1"

  - "Deville-2"

  - "Beaumont-Emond"

## Value

A symmetric matrix whose dimension matches the length of `probs`.

## Deville's Estimators

The "Deville-1" and "Deville-2" approximations have been shown to be
effective for designs that use a fixed sample size with a high-entropy
sampling method. This includes most PPSWOR sampling methods, but
unequal-probability systematic sampling is an important exception.

Deville's variance estimators generally take the following form: \$\$
\hat{v}(\hat{Y}) = \sum\_{i=1}^{n} c_i (\breve{y}\_i -
\frac{1}{\sum\_{i=k}^{n}c_k}\sum\_{k=1}^{n}c_k \breve{y}\_k)^2 \$\$
where \\\breve{y}\_i = y_i/\pi_i\\ is the weighted value of the the
variable of interest, and \\c_i\\ are constants that depend on the
approximation method used.  
  
The matrix of the quadratic form, denoted \\\Sigma\\, has its \\ij\\-th
entry defined as follows: \$\$ \sigma\_{ii} = c_i (1 -
\frac{c_i}{\sum\_{k=1}^{n}c_k}) \textit{ when } i = j \\
\sigma\_{ij}=\frac{-c_i c_j}{\sum\_{k=1}^{n}c_k} \textit{ when } i \neq
j \\ \$\$ When \\\pi\_{i} = 1\\ for every unit, then \\\sigma\_{ij}=0\\
for all \\i,j\\. If there is only one sampling unit, then
\\\sigma\_{11}=0\\; that is, the unit is treated as if it was sampled
with certainty.

The constants \\c_i\\ are defined for each approximation method as
follows, with the names taken directly from Matei and Tillé (2005).

- **"Deville-1"**: \$\$c_i=\left(1-\pi_i\right) \frac{n}{n-1}\$\$

- **"Deville-2"**: \$\$c_i = (1-\pi_i) \left\[1 - \sum\_{k=1}^{n}
  \left(\frac{1-\pi_k}{\sum\_{k=1}^{n}(1-\pi_k)}\right)^2
  \right\]^{-1}\$\$

Both of the approximations **"Deville-1"** and **"Deville-2"** were
shown in the simulation studies of Matei and Tillé (2005) to perform
much better in terms of MSE compared to the strictly-unbiased
Horvitz-Thompson and Yates-Grundy variance estimators. In the case of
simple random sampling without replacement (SRSWOR), these estimators
are identical to the usual Horvitz-Thompson variance estimator.

## Beaumont-Emond Estimator

Beaumont and Emond (2022) proposed a variance estimator for unequal
probability sampling without replacement. This estimator is simply the
Horvitz-Thompson variance estimator with the following approximation for
the joint inclusion probabilities. \$\$ \pi\_{kl} \approx \pi_k \pi_l
\frac{n - 1}{(n-1) + \sqrt{(1-\pi_k)(1-\pi_l)}} \$\$ In the case of
cluster sampling, this approximation should be applied to the clusters
rather than the units within clusters.

## References

Matei, Alina, and Yves Tillé. 2005. "Evaluation of Variance
Approximations and Estimators in Maximum Entropy Sampling with Unequal
Probability and Fixed Sample Size." Journal of Official Statistics
21(4):543-70.
