# Factors for Fay's Generalized Replication Method

Generate a matrix of replication factors using Fay's generalized
replication method. This method yields a fully efficient variance
estimator if a sufficient number of replicates is used.

## Usage

``` r
make_fays_gen_rep_factors(Sigma, max_replicates = Inf, balanced = TRUE)
```

## Arguments

- Sigma:

  A quadratic form matrix corresponding to a target variance estimator.
  Must be positive semidefinite.

- max_replicates:

  The maximum number of replicates to allow. The function will attempt
  to create the minimum number of replicates needed to produce a
  fully-efficient variance estimator. If more replicates are needed than
  `max_replicates`, then the full number of replicates needed will be
  created, but only a random subsample will be retained.

- balanced:

  If `balanced=TRUE`, the replicates will all contribute equally to
  variance estimates, but the number of replicates needed may slightly
  increase.

## Value

A matrix of replicate factors, with the number of rows matching the
number of rows of `Sigma` and the number of columns less than or equal
to `max_replicates`. To calculate variance estimates using these
factors, use the overall scale factor given by calling
`attr(x, "scale")` on the result.

## Statistical Details

See Fay (1989) for a full explanation of Fay's generalized replication
method. This documentation provides a brief overview.

Let \\\boldsymbol{\Sigma}\\ be the quadratic form matrix for a target
variance estimator, which is assumed to be positive semidefinite.
Suppose the rank of \\\boldsymbol{\Sigma}\\ is \\k\\, and so
\\\boldsymbol{\Sigma}\\ can be represented by the spectral decomposition
of \\k\\ eigenvectors and eigenvalues, where the \\r\\-th eigenvector
and eigenvalue are denoted \\\mathbf{v}\_{(r)}\\ and \\\lambda_r\\,
respectively. \$\$ \boldsymbol{\Sigma} = \sum\_{r=1}^k \lambda_r
\mathbf{v}\_{(r)} \mathbf{v^{\prime}}\_{(r)} \$\$ If `balanced = FALSE`,
then we let \\\mathbf{H}\\ denote an identity matrix with \\k' = k\\
rows/columns. If `balanced = TRUE`, then we let \\\mathbf{H}\\ be a
Hadamard matrix (with all entries equal to \\1\\ or \\-1\\), of order
\\k^{\prime} \geq k\\. Let \\\mathbf{H}\_{mr}\\ denote the entry in row
\\m\\ and column \\r\\ of \\\mathbf{H}\\.

Then \\k^{\prime}\\ replicates are formed as follows. Let \\r\\ denote a
given replicate, with \\r = 1, ..., k^{\prime}\\, and let \\c\\ denote
some positive constant (yet to be specified).

The \\r\\-th replicate adjustment factor \\\mathbf{f}\_{r}\\ is formed
as: \$\$ \mathbf{f}\_{r} = 1 + c \sum\_{m=1}^k H\_{m r}
\lambda\_{(m)}^{\frac{1}{2}} \mathbf{v}\_{(m)} \$\$

If `balanced = FALSE`, then \\c = 1\\. If `balanced = TRUE`, then \\c =
\frac{1}{\sqrt{k^{\prime}}}\\.

If any of the replicates are negative, you can use
[`rescale_replicates`](https://bschneidr.github.io/svrep/reference/rescale_replicates.md),
which recalculates the replicate factors with a smaller value of \\c\\.

If all \\k^{\prime}\\ replicates are used, then variance estimates are
calculated as: \$\$ v\_{rep}\left(\hat{T}\_y\right) =
\sum\_{r=1}^{k^{\prime}}\left(\hat{T}\_y^{\*(r)}-\hat{T}\_y\right)^2
\$\$ For population totals, this replication variance estimator will
*exactly* match the target variance estimator if the number of
replicates \\k^{\prime}\\ matches the rank of \\\Sigma\\.

## The Number of Replicates

If `balanced=TRUE`, the number of replicates created may need to
increase slightly. This is due to the fact that a Hadamard matrix of
order \\k^{\prime} \geq k\\ is used to balance the replicates, and it
may be necessary to use order \\k^{\prime} \> k\\.

If the number of replicates \\k^{\prime}\\ is too large for practical
purposes, then one can simply retain only a random subset of \\R\\ of
the \\k^{\prime}\\ replicates. In this case, variances are calculated as
follows: \$\$ v\_{rep}\left(\hat{T}\_y\right) = \frac{k^{\prime}}{R}
\sum\_{r=1}^{R}\left(\hat{T}\_y^{\*(r)}-\hat{T}\_y\right)^2 \$\$ This is
what happens if `max_replicates` is less than the matrix rank of
`Sigma`: only a random subset of the created replicates will be
retained.

Subsampling replicates is only recommended when using `balanced=TRUE`,
since in this case every replicate contributes equally to variance
estimates. If `balanced=FALSE`, then randomly subsampling replicates is
valid but may produce large variation in variance estimates since
replicates in that case may vary greatly in their contribution to
variance estimates.

## Reproducibility

If `balanced=TRUE`, a Hadamard matrix is used as described above. The
Hadamard matrix is deterministically created using the function
[`hadamard()`](https://rdrr.io/pkg/survey/man/hadamard.html) from the
'survey' package. However, the order of rows/columns is randomly
permuted before forming replicates.

In general, column-ordering of the replicate weights is random. To
ensure exact reproducibility, it is recommended to call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before using this
function.

## References

Fay, Robert. 1989. "Theory And Application Of Replicate Weighting For
Variance Calculations." In, 495-500. Alexandria, VA: American
Statistical Association.
http://www.asasrms.org/Proceedings/papers/1989_033.pdf

## See also

Use
[`rescale_replicates`](https://bschneidr.github.io/svrep/reference/rescale_replicates.md)
to eliminate negative adjustment factors.

## Examples

``` r
# \donttest{

# Load an example dataset that uses unequal probability sampling ----
  data('election', package = 'survey')

# Create matrix to represent the Horvitz-Thompson estimator as a quadratic form ----
  n <- nrow(election_pps)
  pi <- election_jointprob
  horvitz_thompson_matrix <- matrix(nrow = n, ncol = n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      horvitz_thompson_matrix[i,j] <- 1 - (pi[i,i] * pi[j,j])/pi[i,j]
    }
  }

  ## Equivalently:

  horvitz_thompson_matrix <- make_quad_form_matrix(
    variance_estimator = "Horvitz-Thompson",
    joint_probs = election_jointprob
  )

# Make generalized replication adjustment factors ----

  adjustment_factors <- make_fays_gen_rep_factors(
    Sigma = horvitz_thompson_matrix,
    max_replicates = 50
  )
  attr(adjustment_factors, 'scale')
#> [1] 1

# Compute the Horvitz-Thompson estimate and the replication estimate

ht_estimate <- svydesign(data = election_pps, ids = ~ 1,
                         prob = diag(election_jointprob),
                         pps = ppsmat(election_jointprob)) |>
  svytotal(x = ~ Kerry)

rep_estimate <- svrepdesign(
  data = election_pps,
  weights = ~ wt,
  repweights = adjustment_factors,
  combined.weights = FALSE,
  scale = attr(adjustment_factors, 'scale'),
  rscales = rep(1, times = ncol(adjustment_factors)),
  type = "other",
  mse = TRUE
) |>
  svytotal(x = ~ Kerry)

SE(rep_estimate)
#> [1] 2523712
SE(ht_estimate)
#>         [,1]
#> [1,] 2523712
SE(rep_estimate) / SE(ht_estimate)
#>      [,1]
#> [1,]    1
# }
```
