# Create a quadratic form's matrix for a Deville-Tillé variance estimator for balanced samples

Creates the quadratic form matrix for a variance estimator for balanced
samples, proposed by Deville and Tillé (2005).

## Usage

``` r
make_deville_tille_matrix(probs, aux_vars)
```

## Arguments

- probs:

  A vector of first-order inclusion probabilities

- aux_vars:

  A matrix of auxiliary variables, with the number of rows matching the
  number of elements of `probs`.

## Value

A symmetric matrix whose dimension matches the length of `probs`.

## Details

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

See Li, Chen, and Krenzke (2014) for an example of this estimator's use
as the basis for a generalized replication estimator. See Breidt and
Chauvet (2011) for a discussion of alternative simulation-based
estimators for the specific application of variance estimation for
balanced samples selected using the cube method.

## References

\- Breidt, F.J. and Chauvet, G. (2011). "Improved variance estimation
for balanced samples drawn via the cube method." Journal of Statistical
Planning and Inference, 141, 411-425.

\- Deville, J.C., and Tillé, Y. (2005). "*Variance approximation under
balanced sampling.*" **Journal of Statistical Planning and Inference**,
128, 569-591.

\- Li, J., Chen, S., and Krenzke, T. (2014). "Replication Variance
Estimation for Balanced Sampling: An Application to the PIAAC Study."
Proceedings of the Survey Research Methods Section, 2014: 985-994.
Alexandria, VA: American Statistical Association.
http://www.asasrms.org/Proceedings/papers/1984_094.pdf.

\- Tillé, Y. (2020). "*Sampling and estimation from finite
populations*." (I. Hekimi, Trans.). Wiley.
