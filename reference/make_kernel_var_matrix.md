# Quadratic Form Matrix of Kernel-based Variance Estimator

Constructs the quadratic form matrix for the kernel-based variance
estimator of Breidt, Opsomer, and Sanchez-Borrego (2016). The bandwidth
is automatically chosen to result in the smallest possible nonempty
kernel window.

## Usage

``` r
make_kernel_var_matrix(x, kernel = "Epanechnikov", bandwidth = "auto")
```

## Arguments

- x:

  A numeric vector, giving the values of an auxiliary variable.

- kernel:

  The name of a kernel function. Currently only "Epanechnikov" is
  supported.

- bandwidth:

  The bandwidth to use for the kernel. The default value is `"auto"`,
  which means that the bandwidth will be chosen automatically to produce
  the smallest window size while ensuring that every unit has a nonempty
  window, as suggested by Breidt, Opsomer, and Sanchez-Borrego (2016).
  Otherwise, the user can supply their own value, which can be a single
  positive number.

## Value

The quadratic form matrix for the variance estimator, with dimension
equal to the length of `x`. The resulting object has an attribute
`bandwidth` that can be retrieved using `attr(Q, 'bandwidth')`

## Details

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

If \\n=2\\, then the estimator is simply the estimator used for simple
random sampling without replacement.

If \\n=1\\, then the matrix simply has an entry equal to 0.

## References

Breidt, F. J., Opsomer, J. D., & Sanchez-Borrego, I. (2016).
"*Nonparametric Variance Estimation Under Fine Stratification: An
Alternative to Collapsed Strata*." **Journal of the American Statistical
Association**, 111(514), 822-833.
https://doi.org/10.1080/01621459.2015.1058264

## Examples

``` r
# The auxiliary variable has the same value for all units
make_kernel_var_matrix(c(1, 1, 1))
#> 3 x 3 Matrix of class "dsyMatrix"
#>      [,1] [,2] [,3]
#> [1,]  1.0 -0.5 -0.5
#> [2,] -0.5  1.0 -0.5
#> [3,] -0.5 -0.5  1.0

# The auxiliary variable differs across units
make_kernel_var_matrix(c(1, 2, 3))
#> 3 x 3 Matrix of class "dsyMatrix"
#>            [,1]       [,2]       [,3]
#> [1,]  0.6440922 -0.8559078  0.2118156
#> [2,] -0.8559078  1.7118156 -0.8559078
#> [3,]  0.2118156 -0.8559078  0.6440922

# View the bandwidth that was automatically selected
Q <- make_kernel_var_matrix(c(1, 2, 4))
attr(Q, 'bandwidth')
#> [1] 3
```
