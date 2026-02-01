# Create a quadratic form's matrix to represent the basic variance estimator for a total under simple random sampling without replacement

The usual variance estimator for simple random sampling without
replacement can be represented as a quadratic form. This function
determines the matrix of the quadratic form.

## Usage

``` r
make_srswor_matrix(n, f = 0)
```

## Arguments

- n:

  Sample size

- f:

  A single number between `0` and `1`, representing the sampling
  fraction. Default value is `0`.

## Value

A symmetric matrix of dimension `n`

## Details

The basic variance estimator of a total for simple random sampling
without replacement is as follows: \$\$ \hat{v}(\hat{Y}) = (1 -
f)\frac{n}{n - 1} \sum\_{i=1}^{n} (y_i - \bar{y})^2 \$\$ where \\f\\ is
the sampling fraction \\\frac{n}{N}\\.  
  
If \\f=0\\, then the matrix of the quadratic form has all non-diagonal
elements equal to \\-(n-1)^{-1}\\, and all diagonal elements equal to
\\1\\. If \\f \> 0\\, then each element is multiplied by \\(1-f)\\.  
  
If \\n=1\\, then this function returns a \\1 \times 1\\ matrix whose
sole element equals \\0\\ (essentially treating the sole sampled unit as
a selection made with probability \\1\\).
