# Create a quadratic form's matrix to represent a successive-difference variance estimator

A successive-difference variance estimator can be represented as a
quadratic form. This function determines the matrix of the quadratic
form.

## Usage

``` r
make_sd_matrix(n, f = 0, type = "SD1")
```

## Arguments

- n:

  Number of rows or columns for the matrix

- f:

  A single number between `0` and `1`, representing the sampling
  fraction. Default value is `0`.

- type:

  Either "SD1" or "SD2". See the "Details" section for definitions.

## Value

A matrix of dimension `n`

## Details

Ash (2014) describes each estimator as follows: \$\$
\hat{v}\_{SD1}(\hat{Y}) = (1-f) \frac{n}{2(n-1)}
\sum\_{k=2}^n\left(\breve{y}\_k-\breve{y}\_{k-1}\right)^2 \$\$ \$\$
\hat{v}\_{SD2}(\hat{Y}) =
\frac{1}{2}(1-f)\left\[\sum\_{k=2}^n\left(\breve{y}\_k-\breve{y}\_{k-1}\right)^2+\left(\breve{y}\_n-\breve{y}\_1\right)^2\right\]
\$\$ where \\\breve{y}\_k\\ is the weighted value \\y_k/\pi_k\\ of unit
\\k\\ with selection probability \\\pi_k\\, and \\f\\ is the sampling
fraction \\\frac{n}{N}\\.

## References

Ash, S. (2014). "*Using successive difference replication for estimating
variances*." **Survey Methodology**, Statistics Canada, 40(1), 47-59.
