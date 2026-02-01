# Row Assignment for Successive Difference Replication

Creates a row assignment matrix: for each row of a dataset of size
\\n\\, assigns two rows of a Hadamard matrix.

## Usage

``` r
assign_hadamard_rows(
  n,
  hadamard_order,
  number_of_cycles = ceiling(n/hadamard_order),
  use_first_row = TRUE,
  circular = TRUE
)
```

## Arguments

- n:

  The sample size of the data.

- hadamard_order:

  The order of the Hadamard matrix (i.e., the number of rows/columns)

- number_of_cycles:

  The number of cycles to use in the row assignment. Must be at least as
  large as `n/hadamard_order`. Only applies when `n` exceeds the number
  of available rows in the Hadamard matrix. The number of available rows
  is `hadamard_order` when `use_first_row = TRUE`, and
  `hadamard_order - 1` when `use_first_row = FALSE`.

- use_first_row:

  Whether to use the first row of the Hadamard matrix. The first row of
  a Hadamard matrix is often all 1's, and so using the first row to
  create replicate factors leads to the creation of a replicate whose
  weights exactly match the full-sample weights. Thus, using the first
  row of the Hadamard matrix may be undesirable for practical purposes,
  even if it is valid for the purpose of variance estimation.

- circular:

  `TRUE` or `FALSE`. Only applies when the number of available rows in
  the Hadamard matrix is at least as large as `n`. Whether to make a
  circular row assignment, so that the resulting successive-difference
  replication variance estimator is equivalent to the SD2 variance
  estimator rather than the SD1 variance estimator (see Ash 2014). The
  number of available rows is `hadamard_order` when
  `use_first_row = TRUE`, and `hadamard_order - 1` when
  `use_first_row = FALSE`.

## Value

A matrix with `n` rows and two columns. Each row gives the assignment of
two rows of a Hadamard matrix to the row of data.

## Details

Implements row-assignment methods described in Ash (2014) and in Fay and
Train (1995). The row-assignment method depends on the number of
available rows of the Hadamard matrix used. The number of available rows
is `hadamard_order` when `use_first_row = TRUE`, and
`hadamard_order - 1` when `use_first_row = FALSE`.

When the number of available Hadamard rows is at least as large as `n`,
then the row assignment is as follows. Let \\i=1,\dots,n\\ be the index
for the data that will receive row assignments, and let \\a_j\\ denote
row \\j\\ of a Hadamard matrix. For \\i \< n\\, the assignments are
entries \\a_i\\ and \\a\_{i+1}\\ when `use_first_row = TRUE`, and when
`use_first_row = FALSE`, the assignments are entries \\a\_{i+1}\\ and
\\a\_{i+2}\\. The assignment for \\i=n\\ depends on whether
`circular = TRUE`. If \\circular=TRUE\\, then the assignment for \\i=n\\
is entries \\a_i\\ and \\a_1\\ if `use_first_row = TRUE`, and entries
\\a\_{i+1}\\ and \\a_2\\ if `use_first_row = FALSE`.

When the number of available Hadamard rows is less than `n`, then the
row assignment method is the method denoted as RA1 in Ash (2014). This
method uses the argument `number_of_cycles` and does *not* use the
argument `circular`.

## References

Ash, S. (2014). "*Using successive difference replication for estimating
variances*." **Survey Methodology**, Statistics Canada, 40(1), 47-59.
