# Find the smallest Hadamard order available in the 'survey' package

Identifies the order of the smallest Hadamard matrix available in the
'survey' package, such that the order is greater than or equal to `n`.
This is useful for identifying the minimum number of replicates that can
be constructed for replication methods such as BRR or Fay's generalized
replication method.

## Usage

``` r
find_minimum_hadamard_order(n)
```

## Arguments

- n:

  A single positive integer.

## Value

A single positive integer.

## Details

To get the Hadamard matrix of this size with the 'survey' package, use
`survey::hadamard(k - 1)`, where `k` is the output of
`find_minimum_hadamard_order()`.
