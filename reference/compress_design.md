# Produce a compressed representation of a survey design object

Produce a compressed representation of a survey design object

## Usage

``` r
compress_design(design, vars_to_keep = NULL)
```

## Arguments

- design:

  A survey design object

- vars_to_keep:

  (Optional) A character vector of variables in the design to keep in
  the compressed design. By default, none of the variables are retained.

## Value

A list with two elements. The `design_subset` element is a a design
object with only the minimal rows needed to represent the survey design.
The `index` element links each row of the original design to a row of
`design_subset`, so that the design can be "uncompressed."
