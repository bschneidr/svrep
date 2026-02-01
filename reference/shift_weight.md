# (Internal function) Shift weight from one set of cases to another

You likely want to use `redistribute_weights` instead. The function
`shift_weight` is internal to this package and is used only
"under-the-hood."

## Usage

``` r
shift_weight(wt_set, is_upweight_case, is_downweight_case)
```

## Arguments

- wt_set:

  A numeric vector of weights

- is_upweight_case:

  A logical vector indicating cases whose weight should be increased

- is_downweight_case:

  A logical vector indicating cases whose weight should be decreased

## Value

A numeric vector of adjusted weights, of the same length as `wt_set`.
