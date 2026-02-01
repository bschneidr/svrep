# Draw a doubled half-sample

Draw a doubled half-sample using the method described in Antal and Tille
(2014)

## Usage

``` r
draw_doubled_half_sample(n)
```

## Arguments

- n:

  An integer greater than 1.

## Value

A vector of length `n`, with elements equal to 0, 1, 2, or 3. This
vector indicate how many times each unit is resampled.

When `n` is even, half the values will be 0 and the other half will be
2. When `n` is odd, the values can be 0, 1, 2, or 3.
