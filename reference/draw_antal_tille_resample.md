# Draw a replicate sample using the method of Antal and Tille (2014)

Draw a replicate sample using the method described in Section 7 of Antal
and Tillé (2014).

## Usage

``` r
draw_antal_tille_resample(sel_probs)
```

## Arguments

- sel_probs:

  A vector of selection probabilities (i.e., inclusion probabilities)
  for every unit in the sample.

## Value

A vector the same length as `sel_probs`, with elements equal to 0, 1, 2,
or 3. This vector indicate how many times each unit is resampled.

## References

Antal, E. and Tillé, Y. (2014). "A new resampling method for sampling
designs without replacement: The doubled half bootstrap."
**Computational Statistics**, *29*(5), 1345-1363.
https://doi.org/10.1007/s00180-014-0495-0
