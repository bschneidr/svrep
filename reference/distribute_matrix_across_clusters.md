# Helper function to turn a cluster-level matrix into an element-level matrix by duplicating rows or columns of the matrix

Turns a cluster-level matrix into an element-level matrix by suitably
duplicating rows or columns of the matrix.

## Usage

``` r
distribute_matrix_across_clusters(
  cluster_level_matrix,
  cluster_ids,
  rows = TRUE,
  cols = TRUE
)
```

## Arguments

- cluster_level_matrix:

  A square matrix, whose number of rows/columns matches the number of
  clusters.

- cluster_ids:

  A vector of cluster identifiers. If `rows=TRUE`, the number of unique
  elements of `cluster_ids` must match the number of rows of
  `cluster_level_matrix`. If `cols=TRUE`, the number of unique elements
  of `cluster_ids` must match the number of columns of
  `cluster_level_matrix`.

- rows:

  Whether to duplicate rows of the `cluster_level_matrix` for elements
  from the same cluster.

- cols:

  Whether to duplicate columns of the `cluster_level_matrix` for
  elements from the same cluster.

## Value

The input `cluster_level_matrix` has its rows/columns duplicated so that
the number of rows (if `rows=TRUE`) or columns (if `cols=TRUE`) equals
the length of `cluster_ids`.
