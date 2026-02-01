# Get variables from a database

A database helper function copied from the 'survey' package

## Usage

``` r
getvars(
  formula,
  dbconnection,
  tables,
  db.only = TRUE,
  updates = NULL,
  subset = NULL
)
```

## Arguments

- formula:

  Either a formula or a character vector giving names of variables

- dbconnection:

  A database connection

- tables:

  Name(s) of table(s) to pull from

- db.only:

  Unclear parameter inherited from the 'survey' package

- updates:

  Updates to potentially make

- subset:

  Optional indices of data to subset when returning result

## Value

A data frame
