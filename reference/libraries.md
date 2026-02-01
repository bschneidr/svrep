# U.S. Public Libraries Survey

Data taken from a complete census of public libraries in the United
States in FY2020 (April 2020 to March 2021). The Public Libraries Survey
(PLS) is an annual census of public libraries in the U.S., including all
public libraries identified by state library administrative agencies in
the 50 states, the District of Columbia, and the outlying territories of
American Samoa, Guam, the Northern Mariana Islands, and the U.S. Virgin
Islands (Puerto Rico did not participate in FY2020).  
  
The primary dataset, `library_census`, represents the full microdata
from the census. The datasets `library_multistage_sample` and
`library_stsys_sample` are samples drawn from `library_census` using
different sampling methods.

## Usage

``` r
data(library_census)

data(library_multistage_sample)

data(library_stsys_sample)
```

## Format

*Library Census (`library_census`):*  
  
The dataset includes 9,245 records (one per library) and 23 variables.
Each column has a variable label, accessible using the function
`var_label()` from the 'labelled' package or simply by calling
`attr(x, 'label')` to a given column. These data include a subset of the
variables included in the public-use data published by PLS, specifically
from the Public Library System Data File. Particularly relevant
variables include:  
  
Identifier variables and survey response status:

- FSCSKEY: A unique identifier for libraries.

- LIBNAME: The name of the library.

- RESPONSE_STATUS: Response status for the Public Library Survey:
  indicates whether the library was a respondent, nonrespondent, or was
  closed.

Numeric summaries:

- TOTCIR: Total circulation

- VISITS: Total visitors

- REGBOR: Total number of registered users

- TOTSTAFF: Total staff (measured in full-time equivalent staff)

- LIBRARIA: Total librarians (measured in full-time equivalent staff)

- TOTOPEXP: Total operating expenses

- TOTINCM: Total income

- BRANLIB: Number of library branches

- CENTLIB: Number of central library locations

Location:

- LONGITUD: Geocoded longitude (in WGS84 CRS)

- LATITUD: Geocoded latitude (in WGS84 CRS)

- STABR: Two-letter state abbreviation

- CBSA: Five-digit identifier for a core-based statistical area (CBSA)

- MICROF: Flag for a metropolitan or micropolitan statistical area

*Library Multistage Sample (`library_multistage_sample`):*  
  
These data represent a two-stage sample (PSUs and SSUs), where the first
stage sample is selected using unequal probability sampling without
replacement (PPSWOR) and the second stage sample is selected using
simple random sampling without replacement (SRSWOR).  
  
Includes the same variables as `library_census`, but with additional
design variables.

- PSU_ID: A unique identifier for primary sampling units

- SSU_ID: A unique identifier for secondary sampling units

- SAMPLING_PROB: Overall inclusion probability

- PSU_SAMPLING_PROB: Inclusion probability for the PSU

- SSU_SAMPLING_PROB: Inclusion probability for the SSU

- PSU_POP_SIZE: The number of PSUs in the population

- SSU_POP_SIZE: The number of population SSUs within the PSU

*Library Stratified Systematic Sample (`library_stsys_sample`):*  
  
These data represent a stratified systematic sample.  
  
Includes the same variables as `library_census`, but with additional
design variables.

- SAMPLING_STRATUM: Unique identifier for sampling strata

- STRATUM_POP_SIZE: The population size in the stratum

- SAMPLING_SORT_ORDER: The sort order used before selecting a random
  systematic sample

- SAMPLING_PROB: Overall inclusion probability

## References

Pelczar, M., Soffronoff, J., Nielsen, E., Li, J., & Mabile, S. (2022).
Data File Documentation: Public Libraries in the United States Fiscal
Year 2020. Institute of Museum and Library Services: Washington, D.C.
