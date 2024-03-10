#' @title Public Libraries Survey (PLS): A Census of U.S. Public Libraries in FY2020
#'
#' @description Data taken from a complete census of public libraries in the United States in FY2020 (April 2020 to March 2021).
#' The Public Libraries Survey (PLS) is an annual census of public libraries in the U.S.,
#' including all public libraries identified by state library administrative
#' agencies in the 50 states, the District of Columbia, and the outlying territories
#' of American Samoa, Guam, the Northern Mariana Islands, and the U.S. Virgin Islands
#' (Puerto Rico did not participate in FY2020). \cr \cr
#' The primary dataset, \code{library_census}, represents the full microdata from the census.
#' The datasets \code{library_multistage_sample} and \code{library_stsys_sample}
#' are samples drawn from \code{library_census} using different sampling methods.
#' @format
#' \emph{Library Census (\code{library_census}): } \cr \cr
#' The dataset includes 9,245 records (one per library) and 23 variables.
#' Each column has a variable label, accessible using the function \code{var_label()} from the 'labelled' package
#' or simply by calling \code{attr(x, 'label')} to a given column.
#' These data include a subset of the variables included in the public-use data published by PLS,
#' specifically from the Public Library System Data File. Particularly relevant variables include:
#' \cr \cr
#' Identifier variables and survey response status:
#' \itemize{
#'   \item FSCSKEY: A unique identifier for libraries.
#'   \item LIBNAME: The name of the library.
#'   \item RESPONSE_STATUS: Response status for the Public Library Survey:
#'     indicates whether the library was a respondent, nonrespondent, or was closed.
#' }
#' Numeric summaries:
#' \itemize{
#'   \item TOTCIR: Total circulation
#'   \item VISITS: Total visitors
#'   \item REGBOR: Total number of registered users
#'   \item TOTSTAFF: Total staff (measured in full-time equivalent staff)
#'   \item LIBRARIA: Total librarians (measured in full-time equivalent staff)
#'   \item TOTOPEXP: Total operating expenses
#'   \item TOTINCM: Total income
#'   \item BRANLIB: Number of library branches
#'   \item CENTLIB: Number of central library locations
#' }
#' Location:
#' \itemize{
#'   \item LONGITUD: Geocoded longitude (in WGS84 CRS)
#'   \item LATITUD: Geocoded latitude (in WGS84 CRS)
#'   \item STABR: Two-letter state abbreviation
#'   \item CBSA: Five-digit identifer for a core-based statistical area (CBSA)
#'   \item MICROF: Flag for a metropolitan or micropolitan statistical area
#'
#' }
#' @references
#' Pelczar, M., Soffronoff, J., Nielsen, E., Li, J., & Mabile, S. (2022). Data File Documentation: Public
#' Libraries in the United States Fiscal Year 2020. Institute of Museum and Library Services: Washington,
#' D.C.
#'
#' @keywords datasets libraries
#' @name libraries
#' @usage data(library_census)
"library_census"

#' @rdname libraries
#' @usage data(library_multistage_sample)
#' @format
#'
#' \emph{Library Multistage Sample (\code{library_multistage_sample}): } \cr \cr
#' These data represent a two-stage sample (PSUs and SSUs),
#' where the first stage sample is selected using unequal probability sampling
#' without replacement (PPSWOR) and the second stage sample is selected
#' using simple random sampling without replacement (SRSWOR). \cr \cr
#' Includes the same variables as \code{library_census},
#' but with additional design variables.
#' \itemize{
#'   \item PSU_ID: A unique identifier for primary sampling units
#'   \item SSU_ID: A unique identifer for secondary sampling units
#'   \item SAMPLING_PROB: Overall inclusion probability
#'   \item PSU_SAMPLING_PROB: Inclusion probability for the PSU
#'   \item SSU_SAMPLING_PROB: Inclusion probability for the SSU
#'   \item PSU_POP_SIZE: The number of PSUs in the population
#'   \item SSU_POP_SIZE: The number of population SSUs within the PSU
#' }
"library_multistage_sample"

#' @rdname libraries
#' @usage data(library_stsys_sample)
#' @format
#'
#' \emph{Library Stratified Systematic Sample (\code{library_stsys_sample}): } \cr \cr
#' These data represent a stratified systematic sample. \cr \cr
#' Includes the same variables as \code{library_census},
#' but with additional design variables.
#' \itemize{
#'   \item SAMPLING_STRATUM: Unique identifier for sampling strata
#'   \item STRATUM_POP_SIZE: The population size in the stratum
#'   \item SAMPLING_SORT_ORDER: The sort order used before selecting a random systematic sample
#'   \item SAMPLING_PROB: Overall inclusion probability
#' }
"library_stsys_sample"
