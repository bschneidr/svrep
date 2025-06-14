% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lou_vax_data.R
\docType{data}
\name{lou_vax_survey}
\alias{lou_vax_survey}
\title{Data of Louisville Vaccination Survey}
\format{
A data frame with 1,000 rows and 6 variables
\describe{
  \item{RESPONSE_STATUS}{Response status to the survey ('Respondent' or 'Nonrespondent')}
  \item{RACE_ETHNICITY}{Race and Hispanic/Latino ethnicity
  derived from RAC1P and HISP variables
  of ACS microdata and collapsed to a smaller number of categories.}
  \item{SEX}{Male or Female}
  \item{EDUC_ATTAINMENT}{Highest level of education attained ('Less than high school' or 'High school or beyond')
  derived from SCHL variable in ACS microdata and collapsed to a smaller number of categories.}
  \item{VAX_STATUS}{Covid-19 vaccination status ('Vaccinated' or 'Unvaccinated')}
  \item{SAMPLING_WEIGHT}{Sampling weight: equal for all cases since data come from a simple random sample}
}
}
\usage{
data(lou_vax_survey)
}
\description{
A survey measuring Covid-19 vaccination status and a handful of demographic variables,
based on a simple random sample of 1,000 residents of Louisville, Kentucky
with an approximately 50\% response rate. \cr

These data were created using simulation.
}
\keyword{datasets}
