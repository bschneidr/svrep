#' @title Louisville Vaccination Survey
#'
#' @description A survey measuring Covid-19 vaccination status and a handful of demographic variables,
#' based on a simple random sample of 1,000 residents of Louisville, Kentucky
#' with an approximately 50\% response rate. \cr
#'
#' These data were created using simulation.
#'
#' @format A data frame with 1,000 rows and 6 variables
#' \describe{
#'   \item{RESPONSE_STATUS}{Response status to the survey ('Respondent' or 'Nonrespondent')}
#'   \item{RACE_ETHNICITY}{Race and Hispanic/Latino ethnicity
#'   derived from RAC1P and HISP variables
#'   of ACS microdata and collapsed to a smaller number of categories.}
#'   \item{SEX}{Male or Female}
#'   \item{EDUC_ATTAINMENT}{Highest level of education attained ('Less than high school' or 'High school or beyond')
#'   derived from SCHL variable in ACS microdata and collapsed to a smaller number of categories.}
#'   \item{VAX_STATUS}{Covid-19 vaccination status ('Vaccinated' or 'Unvaccinated')}
#'   \item{SAMPLING_WEIGHT}{Sampling weight: equal for all cases since data come from a simple random sample}
#' }
#'
#' @keywords datasets
#' @name lou_vax_survey
#' @usage data(lou_vax_survey)
"lou_vax_survey"

#' @title ACS PUMS Data for Louisville
#'
#' @description Person-level microdata from the American Community Survey (ACS) 2015-2019
#' public-use microdata sample (PUMS) data for Louisville, KY. This microdata sample
#' represents all adults (persons aged 18 or over) in Louisville, KY. \cr
#'
#' These data include replicate weights to use for variance estimation.
#'
#'
#' @format A data frame with 80 rows and 85 variables
#' \itemize{
#'   \item UNIQUE_ID: Unique identifier for records
#'   \item AGE: Age in years (copied from the AGEP variable in the ACS microdata)
#'   \item RACE_ETHNICITY: Race and Hispanic/Latino ethnicity
#'     derived from RAC1P and HISP variables
#'     of ACS microdata and collapsed to a smaller number of categories.
#'   \item SEX: Male or Female
#'   \item EDUC_ATTAINMENT: Highest level of education attained ('Less than high school' or 'High school or beyond')
#'     derived from SCHL variable in ACS microdata and collapsed to a smaller number of categories.
#'   \item PWGTP: Weights for the full-sample
#'   \item PWGTP1-PWGTP80: 80 columns of replicate weights
#'     created using the Successive Differences Replication (SDR) method.
#' }
#'
#' @keywords datasets
#' @name lou_pums_microdata
#' @usage data(lou_pums_microdata)
#' @examples
#' \dontrun{
#' data(lou_pums_microdata)
#'
#' # Prepare the data for analysis with the survey package
#'   library(survey)
#'
#'   lou_pums_rep_design <- survey::svrepdesign(
#'     data = lou_pums_microdata,
#'     variables = ~ UNIQUE_ID + AGE + SEX + RACE_ETHNICITY + EDUC_ATTAINMENT,
#'     weights = ~ PWGTP, repweights = "PWGTP\\d{1,2}",
#'     type = "successive-difference",
#'     mse = TRUE
#'   )
#'
#' # Estimate population proportions
#'   svymean(~ SEX, design = lou_pums_rep_design)
#' }
"lou_pums_microdata"

#' @title Control totals for the Louisville Vaccination Survey
#'
#' @description Control totals to use for raking or post-stratification
#' for the Louisville Vaccination Survey data. Control totals are population size estimates
#' from the ACS 2015-2019 5-year Public Use Microdata Sample (PUMS)
#' for specific demographic categories among adults in Jefferson County, KY. \cr
#'
#' These data were created using simulation.
#'
#' @format A nested list object with two lists, \code{poststratification} and \code{raking},
#' each of which contains two elements: \code{estimates} and \code{variance-covariance}.
#' \describe{
#'   \item{poststratification}{Control totals for the combination of
#'   \code{RACE_ETHNICITY}, \code{SEX}, and \code{EDUC_ATTAINMENT}.
#'     \itemize{
#'     \item estimates: A numeric vector of estimated population totals.
#'     \item variance-covariance: A variance-covariance matrix for the estimated population totals.
#'     }
#'   }
#'   \item{raking}{Separate control totals for each of
#'   \code{RACE_ETHNICITY}, \code{SEX}, and \code{EDUC_ATTAINMENT}.
#'     \itemize{
#'     \item estimates: A numeric vector of estimated population totals.
#'     \item variance-covariance: A variance-covariance matrix for the estimated population totals.
#'     }
#'   }
#' }
#'
#' @keywords datasets
#' @name lou_vax_survey_control_totals
#' @usage data(lou_vax_survey_control_totals)
"lou_vax_survey_control_totals"
