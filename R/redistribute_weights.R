#' Shift weight from one set of cases to another, potentially separately for each value of a grouping variable.
#'
#' @param wt_set A numeric vector of weights
#' @param is_upweight_case A logical vector indicating cases whose weight should be increased
#' @param is_downweight_case A logical vector indicating cases whose weight should be decreased
#' @param by_values A vector of values for which adjustments should be conducted separately
#'
#' @return A numeric vector of adjusted weights, of the same length as \code{wt_set}.
#'
#' @examples
#'
#' example_data <- data.frame(weight = c(1.5, 1.25, 1.25, 2, 0.5, 1.5),
#'                            response_status = c(1, 1, 0, 0, 1, 1),
#'                            stratum = c('A', 'B', 'A', 'A', 'B', 'B'))
#'
#' shift_weight_by_group(wt_set = example_data$weight,
#'                       is_upweight_case = example_data$response_status == 1,
#'                       is_downweight_case = example_data$response_status == 0,
#'                       by_values = example_data$stratum)

shift_weight_by_group <- function(wt_set, is_upweight_case, is_downweight_case, by_values) {
  if (is.null(by_values)) {
    by_values <- rep(1, length(wt_set))
  }
  adjusted_wts_by_group <- by(data = data.frame(ordering = 1:length(wt_set),
                                                wts = wt_set,
                                                cases_to_upweight = is_upweight_case,
                                                cases_to_downweight = is_downweight_case),
                              INDICES = by_values,
                              simplify = FALSE, FUN = function(df) {

                                adj_factors <- rep(1, length(df$wts))
                                upweight_sum <- sum(df$wts[df$cases_to_upweight])
                                downweight_sum <- sum(df$wts[df$cases_to_downweight])

                                if (upweight_sum != 0) {
                                  upweight_factor <- 1 + (downweight_sum/upweight_sum)
                                } else {
                                  upweight_factor <- 1
                                }
                                downweight_factor <- 0

                                adj_factors[df$cases_to_upweight] <- upweight_factor
                                adj_factors[df$cases_to_downweight] <- downweight_factor

                                return(data.frame(ordering = df$ordering,
                                                  wts = df$wts * adj_factors))
                              })
  adjusted_wt_set <- do.call(rbind, adjusted_wts_by_group)
  adjusted_wt_set <- adjusted_wt_set[['wts']][order(adjusted_wt_set[['ordering']])]
  return(adjusted_wt_set)
}

#' Redistribute weight from one group to another
#'
#' @description Redistributes weight from one group to another: for example, from non-respondents to respondents.
#' Redistribution can be done separately by values of a categorical variable, for example to implement a nonresponse weighting class adjustment.
#' Redistribution is conducted for the full-sample weights as well as each set of replicate weights.
#'
#' @param design A survey design object, created with either the \code{survey} or \code{srvyr} packages.
#' @param cat_var A string (e.g. \code{"response_status"}) giving the name of a variable dividing the data
#'                into groups which should be upweighted or downweighted.
#' @param from_values The values of \code{cat_var} for cases which should have their weights increased.
#' @param to_values The values of \code{cat_var} for cases which should have their weights decreased.
#' @param by_var (Optional) The name of a variable used to group the redistribution of weights.
#'               For example, if the data include a variable named \code{"wt_class"}, one could specify \code{by_var = "wt_class"}.
#'
#' @return The survey design object, but with updated full-sample weights and updated replicate weights.
#' @export
#'
#' @examples
#'
#' # Load example data
#' suppressPackageStartupMessages(library(survey))
#' data(api)
#'
#' dclus1 <- svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
#' dclus1$variables$response_status <- sample(x = c("Respondent", "Nonrespondent", "Ineligible", "Unknown eligibility"),
#'                                            size = nrow(dclus1),
#'                                            replace = TRUE)
#' rep_design <- as.svrepdesign(dclus1)
#'
#' # Adjust weights for cases with unknown eligibility
#' ue_adjusted_design <- redistribute_weights(design = rep_design,
#'                                            cat_var = "response_status",
#'                                            from_values = c("Unknown eligibility"),
#'                                            to_values = c("Respondent", "Nonrespondent", "Ineligible"),
#'                                            by_var = "stype")
#'
#' # Adjust weights for nonresponse
#' nr_adjusted_design <- redistribute_weights(design = ue_adjusted_design,
#'                                            cat_var = "response_status",
#'                                            from_values = "Nonrespondent",
#'                                            to_values = "Respondent",
#'                                            by_var = "stype")
#
redistribute_weights <- function(design, cat_var,
                                 from_values = c("EN"), to_values = c("ER"),
                                 by_var = NULL) {
  UseMethod("redistribute_weights", design)
}

redistribute_weights.svyrep.design <- function(design, cat_var,
                                               from_values = c("EN"), to_values = c("ER"),
                                               by_var = NULL) {

  # Check validity of inputs
  if (is.null(from_values) || is.null(to_values)) {
    stop("Must supply values to the arguments `from_values` and `to_values`.")
  }

  cat_var_type <- typeof(design[[cat_var]])
  if (! typeof(from_values) %in% c("character", cat_var_type) | !typeof(to_values) %in% c("character", cat_var_type)) {
    stop(sprintf("The values supplied to `from_values` and `to_values` must be either character or the same type as `%s`.",
                 cat_var))
  }

  # Determine whether replicate weights in input are compressed for storage
  is_compressed <- "repweights_compressed" %in% class(design$repweights)

  # Extract the matrix of replicate weights
  rep_wts <- as.matrix(design[['repweights']])

  # Extract the full-sample weights
  full_sample_wts <- design[['pweights']]

  # Determine row indices of cases to upweight or downweight
  cases_to_upweight <- design[['variables']][[cat_var]] %in% from_values
  cases_to_downweight <- design[['variables']][[cat_var]] %in% to_values

  # Group data according to the by variable

  if (is.null(by_var)) {
    by_values <- rep(1, nrow(design))
  } else {
    by_values <- design[['variables']][[by_var]]
  }

  # Adjust the full-sample weights

  adjusted_full_sample_weights <- shift_weight_by_group(
    wt_set = full_sample_wts,
    is_upweight_case = cases_to_upweight,
    is_downweight_case = cases_to_downweight,
    by_values = by_values
  )


  # Adjust the replicate weights
  adjusted_rep_wts <- apply(X = rep_wts, MARGIN = 2, FUN = function(rep_wt_set) {
    shift_weight_by_group(wt_set = rep_wt_set,
                          is_upweight_case =  cases_to_upweight,
                          is_downweight_case = cases_to_downweight,
                          by_values = by_values)
  })


  # Update the survey design object to use the adjusted weights
  result <- design
  result[['pweights']] <- adjusted_full_sample_weights
  if (is_compressed) {
    result[['repweights']] <- survey::compressWeights(adjusted_rep_wts)
  } else {
    result[['repweights']][['weights']] <- adjusted_rep_wts
  }
  return(result)
}
