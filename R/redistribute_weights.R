#' Shift weight from one set of cases to another, potentially separately for each value of a grouping variable.
#'
#' @param wt_set A numeric vector of weights
#' @param is_upweight_case A logical vector indicating cases whose weight should be increased
#' @param is_downweight_case A logical vector indicating cases whose weight should be decreased
#'
#' @return A numeric vector of adjusted weights, of the same length as \code{wt_set}.
#'
#' @examples
#'
#' example_data <- data.frame(weight = c(1.5, 1.25, 1.25, 2, 0.5, 1.5),
#'                            response_status = c(1, 1, 0, 0, 1, 1),
#'                            stratum = c('A', 'B', 'A', 'A', 'B', 'B'))
#'
#' shift_weight(wt_set = example_data$weight,
#'              is_upweight_case = example_data$response_status == 1,
#'              is_downweight_case = example_data$response_status == 0)

shift_weight <- function(wt_set, is_upweight_case, is_downweight_case) {

                                adj_factors <- rep(1, length(wt_set))
                                upweight_sum <- sum(wt_set[is_upweight_case])
                                downweight_sum <- sum(wt_set[is_downweight_case])

                                if (upweight_sum != 0) {
                                  upweight_factor <- 1 + (downweight_sum/upweight_sum)
                                } else {
                                  upweight_factor <- 1
                                }
                                downweight_factor <- 0

                                adj_factors[is_upweight_case] <- upweight_factor
                                adj_factors[is_downweight_case] <- downweight_factor

  return(wt_set * adj_factors)
}

#' Redistribute weight from one group to another
#'
#' @description Redistributes weight from one group to another: for example, from non-respondents to respondents.
#' Redistribution is conducted for the full-sample weights as well as each set of replicate weights.
#' This can be done separately for each combination of a set of grouping variables, for example to implement a nonresponse weighting class adjustment.
#'
#' @param design A survey design object, created with either the \code{survey} or \code{srvyr} packages.
#' @param reduce_if An expression indicating which cases should have their weights set to zero.
#' Must evaluate to a logical vector with only values of TRUE or FALSE.
#' @param increase_if An expression indicating which cases should have their weights increased.
#' Must evaluate to a logical vector with only values of TRUE or FALSE.
#' @param by (Optional) A character vector with the names of variables used to group the redistribution of weights.
#' For example, if the data include variables named \code{"stratum"} and \code{"wt_class"}, one could specify \code{by_var = c("stratum", "wt_class")}.
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
#'                                            reduce_if = response_status %in% c("Unknown eligibility"),
#'                                            increase_if = !response_status %in% c("Unknown eligibility"),
#'                                            by = c("stype", "cname"))
#'
#' # Adjust weights for nonresponse
#' nr_adjusted_design <- redistribute_weights(design = ue_adjusted_design,
#'                                            reduce_if = response_status %in% c("Nonrespondent"),
#'                                            increase_if = response_status == "Respondent",
#'                                            by = c("stype", "cname"))
#
redistribute_weights <- function(design,
                                 reduce_if, increase_if,
                                 by) {
  if (!'svyrep.design' %in% class(design)) {
    stop("`design` must be a replicate design object.")
  }
  UseMethod("redistribute_weights", design)
}

redistribute_weights.svyrep.design <- function(design, reduce_if, increase_if, by) {

  # Check validity of inputs

  if (missing(by) || is.null(by)) {
    by <- NULL
  } else {
    if (!is.character(by)) {
      stop("`by` must be a character vector with at least one column name from `design`.")
    }
    if (!all(by %in% colnames(design[['variables']]))) {
      missing_variables <- setdiff(by, colnames(design[['variables']]))
      error_msg <- sprintf("The following `by` variables are missing from the data: %s",
                           paste(missing_variables, collapse = ", "))
      stop(error_msg)
    }
  }

  if (missing(reduce_if) || missing(increase_if)) {
    stop("Must supply expressions to the arguments `reduce_if` and `increase_if`.")
  }

  case_groupings <- dplyr::transmute(design[['variables']],
                                     `_IS_DOWNWT_CASE_` = {{ reduce_if }},
                                     `_IS_UPWT_CASE_` = {{ increase_if }})

  if (!is.logical(case_groupings[["_IS_DOWNWT_CASE_"]]) || !is.logical(case_groupings[['_IS_UPWT_CASE_']])) {
    stop("The expressions supplied to `reduce_if` and `increase_if` must result in logical values of TRUE or FALSE.")
  }
  if (any(is.na(case_groupings[["_IS_DOWNWT_CASE_"]])) || any(is.na(case_groupings[["_IS_UPWT_CASE_"]]))) {
    stop("The result of the expressions supplied to `reduce_if` and `increase_if` must be TRUE or FALSE, not NA.")
  }

  # Determine whether replicate weights in input are compressed for storage
  is_compressed <- "repweights_compressed" %in% class(design$repweights)

  # Extract the matrix of replicate weights
  rep_wts <- as.matrix(design[['repweights']])


  # Divide data into groups, and obtain list of row indices for each group

  if (!is.null(by)) {
    grouping_df <- dplyr::group_by(design[['variables']],
                                   across(one_of({{by}})))
    grouping_df <- attributes(grouping_df)[['groups']]
    group_row_indices <- grouping_df[['.rows']]
  } else {
    group_row_indices <- list(seq_len(nrow(design)))
  }


  # Extract the full-sample weights
  full_sample_wts <- design[['pweights']]

  # Adjust the full-sample weights

  adjusted_full_sample_weights <- full_sample_wts

  for (set_of_indices in group_row_indices) {
    adjusted_full_sample_weights[set_of_indices] <- shift_weight(
      wt_set = adjusted_full_sample_weights[set_of_indices],
      is_upweight_case = case_groupings[['_IS_UPWT_CASE_']][set_of_indices],
      is_downweight_case = case_groupings[['_IS_DOWNWT_CASE_']][set_of_indices]
    )
  }


  # Adjust the replicate weights
  adjusted_rep_wts <- apply(X = rep_wts, MARGIN = 2, FUN = function(rep_wt_set) {

    adjusted_wt_set <- rep_wt_set

    for (set_of_indices in group_row_indices) {
      adjusted_wt_set[set_of_indices] <- shift_weight(
        wt_set = rep_wt_set[set_of_indices],
        is_upweight_case = case_groupings[['_IS_UPWT_CASE_']][set_of_indices],
        is_downweight_case = case_groupings[['_IS_DOWNWT_CASE_']][set_of_indices]
      )
    }
    return(adjusted_wt_set)
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
