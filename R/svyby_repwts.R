#' Compare survey statistics calculated separately from different sets of replicate weights
#'
#' @description A modified version of the \code{svyby()} function from the \code{survey} package.
#' Whereas \code{svyby()} calculates statistics separately for each subset formed by a specified grouping variable,
#' \code{svyby_repwts()} calculates statistics separately for each replicate design, in addition to any additional user-specified grouping variables.
#'
#' @param rep_designs The replicate-weights survey designs to be compared. Supplied either as:
#' \itemize{
#'   \item A named list of replicate-weights survey design objects, for example \code{list('nr' = nr_adjusted_design, 'ue' = ue_adjusted_design)}.
#'   \item A 'stacked' replicate-weights survey design object created by \code{stack_replicate_designs()}.
#' }
#' The designs must all have the same number of columns of replicate weights, of the same type (bootstrap, JKn, etc.)
#' @param formula A formula specifying the variables to pass to \code{FUN}
#' @param by A formula specifying factors that define subsets
#' @param FUN A function taking a formula and survey design object as its first two arguments.
#' Usually a function from the \code{survey} package, such as \code{svytotal} or \code{svymean}.
#' @param ... Other arguments to \code{FUN}
#' @param deff A value of \code{TRUE} or \code{FALSE}, indicating whether design effects should be estimated if possible.
#' @param keep.var A value of \code{TRUE} or \code{FALSE}. If \code{FUN} returns a \code{svystat} object, indicates whether to extract standard errors from it.
#' @param keep.names Define row names based on the subsets
#' @param verbose If \code{TRUE}, print a label for each subset as it is processed.
#' @param vartype Report variability as one or more of standard error, confidence interval, coefficient of variation,  percent coefficient of variation, or variance
#' @param drop.empty.groups If \code{FALSE}, report \code{NA} for empty groups, if \code{TRUE} drop them from the output
#' @param return.replicates If \code{TRUE}, return all the replicates as the "replicates" attribute of the result.
#' This can be useful if you want to produce custom summaries of the estimates from each replicate.
#' @param na.rm.by If true, omit groups defined by \code{NA} values of the \code{by} variables
#' @param na.rm.all If true, check for groups with no non-missing observations for variables defined by \code{formula} and treat these groups as empty
#' @param multicore Use \code{multicore} package to distribute subsets over multiple processors?
#' @return An object of class \code{"svyby"}: a data frame showing the grouping factors and results of \code{FUN} for each combination of the grouping factors.
#' The first grouping factor always consists of indicators for which replicate design was used for an estimate.
#' @export
#'
#' @examples
#' \dontrun{
#' suppressPackageStartupMessages(library(survey))
#' data(api)
#'
#' dclus1 <- svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
#' dclus1$variables$response_status <- sample(x = c("Respondent", "Nonrespondent",
#'                                                  "Ineligible", "Unknown eligibility"),
#'                                            size = nrow(dclus1),
#'                                            replace = TRUE)
#' orig_rep_design <- as.svrepdesign(dclus1)
#'
#' # Adjust weights for cases with unknown eligibility
#' ue_adjusted_design <- redistribute_weights(
#'     design = orig_rep_design,
#'     reduce_if = response_status %in% c("Unknown eligibility"),
#'     increase_if = !response_status %in% c("Unknown eligibility"),
#'     by = c("stype")
#' )
#'
#' # Adjust weights for nonresponse
#' nr_adjusted_design <- redistribute_weights(
#'     design = ue_adjusted_design,
#'     reduce_if = response_status %in% c("Nonrespondent"),
#'     increase_if = response_status == "Respondent",
#'     by = c("stype")
#' )
#'
#' # Compare estimates from the three sets of replicate weights
#'
#'   list_of_designs <- list('original' = orig_rep_design,
#'                           'unknown eligibility adjusted' = ue_adjusted_design,
#'                           'nonresponse adjusted' = nr_adjusted_design)
#'
#'   ##_ First compare overall means for two variables
#'   means_by_design <- svyby_repwts(formula = ~ api00 + api99,
#'                                   FUN = svymean,
#'                                   rep_design = list_of_designs)
#'
#'   print(means_by_design)
#'
#'   ##_ Next compare domain means for two variables
#'   domain_means_by_design <- svyby_repwts(formula = ~ api00 + api99,
#'                                          by = ~ stype,
#'                                          FUN = svymean,
#'                                          rep_design = list_of_designs)
#'
#'   print(domain_means_by_design)
#'
#' # Calculate confidence interval for difference between estimates
#'
#' ests_by_design <- svyby_repwts(rep_designs = list('NR-adjusted' = nr_adjusted_design,
#'                                                   'Original' = orig_rep_design),
#'                                FUN = svymean, formula = ~ api00 + api99)
#'
#' differences_in_estimates <- svycontrast(stat = ests_by_design, contrasts = list(
#'   'Mean of api00: NR-adjusted vs. Original' = c(1,-1,0,0),
#'   'Mean of api99: NR-adjusted vs. Original' = c(0,0,1,-1)
#' ))
#'
#' print(differences_in_estimates)
#'
#' confint(differences_in_estimates, level = 0.95)
#' }
svyby_repwts <- function(rep_designs,
                         formula, by, FUN, ...,
                         deff = FALSE,
                         keep.var = TRUE,
                         keep.names = TRUE,
                         verbose = FALSE,
                         vartype = c("se","ci","ci","cv","cvpct","var"),
                         drop.empty.groups = TRUE,
                         return.replicates = FALSE, na.rm.by=FALSE, na.rm.all=FALSE,
                         multicore = getOption("survey.multicore")) {


  if (!'svyrep.stacked' %in% class(rep_designs)) {
    if (!is.list(rep_designs) || !all(sapply(rep_designs, function(x) 'svyrep.design' %in% class(x)))) {
      stop("`design_list` must be either a list of replicate-weights survey designs, or the output of `stack_replicate_designs()`.")
    }
  }

  if ('svyrep.stacked' %in% class(rep_designs)) {
    variable_for_source_design <- attributes(rep_designs)[['variable_for_source_design']]
    stacked_design <- rep_designs
  }
  if (!'svyrep.stacked' %in% class(rep_designs)) {
    stacked_design <- stack_replicate_designs(rep_designs, .id = "Design_Name")
    variable_for_source_design <- "Design_Name"
  }

  if (missing(by) || is.null(by)) {
    by_input <- as.formula(sprintf("~ %s", variable_for_source_design))
  } else if ('formula' %in% class(by)) {
    by_input <- stats::update(by, as.formula(sprintf("~ %s + .", variable_for_source_design)))
  } else if (is.list(by)) {
    by_input <- append(x = by,
                       values = setNames(list(stacked_design[['variables']][[variable_for_source_design]]),
                                         variable_for_source_design),
                       after = 0)
  }

  if (missing(vartype)) {
    vartype <- "se"
  }

  svyby_result <- survey::svyby(formula = formula, by = by_input,
                                design = stacked_design,
                                FUN = FUN,
                                ... = ...,
                                deff = deff,
                                keep.var = keep.var,
                                keep.names = keep.names,
                                verbose = verbose,
                                vartype = vartype,
                                drop.empty.groups = drop.empty.groups, covmat = TRUE,
                                return.replicates = return.replicates, na.rm.by=na.rm.by, na.rm.all=na.rm.all,
                                multicore = getOption("survey.multicore"))

  value_for_statistic_attribute <- deparse(substitute(FUN))
  attr(svyby_result, 'svyby')$statistic <- value_for_statistic_attribute

  return(svyby_result)
}

