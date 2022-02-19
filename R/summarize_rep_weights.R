#' Summarize the replicate weights
#'
#' @description Summarize the replicate weights of a design
#'
#' @param rep_design A replicate design object, created with either the \code{survey} or \code{srvyr} packages.
#' @param type Default is \code{"both"}. Use \code{type = "overall"}, for an overall summary of the replicate weights.
#' Use \code{type = "specific"} for a summary of each column of replicate weights,
#' with each column of replicate weights summarized in a given row of the summary. \cr
#' \cr
#' Use \code{type = "both"} for a list containing both summaries,
#' with the list containing the names \code{"overall"} and \code{"both"}.
#'
#' @return If \code{type = "both"} (the default), the result is a list of data frames
#' with names \code{"overall"} and \code{"specific"}. If \code{type = "overall"}, the result is
#' a data frame providing an overall summary of the replicate weights. \cr
#' \cr
#' The contents of the \code{"overall"} summary are the following:
#' \itemize{
#'  \item{"nrows"}{: Number of rows for the weights}
#'  \item{"ncols"}{: Number of columns of replicate weights}
#'  \item{"degf_svy_pkg"}{: The degrees of freedom according to the survey package in R}
#'  \item{"rank"}{: The matrix rank as determined by a QR decomposition}
#'  \item{"avg_wgt_sum"}{: The average column sum}
#'  \item{"sd_wgt_sums"}{: The standard deviation of the column sums}
#'  \item{"min_rep_wgt"}{: The minimum value of any replicate weight}
#'  \item{"max_rep_wgt"}{: The maximum value of any replicate weight}
#' }
#'
#' If \code{type = "specific"}, the result is a data frame providing a
#' summary of each column of replicate weights, with each column of replicate weights
#' described in a given row of the data frame.
#' The contents of the \code{"specific"} summary are the following:
#' \itemize{
#'  \item{"Rep_Column"}{: The name of a given column of replicate weights.
#'   If columns are unnamed, the column number is used instead}
#'  \item{"N"}{: The number of entries}
#'  \item{"N_NONZERO"}{: The number of nonzero entries}
#'  \item{"SUM"}{: The sum of the weights}
#'  \item{"MEAN"}{: The average of the weights}
#'  \item{"CV"}{: The coefficient of variation of the weights (standard deviation divided by mean)}
#'  \item{"MIN"}{: The minimum weight}
#'  \item{"MAX"}{: The maximum weight}
#' }
#'
#' @export
#'
#' @examples
#'
#' # Load example data
#' suppressPackageStartupMessages(library(survey))
#' data(api)
#'
#' dclus1 <- svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
#' dclus1$variables$response_status <- sample(x = c("Respondent", "Nonrespondent",
#'                                                  "Ineligible", "Unknown eligibility"),
#'                                            size = nrow(dclus1),
#'                                            replace = TRUE)
#' rep_design <- as.svrepdesign(dclus1)
#'
#' # Adjust weights for cases with unknown eligibility
#' ue_adjusted_design <- redistribute_weights(
#'     design = rep_design,
#'     reduce_if = response_status %in% c("Unknown eligibility"),
#'     increase_if = !response_status %in% c("Unknown eligibility"),
#'     by = c("stype")
#' )
#'
#' # Summarize replicate weights
#'
#' summarize_rep_weights(rep_design, type = "both")
#'
#' # Compare replicate weights
#'
#' rep_wt_summaries <- lapply(list('original' = rep_design,
#'                                 'adjusted' = ue_adjusted_design),
#'                            summarize_rep_weights,
#'                            type = "overall")
#' print(rep_wt_summaries)
#'
summarize_rep_weights <- function(rep_design, type = 'both') {

  if (!'svyrep.design' %in% class(rep_design)) {
    stop("`rep_design` must be a replicate design object, with class 'svyrep.design'.")
  }

  if (is.null(type) || !is.character(type) || length(type) != 1) {
    stop("`type` must be one of 'overall', 'specific', or 'both'.")
  } else {
    if (!tolower(type) %in% c("both", "overall", "specific")) {
      stop("`type` must be one of 'overall', 'specific', or 'both'.")
    } else {
      type <- tolower(type)
    }
  }

  # Make sure replicate weights
  # are the "combined" weights
  # (i.e. incorporating base sampling weights)
  if (rep_design$combined.weights) {
    rep_wts_matrix <- as.matrix(rep_design$repweights)
  } else {
    rep_wts_matrix <- apply(as.matrix(rep_design$repweights),
                            MARGIN = 2, FUN = function(rep_wt_set) {
                              rep_wt_set * rep_design$pweights
                            })
  }

  # Summarize each set of replicate weights
  wt_set_summaries <- apply(rep_wts_matrix,
                            simplify = FALSE, MARGIN = 2,
                            FUN = function(wt_set) {

                              SD <- stats::sd(wt_set)
                              MEAN <- mean(wt_set)
                              CV <- SD/MEAN

                              c(N = nrow(rep_wts_matrix),
                                N_NONZERO = sum(wt_set != 0),
                                SUM = sum(wt_set),
                                MEAN = MEAN,
                                CV = CV,
                                MIN = min(wt_set),
                                MAX = max(wt_set))
                            })

  wt_set_summaries <- do.call(what = rbind,
                              args = wt_set_summaries)
  wt_set_summaries <- as.data.frame(wt_set_summaries)

  if (!is.null(colnames(rep_wts_matrix))) {
    wt_set_summaries <- cbind(data.frame("Rep_Column" = colnames(rep_wts_matrix)),
                              wt_set_summaries)
  } else {
    wt_set_summaries <- cbind(data.frame("Rep_Column" = seq_len(nrow(wt_set_summaries))),
                              wt_set_summaries)
  }
  rownames(wt_set_summaries) <- NULL

  # Produce an overall summary
  max_rep_wgt <- max(wt_set_summaries$MAX)
  min_rep_wgt <- min(wt_set_summaries$MIN)
  sd_wgt_sums <- stats::sd(wt_set_summaries$SUM)
  avg_wgt_sum <- mean(wt_set_summaries$SUM)

  matrix_rank <- qr(x = rep_wts_matrix, tol = 1e-05)[['rank']]

  rep_wt_dims <- setNames(dim(rep_wts_matrix),
                          nm = c("rows", "cols"))

  degf_svy_pkg <- survey::degf(rep_design)

  summary_df <- data.frame(
    'nrows' = rep_wt_dims['rows'],
    'ncols' = rep_wt_dims['cols'],
    'degf_svy_pkg' = degf_svy_pkg,
    'rank' = matrix_rank,
    'avg_wgt_sum' = avg_wgt_sum,
    'sd_wgt_sums' = sd_wgt_sums,
    'min_rep_wgt' = min_rep_wgt,
    'max_rep_wgt' = max_rep_wgt,
    row.names = NULL
  )

  # Combine into a list
  # to return as output
  if (type == "overall") {
    result <- summary_df
  }
  if (type == "specific") {
    result <- wt_set_summaries
  }
  if (type == "both") {
    result <- list(
      'overall' = summary_df,
      'specific' = wt_set_summaries
    )
  }
  return(result)

}
