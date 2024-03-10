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
#' @param by (Optional) A character vector with the names of variables used to group the summaries.
#'
#' @return If \code{type = "both"} (the default), the result is a list of data frames
#' with names \code{"overall"} and \code{"specific"}. If \code{type = "overall"}, the result is
#' a data frame providing an overall summary of the replicate weights. \cr
#' \cr
#' The contents of the \code{"overall"} summary are the following:
#' \itemize{
#'  \item "nrows": Number of rows for the weights
#'  \item "ncols": Number of columns of replicate weights
#'  \item "degf_svy_pkg": The degrees of freedom according to the survey package in R
#'  \item "rank": The matrix rank as determined by a QR decomposition
#'  \item "avg_wgt_sum": The average column sum
#'  \item "sd_wgt_sums": The standard deviation of the column sums
#'  \item "min_rep_wgt": The minimum value of any replicate weight
#'  \item "max_rep_wgt": The maximum value of any replicate weight
#' }
#'
#' If \code{type = "specific"}, the result is a data frame providing a
#' summary of each column of replicate weights, with each column of replicate weights
#' described in a given row of the data frame.
#' The contents of the \code{"specific"} summary are the following:
#' \itemize{
#'  \item "Rep_Column": The name of a given column of replicate weights.
#'    If columns are unnamed, the column number is used instead
#'  \item "N": The number of entries
#'  \item "N_NONZERO": The number of nonzero entries
#'  \item "SUM": The sum of the weights
#'  \item "MEAN": The average of the weights
#'  \item "CV": The coefficient of variation of the weights (standard deviation divided by mean)
#'  \item "MIN": The minimum weight
#'  \item "MAX": The maximum weight
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
#' # Summarize replicate weights by grouping variables
#'
#' summarize_rep_weights(ue_adjusted_design, type = 'overall',
#'                       by = c("response_status"))
#'
#' summarize_rep_weights(ue_adjusted_design, type = 'overall',
#'                       by = c("stype", "response_status"))
#'
#' # Compare replicate weights
#'
#' rep_wt_summaries <- lapply(list('original' = rep_design,
#'                                 'adjusted' = ue_adjusted_design),
#'                            summarize_rep_weights,
#'                            type = "overall")
#' print(rep_wt_summaries)
#'
summarize_rep_weights <- function(rep_design, type = 'both', by) {

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

  if (missing(by) || is.null(by)) {
    by <- NULL
  } else {
    if (!is.character(by)) {
      stop("`by` must be a character vector with at least one column name from `design`.")
    }
    if (!all(by %in% colnames(rep_design[['variables']]))) {
      missing_variables <- setdiff(by, colnames(rep_design[['variables']]))
      error_msg <- sprintf("The following `by` variables are missing from the data: %s",
                           paste(missing_variables, collapse = ", "))
      stop(error_msg)
    }
  }

  # Make sure replicate weights
  # are the "combined" weights
  # (i.e. incorporating base sampling weights)
  rep_wts_matrix <- weights(object = rep_design, type = 'analysis')

  # If summaries should be grouped, divide data into groups, and obtain list of row indices for each group

  if (!is.null(by)) {
    is_grouped <- TRUE
    ## Create factor whose levels are observed combinations of grouping variables
    unique_values_factor <- interaction(rep_design$variables[,by],
                                        sep = "_<>_", drop = TRUE)
    ## Obtain a dataframe giving the unique values of the combinations of grouping variables
    first_index_per_unique_combo <- !duplicated(unique_values_factor)
    groups_df <- rep_design$variables[first_index_per_unique_combo, by, drop=FALSE]
    ## Reorder rows of the dataframe to match the order of the levels of the unique values factor
    desired_row_order <- unique_values_factor[first_index_per_unique_combo] |>
      as.numeric() |> order()
    groups_df <- groups_df[desired_row_order, , drop = FALSE]

    ## Get row indices for each group
    group_row_indices <- lapply(levels(unique_values_factor),
                                function(factor_level) which(unique_values_factor == factor_level))
    names(group_row_indices) <- levels(unique_values_factor)
  } else {
    is_grouped <- FALSE
    group_row_indices <- list(seq_len(nrow(rep_design)))
  }

  if (nrow(rep_design) == 0) {
    group_row_indices <- lapply(group_row_indices, FUN = function(x) 0L)
  }

  # Special handling if the replicate design has zero rows of data
  if (nrow(rep_design) == 0) {

    rep_wt_dims <- setNames(dim(rep_wts_matrix),
                            nm = c("rows", "cols"))

    ## Prepare overall summary
    overall_summary <- data.frame(
      'nrows' = rep_wt_dims['rows'],
      'ncols' = rep_wt_dims['cols'],
      'degf_svy_pkg' = NA_integer_,
      'rank' = NA_integer_,
      'avg_wgt_sum' = NA_real_,
      'sd_wgt_sums' = NA_real_,
      'min_rep_wgt' = NA_real_,
      'max_rep_wgt' = NA_real_,
      row.names = NULL
    )

    ## Prepare replicate-specific summaries
    wt_set_summary <- c(N = 0,
                        N_NONZERO = 0,
                        SUM = NA_real_,
                        MEAN = NA_real_,
                        CV = NA_real_,
                        MIN = NA_real_,
                        MAX = NA_real_)

    wt_set_summaries <- lapply(X = seq_len(rep_wt_dims['cols']), FUN = function(rep_column) {
      wt_set_summary
    }) |> do.call(what = rbind)

    if (!is.null(colnames(rep_wts_matrix))) {
      wt_set_summaries <- cbind(data.frame("Rep_Column" = colnames(rep_wts_matrix)),
                                wt_set_summaries,
                                row.names = NULL)
    } else {
      wt_set_summaries <- cbind(data.frame("Rep_Column" = seq_len(nrow(wt_set_summaries))),
                                wt_set_summaries,
                                row.names = NULL)
    }

    if (type == "overall") {
      result <- overall_summary
    }
    if (type == "specific") {
      result <- wt_set_summaries
    }
    if (type == "both") {
      result <- list(
        'overall' = overall_summary,
        'specific' = wt_set_summaries
      )
    }
    return(result)
  }

  # Summarize each set of replicate weights
  wt_set_summaries_by_group <- lapply(
    X = group_row_indices, FUN = function(row_indices) {

      rep_wts_matrix <- rep_wts_matrix[row_indices,,drop=FALSE]

      wt_set_summaries <- apply(rep_wts_matrix,
                                simplify = FALSE, MARGIN = 2,
                                FUN = function(wt_set) {

                                  SD <- stats::sd(wt_set)
                                  MEAN <- mean(wt_set)
                                  CV <- SD/abs(MEAN)

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
                                  wt_set_summaries,
                                  row.names = NULL)
      } else {
        wt_set_summaries <- cbind(data.frame("Rep_Column" = seq_len(nrow(wt_set_summaries))),
                                  wt_set_summaries,
                                  row.names = NULL)
      }
      rownames(wt_set_summaries) <- NULL

      return(wt_set_summaries)
    }
  )

  # Produce an overall summary
  if (type %in% c("overall", "both")) {

    overall_summaries_by_group <- lapply(
      X = seq_along(group_row_indices), FUN = function(group_index) {

        wt_set_summaries <- wt_set_summaries_by_group[[group_index]]
        rep_wts_matrix <- rep_wts_matrix[group_row_indices[[group_index]],,drop=FALSE]

        max_rep_wgt <- max(wt_set_summaries$MAX)
        min_rep_wgt <- min(wt_set_summaries$MIN)
        sd_wgt_sums <- stats::sd(wt_set_summaries$SUM)
        avg_wgt_sum <- mean(wt_set_summaries$SUM)

        matrix_rank <- qr(x = rep_wts_matrix, tol = 1e-05)[['rank']]

        rep_wt_dims <- setNames(dim(rep_wts_matrix),
                                nm = c("rows", "cols"))

        degf_svy_pkg <- survey::degf(rep_design[group_row_indices[[group_index]],])

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
        return(summary_df)
      }
    )
  }

  # Combine grouped results into single dataframe, if applicable

  if (is_grouped) {
    wt_set_summaries <- lapply(X = seq_along(group_row_indices), FUN = function(group_index) {
      cbind(groups_df[group_index,,drop=FALSE],
            wt_set_summaries_by_group[[group_index]],
            row.names = NULL)
    }) |> do.call(what = rbind)

    if (type %in% c("overall", "both")) {
      overall_summary <- lapply(X = seq_along(group_row_indices), FUN = function(group_index) {
        cbind(groups_df[group_index,,drop=FALSE],
              overall_summaries_by_group[[group_index]],
              row.names = NULL)
      }) |> do.call(what = rbind)
    }
  }

  if (!is_grouped) {
    wt_set_summaries <- wt_set_summaries_by_group[[1]]

    if (type %in% c("overall", "both")) {
      overall_summary <- overall_summaries_by_group[[1]]
    }
  }

  # Combine into a list
  # to return as output
  if (type == "overall") {
    result <- overall_summary
  }
  if (type == "specific") {
    result <- wt_set_summaries
  }
  if (type == "both") {
    result <- list(
      'overall' = overall_summary,
      'specific' = wt_set_summaries
    )
  }
  return(result)

}
