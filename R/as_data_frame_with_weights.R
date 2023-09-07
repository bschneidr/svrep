#' Convert a survey design object to a data frame with weights stored as columns
#'
#' @param design A survey design object, created with either the \code{survey} or \code{srvyr} packages.
#' @param full_wgt_name The column name to use for the full-sample weights
#' @param rep_wgt_prefix For replicate design objects, a prefix to use for the column names
#' of the replicate weights. The column names will be created by appending
#' the replicate number after the prefix.
#' @param vars_to_keep By default, all variables in the data will be kept.
#' To select only a subset of the non-weight variables,
#' you can supply a character vector of variable names to keep.
#'
#' @return A data frame, with new columns containing the weights from the survey design object
#' @export
#'
#' @examples
#'
#' data("lou_vax_survey", package = 'svrep')
#'
#' library(survey)
#' # Create a survey design object
#' survey_design <- svydesign(data = lou_vax_survey,
#'                            weights = ~ SAMPLING_WEIGHT,
#'                            ids = ~ 1)
#'
#' rep_survey_design <- as.svrepdesign(survey_design,
#'                                     type = "boot",
#'                                     replicates = 10)
#'
#' # Adjust the weights for nonresponse
#' nr_adjusted_design <- redistribute_weights(
#'   design = rep_survey_design,
#'   reduce_if = RESPONSE_STATUS == "Nonrespondent",
#'   increase_if = RESPONSE_STATUS == "Respondent",
#'   by = c("RACE_ETHNICITY", "EDUC_ATTAINMENT")
#' )
#'
#' # Save the survey design object as a data frame
#' nr_adjusted_data <- as_data_frame_with_weights(
#'   nr_adjusted_design,
#'   full_wgt_name = "NR_ADJUSTED_WGT",
#'   rep_wgt_prefix = "NR_ADJUSTED_REP_WGT_"
#' )
#' head(nr_adjusted_data)
#'
#' # Check the column names of the result
#' colnames(nr_adjusted_data)
#'
as_data_frame_with_weights <- function(design, full_wgt_name = "FULL_SAMPLE_WGT",
                                       rep_wgt_prefix = "REP_WGT_",
                                       vars_to_keep = NULL) {

  # Extract weights from the design object
  if (inherits(design, 'svyrep.design')) {
    replicate_weights <- stats::weights(design, type = "analysis")
    n_replicates <- ncol(replicate_weights)
  }

  sampling_weights <- stats::weights(design, type = "sampling")

  # Extract variables from the design
  if (inherits(design, "DBIsvydesign")) {
    if (is.null(vars_to_keep)) {
      vars_to_keep <- dimnames(design)[[2]]
    }
    survey_data <- getvars(
      formula = vars_to_keep,
      dbconnection = design$db$connection,
      tables = design$db$tablename,
      updates = design$updates,
      subset = design$subset
    )
  }

  if (!inherits(design, "DBIsvydesign")) {
    if (!is.null(vars_to_keep)) {
      survey_data <- design[['variables']][, vars_to_keep, drop = FALSE]
    } else {
      survey_data <- design[['variables']]
    }
  }

  # Add full-sample weights to the data frame

  if (missing(full_wgt_name) || is.null(full_wgt_name)) {
    message("Column name for full-sample weights not supplied. Using default of 'FULL_SAMPLE_WGT'.")
    full_wgt_name <- 'FULL_SAMPLE_WGT'
  }
  if (length(full_wgt_name) != 1 || !is.character(full_wgt_name) || is.na(full_wgt_name)) {
    error_msg <- sprintf("If supplied, `full_wgt_name` must be a single character string.",
                         full_wgt_name)
    stop(error_msg)
  }
  if (full_wgt_name %in% colnames(survey_data)) {
    warning_message <- sprintf("There is already a column named `%s` in the survey data. It will be replaced.",
                               full_wgt_name)
    warning(warning_message)
  }
  survey_data[[full_wgt_name]] <- sampling_weights

  # Add replicate weights to the data frame

  if (inherits(design, 'svyrep.design')) {
    if (missing(rep_wgt_prefix) || is.null(rep_wgt_prefix)) {
      message("Prefix for replicate weight column names not supplied. Using default of 'REP_WGT_'.")
      rep_wgt_prefix <- 'REP_WGT_'
    }
    if (length(rep_wgt_prefix) != 1 || !is.character(rep_wgt_prefix) || is.na(rep_wgt_prefix)) {
      error_msg <- sprintf("If supplied, `rep_wgt_prefix` must be a single character string.",
                           rep_wgt_prefix)
      stop(error_msg)
    }
    rep_col_names <- paste0(rep_wgt_prefix, seq_len(n_replicates))

    if (any(rep_col_names %in% colnames(survey_data))) {
      col_names_string <- paste(sprintf("'%s'", intersect(rep_col_names, colnames(survey_data))),
                                collapse = " ")
      error_msg <- paste0("The survey data already contain columns named ",
                          col_names_string)
      stop(error_msg)
    }

    replicate_weights <- as.data.frame(replicate_weights)
    colnames(replicate_weights) <- rep_col_names
    survey_data <- cbind(survey_data, replicate_weights)
  }

  # Drop rownames
  rownames(survey_data) <- NULL

  # Return the result
  return(survey_data)
}
