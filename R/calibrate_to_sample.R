#' @title Calibrate weights from a primary survey to estimated totals from a control survey,
#' with replicate-weight adjustments that account for variance of the control totals
#' @description Calibrate the weights of a primary survey to match estimated totals from a control survey,
#' using adjustments to the replicate weights to account for the variance of the estimated control totals.
#' The adjustments to replicate weights are conducted using the method proposed by Opsomer and Erciulescu (2021).
#' This method can be used to implement general calibration as well as post-stratification or raking specifically
#' (see the details for the \code{calfun} parameter).
#' @details With the Opsomer-Erciulescu method, each column of replicate weights from the control survey
#' is randomly matched to a column of replicate weights from the primary survey,
#' and then the column from the primary survey is calibrated to control totals estimated by
#' perturbing the control sample's full-sample estimates using the estimates from the
#' matched column of replicate weights from the control survey.
#' \cr \cr
#' If there are fewer columns of replicate weights in the control survey than in the primary survey,
#' then not all primary replicate columns will be matched to a replicate column from the control survey. \cr
#'
#' If there are more columns of replicate weights in the control survey than in the primary survey,
#' then the columns of replicate weights in the primary survey will be duplicated \code{k} times, where \code{k} is the smallest
#' positive integer such that the resulting number of columns of replicate weights for the primary survey is greater than or equal
#' to the number of columns of replicate weights in the control survey. \cr
#'
#' Because replicate columns of the control survey are matched \emph{at random} to primary survey replicate columns,
#' there are multiple ways to ensure that this matching is reproducible.
#' The user can either call \link[base]{set.seed} before using the function,
#' or supply a mapping to the argument \code{control_col_matches}.
#'
#' @param primary_rep_design A replicate design object for the primary survey, created with either the \code{survey} or \code{srvyr} packages.
#' @param control_rep_design A replicate design object for the control survey.
#' @param cal_formula A formula listing the variables to use for calibration.
#' All of these variables must be included in both \code{primary_rep_design} and \code{control_rep_design}.
#' @param calfun A calibration function from the \code{survey} package,
#' such as \link[survey]{cal.linear}, \link[survey]{cal.raking}, or \link[survey]{cal.logit}.
#' Use \code{cal.linear} for ordinary post-stratification, and \code{cal.raking} for raking.
#' See \link[survey]{calibrate} for additional details.
#' @param bounds Parameter passed to \link[survey]{grake} for calibration. See \link[survey]{calibrate} for details.
#' @param verbose Parameter passed to \link[survey]{grake} for calibration. See \link[survey]{calibrate} for details.
#' @param maxit Parameter passed to \link[survey]{grake} for calibration. See \link[survey]{calibrate} for details.
#' @param epsilon Parameter passed to \link[survey]{grake} for calibration. \cr
#' After calibration, the absolute difference between each calibration target and the calibrated estimate
#' will be no larger than \code{epsilon} times (1 plus the absolute value of the target).
#' See \link[survey]{calibrate} for details.
#' @param variance Parameter passed to \link[survey]{grake} for calibration. See \link[survey]{calibrate} for details.
#' @param control_col_matches Optional parameter to specify which control survey replicate
#' is matched to each primary survey replicate. If the \eqn{i-th} entry of \code{control_col_matches}
#' equals \eqn{k}, then replicate \eqn{i} in \code{primary_rep_design} is matched
#' to replicate \eqn{k} in \code{control_rep_design.}
#' Entries of \code{NA} denote a primary survey replicate not matched to any control survey replicate.
#' If this parameter is not used, matching is done at random.
#' @return A replicate design object, with full-sample weights calibrated to totals from \code{control_rep_design},
#' and replicate weights adjusted to account for variance of the control totals.
#' If \code{primary_rep_design} had fewer columns of replicate weights than \code{control_rep_design},
#' then the number of replicate columns and the length of \code{rscales} will be increased by a multiple \code{k},
#' and the \code{scale} will be updated by dividing by \code{k}. \cr \cr
#' The element \code{control_column_matches} indicates, for each replicate column of the calibrated primary survey,
#' which column of replicate weights it was matched to from the control survey.
#' Columns which were not matched to control survey replicate column are indicated by \code{NA}. \cr \cr
#' The element \code{degf} will be set to match that of the primary survey
#' to ensure that the degrees of freedom are not erroneously inflated by
#' potential increases in the number of columns of replicate weights.
#' @section Syntax for Common Types of Calibration:
#' For ratio estimation with an auxiliary variable \code{X},
#' use the following options: \cr
#'   - \code{cal_formula = ~ -1 + X} \cr
#'   - \code{variance = 1}, \cr
#'   - \code{cal.fun = survey::cal.linear}
#'
#' For post-stratification, use the following option:
#'
#'   - \code{cal.fun = survey::cal.linear}
#'
#' For raking, use the following option:
#'
#'   - \code{cal.fun = survey::cal.raking}
#' @references
#' Opsomer, J.D. and A. Erciulescu (2021).
#' "Replication variance estimation after sample-based calibration."
#' \strong{Survey Methodology}, \emph{47}: 265-277.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Load example data for primary survey ----
#'
#'   suppressPackageStartupMessages(library(survey))
#'   data(api)
#'
#'   primary_survey <- svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc) |>
#'     as.svrepdesign(type = "JK1")
#'
#' # Load example data for control survey ----
#'
#'   control_survey <- svydesign(id = ~ 1, fpc = ~fpc, data = apisrs) |>
#'     as.svrepdesign(type = "JK1")
#'
#' # Calibrate totals for one categorical variable and one numeric ----
#'
#'   calibrated_rep_design <- calibrate_to_sample(
#'     primary_rep_design = primary_survey,
#'     control_rep_design = control_survey,
#'     cal_formula = ~ stype + enroll,
#'   )
#'
#' # Inspect estimates before and after calibration ----
#'
#'   ##_ For the calibration variables, estimates and standard errors
#'   ##_ from calibrated design will match those of the control survey
#'
#'     svytotal(x = ~ stype + enroll, design = primary_survey)
#'     svytotal(x = ~ stype + enroll, design = control_survey)
#'     svytotal(x = ~ stype + enroll, design = calibrated_rep_design)
#'
#'   ##_ Estimates from other variables will be changed as well
#'
#'     svymean(x = ~ api00 + api99, design = primary_survey)
#'     svymean(x = ~ api00 + api99, design = control_survey)
#'     svymean(x = ~ api00 + api99, design = calibrated_rep_design)
#'
#' # Inspect weights before and after calibration ----
#'
#'   summarize_rep_weights(primary_survey, type = 'overall')
#'   summarize_rep_weights(calibrated_rep_design, type = 'overall')
#'
#' # For reproducibility, specify how to match replicates between surveys ----
#'
#'   column_matching <- calibrated_rep_design$control_col_matches
#'   print(column_matching)
#'
#'   calibrated_rep_design <- calibrate_to_sample(
#'     primary_rep_design = primary_survey,
#'     control_rep_design = control_survey,
#'     cal_formula = ~ stype + enroll,
#'     control_col_matches = column_matching
#'   )
#' }

calibrate_to_sample <- function(primary_rep_design, control_rep_design,
                                cal_formula,
                                calfun = survey::cal.linear,
                                bounds = list('lower' = -Inf, 'upper' = Inf),
                                verbose = FALSE, maxit = 50,
                                epsilon = 1e-7, variance = NULL,
                                control_col_matches = NULL) {

  if (!inherits(primary_rep_design, "svyrep.design")) {
    stop("`primary_rep_design` must be a replicate survey design object, with class `svyrep.design`")
  }
  if (!inherits(control_rep_design, "svyrep.design")) {
    stop("`control_rep_design` must be a replicate survey design object, with class `svyrep.design`")
  }

  is_tbl_svy <- inherits(primary_rep_design, 'tbl_svy')

  # For database-backed design, obtain calibration variables ----
  is_db_backed_design <- inherits(primary_rep_design, 'DBIrepdesign')

  if (is_db_backed_design) {
    primary_rep_design$variables <- getvars(
      formula = cal_formula,
      dbconnection = primary_rep_design$db$connection,
      tables = primary_rep_design$db$tablename,
      updates = primary_rep_design$updates,
      subset = primary_rep_design$subset
    )
    db_info <- primary_rep_design$db
  } else {
    db_info <- NULL
  }

  # Determine parameters describing replicate designs ----
  R_control <- ncol(control_rep_design$repweights)
  R_primary <- ncol(primary_rep_design$repweights)

  A_primary <- primary_rep_design$scale
  A_control <- control_rep_design$scale

  rscales_primary <- primary_rep_design$rscales
  rscales_control <- control_rep_design$rscales

  # If the number of control replicates exceeds number of primary replicates, duplicate primary replicates ----

  R_ratio <- R_control/R_primary

  if (R_ratio > 1) {

    k <- ceiling(R_ratio)

    col_indices_duped <- rep(seq_len(R_primary), each = k)
    A_primary <- A_primary / k

    R_primary <- k * R_primary
    rscales_primary <- rep(primary_rep_design$rscales, each = k)


    primary_rep_design <- survey::svrepdesign(
      data = primary_rep_design$variables,
      repweights = weights(primary_rep_design, type = 'replication')[,col_indices_duped],
      weights = weights(primary_rep_design, type = 'sampling'),
      type = primary_rep_design$type,
      combined.weights = primary_rep_design$combined.weights,
      rho = if (primary_rep_design$type %in% c("JK1", "JKn")) NULL else primary_rep_design$rho,
      scale = A_primary,
      rscales = rscales_primary,
      fpc = primary_rep_design$fpc,
      fpctype = primary_rep_design$fpctype
    )

    duplication_message <- paste(
      "The primary survey has fewer replicates than the control survey,",
      "so columns in the primary survey will be duplicated %s times,",
      "with suitable adjustments made to `scale` and `rscales`."
    )
    duplication_message <- sprintf(fmt = duplication_message, k)
    message(duplication_message)
  }

  # Match control replicate columns to primary replicate columns ----

  ##_ Use user-supplied matching if supplied, otherwise match at random
  if (!is.null(control_col_matches)) {
    if (length(control_col_matches) != R_primary) {
      stop(sprintf("`control_col_matches` must have %s entries.", R_primary))
    }
    if (length(setdiff(seq_len(R_control), control_col_matches)) > 0) {
      stop(sprintf("All elements of the sequence 1,...,%s must be in `control_col_matches`", R_control))
    }
    if (length(setdiff(control_col_matches, c(seq_len(R_control), NA)))) {
      stop(sprintf("`control_col_matches` should only contain values of NA or the sequence 1,...,%s",
                   R_control))
    }

    matched_control_cols <- control_col_matches
    matched_primary_cols <- sapply(
      seq_len(R_control), function(i) {
      result <- which(matched_control_cols == i)
      return(result)
    })

  } else {
    matched_primary_cols <- sample(x = seq_len(R_primary),
                                   size = R_control,
                                   replace = FALSE)
    matched_control_cols <- sapply(seq_len(R_primary), function(i) {
      result <- which(matched_primary_cols == i)
      if (length(result) == 0) {
        NA_integer_
      } else {
        result
      }
    })

    matching_msg <- paste("Matching between primary and control replicates will be done at random.",
                          "For tips on reproducible matching, see `help('calibrate_to_sample')`",
                          sep = "\n")
    message(matching_msg)

  }



  # Generate replicate factors to account for difference in methods ----

  a_r <- rep(0, times = R_primary)

  A_updated <- A_control/A_primary
  rscales_updated <- sapply(seq_len(R_control), function(i) {
    rscales_control[i]/rscales_primary[matched_primary_cols[i]]
  })

  a_r[matched_primary_cols] <- sqrt(A_updated * rscales_updated)

  # Create needed data matrices for primary design ----

    ##_ Get dataframe with all of the variables mentioned in the formula
    mf <- model.frame(cal_formula, primary_rep_design$variables, na.action=na.pass)

    ##__For each factor/character variable, obtain a matrix of dummy variables, one dummy per category.
    ##_ For numeric variables, simply return a matrix with the variable itself
    xx <- lapply(
      attr(terms(cal_formula),"variables")[-1],
      function(tt) model.matrix( eval(bquote( ~0 + .(tt))), mf)
    )

    ##_ Determine dimensions of matrix of all variables
    cols <- sapply(xx, NCOL)

    x <- matrix(nrow = NROW(xx[[1]]), ncol = sum(cols))

    scols <- c(0, cumsum(cols))

    ##_ Combine all of the separate matrices into a single matrix
    for (i in 1:length(xx)){
      x[,scols[i]+1:cols[i]] <- xx[[i]]
    }
    colnames(x) <- do.call("c",lapply(xx,colnames))

  # Extract replicate weights matrix ----

  primary_replicate_weights <- weights(primary_rep_design, type = 'analysis')

  # Generate replicate-specific control totals ----

  unadjusted_control_totals <- survey::svytotal(x = cal_formula,
                                                design = control_rep_design,
                                                return.replicates = TRUE)

  unadjusted_control_totals <- list(
    'full-sample' = coef(unadjusted_control_totals),
    'replicate-specific' = as.matrix(unadjusted_control_totals$replicates)
  )

  replicate_control_totals <- matrix(data = unadjusted_control_totals[['full-sample']],
                                     nrow = R_primary,
                                     ncol = length(unadjusted_control_totals[['full-sample']]),
                                     byrow = TRUE)

  for (i in seq_len(R_control)) {
    i_star <- matched_primary_cols[i]
    replicate_control_totals[i_star,] <- unadjusted_control_totals[['full-sample']] +
      a_r[i_star] * (unadjusted_control_totals[['replicate-specific']][i,] - unadjusted_control_totals[['full-sample']])
  }

  # Ensure that order of control totals matches order of data variables ----

  primary_calib_variables <- colnames(x)
  control_calib_variables <- names(unadjusted_control_totals[['full-sample']])

  differing_variables <- union(setdiff(primary_calib_variables, control_calib_variables),
                               setdiff(control_calib_variables, primary_calib_variables))
  if (length(differing_variables) > 0) {
    error_msg <- paste(
      "There are differences between `primary_rep_design` and `control_rep_design`",
      "in the type or categories for the calibration variables."
    )
    stop(error_msg)
  }
  x <- x[,names(unadjusted_control_totals[['full-sample']]), drop = FALSE]

  # Calibrate the replicate weights ----

  adjusted_replicate_weights <- matrix(nrow = nrow(primary_replicate_weights),
                                       ncol = ncol(primary_replicate_weights))
  for (i in seq_len(R_primary)) {
    g_weights <- survey::grake(mm = x, ww = primary_replicate_weights[,i, drop = TRUE],
                               population = replicate_control_totals[i, ,drop = TRUE],
                               calfun = calfun,
                               bounds = bounds,
                               verbose = verbose, maxit = maxit,
                               epsilon = epsilon, variance = variance)
    if (is.null(attr(g_weights, 'failed'))) {
      convergence_achieved <- TRUE
    } else {
      convergence_achieved <- FALSE
    }
    if (!convergence_achieved) {
      error_msg <- sprintf("Convergence was not achieved for replicate %s. Consider increasing `maxit` or relaxing `epsilon`.", i)
      stop(error_msg)
    }
    adjusted_replicate_weights[,i] <- as.vector(primary_replicate_weights[,i]) * g_weights
  }

  # Calibrate the full-sample weights ----

  g_weights <- survey::grake(mm = x, ww = as.vector(primary_rep_design$pweights),
                             population = unadjusted_control_totals[['full-sample']],
                             calfun = calfun,
                             bounds = bounds,
                             verbose = verbose, maxit = maxit,
                             epsilon = epsilon, variance = variance)

  if (is.null(attr(g_weights, 'failed'))) {
    convergence_achieved <- TRUE
  } else {
    convergence_achieved <- FALSE
  }
  if (!convergence_achieved) {
    error_msg <- "Convergence was not achieved for calibration of full-sample weights. Consider increasing `maxit` or relaxing `epsilon`."
    stop(error_msg)
  }

  adjusted_fullsample_weights <- as.vector(primary_rep_design$pweights) * g_weights
  attr(adjusted_fullsample_weights, 'eta') <- NULL

  # Assemble the updated replicate design object ----

  calibrated_rep_design <- primary_rep_design

  calibrated_rep_design$pweights <- adjusted_fullsample_weights
  names(calibrated_rep_design$pweights) <- names(primary_rep_design)

  calibrated_rep_design$repweights <- adjusted_replicate_weights
  class(calibrated_rep_design$repweights) <- 'repweights'
  calibrated_rep_design$combined.weights <- TRUE

  if (primary_rep_design$type %in% c("JK1", "JKn", "JK2", "ACS", "successive-difference")) {
    primary_rep_design_rho <- NULL
  } else {
    primary_rep_design_rho <- primary_rep_design$rho
  }

  primary_rep_design_type <- ifelse(
    primary_rep_design$type %in% c("bootstrap", "subbootstrap", "mrbbootstrap"),
    "bootstrap", primary_rep_design$type
  )

  calibrated_rep_design <- survey::svrepdesign(
    data = primary_rep_design$variables,
    repweights = adjusted_replicate_weights,
    weights = adjusted_fullsample_weights,
    type = primary_rep_design_type,
    combined.weights = TRUE,
    rho = primary_rep_design_rho,
    scale = primary_rep_design$scale,
    rscales = primary_rep_design$rscales,
    fpc = primary_rep_design$fpc,
    fpctype = primary_rep_design$fpctype,
    mse = TRUE
  )

  if (is_db_backed_design) {
    # Replace 'variables' with a database connection
    # and make the object have the appropriate class
    calibrated_rep_design$variables <- NULL
    if (db_info$dbtype == "ODBC") {
      stop("'RODBC' no longer supported. Use the odbc package")
    } else {
      db <- DBI::dbDriver(db_info$dbtype)
      dbconn <- DBI::dbConnect(db, db_info$dbname)
    }
    calibrated_rep_design$db <- list(
      dbname = db_info$dbname, tablename = db_info$tablename,
      connection = dbconn,
      dbtype = db_info$dbtype
    )
    class(calibrated_rep_design) <- c(
      "DBIrepdesign", "DBIsvydesign", class(calibrated_rep_design)
    )
  }

  if (is_tbl_svy && ('package:srvyr' %in% search())) {
    calibrated_rep_design <- srvyr::as_survey_rep(
      calibrated_rep_design
    )
  }

  calibrated_rep_design$rho <- primary_rep_design$rho

  if (!primary_rep_design$mse) {
    warning("Setting `mse` to TRUE; variance estimates will be centered around full-sample estimate, not mean of replicates.")
  }

  # Indicate which replicate columns correspond ----
  # to which replicate columns of the control survey ----
  calibrated_rep_design$control_column_matches <- matched_control_cols

  # Set degrees of freedom to match that of the primary survey ----
  calibrated_rep_design$degf <- survey::degf(primary_rep_design)

  # Return the result ----
  return(calibrated_rep_design)
}
