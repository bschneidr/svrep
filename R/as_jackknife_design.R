#' @title Convert a survey design object to a random-groups jackknife design
#' @description
#' Forms a specified number of jackknife replicates
#' based on grouping primary sampling units (PSUs)
#' into random, (approximately) equal-sized groups.
#' @section Formation of Random Groups:
#' Within each value of \code{VAR_STRAT},
#' the data are sorted by first-stage sampling strata,
#' and then the PSUs in each stratum are randomly arranged.
#' Groups are then formed by serially placing PSUs
#' into each group.
#' The first PSU in the \code{VAR_STRAT} is placed into the first group,
#' the second PSU into the second group, and so on.
#' Once a PSU has been assigned to the last group,
#' the process begins again by assigning the next PSU to the first group,
#' the PSU after that to the second group, and so on.
#'
#'
#' The random group that each observation is assigned to
#' can be saved as a variable in the data
#' by using the function argument \code{group_var_name}.
#' @param design A survey design object created using the 'survey' (or 'srvyr') package,
#' with class \code{'survey.design'} or \code{'svyimputationList'}.
#' @param replicates The number of replicates to create
#' for each variance stratum. The total number of replicates
#' created is the number of variance strata times \code{replicates}.
#' Every design stratum must have at least as many primary sampling units (PSUs),
#' as \code{replicates}.
#' @param adj_method Specifies how to calculate the
#' replicate weight adjustment factor.
#' Available options for \code{adj_method} include:
#' \itemize{
#'   \item \code{"variance-stratum-psus"} (the default) \cr
#'     The replicate weight adjustment for a unit
#'     is based on the number of PSUs in its variance stratum.
#'   \item \code{"variance-units"} \cr
#'     The replicate weight adjustment for a unit
#'     is based on the number of variance units
#'     in its variance stratum.
#' }
#' See the section "Adjustment and Scale Methods" for details.
#' @param scale_method Specifies how to calculate the
#' scale factor for each replicate.
#' Available options for \code{scale_method} include:
#' \itemize{
#'   \item \code{"variance-stratum-psus"} \cr
#'     The scale factor for a variance unit
#'     is based on its number of PSUs compared
#'     to the number of PSUs in its variance stratum.
#'   \item \code{"variance-units"} \cr
#'     The scale factor for a variance unit is
#'     based on the number of variance units in
#'     its variance stratum.
#' }
#' See the section "Adjustment and Scale Methods" for details.
#' @param var_strat Specifies the name of a variable
#' in the data that defines variance strata to use
#' for the grouped jackknife. If \code{var_strat = NULL},
#' then there is effectively only one variance stratum.
#' @param var_strat_frac Specifies the sampling fraction
#' to use for finite population corrections in each
#' value of \code{var_strat}. Can use either a single number
#' or a variable in the data corresponding to \code{var_strat}.
#' @param sort_var (Optional) Specifies the name of a variable
#' in the data which should be used to sort the data before
#' assigning random groups. If a variable is specified for
#' \code{var_strat}, the sorting will happen within values of
#' that variable.
#' @param group_var_name (Optional) The name of a new variable created to save
#' identifiers for which random group each PSU was grouped into
#' for the purpose of forming replicates.
#' Specify \code{group_var_name = NULL} to avoid creating the variable in the data.
#' @param compress Use a compressed representation of the replicate weights matrix.
#' This reduces the computer memory required to represent the replicate weights and has no
#' impact on estimates.
#' @param mse If \code{TRUE}, compute variances from sums of squares around the point estimate from the full-sample weights.
#' If \code{FALSE}, compute variances from sums of squares around the mean estimate from the replicate weights.
#'
#' @return
#' A replicate design object, with class \code{svyrep.design}, which can be used with the usual functions,
#' such as \code{svymean()} or \code{svyglm()}.
#'
#' Use \code{weights(..., type = 'analysis')} to extract the matrix of replicate weights. \cr
#' Use \code{as_data_frame_with_weights()} to convert the design object to a data frame with columns
#' for the full-sample and replicate weights.
#' @section Adjustment and Scale Methods:
#'
#' The jackknife replication variance estimator based on \eqn{R} replicates takes the following form:
#' \deqn{
#'   v(\hat{\theta}) = \sum_{r=1}^{R} (1 - f_r) \times c_r \times \left(\hat{\theta}_r - \hat{\theta}\right)^2
#' }
#' where \eqn{r} indexes one of the \eqn{R} sets of replicate weights,
#' \eqn{c_r} is a corresponding scale factor for the \eqn{r}-th replicate,
#' and \eqn{1 - f_r} is an optional finite population correction factor
#' that can potentially differ across variance strata.
#'
#' To form the replicate weights, the PSUs are divided into \eqn{\tilde{H}} variance strata,
#' and the \eqn{\tilde{h}}-th variance stratum contains \eqn{G_{\tilde{h}}}
#' random groups. The number of replicates \eqn{R} equals the total number
#' of random groups across all variance strata:
#' \eqn{R = \sum_{\tilde{h}}^{\tilde{H}} G_{\tilde{h}}}. In other words,
#' each replicate corresponds to one of the random groups from one of the variance strata.
#'
#' The weights for replicate \eqn{r} corresponding to random group \eqn{g} within
#' variance stratum \eqn{\tilde{h}} is defined as follows.
#'
#' If case \eqn{i}
#' is not in variance stratum \eqn{\tilde{h}}, then \eqn{w_{i}^{(r)} = w_i}.
#'
#' If case \eqn{i} is in variance stratum \eqn{\tilde{h}} and not in random group \eqn{g},
#' then \eqn{w_{i}^{(r)} = a_{\tilde{h}g} w_i}.
#'
#' Otherwise, if case \eqn{i} is in random group \eqn{g}
#' of variance stratum \eqn{\tilde{h}}, then \eqn{w_{i}^{(r)} = 0}.
#'
#' The R function argument \code{adj_method} determines how
#' the adjustment factor \eqn{a_{\tilde{h} g}} is calculated.
#' When \code{adj_method = "variance-units"}, then
#' \eqn{a_{\tilde{h} g}} is calculated based on \eqn{G_{\tilde{h}}},
#' which is the number of random groups in variance stratum \eqn{\tilde{h}}.
#' When \code{adj_method = "variance-stratum-psus"}, then
#' \eqn{a_{\tilde{h} g}} is calculated based on \eqn{n_{\tilde{h}g}},
#' which is the number of PSUs in random group \eqn{g} in variance stratum \eqn{\tilde{h}},
#' as well as \eqn{n_{\tilde{h}}}, the total number of PSUs in variance stratum \eqn{\tilde{h}}.
#'
#' If \code{adj_method = "variance-units"}, then: \deqn{a_{\tilde{h}g} = \frac{G_{\tilde{h}}}{G_{\tilde{h}} - 1}}
#'
#' If \code{adj_method = "variance-stratum-psus"}, then: \deqn{a_{\tilde{h}g} = \frac{n_{\tilde{h}}}{n_{\tilde{h}} - n_{\tilde{h}g}}}
#'
#' The scale factor \eqn{c_r} for replicate \eqn{r}
#' corresponding to random group \eqn{g} within variance stratum \eqn{\tilde{h}} is
#' calculated according to the function argument \code{scale_method}.
#'
#' If \code{scale_method = "variance-units"}, then: \deqn{c_r = \frac{G_{\tilde{h}} - 1}{G_{\tilde{h}}}}
#'
#' If \code{scale_method = "variance-stratum-psus"}, then: \deqn{c_r = \frac{n_{\tilde{h}} - n_{\tilde{h}g}}{n_{\tilde{h}}}}
#'
#' The sampling fraction \eqn{f_r} used for finite population correction \eqn{1 - f_r}
#' is by default assumed to equal 0. However, the user can supply a sampling fraction
#' for each variance stratum using the argument \code{var_strat_frac}.
#'
#' When variance units in a variance stratum
#' have differing numbers of PSUs,
#' the combination \code{adj_method = "variance-stratum-psus"}
#' and \code{scale_method = "variance-units"} is
#' recommended by Valliant, Brick, and Dever (2008),
#' corresponding to their method \code{"GJ2"}.
#'
#' The random-groups jackknife method often referred to as "DAGJK"
#' corresponds to the options \code{var_strat = NULL},
#' \code{adj_method = "variance-units"}, and \code{scale_method = "variance-units"}.
#' The DAGJK method will yield upwardly-biased variance estimates for totals
#' if the total number of PSUs is not a multiple of the total number of replicates (Valliant, Brick, and Dever 2008).
#' @references
#' See Section 15.5 of Valliant, Dever, and Kreuter (2018)
#' for an introduction to the grouped jackknife and
#' guidelines for creating the random groups.
#'
#' - Valliant, R., Dever, J., Kreuter, F. (2018).
#' "Practical Tools for Designing and Weighting Survey Samples, 2nd edition." New York: Springer.
#'
#' See Valliant, Brick, and Dever (2008)
#' for statistical details related to the
#' \code{adj_method} and \code{scale_method} arguments.
#'
#' - Valliant, Richard, Michael Brick, and Jill Dever. 2008.
#' "Weight Adjustments for the Grouped Jackknife Variance Estimator."
#' \emph{Journal of Official Statistics}. 24: 469-88.
#'
#' See Chapter 4 of Wolter (2007) for additional details of the jackknife,
#' including the method based on random groups.
#'
#' - Wolter, Kirk. 2007. "Introduction to Variance Estimation." New York, NY: Springer New York. https://doi.org/10.1007/978-0-387-35099-8.
#' @export
#'
#' @examples
#' library(survey)
#'
#' # Load example data
#'
#'  data('api', package = 'survey')
#'
#'  api_strat_design <- svydesign(
#'    data = apistrat,
#'    id = ~ 1,
#'    strata = ~stype,
#'    weights = ~pw
#'  )
#'
#' # Create a random-groups jackknife design
#'
#'  jk_design <- as_random_group_jackknife_design(
#'    api_strat_design,
#'    replicates = 15
#'  )
#'  print(jk_design)
as_random_group_jackknife_design <- function(
    design,
    replicates = 50,
    var_strat = NULL,
    var_strat_frac = NULL,
    sort_var = NULL,
    adj_method = "variance-stratum-psus",
    scale_method = "variance-stratum-psus",
    group_var_name = ".random_group",
    compress = TRUE,
    mse = getOption("survey.replicates.mse")
) {
  UseMethod("as_random_group_jackknife_design", design)
}


#' @export
as_random_group_jackknife_design.survey.design <- function(
    design,
    replicates = 50,
    var_strat = NULL,
    var_strat_frac = NULL,
    sort_var = NULL,
    adj_method = "variance-stratum-psus",
    scale_method = "variance-stratum-psus",
    group_var_name = ".random_group",
    compress = TRUE,
    mse = getOption("survey.replicates.mse")
  ) {

  # Check that valid parameters specified
  if (!adj_method %in% c("variance-stratum-psus", "variance-units")) {
    stop("`adj_method` must be either 'variance-stratum-psus' or 'variance-units'")
  }
  if (!scale_method %in% c("variance-stratum-psus", "variance-units")) {
    stop("`scale_method` must be either 'variance-stratum-psus' or 'variance-units'")
  }

  # Begin extracting necessary design information
  n <- nrow(design)
  design_vars <- data.frame(
    'ROW_ID' = seq_len(n),
    'STRATUM' = design$strata[,1,drop=TRUE]
  )

  # Handle VAR_STRAT variable
  if (!is.null(var_strat)) {
    if (!(var_strat %in% colnames(design$variables))) {
      stop("`var_strat` must be either NULL or the name of a variable in the data.")
    }
    design_vars[['VAR_STRAT']] <- design$variables[[var_strat]]
    if (any(is.na(design_vars[['VAR_STRAT']]))) {
      stop("The `var_strat` variable cannot have any missing values.")
    }
  }
  if (is.null(var_strat)) {
    design_vars[['VAR_STRAT']] <- 1
  }
  if (!is.null(sort_var)) {
    if (!(sort_var %in% colnames(design$variables))) {
      stop("`sort_var` must be either NULL or the name of a variable in the data.")
    }
    design_vars[['SORT_VAR']] <- design$variables[[sort_var]]
    if (any(is.na(design_vars[['SORT_VAR']]))) {
      stop("The `sort_var` variable cannot have any missing values.")
    }
  }
  if (is.null(sort_var)) {
    design_vars[['SORT_VAR']] <- 1
  }

  # Create unique PSU IDs within
  # combinations of variance strata and design strata
  design_vars[['PSU']] <- interaction(design_vars[['VAR_STRAT']],
                                      design$strata[,1,drop=TRUE],
                                      design$cluster[,1,drop=TRUE],
                                      drop = TRUE) |> as.numeric()

  # Warn the user about FPCs and check that `var_strat_frac` is valid

  if (!is.null(design$fpc$popsize)) {
    warning_msg <- paste(
      "Ignoring finite population corrections in `design`.",
      sprintf("Instead using `var_strat_frac`: %s.",
              var_strat_frac), paste = "\n"
    ) |> warning()
  }

  if (is.character(var_strat_frac)) {
    if (!var_strat_frac %in% colnames(design$variables)) {
      stop("`var_strat_frac` must be either a single number or the name of a variable in the data.")
    }
    design_vars[['VAR_STRAT_FRAC']] <- design$variables[[var_strat_frac]]
    if (any(is.na(design_vars[[var_strat_frac]]))) {
      stop("The `var_strat_frac` variable cannot have any missing values.")
    }

    n_unique_combos <- nrow(unique(design_vars[,c("VAR_STRAT_FRAC", "VAR_STRAT_FRAC")]))
    n_unique_var_strat <- length(unique(design_vars[['VAR_STRAT_FRAC']]))
    if (n_unique_combos != n_unique_var_strat) {
      stop("Each value of `var_strat` must correspond to only one value of `var_strat_frac`.")
    }
  }

  if (is.numeric(var_strat_frac)) {
    if ((length(var_strat_frac) != 1) || is.na(var_strat_frac)) {
      stop("`var_strat_frac` must be either a single number or the name of a variable in the data.")
    }
    if ((var_strat_frac < 0) | (var_strat_frac > 1)) {
      stop("`var_strat_frac` must be between 0 and 1.")
    }
    design_vars[['VAR_STRAT_FRAC']] <- var_strat_frac
  }

  if (is.null(var_strat_frac)) {
    design_vars[['VAR_STRAT_FRAC']] <- 0
  }

  # Check that the number of replicates is valid

  min_psus_in_any_strata <- min(design$fpc$sampsize[,1,drop=TRUE])

  if (min_psus_in_any_strata < replicates) {
    stop("There is at least one stratum with fewer PSUs than the desired number of replicates")
  }

  # Randomly shuffle the order of the PSUs in the data, within strata
  # (implemented by randomly relabeling the PSUs,
  #  then sorting the data by the random PSU labels)
  n_psus <- design_vars$PSU |> unique() |> length()
  psu_random_labels <- sample(x = n_psus, size = n_psus, replace = FALSE)
  design_vars[['RAND_PSU_ID']] <- interaction(
    design_vars[['STRATUM']],
    psu_random_labels[design_vars[['PSU']]],
    drop = TRUE, lex.order = TRUE
  ) |> as.numeric()

  # Order the data by varstrat, then stratum, then sort variable, then by PSU
  if (packageVersion("base") <= "4.4.0") {
    sort_by <- function(x, y, ...) {
        if (inherits(y, "formula")) 
            y <- .formula2varlist(y, x)
        if (!is.list(y)) 
            y <- list(y)
        o <- do.call(order, c(unname(y), list(...)))
        x[o, , drop = FALSE]
    }
  }
  design_vars <- design_vars |> sort_by(
    ~ VAR_STRAT + STRATUM + SORT_VAR + RAND_PSU_ID
  )

  # Create random groups separately by VAR_STRAT
  design_vars[['RANDOM_GROUP_VAR_UNIT']] <- NA_real_
  unique_var_strat <- design_vars[['VAR_STRAT']] |> unique() |> sort()

  for (i in seq_along(unique_var_strat)) {
    variance_stratum <- unique_var_strat[i]

    row_indices <- (design_vars[['VAR_STRAT']] == variance_stratum)
    n_psus <- design_vars[['PSU']][row_indices] |> unique() |> length()

    group_size <- n_psus / replicates

    random_group_assignments <- rep(
      seq_len(replicates), times = ceiling(group_size)
    )[seq_len(n_psus)]

    design_vars[['RANDOM_GROUP_VAR_UNIT']][row_indices] <- (
      random_group_assignments[
        design_vars[['RAND_PSU_ID']][row_indices] |>
          factor() |> as.numeric()
      ]
    )
  }

  design_vars[['RANDOM_GROUP_VAR_UNIT']] <- interaction(
    design_vars$VAR_STRAT,
    design_vars$RANDOM_GROUP_VAR_UNIT,
    drop = TRUE, lex.order = TRUE
  ) |> as.numeric()




  # Create the jackknife replicate weights

  if ((adj_method == "variance-units") & (scale_method == "variance-units")) {

    # Reorder the data back to its original order
    design_vars <- design_vars[order(design_vars[['ROW_ID']]),]

    jk_design <- survey::svydesign(
      data = design$variables,
      ids = design_vars[['RANDOM_GROUP_VAR_UNIT']],
      strata = design_vars[['VAR_STRAT']],
      weights = weights(design, type = 'sampling')
    ) |> survey::as.svrepdesign(
      type = "JKn",
      compress = compress,
      mse = mse
    )
  }

  if ((adj_method != "variance-units") || (scale_method != "variance-units")) {

    # Obtain counts of PSUs in each VAR_UNIT and VAR_STRAT
    var_unit_psu_counts <- by(
      design_vars, INDICES = design_vars$RANDOM_GROUP_VAR_UNIT,
      FUN = function(df) {
        data.frame(VAR_STRAT = df$VAR_STRAT[1],
                   VAR_STRAT_FRAC = df$VAR_STRAT_FRAC[1],
                   RANDOM_GROUP_VAR_UNIT = df$RANDOM_GROUP_VAR_UNIT[1],
                   N_PSUS_IN_VAR_UNIT = length(unique(df$RAND_PSU_ID)))
      }) |> do.call(what = rbind)

    var_unit_psu_counts <- by(
      var_unit_psu_counts, INDICES = var_unit_psu_counts$VAR_STRAT,
      FUN = function(df) {
        data.frame(VAR_STRAT = df$VAR_STRAT,
                   VAR_STRAT_FRAC = df$VAR_STRAT_FRAC,
                   RANDOM_GROUP_VAR_UNIT = df$RANDOM_GROUP_VAR_UNIT,
                   N_VAR_UNITS_IN_VAR_STRAT = length(df$RANDOM_GROUP_VAR_UNIT),
                   N_PSUS_IN_VAR_UNIT = df$N_PSUS_IN_VAR_UNIT,
                   N_PSUS_IN_VAR_STRAT = sum(df$N_PSUS_IN_VAR_UNIT))
      }
    ) |> do.call(what = rbind)  |> `rownames<-`(NULL)

    var_unit_psu_counts <- (
      var_unit_psu_counts[sort(var_unit_psu_counts$RANDOM_GROUP_VAR_UNIT),]
    )

    # Calculate adjustment factors
    adj_factors <- switch(
      adj_method,
      "variance-stratum-psus" = (
        var_unit_psu_counts$N_PSUS_IN_VAR_STRAT / (var_unit_psu_counts$N_PSUS_IN_VAR_STRAT - var_unit_psu_counts$N_PSUS_IN_VAR_UNIT)
      ),
      "variance-units" = (
        var_unit_psu_counts$N_VAR_UNITS_IN_VAR_STRAT / (var_unit_psu_counts$N_VAR_UNITS_IN_VAR_STRAT - 1)
      )
    )

    rep_factors <- matrix(
      data = 1,
      nrow = nrow(var_unit_psu_counts),
      ncol = nrow(var_unit_psu_counts)
    )

    for (variance_stratum in unique(var_unit_psu_counts$VAR_STRAT)) {
      indices <- which(var_unit_psu_counts$VAR_STRAT == variance_stratum)
      for (index in indices) {
        rep_factors[indices, index] <- adj_factors[index]
      }

    }

    diag(rep_factors) <- 0

    colnames(rep_factors) <- paste0("REP_", seq_len(ncol(rep_factors)))

    # Calculate scale factor for each replicate
    scale_factors <- switch(
      scale_method,
      "variance-stratum-psus" = (
        var_unit_psu_counts$N_PSUS_IN_VAR_STRAT / (var_unit_psu_counts$N_PSUS_IN_VAR_STRAT - var_unit_psu_counts$N_PSUS_IN_VAR_UNIT)
      )^(-1),
      "variance-units" = (
        var_unit_psu_counts$N_VAR_UNITS_IN_VAR_STRAT / (var_unit_psu_counts$N_VAR_UNITS_IN_VAR_STRAT - 1)
      )^(-1)
    )

    scale_factors <- scale_factors * (1 - var_unit_psu_counts$VAR_STRAT_FRAC)

    # Combine replicate factors and the design variables
    design_vars <- merge(
      x = design_vars,
      y = cbind(
        RANDOM_GROUP_VAR_UNIT = var_unit_psu_counts[['RANDOM_GROUP_VAR_UNIT']],
        rep_factors
      ),
      by = "RANDOM_GROUP_VAR_UNIT"
    )

    if (length(unique(scale_factors)) == 1) {
      overall_scale <- scale_factors[1]
      rscales <- rep(1, times = length(scale_factors))
    } else {
      overall_scale <- 1
      rscales <- scale_factors
    }

    # Reorder the data back to its original order
    design_vars <- design_vars[order(design_vars[['ROW_ID']]),]

    # Combine data, replicate factors, and replication coefficients
    jk_design <- survey::svrepdesign(
      variables = design$variables,
      repweights = design_vars[,colnames(rep_factors)],
      weights = weights(design, type = "sampling"),
      scale = overall_scale,
      rscales = rscales,
      type = "other",
      mse = mse,
      combined.weights = FALSE
    )

    if (compress) {
      jk_design <- survey::compressWeights(jk_design)
    }
  }

  jk_design$type <- "Random-groups jackknife"

  # Add the groups variable to the data
  if (!is.null(group_var_name)) {
    if (group_var_name %in% colnames(design$variables)) {
      sprintf("Overwriting an existing variable named `%s`",
              group_var_name) |> warning()
    }
    jk_design$variables[[group_var_name]] <- design_vars[['RANDOM_GROUP_VAR_UNIT']]
  }

  if (inherits(design, 'tbl_svy') && ('package:srvyr' %in% search())) {
    jk_design <- srvyr::as_survey_rep(
      jk_design
    )
  }

  # Return the result
  jk_design$call <- sys.call(which = -1)
  return(jk_design)
}

#' @export
as_random_group_jackknife_design.DBIsvydesign <- function(
    design,
    replicates = 50,
    var_strat = NULL,
    var_strat_frac = NULL,
    sort_var = NULL,
    adj_method = "variance-stratum-psus",
    scale_method = "variance-stratum-psus",
    group_var_name = ".random_group",
    compress = TRUE,
    mse = getOption("survey.replicates.mse")
) {

  rep_design <- NextMethod(design)

  # Replace 'variables' with a database connection
  # and make the object have the appropriate class
  rep_design$variables <- NULL
  if (design$db$dbtype == "ODBC") {
    stop("'RODBC' no longer supported. Use the odbc package")
  } else {
    db <- DBI::dbDriver(design$db$dbtype)
    dbconn <- DBI::dbConnect(db, design$db$dbname)
  }
  rep_design$db <- list(
    dbname = design$db$dbname, tablename = design$db$tablename,
    connection = dbconn,
    dbtype = design$db$dbtype
  )
  class(rep_design) <- c(
    "DBIrepdesign", "DBIsvydesign",
    setdiff(class(rep_design), c("DBIrepdesign", "DBIsvydesign"))
  )

  rep_design$call <- sys.call(which = -1)

  return(rep_design)
}
