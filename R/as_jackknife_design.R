#' @title Convert a survey design object to a random-groups jackknife design
#' @description
#' Forms a specified number of jackknife replicates
#' based on grouping primary sampling units (PSUs)
#' into random, (approximately) equal-sized groups.
#' @details Within each value of \code{VAR_STRAT},
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
#' @param design A survey design object created using the 'survey' (or 'srvyr') package,
#' with class \code{'survey.design'} or \code{'svyimputationList'}.
#' @param replicates The number of replicates to create.
#' Every stratum must have at least this many primary sampling units (PSUs),
#' or else an error will occur.
#' @param adj_method Specifies how to calculate the
#' replicate weight adjustment factor. These adjustment factors
#' generally take the following form:
#' \deqn{
#'  w_{h i(\tilde{h} g)} =
#'  \begin{cases} w_{h i} & (h i) \notin S_{\tilde{h}} \\
#'  a_{h(\tilde{h} g)} w_{h i} & (h i) \in S_{\tilde{h}},(h i) \notin S_{\tilde{h} h g} \\
#'  0 & (h i) \in S_{\tilde{h} h g}
#'  \end{cases}
#' }
#' The options for \code{adj_method} differ in how
#' \eqn{a_{h(\tilde{h} g)}} is calculated.
#' Available options include:
#' \itemize{
#'   \item{"variance-stratum-psus" (the default)}{
#'     The replicate weight adjustment for a unit
#'     is based on the number of PSUs in its variance stratum.
#'   }
#'   \item{"variance-units"}{
#'     The replicate weight adjustment for a unit
#'     is based on the number of variance units
#'     in its variance stratum.
#'   }
#' }
#' @param scale_method Specifies how to calculate the
#' scale factor for each replicate.
#' Available options include:
#' \itemize{
#'   \item{"variance-stratum-psus"}{
#'     The scale factor for a variance unit
#'     is based on its number of PSUs compared
#'     to the number of PSUs in its variance stratum.
#'   }
#'   \item{"variance-units"}{
#'     The scale factor for a variance unit is
#'     based on the number of variance units in
#'     its variance stratum.
#'   }
#' }
#' When variance units in a variance stratum
#' have differing numbers of PSUs,
#' the combination \code{adj_method = "variance-stratum-psus"}
#' and \code{scale_method = "variance-units"} is
#' recommended by Valliant, Brick, and Dever (2008),
#' corresponding to their method \code{"GJ2"}.
#' @param var_strat Specifies the name of a variable
#' in the data that defines variance strata to use
#' for the grouped jackknife. If \code{var_strat = NULL},
#' then there is effectively only one variance stratum.
#' @param var_strat_frac Specifies the sampling fraction
#' to use for finite population corrections in each
#' value of \code{var_strat}. Can use either a single number
#' or a variable in the data corresponding to \code{var_strat}.
#' @param group_var_name (Optional) The name of a new variable created to save
#' identifiers for which random group each PSU was grouped into
#' for the purpose of forming replicates.
#' Specify \code{group_var_name = NULL} to avoid creating the variable in the data.
#' @param compress Use a compressed representation of the replicate weights matrix.
#' This reduces the computer memory required to represent the replicate weights and has no
#' impact on estimates.
#' @param mse If \code{TRUE}, compute variances from sums of squares around the point estimate from the full-sample weights,
#' If \code{FALSE}, compute variances from sums of squares around the mean estimate from the replicate weights.
#'
#' @return
#' A replicate design object, with class \code{svyrep.design}, which can be used with the usual functions,
#' such as \code{svymean()} or \code{svyglm()}.
#'
#' Use \code{weights(..., type = 'analysis')} to extract the matrix of replicate weights. \cr
#' Use \code{as_data_frame_with_weights()} to convert the design object to a data frame with columns
#' for the full-sample and replicate weights.
#' @details
#' The random-groups jackknife method often referred to as "DAGJK"
#' corresponds to the options \code{var_strat = NULL},
#' \code{adj_method = "variance-units"}, and \code{scale_method = "variance-units"}.
#' The DAGJK method will yield upwardly-biased variance estimates for totals
#' if the total number of PSUs is not a multiple of the total number of replicates (Valliant, Brick, and Dever 2008).
#'
#'
#' @references
#' See Valliant, Brick, and Dever (2008) for an overview of the
#' grouped jackknife and statistical details related to the
#' \code{adj_method} and \code{scale_method} arguments.
#'
#' - Valliant, Richard, Michael Brick, and Jill Dever. 2008.
#' "Weight Adjustments for the Grouped Jackknife Variance Estimator."
#' \emph{Journal of Official Statistics}. 24: 469–88.
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

  # Extract necessary design information
  n <- nrow(design)

  design_vars <- data.frame(
    'ROW_ID' = seq_len(n),
    'STRATUM' = design$strata[,1,drop=TRUE],
    'PSU' = interaction(design$strata[,1,drop=TRUE],
                        design$cluster[,1,drop=TRUE],
                        drop = TRUE) |> as.numeric(),
    stringsAsFactors = FALSE
  )

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

  # Order the data by stratum, then by PSU
  design_vars <- design_vars[order(design_vars[['RAND_PSU_ID']]),,drop=FALSE]
  design_vars <- design_vars[order(design_vars[['STRATUM']]),]
  design_vars <- design_vars[order(design_vars[['VAR_STRAT']]),]

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
        design_vars[['RAND_PSU_ID']][row_indices]
      ]
    )
  }

  design_vars[['RANDOM_GROUP_VAR_UNIT']] <- interaction(
    design_vars$VAR_STRAT,
    design_vars$RANDOM_GROUP_VAR_UNIT,
    drop = TRUE
  ) |> as.numeric()

  # Create the jackknife replicate weights

  if ((adj_method == "variance-units") & (scale_method == "variance-units")) {
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
        var_unit_psu_counts$N_VAR_UNITS_IN_VAR_STRAT / (var_unit_psu_counts$N_VAR_UNITS_IN_VAR_STRAT - var_unit_psu_counts$N_VAR_UNITS_IN_VAR_STRAT)
      )
    )

    rep_factors <- matrix(
      data = 1,
      nrow = nrow(var_unit_psu_counts),
      ncol = nrow(var_unit_psu_counts)
    )

    for (variance_stratum in unique(var_unit_psu_counts$VAR_STRAT)) {
      indices <- (var_unit_psu_counts$VAR_STRAT == variance_stratum)
      rep_factors[indices, indices] <- adj_factors[indices]
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
        var_unit_psu_counts$N_VAR_UNITS_IN_VAR_STRAT / (var_unit_psu_counts$N_VAR_UNITS_IN_VAR_STRAT - var_unit_psu_counts$N_VAR_UNITS_IN_VAR_STRAT)
      )^(-1)
    )

    scale_factors <- scale_factors * (1 - var_unit_psu_counts$VAR_STRAT_FRAC)

    # Reorder the data by the original order

    design_vars <- merge(
      x = design_vars,
      y = cbind(
        RANDOM_GROUP_VAR_UNIT = var_unit_psu_counts[['RANDOM_GROUP_VAR_UNIT']],
        rep_factors
      ),
      by = "RANDOM_GROUP_VAR_UNIT"
    )

    design_vars <- design_vars[order(design_vars[['ROW_ID']]),]

    if (length(unique(scale_factors)) == 1) {
      overall_scale <- scale_factors[1]
      rscales <- rep(1, times = length(scale_factors))
    } else{
      overall_scale <- 1
      rscales <- scale_factors
    }

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

  # Return the result
  jk_design$call <- sys.call(which = -1)
  return(jk_design)
}