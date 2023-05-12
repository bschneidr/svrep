#' @title Convert a survey design object to a random-groups jackknife design
#' @description
#' Forms a specified number of jackknife replicates
#' based on grouping primary sampling units (PSUs)
#' into random, equal-sized groups. The random-groups jackknife
#' is only valid for designs where the first-stage sampling fraction
#' is ignorable or is identical in every stratum.
#' @details The data are sorted by first-stage sampling strata,
#' and then the PSUs in each stratum are randomly arranged.
#' Groups are then formed by sequentially placing PSUs
#' into each group.
#' The first PSU in the dataset is placed into the first group,
#' the second PSU into the second group, and so on.
#' Once a PSU has been assigned to the last group,
#' the process begins again by assigning the next PSU to the first group,
#' the next PSU to the second group, and so on.
#'
#' @param design A survey design object created using the 'survey' (or 'srvyr') package,
#' with class \code{'survey.design'} or \code{'svyimputationList'}.
#' @param replicates The number of replicates to create.
#' Every stratum must have at least this many primary sampling units (PSUs),
#' or else an error will occur.
#' @param group_var_name (Optional) The name of a new variable created to save
#' identifiers for which random group each PSU was grouped into
#' for the purpose of forming replicates.
#' Specify \code{group_var_name = NULL} to avoid creating the variable in the data.
#' @param sampling_fraction By default, the random groups jackknife
#' assumes that sampling is with-replacement (i.e., the sampling fraction is ignorable).
#' However, if PSUs were not stratified or if every stratum had the same sampling fraction,
#' you can specify a single value of \code{sampling_fraction}, which will
#' be used to adjust the overall replication scale factor through a finite population correction.
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
#' @references
#' See Chapter 4 of Wolter (2007) for an overview of the jackknife,
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
    sampling_fraction = 0,
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
    sampling_fraction = 0,
    group_var_name = ".random_group",
    compress = TRUE,
    mse = getOption("survey.replicates.mse")
  ) {

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

  n_psus <- design_vars[['PSU']] |> unique() |> length()

  # Warn the user about FPCs

  if (!is.null(design$fpc$popsize)) {
    warning_msg <- paste(
      "Ignoring finite population corrections in `design`.",
      sprintf("Instead using `sampling_fraction`: %s.",
              sampling_fraction), paste = "\n"
    ) |> warning()
  }

  if (!is.numeric(sampling_fraction) || (length(sampling_fraction) != 1) || is.na(sampling_fraction)) {
    stop("`sampling_fraction` must be a single number.")
  }

  if ((sampling_fraction < 0) | (sampling_fraction > 1)) {
    stop("`samping_fraction` must be between 0 and 1.")
  }

  # Check that the number of replicates is valid

  min_psus_in_any_strata <- min(design$fpc$sampsize[,1,drop=TRUE])

  if (min_psus_in_any_strata < replicates) {
    stop("There is at least one stratum with fewer PSUs than the desired number of replicates")
  }

  # Randomly shuffle the order of the PSUs in the data, within strata
  # (implemented by randomly relabeling the PSUs,
  #  then sorting the data by the random PSU labels)
  psu_random_labels <- sample(x = n_psus, size = n_psus, replace = FALSE)
  design_vars[['RAND_PSU_ID']] <- interaction(
    design_vars[['STRATUM']],
    psu_random_labels[design_vars[['PSU']]],
    drop = TRUE, lex.order = TRUE
  )

  # Order the data by stratum, then by PSU
  design_vars <- design_vars[order(design_vars[['RAND_PSU_ID']]),,drop=FALSE]
  design_vars <- design_vars[order(design_vars[['STRATUM']]),]

  # Create random groups
  group_size <- n_psus / replicates

  random_group_assignments <- rep(
    seq_len(replicates), times = ceiling(group_size)
  )[seq_len(n_psus)]

  random_group_assignments[design_vars[['RAND_PSU_ID']]]

  design_vars[['RANDOM_GROUP']] <- random_group_assignments[design_vars[['RAND_PSU_ID']]]

  # Reorder the data by the original order

  design_vars <- design_vars[order(design_vars[['ROW_ID']]),]

  # Create the jackknife replicate weights
  jk_design <- survey::svydesign(
    data = design$variables,
    ids = design_vars[['RANDOM_GROUP']],
    strata = NULL,
    weights = weights(design, type = 'sampling')
  ) |> survey::as.svrepdesign(
    type = "JK1",
    compress = compress,
    mse = mse
  )

  jk_design$type <- "Random-groups jackknife"

  # Adjust the scale to account for sampling fraction
  jk_design$scale <- jk_design$scale * (1 - sampling_fraction)

  # Add the groups variable to the data
  if (!is.null(group_var_name)) {
    if (group_var_name %in% colnames(design$variables)) {
      sprintf("Overwriting an existing variable named `%s`",
              group_var_name) |> warning()
    }
    jk_design$variables[[group_var_name]] <- design_vars[['RANDOM_GROUP']]
  }

  # Return the result
  jk_design$call <- sys.call(which = -1)
  return(jk_design)
}
