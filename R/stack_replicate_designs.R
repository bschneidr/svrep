#' Stack replicate designs, combining data and weights into a single object
#'
#' @description Stack replicate designs: combine rows of data, rows of replicate weights, and the respective full-sample weights.
#' This can be useful when comparing estimates before and after a set of adjustments made to the weights.
#' Another more delicate application is when combining sets of replicate weights from multiple years of data for a survey, although this must be done carefully based on guidance from a data provider.
#' @param ... Replicate-weights survey design objects to combine. These can be supplied in one of two ways. \cr
#' \itemize{
#'   \item Option 1 - A series of design objects, for example \code{'adjusted' = adjusted_design, 'orig' = orig_design}.
#'   \item Option 2 - A list object containing design objects, for example  \code{list('nr' = nr_adjusted_design, 'ue' = ue_adjusted_design)}.
#' }
#' All objects must have the same specifications for \code{type}, \code{rho},
#' \code{mse}, \code{scales}, and \code{rscales}.
#' @param .id A single character value, which becomes the name of a new column of identifiers
#' created in the output data to link each row to the design from which it was taken. \cr
#' The labels used for the identifiers are taken from named arguments.
#'
#' @return A replicate-weights survey design object, with class \code{svyrep.design} and \code{svyrep.stacked}.
#' The resulting survey design object always has its value of \code{combined.weights} set to \code{TRUE}.
#' @export
#'
#' @examples
#' # Load example data, creating a replicate design object
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
#' # Stack the three designs, using any of the following syntax options
#' stacked_design <- stack_replicate_designs(orig_rep_design, ue_adjusted_design, nr_adjusted_design,
#'                                           .id = "which_design")
#' stacked_design <- stack_replicate_designs('original' = orig_rep_design,
#'                                           'unknown eligibility adjusted' = ue_adjusted_design,
#'                                           'nonresponse adjusted' = nr_adjusted_design,
#'                                           .id = "which_design")
#' list_of_designs <- list('original' = orig_rep_design,
#'                         'unknown eligibility adjusted' = ue_adjusted_design,
#'                         'nonresponse adjusted' = nr_adjusted_design)
#' stacked_design <- stack_replicate_designs(list_of_designs, .id = "which_design")

stack_replicate_designs <- function(..., .id = "Design_Name") {

  if (is.null(.id)) {
    stop("`.id` argument must be a single character value.")
  } else {
    if (!is.character(.id) || is.na(.id) || length(.id) != 1) {
      stop("`.id` argument must be a single character value.")
    }
  }

  # If the input is a series of arguments, ensure it's handled appropriately
  design_list <- list(...)

  if (length(design_list) == 0) {
    stop("No design objects were supplied.")
  }

  # If the input is a single list, ensure it's handled as such
  if (length(design_list) == 1 && is.list(design_list[[1]]) && !"svyrep.design" %in% class(design_list)) {
    input_type <- "list"
    design_list <- design_list[[1]]
    if (is.null(names(design_list))) {
      names(design_list) <- rep("", length(design_list))
    }
  } else {
    input_type <- "arguments"
  }

  is_svyrep <- sapply(design_list, function(x) "svyrep.design" %in% class(x))

  if (!all(is_svyrep)) {
    stop("All objects must be replicate-weights survey design objects, of class `svyrep.design`.")
  }

  # Get verbatim R code supplied to the argument
  # and use it to fill in missing names
  calls <- match.call(expand.dots = FALSE)$`...`
  calls <- as.character(lapply(calls, deparse))
  if (!is.null(names(design_list))) {
    design_names <- names(design_list)
  } else {
    design_names <- rep("", length(design_list))
  }
  missing_names <- which(is.na(design_names) | design_names == "")
  design_names <- names(design_list)
  if (input_type != "list") {
    design_names[missing_names] <- calls[missing_names]
  } else {
    design_names[missing_names] <- seq_along(missing_names)
  }
  names(design_list) <- design_names

  # Determine the replicate specifications for each design
  design_specs <- lapply(design_list, function(design_obj) {
    list(type = design_obj$type,
         rho = design_obj$rho,
         mse = design_obj$mse,
         scale = design_obj$scale,
         rscales = design_obj$rscales)
  })

  design_specs <- lapply(setNames(c("type", "rho", "mse", "scale", "rscales", "fpc", "fpctype"),
                                  c("type", "rho", "mse", "scale", "rscales", "fpc", "fpctype")),
                         function(spec_element) {
                           lapply(design_list, function(design_obj) {
                             design_obj[[spec_element]]
                           })
                         })

  design_specs_match <- sapply(design_specs, function(x) length(unique(x)) == 1)
  rscales_lengths <- sapply(design_list, function(x) length(x$rscales))
  diff_ncols <- length(unique(rscales_lengths)) > 1

  if (diff_ncols) {
    stop("The designs must all have the same number of columns of replicate weights.")
  }
  if (!all(design_specs_match)) {

    stop(sprintf("The following specifications differ across designs: %s",
                 paste(names(design_specs_match)[!design_specs_match],
                       collapse = ", ")))
  }

  design_specs <- lapply(design_specs, function(x) x[[1]])

  if (design_specs[['type']] %in% c("JK1", "JKn", "ACS", "successive-difference", "JK2")) {
    design_specs[['rho']] <- NULL
  }

  if (design_specs[['type']] %in% c("subbootstrap", "mrbbootstrap")) {
    design_specs[['type']] <- 'bootstrap'
  }

  known_types <- c("bootstrap", "JK1", "JK2", "JKn", "BRR", "Fay",
                   "ACS", "successive-difference")
  if (!design_specs[['type']] %in% known_types) {
    design_specs[['type']] <- "other"
  }

  # Extract the matrices of replicate weights
  # and determine which were compressed and which are combined with full-sample weights

  which_rep_wts_combined <- sapply(design_list, function(design_obj) design_obj[['combined.weights']])
  rep_wts_compressed <- sapply(design_list,
                               function(design_obj) {
                                 "repweights_compressed" %in% class(design_obj$repweights)
                               })
  any_rep_wts_compressed <- any(rep_wts_compressed)
  rep_wts_matrices <- lapply(design_list, function(design_obj) {
    if (design_obj$combined.weights) {
      as.matrix(design_obj$repweights)
    } else {
      as.matrix(design_obj$repweights) * design_obj$pweights
    }
  })


  rep_wts_datasets <- lapply(seq_along(design_list), function(design_index) {
    df <- design_list[[design_index]][['variables']]
    if (!is.null(.id)) {
      df[[.id]] <- names(design_list)[design_index]
    }
    return(df)
  })
  list_of_fullsamp_wts <- lapply(design_list, function(design_obj) unname(design_obj$pweights))

  # Determine the dimensions of each matrix of replicate weights
  rep_wts_ncols <- sapply(rep_wts_matrices, ncol, simplify = TRUE)
  rep_wts_nrows <- sapply(rep_wts_matrices, ncol, simplify = TRUE)
  if (length(unique(rep_wts_ncols)) > 1) {
    stop("The design objects must have the same number of columns of replicate weights.")
  }
  if (any(rep_wts_nrows) == 0) {
    stop("One of the designs lacks any rows of replicate weights.")
  }

  # Stack the matrices of replicate weights
  stacked_rep_wts_matrix <- Reduce(x = rep_wts_matrices, f = rbind)

  # Stack the datasets
  rep_wts_datasets_col_names <- rep_wts_datasets |>
    lapply(colnames) |> 
    unlist() |> 
    unique()
  
  stacked_datasets <- lapply(
    rep_wts_datasets, function(df) {
      data.frame(
        c(df, sapply(setdiff(rep_wts_datasets_col_names, colnames(df)),
                     function(y) NA)),
        check.names = FALSE
      )
    }
  ) |> Reduce(f = rbind)
  
  # Stack the full-sample weights
  stacked_fs_weights <- Reduce(x = list_of_fullsamp_wts, f = c)

  # Create a survey design object

  combined_svrepdesign <- survey::svrepdesign(
    data = stacked_datasets,
    repweights = stacked_rep_wts_matrix,
    weights = stacked_fs_weights,
    type = design_specs$type,
    combined.weights = TRUE,
    rho = design_specs$rho,
    scale = design_specs$scale,
    rscales = design_specs$rscales,
    fpc = design_specs$fpc, fpctype = design_specs$fpctype,
    mse = design_specs$mse
  )

  class(combined_svrepdesign) <- append(x = class(combined_svrepdesign), "svyrep.stacked")

  combined_svrepdesign <- `attr<-`(combined_svrepdesign, 'variable_for_source_design', .id)

  # Check whether the object should be a `tbl_svy` from the 'srvyr' package
  if (any(sapply(design_list, function(design) inherits(design, 'tbl_svy')))) {
    if ('package:srvyr' %in% search()) {
      combined_svrepdesign <- srvyr::as_survey_rep(
        combined_svrepdesign
      )
    }
  }

  return(combined_svrepdesign)
}
