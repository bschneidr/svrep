#' @title Control Bootstrap Simulation Error
#'
#' @description This function estimates the number of bootstrap replicates
#' needed to reduce the simulation error of a bootstrap variance estimator to a target level,
#' where "simulation error" is defined as error caused by using only a finite number of bootstrap replicates
#' and this simulation error is measured as a simulation coefficient of variation ("simulation CV").
#'
#' @param svrepstat An estimate obtained from a bootstrap replicate survey design object,
#' with a function such as \code{svymean(..., return.replicates = TRUE)} or \code{withReplicates(..., return.replicates = TRUE)}.
#' @param target_cv A numeric value (or vector of numeric values) between 0 and 1.
#'  This is the target simulation CV for the bootstrap variance estimator.
#'
#' @section Suggested Usage:
#'  - \strong{Step 1}: Determine the largest acceptable level of simulation error for key survey estimates,
#'  where the level of simulation error is measured in terms of the simulation CV. We refer to this as the "target CV."
#'  A conventional value for the target CV is 5\%.
#'
#'  - \strong{Step 2}: Estimate key statistics of interest using a large number of bootstrap replicates (such as 5,000)
#'  and save the estimates from each bootstrap replicate. This can be conveniently done using a function
#'  from the survey package such as \code{svymean(..., return.replicates = TRUE)} or \code{withReplicates(..., return.replicates = TRUE)}.
#'
#'  - \strong{Step 3}: Use the function \code{estimate_boot_reps_for_target_cv()} to estimate the minimum number of bootstrap
#'  replicates needed to attain the target CV.
#'
#' @section Statistical Details:
#'  Unlike other replication methods such as the jackknife or balanced repeated replication,
#'  the bootstrap variance estimator's precision can always be improved by using a
#'  larger number of replicates, as the use of only a finite number of bootstrap replicates
#'  introduces simulation error to the variance estimation process.
#'  Simulation error can be measured as a "simulation coefficient of variation" (CV), which is
#'  the ratio of the standard error of a bootstrap estimator to the expectation of that bootstrap estimator,
#'  where the expectation and standard error are evaluated with respect to the bootstrapping process
#'  given the selected sample.
#'
#'  For a statistic \eqn{\hat{\theta}}, the simulation CV of the bootstrap variance estimator
#'  \eqn{v_{B}(\hat{\theta})} based on \eqn{B} replicate estimates \eqn{\hat{\theta}^{\star}_1,\dots,\hat{\theta}^{\star}_B} is defined as follows:
#'  \deqn{
#'    CV_{\star}(v_{B}(\hat{\theta})) = \frac{\sqrt{var_{\star}(v_B(\hat{\theta}))}}{E_{\star}(v_B(\hat{\theta}))} = \frac{CV_{\star}(E_2)}{\sqrt{B}}
#'  }
#'  where
#'  \deqn{
#'    E_2 = (\hat{\theta}^{\star} - \hat{\theta})^2
#'  }
#'  \deqn{
#'    CV_{\star}(E_2) = \frac{\sqrt{var_{\star}(E_2)}}{E_{\star}(E_2)}
#'  }
#'  and \eqn{var_{\star}} and \eqn{E_{\star}} are evaluated with respect to
#'  the bootstrapping process, given the selected sample.
#' \cr \cr
#'  The simulation CV, denoted \eqn{CV_{\star}(v_{B}(\hat{\theta}))}, is estimated for a given number of replicates \eqn{B}
#'  by estimating \eqn{CV_{\star}(E_2)} using observed values and dividing this by \eqn{\sqrt{B}}. If the bootstrap errors
#'  are assumed to be normally distributed, then \eqn{CV_{\star}(E_2)=\sqrt{2}} and so \eqn{CV_{\star}(v_{B}(\hat{\theta}))} would not need to be estimated.
#'  Using observed replicate estimates to estimate the simulation CV instead of assuming normality allows simulation CV to be
#'  used for a a wide array of bootstrap methods.
#'
#' @references
#'  See Section 3.3 and Section 8 of Beaumont and Patak (2012) for details and an example where the simulation CV is used
#'  to determine the number of bootstrap replicates needed for various alternative bootstrap methods in an empirical illustration.
#'
#'  Beaumont, J.-F. and Z. Patak. (2012),
#'  "On the Generalized Bootstrap for Sample Surveys with Special Attention to Poisson Sampling."
#'  \strong{International Statistical Review}, \emph{80}: 127-148. \doi{https://doi.org/10.1111/j.1751-5823.2011.00166.x}.
#'
#' @return
#' An object of class \code{sim_cv_curve}.
#' The functions \code{summary()}, \code{plot()}, and \code{print()} can be used on this object.
#' Call \code{summary()} to obtain a summary data frame.
#' The summary data frame has one row for each value of \code{target_cv}.
#' The column \code{TARGET_CV} gives the target coefficient of variation.
#' The column \code{MAX_REPS} gives the maximum number of replicates needed
#' for all of the statistics included in \code{svrepstat}. The remaining columns
#' give the number of replicates needed for each statistic.
#' 
#' Calling \code{plot()} produces a plot with 'ggplot2'.
#' @export
#' @seealso Use \code{\link[svrep]{estimate_boot_sim_cv}} to estimate the simulation CV for the number of bootstrap replicates actually used.
#' @examples
#' \donttest{
#' set.seed(2022)
#'
#' # Create an example bootstrap survey design object ----
#' data('api', package = 'survey')
#'
#' boot_design <- svydesign(id=~1,strata=~stype, weights=~pw,
#'                          data=apistrat, fpc=~fpc) |>
#'  svrep::as_bootstrap_design(replicates = 5000)
#'
#' # Calculate estimates of interest and retain estimates from each replicate ----
#'
#' estimated_means_and_proportions <- svymean(x = ~ api00 + api99 + stype, design = boot_design,
#'                                            return.replicates = TRUE)
#' custom_statistic <- withReplicates(design = boot_design,
#'                                    return.replicates = TRUE,
#'                                    theta = function(wts, data) {
#'                                       numerator   <- sum(data$api00 * wts)
#'                                       denominator <- sum(data$api99 * wts)
#'                                       statistic   <- numerator/denominator
#'                                       return(statistic)
#'                                    })
#' # Determine minimum number of bootstrap replicates needed to obtain given simulation CVs ----
#'
#'   ## For means and proportions
#'   sim_cv_curve <- estimate_boot_reps_for_target_cv(
#'     svrepstat = estimated_means_and_proportions,
#'     target_cv = c(0.01, 0.05, 0.10)
#'   )
#' 
#'   sim_cv_curve
#' 
#'   if (require('ggplot2')) {
#'     plot(sim_cv_curve)
#'   }
#'
#'   ## For custom statistic
#'   estimate_boot_reps_for_target_cv(
#'     svrepstat = custom_statistic,
#'     target_cv = c(0.01, 0.05, 0.10)
#'   )
#' }
estimate_boot_reps_for_target_cv <- function(svrepstat, target_cv = 0.05) {

  if (is.null(svrepstat$replicates)) {
    stop("`svrepstat` must have an element named 'replicates', obtained by using 'return.replicates=TRUE' in the estimation function.")
  }

  # Extract matrix of replicate estimates
  replicate_estimates <- as.matrix(svrepstat$replicates)
  colnames(replicate_estimates) <- names(coef(svrepstat))

  if (is.null(colnames(replicate_estimates))) {
    warning("The elements of `svrepstat` are unnamed. Placeholder names (STATISTIC_1, etc.) will be used instead.")
    colnames(replicate_estimates) <- sprintf("STATISTIC_%s", seq_len(ncol(replicate_estimates)))
  }

  # Estimate CV of replicate estimates around their mean
  squared_residuals <- (replicate_estimates - colMeans(replicate_estimates))^2
  mean_squared_residuals <- colMeans(squared_residuals)
  std_dev_squared_residuals <- apply(
    X = squared_residuals,
    MARGIN = 2,
    FUN = stats::sd
  )
  cv_squared_residuals <- std_dev_squared_residuals / mean_squared_residuals

  # Use Beaumont-Patak (2012) formula to estimate number of replicates
  # required to obtain a given target CV
  summary_fn <- function(target_cv_value) {

    if (any(is.na(target_cv_value)) || any(target_cv_value <= 0)) {
      stop('`target_cv_value` must be a numeric vector of positive values.')
    }

    B_values <- sapply(
      X = target_cv_value, 
      simplify = 'array',
      FUN = function(target_cv_value) {
        B <- (cv_squared_residuals / target_cv_value)^2
        B <- ceiling(B)
        return(B)
      }
    )
    if (is.matrix(B_values)) {
      B_values <- t(B_values)
    } else {
      B_values <- as.matrix(B_values)
    }
    B_values <- unname(B_values)

    # Add statistics' names to the columns of the result
    colnames(B_values) <- colnames(replicate_estimates)
    return(B_values)
  }

  B_values <- summary_fn(target_cv)


  # For a given target CV, calculate maximum number of replicates
  # needed across all the statistics
  max_B_values <- apply(X = B_values,
                        MARGIN = 1,
                        FUN = max)

  # Combine the results into a dataframe
  summary_df <- cbind(
    'TARGET_CV' = target_cv,
    'MAX_REPS' = max_B_values,
    B_values
  ) |> as.data.frame()
  
  result <- list(
    'summary' = summary_df,
    'function' = summary_fn
  )
  # Return the result
  class(result) <- union('sim_cv_curve', class(result))
  return(result)
}

#' @export
summary.sim_cv_curve <- function(object, ...) {
  object[['summary']]
}

#' @export
print.sim_cv_curve <- function(x, ...) {
  print(x[['summary']])
}

#' @export
as.data.frame.sim_cv_curve <- function(x, ...) {
  x[['summary']]
}

#' @export
plot.sim_cv_curve <- function(x, ...) {

  rlang::check_installed(c("ggplot2", "scales"), "plotting this object.")

  # Get a long-format data frame useful for plotting
  summary_df <- x[['summary']]
  stat_names <- setdiff(colnames(summary_df), c("TARGET_CV", "MAX_REPS"))
  long_data <- data.frame(
    'STATISTIC' = character(0),
    'TARGET_CV' = numeric(0),
    'REQUIRED_REPS' = numeric(0)
  )
  for (col_name in stat_names) {
    long_data <- rbind(
      long_data,
      data.frame(
        'STATISTIC' = rep(col_name, nrow(summary_df)),
        'TARGET_CV' = summary_df[['TARGET_CV']],
        'REQUIRED_REPS' = summary_df[[col_name]]
      )
    )
  }
  # Create a plot using ggplot2
  long_data |>
    ggplot2::ggplot(
      ggplot2::aes(x = .data[['TARGET_CV']],
                   y = .data[['REQUIRED_REPS']])
    ) +
    ggplot2::theme_linedraw() +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(reformulate('STATISTIC')) +
    ggplot2::geom_line(linetype = 2) +
    ggplot2::labs(
      x = 'Target Simulation CV',
      y = 'Required Bootstrap Replicates',
      title = 'Required Bootstrap Replicates for Target Simulation CV'
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_percent()
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_comma()
    )
}

#' @title Estimate Bootstrap Simulation Error
#' @description Estimates the bootstrap simulation error, expressed as a "simulation coefficient of variation" (CV).
#' @param svrepstat An estimate obtained from a bootstrap replicate survey design object,
#' with a function such as \code{svymean(..., return.replicates = TRUE)} or \code{withReplicates(..., return.replicates = TRUE)}.
#' @section Statistical Details:
#'  Unlike other replication methods such as the jackknife or balanced repeated replication,
#'  the bootstrap variance estimator's precision can always be improved by using a
#'  larger number of replicates, as the use of only a finite number of bootstrap replicates
#'  introduces simulation error to the variance estimation process.
#'  Simulation error can be measured as a "simulation coefficient of variation" (CV), which is
#'  the ratio of the standard error of a bootstrap estimator to the expectation of that bootstrap estimator,
#'  where the expectation and standard error are evaluated with respect to the bootstrapping process
#'  given the selected sample.
#'
#'  For a statistic \eqn{\hat{\theta}}, the simulation CV of the bootstrap variance estimator
#'  \eqn{v_{B}(\hat{\theta})} based on \eqn{B} replicate estimates \eqn{\hat{\theta}^{\star}_1,\dots,\hat{\theta}^{\star}_B} is defined as follows:
#'  \deqn{
#'    CV_{\star}(v_{B}(\hat{\theta})) = \frac{\sqrt{var_{\star}(v_B(\hat{\theta}))}}{E_{\star}(v_B(\hat{\theta}))} = \frac{CV_{\star}(E_2)}{\sqrt{B}}
#'  }
#'  where
#'  \deqn{
#'    E_2 = (\hat{\theta}^{\star} - \hat{\theta})^2
#'  }
#'  \deqn{
#'    CV_{\star}(E_2) = \frac{\sqrt{var_{\star}(E_2)}}{E_{\star}(E_2)}
#'  }
#'  and \eqn{var_{\star}} and \eqn{E_{\star}} are evaluated with respect to
#'  the bootstrapping process, given the selected sample.
#' \cr \cr
#'  The simulation CV, denoted \eqn{CV_{\star}(v_{B}(\hat{\theta}))}, is estimated for a given number of replicates \eqn{B}
#'  by estimating \eqn{CV_{\star}(E_2)} using observed values and dividing this by \eqn{\sqrt{B}}. If the bootstrap errors
#'  are assumed to be normally distributed, then \eqn{CV_{\star}(E_2)=\sqrt{2}} and so \eqn{CV_{\star}(v_{B}(\hat{\theta}))} would not need to be estimated.
#'  Using observed replicate estimates to estimate the simulation CV instead of assuming normality allows simulation CV to be
#'  used for a a wide array of bootstrap methods.
#'
#' @references
#'  See Section 3.3 and Section 8 of Beaumont and Patak (2012) for details and an example where the simulation CV is used
#'  to determine the number of bootstrap replicates needed for various alternative bootstrap methods in an empirical illustration.
#'
#'  Beaumont, J.-F. and Z. Patak. (2012),
#'  "On the Generalized Bootstrap for Sample Surveys with Special Attention to Poisson Sampling."
#'  \strong{International Statistical Review}, \emph{80}: 127-148. \doi{https://doi.org/10.1111/j.1751-5823.2011.00166.x}.
#'
#' @return
#' An object of class \code{sim_cv_estimate},
#' with methods \code{print()} and \code{summary()}.
#' Call \code{summary()} to get a data frame with one row for each statistic.
#' The column \code{STATISTIC} gives the name of the statistic.
#' The column \code{SIMULATION_CV} gives the estimated simulation CV of the statistic.
#' The column \code{N_REPLICATES} gives the number of bootstrap replicates.
#' @export
#' @seealso Use \code{\link[svrep]{estimate_boot_reps_for_target_cv}} to help choose the number of bootstrap replicates.
#' @examples
#' \donttest{
#' 
#' set.seed(2022)
#'
#' # Create an example bootstrap survey design object ----
#' data('api', package = 'survey')
#'
#' boot_design <- svydesign(id=~1,strata=~stype, weights=~pw,
#'                          data=apistrat, fpc=~fpc) |>
#'   as_bootstrap_design(replicates = 5000)
#'
#' # Calculate estimates of interest and retain estimates from each replicate ----
#'
#' estimated_means_and_proportions <- svymean(
#'   x = ~ api00 + api99 + stype, 
#'   design = boot_design,
#'   return.replicates = TRUE
#' )
#' 
#' custom_statistic <- withReplicates(
#'   design = boot_design,
#'   return.replicates = TRUE,
#'   theta = function(wts, data) {
#'      numerator   <- sum(data$api00 * wts)
#'      denominator <- sum(data$api99 * wts)
#'      statistic   <- numerator/denominator
#'      return(statistic)
#'   }
#' )
#' # Estimate simulation CV of bootstrap estimates ----
#'
#'   estimate_boot_sim_cv(
#'     svrepstat = estimated_means_and_proportions
#'   )
#'
#'   estimate_boot_sim_cv(
#'     svrepstat = custom_statistic
#'   )
#' }

estimate_boot_sim_cv <- function(svrepstat) {

  if (is.null(svrepstat$replicates)) {
    stop("`svrepstat` must have an element named 'replicates', obtained by using 'return.replicates=TRUE' in the estimation function.")
  }

  # Extract matrix of replicate estimates
  replicate_estimates <- as.matrix(svrepstat$replicates)
  colnames(replicate_estimates) <- names(coef(svrepstat))

  B <- nrow(replicate_estimates)

  if (is.null(colnames(replicate_estimates))) {
    warning("The elements of `svrepstat` are unnamed. Placeholder names (STATISTIC_1, etc.) will be used instead.")
    colnames(replicate_estimates) <- sprintf("STATISTIC_%s", seq_len(ncol(replicate_estimates)))
  }

  # Estimate CV of replicate estimates around their mean
  squared_residuals <- (replicate_estimates - colMeans(replicate_estimates))^2
  mean_squared_residuals <- colMeans(squared_residuals)
  std_dev_squared_residuals <- apply(
    X = squared_residuals,
    MARGIN = 2,
    FUN = stats::sd
  )

  cv_squared_residuals <- std_dev_squared_residuals / mean_squared_residuals

  simulation_cv <- cv_squared_residuals / sqrt(B)

  # Add statistics' names to the columns of the result
  names(simulation_cv) <- colnames(replicate_estimates)

  # Create a summary data frame
  summary_df <- data.frame(
    'STATISTIC' = colnames(replicate_estimates),
    'SIMULATION_CV' = simulation_cv,
    'N_REPLICATES' = B,
    row.names = NULL
  )

  result <- list(
    'summary'             = summary_df,
    'replicate_estimates' = replicate_estimates
  )
  class(result) <- union("sim_cv_estimate", class(result))

  # Return the result
  return(result)
}

#' @export
summary.sim_cv_estimate <- function(object, ...) {
  object[['summary']]
}

#' @export
print.sim_cv_estimate <- function(x, ...) {
  print(x[['summary']])
}

#' @export
as.data.frame.sim_cv_estimate <- function(x, ...) {
  x[['summary']]
}
