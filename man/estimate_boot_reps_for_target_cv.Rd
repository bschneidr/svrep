% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap_helpers.R
\name{estimate_boot_reps_for_target_cv}
\alias{estimate_boot_reps_for_target_cv}
\title{Control Bootstrap Simulation Error}
\usage{
estimate_boot_reps_for_target_cv(svrepstat, target_cv = 0.05)
}
\arguments{
\item{svrepstat}{An estimate obtained from a bootstrap replicate survey design object,
with a function such as \code{svymean(..., return.replicates = TRUE)} or \code{withReplicates(..., return.replicates = TRUE)}.}

\item{target_cv}{A numeric value (or vector of numeric values) between 0 and 1.
This is the target simulation CV for the bootstrap variance estimator.}
}
\value{
A data frame with one row for each value of \code{target_cv}.
The column \code{TARGET_CV} gives the target coefficient of variation.
The column \code{MAX_REPS} gives the maximum number of replicates needed
for all of the statistics included in \code{svrepstat}. The remaining columns
give the number of replicates needed for each statistic.
}
\description{
This function estimates the number of bootstrap replicates
needed to reduce the simulation error of a bootstrap variance estimator to a target level,
where "simulation error" is defined as error caused by using only a finite number of bootstrap replicates
and this simulation error is measured as a simulation coefficient of variation ("simulation CV").
}
\section{Suggested Usage}{

 - \strong{Step 1}: Determine the largest acceptable level of simulation error for key survey estimates,
 where the level of simulation error is measured in terms of the simulation CV. We refer to this as the "target CV."
 A conventional value for the target CV is 5\%.

 - \strong{Step 2}: Estimate key statistics of interest using a large number of bootstrap replicates (such as 5,000)
 and save the estimates from each bootstrap replicate. This can be conveniently done using a function
 from the survey package such as \code{svymean(..., return.replicates = TRUE)} or \code{withReplicates(..., return.replicates = TRUE)}.

 - \strong{Step 3}: Use the function \code{estimate_boot_reps_for_target_cv()} to estimate the minimum number of bootstrap
 replicates needed to attain the target CV.
}

\section{Statistical Details}{

 Unlike other replication methods such as the jackknife or balanced repeated replication,
 the bootstrap variance estimator's precision can always be improved by using a
 larger number of replicates, as the use of only a finite number of bootstrap replicates
 introduces simulation error to the variance estimation process.
 Simulation error can be measured as a "simulation coefficient of variation" (CV), which is
 the ratio of the standard error of a bootstrap estimator to the expectation of that bootstrap estimator,
 where the expectation and standard error are evaluated with respect to the bootstrapping process
 given the selected sample.

 For a statistic \eqn{\hat{\theta}}, the simulation CV of the bootstrap variance estimator
 \eqn{v_{B}(\hat{\theta})} based on \eqn{B} replicate estimates \eqn{\hat{\theta}^{\star}_1,\dots,\hat{\theta}^{\star}_B} is defined as follows:
 \deqn{
   CV_{\star}(v_{B}(\hat{\theta})) = \frac{\sqrt{var_{\star}(v_B(\hat{\theta}))}}{E_{\star}(v_B(\hat{\theta}))} = \frac{CV_{\star}(E_2)}{\sqrt{B}}
 }
 where
 \deqn{
   E_2 = (\hat{\theta}^{\star} - \hat{\theta})^2
 }
 \deqn{
   CV_{\star}(E_2) = \frac{\sqrt{var_{\star}(E_2)}}{E_{\star}(E_2)}
 }
 and \eqn{var_{\star}} and \eqn{E_{\star}} are evaluated with respect to
 the bootstrapping process, given the selected sample.
\cr \cr
 The simulation CV, denoted \eqn{CV_{\star}(v_{B}(\hat{\theta}))}, is estimated for a given number of replicates \eqn{B}
 by estimating \eqn{CV_{\star}(E_2)} using observed values and dividing this by \eqn{\sqrt{B}}. If the bootstrap errors
 are assumed to be normally distributed, then \eqn{CV_{\star}(E_2)=\sqrt{2}} and so \eqn{CV_{\star}(v_{B}(\hat{\theta}))} would not need to be estimated.
 Using observed replicate estimates to estimate the simulation CV instead of assuming normality allows simulation CV to be
 used for a a wide array of bootstrap methods.
}

\examples{
\donttest{
set.seed(2022)

# Create an example bootstrap survey design object ----
library(survey)
data('api', package = 'survey')

boot_design <- svydesign(id=~1,strata=~stype, weights=~pw,
                         data=apistrat, fpc=~fpc) |>
 svrep::as_bootstrap_design(replicates = 5000)

# Calculate estimates of interest and retain estimates from each replicate ----

estimated_means_and_proportions <- svymean(x = ~ api00 + api99 + stype, design = boot_design,
                                           return.replicates = TRUE)
custom_statistic <- withReplicates(design = boot_design,
                                   return.replicates = TRUE,
                                   theta = function(wts, data) {
                                      numerator <- sum(data$api00 * wts)
                                      denominator <- sum(data$api99 * wts)
                                      statistic <- numerator/denominator
                                      return(statistic)
                                   })
# Determine minimum number of bootstrap replicates needed to obtain given simulation CVs ----

  estimate_boot_reps_for_target_cv(
    svrepstat = estimated_means_and_proportions,
    target_cv = c(0.01, 0.05, 0.10)
  )

  estimate_boot_reps_for_target_cv(
    svrepstat = custom_statistic,
    target_cv = c(0.01, 0.05, 0.10)
  )
}
}
\references{
See Section 3.3 and Section 8 of Beaumont and Patak (2012) for details and an example where the simulation CV is used
 to determine the number of bootstrap replicates needed for various alternative bootstrap methods in an empirical illustration.

 Beaumont, J.-F. and Z. Patak. (2012),
 "On the Generalized Bootstrap for Sample Surveys with Special Attention to Poisson Sampling."
 \strong{International Statistical Review}, \emph{80}: 127-148. \doi{https://doi.org/10.1111/j.1751-5823.2011.00166.x}.
}
\seealso{
Use \code{\link[svrep]{estimate_boot_sim_cv}} to estimate the simulation CV for the number of bootstrap replicates actually used.
}
