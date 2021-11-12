suppressPackageStartupMessages(library(survey))

# Create example data ----
set.seed(1999)
data(api)

dclus1 <- svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
dclus1$variables$response_status <- sample(x = c("Respondent", "Nonrespondent",
                                                 "Ineligible", "Unknown eligibility"),
                                           size = nrow(dclus1),
                                           replace = TRUE)
orig_rep_design <- as.svrepdesign(dclus1)

boot_design <- as.svrepdesign(dclus1, 'bootstrap', replicates = ncol(orig_rep_design$repweights$weights))
boot_design_more_cols <- as.svrepdesign(dclus1, "bootstrap",
                                        replicates = ncol(boot_design$repweights$weights) + 1)

# Adjust weights for cases with unknown eligibility
ue_adjusted_design <- redistribute_weights(design = orig_rep_design,
                                           reduce_if = response_status %in% c("Unknown eligibility"),
                                           increase_if = !response_status %in% c("Unknown eligibility"),
                                           by = c("stype", "cname"))

# Adjust weights for nonresponse
nr_adjusted_design <- redistribute_weights(design = ue_adjusted_design,
                                           reduce_if = response_status %in% c("Nonrespondent"),
                                           increase_if = response_status == "Respondent",
                                           by = c("stype", "cname"))

# Test whether different design inputs accepted ----

test_that(desc = "Able to supply list of designs in multiple formats", code = {
  expect_equal(object = svyby_repwts(rep_designs = list('orig' = orig_rep_design,
                                                        'nr-adjusted' = nr_adjusted_design),
                                     formula = ~ api00 + api99, FUN = svymean),
               expected = svyby_repwts(rep_designs = stack_replicate_designs(
                 'orig' = orig_rep_design,
                 'nr-adjusted' = nr_adjusted_design),
                 formula = ~ api00 + api99, FUN = svymean)
  )
})

# Test that point estimates and std errors match results from separate designs ----

sep_estimates <- list(
  'orig' = svymean(x = ~ api00, design = orig_rep_design),
  'nr-adjusted' = svymean(x = ~ api00, design = nr_adjusted_design)
)

stacked_design <- stack_replicate_designs('orig' = orig_rep_design,
                                          'nr-adjusted' = nr_adjusted_design,
                                          .id = "Design_Name")
combined_estimates <- svyby_repwts(formula = ~ api00, by = ~ Design_Name,
                                   FUN = svymean,
                                   rep_designs = stacked_design)

test_that("Estimates from separate and stacked designs match", code = {
  expect_equal(expected = coef(combined_estimates)[c("nr-adjusted", "orig")],
               object = c('nr-adjusted' = unname(coef(sep_estimates[['nr-adjusted']])),
                          'orig' = unname(coef(sep_estimates[['orig']]))))
  expect_equal(expected = SE(combined_estimates),
               object = c(unname(SE(sep_estimates[['nr-adjusted']])),
                          unname(SE(sep_estimates[['orig']]))))
})

# Test that covariances are estimated ----

test_that("Covariances are estimated", code = {
  expect_true(object = !is.na(vcov(combined_estimates)[1,2]) && vcov(combined_estimates)[1,2] > 0)
})

# Test that works with svycontrast() ----

test_that("Works with svycontrast()", code = {
  expect_true(object = "svystat" %in% class(svycontrast(stat = combined_estimates, c(1,-1))))
  expect_equal(
    expected = unname(coef(sep_estimates[['orig']]) - coef(sep_estimates[['nr-adjusted']])),
    object = unname(coef(svycontrast(stat = combined_estimates, c(-1,1))))
  )
})

# Matches results from svyby() ----

test_that("Matches results from svyby()", code = {
  expect_equal(
    expected = vcov(svyby(formula = ~ api00, FUN = svymean,
                          design = stacked_design, by = ~ Design_Name + stype,
                          covmat = TRUE)),
    object = vcov(svyby_repwts(formula = ~ api00, FUN = svymean,
                               rep_designs = stacked_design, by = ~ stype))
  )
})

test_that("Matches results from svyby()", code = {
  expect_equal(
    expected = vcov(svyby(formula = ~ api00, FUN = svymean,
                          design = stacked_design, by = list('Design_Name'= stacked_design$variables$Design_Name,
                                                             'stype' = stacked_design$variables$stype),
                          covmat = TRUE)),
    object = vcov(svyby_repwts(formula = ~ api00, FUN = svymean,
                               rep_designs = stacked_design,
                               by = list('stype' = stacked_design$variables$stype)))
  )
})

# Test that svyby() arguments are passed correctly ----

test_that("Arguments are correctly passed to svyby()", code = {
  expect_equal(
    expected = `attr<-`(svyby(formula = ~ api00, FUN = svymean,
                          design = stacked_design, by = ~ Design_Name + stype,
                          covmat = TRUE, deff = TRUE, vartype = c("ci", "cv"),
                          keep.names = FALSE, na.rm.by = TRUE, na.rm.all = TRUE),
                        'call', NULL),
    object = `attr<-`(svyby_repwts(formula = ~ api00, FUN = svymean,
                               rep_designs = stacked_design, by = ~ stype,
                               deff = TRUE, vartype = c("ci", "cv"),
                               keep.names = FALSE, na.rm.by = TRUE, na.rm.all = TRUE),
                      'call', NULL)
  )
})

# Check for error messages ----

test_that("Error message for invalid inputs", code = {
  expect_error(
    object = {
      svyby_repwts(rep_designs = orig_rep_design,
                   formula = ~ api00, FUN = svymean)
    }, regexp = "must be.+list of.+or.+stack"
  )
})
