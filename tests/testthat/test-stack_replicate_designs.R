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

# Test for informative error that designs are conformable ----

test_that("Informative error message for different types of designs", code = {
 expect_error(stack_replicate_designs(orig_rep_design, boot_design),
              regexp = "specifications differ")
})
test_that("Informative error message for non-comformable designs", {
  expect_error(stack_replicate_designs(boot_design, boot_design_more_cols),
               regexp = "must all have the same number of columns ")
})

# Test that user can supply multiple types of arguments ----

test_that("Can supply list of designs in multiple formats", code = {
  expect_equal(stack_replicate_designs(orig_rep_design, ue_adjusted_design),
               stack_replicate_designs(list('orig_rep_design' = orig_rep_design,
                                            'ue_adjusted_design' = ue_adjusted_design)))
  expect_equal(stack_replicate_designs(orig = orig_rep_design, adjusted = ue_adjusted_design),
               stack_replicate_designs(list('orig' = orig_rep_design,
                                            'adjusted' = ue_adjusted_design)))
})
