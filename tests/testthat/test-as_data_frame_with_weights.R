suppressWarnings({
  suppressPackageStartupMessages({
    library(survey)
  })
})

data("lou_vax_survey", package = 'svrep')

# Create a survey design object
  survey_design <- svydesign(data = lou_vax_survey,
                             weights = ~ SAMPLING_WEIGHT,
                             ids = ~ 1)

  rep_survey_design <- as.svrepdesign(survey_design,
                                      type = "boot",
                                      replicates = 10)

  # Adjust the weights for nonresponse
  nr_adjusted_design <- redistribute_weights(
    design = rep_survey_design,
    reduce_if = RESPONSE_STATUS == "Nonrespondent",
    increase_if = RESPONSE_STATUS == "Respondent",
    by = c("RACE_ETHNICITY", "EDUC_ATTAINMENT")
  )

  # Save the survey design object as a data frame
  nr_adjusted_data <- as_data_frame_with_weights(
    nr_adjusted_design,
    full_wgt_name = "NR_ADJUSTED_WGT",
    rep_wgt_prefix = "NR_ADJUSTED_REP_WGT_"
  )

# Test that works for non-replicate design ----

  test_that("Extracts full-sample weights for non-replicate design", {
    expect_named(
      object = as_data_frame_with_weights(survey_design, full_wgt_name = "FULL_SAMPLE_WGT"),
      expected = c(colnames(survey_design), "FULL_SAMPLE_WGT")
    )
    expect_named(
      object = suppressMessages(
        as_data_frame_with_weights(survey_design, full_wgt_name = NULL)
      ),
      expected = c(colnames(survey_design), "FULL_SAMPLE_WGT")
    )
  })

# Test that works for replicate design ----

  test_that("Extracts full-sample weights for replicate designs", {

    suppressMessages(
      rep_df_names <- as_data_frame_with_weights(rep_survey_design,
                                                 full_wgt_name = "FS_WGT") |>
        colnames()
    )

    expect_equal(
      object = rep_df_names[rep_df_names == "FS_WGT"],
      expected = "FS_WGT"
    )

    suppressMessages(
      rep_df_names <- as_data_frame_with_weights(rep_survey_design,
                                                 full_wgt_name = NULL) |>
        colnames()
    )

    expect_equal(
      object = rep_df_names[rep_df_names == "FULL_SAMPLE_WGT"],
      expected = "FULL_SAMPLE_WGT"
    )

  })

  test_that("Extracts replicate weights from replicate designs", {

    suppressMessages(
      rep_df_names <- as_data_frame_with_weights(rep_survey_design,
                                                 full_wgt_name = "FS_WGT",
                                                 rep_wgt_prefix = "REPLICATE_WGT_") |>
        colnames()
    )

    expect_equal(
      object = rep_df_names[grepl(x = rep_df_names, pattern = "REPLICATE_WGT")],
      expected = paste0("REPLICATE_WGT_", 1:10)
    )

    suppressMessages(
      rep_df_names <- as_data_frame_with_weights(rep_survey_design,
                                                 full_wgt_name = NULL,
                                                 rep_wgt_prefix = "RW_") |>
        colnames()
    )

    expect_equal(
      object = rep_df_names[grepl(x = rep_df_names, pattern = "RW_")],
      expected = paste0("RW_", 1:10)
    )

  })

# Handles name conflicts ----

  test_that("Informative warnings or errors about name conflicts", {

    rep_survey_design2 <- rep_survey_design
    rep_survey_design2$variables[['FS_WGT']] <- rep_survey_design$variables$SAMPLING_WEIGHT

    expect_warning(
      object = {
        as_data_frame_with_weights(rep_survey_design2,
                                   full_wgt_name = "FS_WGT",
                                   rep_wgt_prefix = "REP_WGT_")
      },
      regexp = "already a column named `FS_WGT`"
    )

    rep_survey_design2$variables <- cbind(
      rep_survey_design2$variables,
      weights(rep_survey_design2, type = "analysis") |>
        as.data.frame() |>
        setNames(paste0("REPWT_", 1:10))
    )

    expect_error(
      object = {
        as_data_frame_with_weights(rep_survey_design2,
                                   full_wgt_name = "FULL_WEIGHT",
                                   rep_wgt_prefix = "REPWT_")
      },
      regexp = "already contain columns.+REPWT_2"
    )

  })
