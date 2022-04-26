library(tidycensus)
library(dplyr)

library(survey)
library(srvyr)
library(svrep)

set.seed(2014)

# Load the Census API key ----
  census_api_key(
    readLines("census-api-key")
  )

# Download PUMS data for Louisville ----
  lou_pums_data <- get_pums(
    variables = c("SEX", "RAC1P", "HISP", "AGEP", "SCHL"),
    survey = "acs5",
    year = 2019,
    state = "KY",
    puma = paste0("0170", 1:6),
    recode = TRUE,
    rep_weights = "person"
  )

# Add derived variables to use for calibration ----
  lou_pums_data <- lou_pums_data |>
    mutate(
      RACE_ETHNICITY = case_when(
        HISP_label != "Not Spanish/Hispanic/Latino" ~ "Hispanic or Latino",
        !as.character(RAC1P_label) %in% c("White alone", "Black or African American alone",
                            "Hispanic or Latino") ~ "Other Race, not Hispanic or Latino",
        TRUE ~ paste0(as.character(RAC1P_label), ", not Hispanic or Latino")
      ),
      EDUC_ATTAINMENT = case_when(
        SCHL_label %in% c("Associate's degree", "Bachelor's degree", "Master's degree",
                          "Professional degree beyond a bachelor's degree",
                          "Doctorate degree") ~ "High school or beyond",
        TRUE ~ "Less than high school"
      )
    )

  lou_pums_data <- lou_pums_data |>
    mutate(CONTROL_CATEGORY = interaction(RACE_ETHNICITY, SEX_label,
                                          EDUC_ATTAINMENT, sep = "|"))

# Convert to survey design object ----
  lou_rep_design <- tidycensus::to_survey(
    df = lou_pums_data,
    type = "person"
  )

# Generate population vaccination rates ----

  population_counts <- lou_rep_design |>
    filter(AGEP >= 18) |>
    survey_count(CONTROL_CATEGORY, name = "Population_Size") |>
    mutate(categories = strsplit(x = as.character(CONTROL_CATEGORY),
                                 split = "|", fixed = TRUE),
           RACE_ETHNICITY = sapply(categories, function(x) x[1]),
           SEX = sapply(categories, function(x) x[2]),
           EDUC_ATTAINMENT = sapply(categories, function(x) x[3])) |>
    select(-categories)

  pop_vax_rates <- population_counts |>
    mutate(
      VAX_RATE = case_when(
        grepl(x = RACE_ETHNICITY, "Black or African American alone") ~ 0.55,
        grepl(x = RACE_ETHNICITY, "White alone") ~ 0.58,
        grepl(x = RACE_ETHNICITY, "^Hispanic or Latino") ~ 0.48,
        TRUE ~ 0.50
      )
    ) |>
    mutate(
      VAX_ODDS = case_when(
        SEX == "Female" ~ VAX_RATE/(1-VAX_RATE) * 1.25,
        SEX == "Male" ~ VAX_RATE/(1-VAX_RATE) * 0.75
      )
    ) |>
    mutate(
      VAX_ODDS = case_when(
        EDUC_ATTAINMENT == "High school or beyond" ~ VAX_ODDS * 1.33,
        EDUC_ATTAINMENT == "Less than high school" ~ VAX_ODDS * 0.67
      )
    ) |>
    mutate(VAX_RATE = VAX_ODDS / (1 + VAX_ODDS))

# Generate response propensity variable ----

  pop_vax_rates <- pop_vax_rates |>
    mutate(
      RESP_PROPENSITY = case_when(
        grepl(x = RACE_ETHNICITY, "Black or African American alone") ~ 0.45,
        grepl(x = RACE_ETHNICITY, "White alone") ~ 0.48,
        grepl(x = RACE_ETHNICITY, "^Hispanic or Latino") ~ 0.4,
        TRUE ~ 0.45
      )
    ) |>
    mutate(
      RESP_ODDS = case_when(
        SEX == "Female" ~ RESP_PROPENSITY/(1-RESP_PROPENSITY) * 1.07,
        SEX == "Male" ~ RESP_PROPENSITY/(1-RESP_PROPENSITY) * 0.93
      )
    ) |>
    mutate(
      RESP_ODDS = case_when(
        EDUC_ATTAINMENT == "High school or beyond" ~ RESP_ODDS * 1.5,
        EDUC_ATTAINMENT == "Less than high school" ~ RESP_ODDS * 1
      )
    ) |>
    mutate(RESP_PROPENSITY = RESP_ODDS / (1 + RESP_ODDS))

# Draw simple random sample for vaccination survey ----

  lou_vax_survey <- pop_vax_rates |>
    mutate(SAMPLING_WEIGHT = sum(Population_Size) / 1000) |>
    select(RACE_ETHNICITY, SEX, EDUC_ATTAINMENT,
           Population_Size, VAX_RATE,
           RESP_PROPENSITY, SAMPLING_WEIGHT) |>
    sample_n(size = 1000, weight = Population_Size, replace = TRUE) |>
    select(-Population_Size)

  ##_ Generate vaccination status and response status
  lou_vax_survey <- lou_vax_survey |>
    mutate(VAX_STATUS = sapply(VAX_RATE, FUN = function(vax_prob) {
      ifelse(vax_prob > runif(n = 1), "Vaccinated", "Unvaccinated")
    })) |>
    mutate(VAX_STATUS = sapply(VAX_RATE, FUN = function(vax_prob) {
      ifelse(vax_prob > runif(n = 1), "Vaccinated", "Unvaccinated")
    })) |>
    mutate(RESPONSE_STATUS = sapply(RESP_PROPENSITY, FUN = function(resp_prob) {
      ifelse(resp_prob > runif(n = 1), "Respondent", "Nonrespondent")
    })) |>
    select(-VAX_RATE, -RESP_PROPENSITY)

  lou_vax_survey[['VAX_STATUS']] <- ifelse(
    lou_vax_survey[['RESPONSE_STATUS']] == "Nonrespondent", NA_character_,
    lou_vax_survey[['VAX_STATUS']]
  )

  ##_ Rearrange columns and shuffles rows
  lou_vax_survey <- lou_vax_survey |>
    select(RESPONSE_STATUS, SAMPLING_WEIGHT, everything()) |>
    sample_n(size = 1000, replace = FALSE)

  lou_vax_survey <- lou_vax_survey[,c(setdiff(colnames(lou_vax_survey),
                                              "SAMPLING_WEIGHT"),
                                      "SAMPLING_WEIGHT")]

# Estimate control totals ----

  ##_ For post-stratification

  poststratification_totals <- lou_rep_design |>
    filter(AGEP >= 18) |>
    svytotal(x = ~ CONTROL_CATEGORY)

  vcov_poststratification_totals <- vcov(poststratification_totals) |> as.matrix()

  poststratification_totals <- coef(poststratification_totals)
  names(poststratification_totals) <- gsub(
    x = names(poststratification_totals),
    pattern = "CONTROL_CATEGORY", replacement = ""
  )
  colnames(vcov_poststratification_totals) <- rownames(vcov_poststratification_totals) <- names(
    poststratification_totals
  )

  attributes(vcov_poststratification_totals)$means <- NULL

  lou_vax_survey_poststrat_totals <- list(
    'estimates' = poststratification_totals,
    'variance-covariance' = vcov_poststratification_totals
  )

  ##_ For raking

  raking_totals <- lou_rep_design |>
    filter(AGEP >= 18) |>
    svytotal(x = ~ RACE_ETHNICITY + SEX_label + EDUC_ATTAINMENT)

  vcov_raking_totals <- vcov(raking_totals) |> as.matrix()

  raking_totals <- coef(raking_totals)
  names(raking_totals) <- gsub(
    x = names(raking_totals),
    pattern = "(RACE_ETHNICITY|SEX_label|EDUC_ATTAINMENT)", replacement = ""
  )
  colnames(vcov_raking_totals) <- rownames(vcov_raking_totals) <- names(
    raking_totals
  )

  attributes(vcov_raking_totals)$means <- NULL

  lou_vax_survey_raking_totals <- list(
    'estimates' = raking_totals,
    'variance-covariance' = vcov_raking_totals
  )

  lou_vax_survey_control_totals <- list(
    'poststratification' = lou_vax_survey_poststrat_totals,
    'raking' = lou_vax_survey_raking_totals
  )

# Save the dataset(s) of interest ----

  usethis::use_data(lou_vax_survey,
                    lou_vax_survey_control_totals,
                    overwrite = TRUE)
