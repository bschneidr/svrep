library(tidycensus)
library(dplyr)

library(survey)
library(srvyr)
library(svrep)

set.seed(2014)

# Documentation: https://www.imls.gov/sites/default/files/2022-07/2020_pls_data_file_documentation.pdf

# Download the ZIP file from "imls.gov" ----

  pls_fy2020_spss_url <- "https://www.imls.gov/sites/default/files/2022-07/pls_fy2020_spss.zip"
  temp_zip_file <- tempfile()
  download.file(url = pls_fy2020_spss_url, destfile = temp_zip_file)
  temp_unzipped_folder <- tempdir(check = TRUE)
  unzip(zipfile = temp_zip_file, exdir = temp_unzipped_folder)

# Clean up the outlets-level file ----
  pls_fy20_outlets <- file.path(temp_unzipped_folder, "PLS_FY20_Outlet_pud20i.sav") |>
    haven::read_sav(encoding = "latin1")

  ##_ Fix issue with character encoding
  pls_fy20_outlets <- pls_fy20_outlets |>
    mutate(
      across(where(is.character), function(x) {
        iconv(x, from = "UTF-8", to = "latin1") |>
          `Encoding<-`("UTF-8")
      })
    )

  for (col_name in colnames(pls_fy20_outlets)) {
    if (!is.null(labelled::val_labels(pls_fy20_outlets[[col_name]])))
    names(labelled::val_labels(pls_fy20_outlets[[col_name]])) <- pls_fy20_outlets[[col_name]] |>
      labelled::val_labels() |>
      names() |>
      iconv(from = "UTF-8", to = "latin1") |>
      `Encoding<-`("UTF-8")
  }

# Clean up the system-level file ----

  pls_fy20_systems <- file.path(temp_unzipped_folder, "PLS_FY20_AE_pud20i.sav") |>
    haven::read_sav(encoding = "latin1")

  ##_ Fix issue with character encoding
  pls_fy20_systems <- pls_fy20_systems |>
    mutate(
      across(where(is.character), function(x) {
        iconv(x, from = "UTF-8", to = "latin1") |>
          `Encoding<-`("UTF-8")
      })
    )

  for (col_name in colnames(pls_fy20_systems)) {
    if (!is.null(labelled::val_labels(pls_fy20_systems[[col_name]])))
      names(labelled::val_labels(pls_fy20_systems[[col_name]])) <- pls_fy20_systems[[col_name]] |>
        labelled::val_labels() |>
        names() |>
        iconv(from = "UTF-8", to = "latin1") |>
        `Encoding<-`("UTF-8")
  }

# Prepare data files ----

  library_outlets_census <- pls_fy20_outlets |>
    select(STABR, FSCSKEY, FSCS_SEQ, C_FSCS, LIBID, LIBNAME,
           CITY, ZIP, CNTY, LONGITUD, LATITUDE,
           CENTRACT, CENBLOCK, CBSA, MICROF,
           HOURS, WKS_OPEN,
           C19WKSCL, C19WKSLO)

  library_census <-  pls_fy20_systems |>
    select(
      # Unique identifiers
      FSCSKEY, C_FSCS, LIBID, LIBNAME,
      # Whether closed or temporarily closed
      STATSTRU,
      # Reporting period information
      STARTDAT, ENDDATE,
      # Number of branches and central libraries
      BRANLIB, CENTLIB,
      # Location
      STABR,
       #CITY, ZIP, CNTY,
      LONGITUD, LATITUDE,
       #CENTRACT, CENBLOCK,
      CBSA, MICROF,
      # Circulation
      TOTCIR, PHYSCIR, ELMATCIR,
      # Visits
      VISITS,
      # Number of registered users
      REGBOR,
      # Total paid employees and librarians
      TOTSTAFF, LIBRARIA,
      # Total operating expenses/revenue
      TOTOPEXP, TOTINCM
    )

  library_census <- library_census |>
    mutate(
      RESPONSE_STATUS = case_when(
        STATSTRU == "25" ~ "Survey Nonrespondent",
        STATSTRU == "03" ~ "Closed",
        STATSTRU == "23" ~ "Temporarily Closed",
        TRUE ~ "Survey Respondent"
      ) |> labelled::`var_label<-`(
        "Public Library Survey response status derived from `STATSTRU` in public-use data file."
      )
    ) |> select(-STATSTRU)

  ##_ Remove `haven_labelled` class from vectors

  library_census <- library_census |>
    mutate(across(one_of(c("STARTDAT", "ENDDATE",
                           "BRANLIB", "CENTLIB",
                           "TOTCIR", "PHYSCIR", "ELMATCIR",
                           "VISITS", "REGBOR",
                           "TOTSTAFF", "LIBRARIA",
                           "TOTOPEXP", "TOTINCM")),
                  haven::zap_labels)) |>
    mutate(C_FSCS = factor(C_FSCS, levels = c("Y", "N"),
                           labels = c("Yes", "No")) |>
             labelled::`var_label<-`(
               labelled::var_label(C_FSCS)
             ))

# Draw stratified systematic sample ----

  library_stsys_sample <- library_census |>
    arrange(STABR, TOTCIR) |>
    group_by(STABR) |>
    mutate(
      SAMPLING_STRATUM = stringr::str_pad(cur_group_id(), width = 3,
                                          side = "left", pad = "0"),
      SAMPLING_SORT_ORDER = paste0(SAMPLING_STRATUM, "_",
                                   stringr::str_pad(row_number(),
                                                    width = 4,
                                                    side = "left", pad = "0")),
      STRATUM_POP_SIZE = n(),
      SAMPLING_PROB = case_when(
        n() == 2 ~ 1,
        n() >= 100 ~ sampling::inclusionprobabilities(rep(1, times = n()),
                                                      n = ceiling(n()/50)),
        TRUE ~ sampling::inclusionprobabilities(rep(1, times = n()),
                                                n = 2)
      ),
      SAMPLING_INDICATOR = sampling::UPsystematic(SAMPLING_PROB)
    ) |>
    ungroup() |>
    filter(SAMPLING_INDICATOR == 1) |>
    select(-SAMPLING_INDICATOR)

# Draw two-stage sample, first stage PPS, second stage SRSWOR ----

  first_stage_pps_sample <- library_census |>
    group_by(STABR, CBSA) |>
    summarize(PSU_ID = cur_group_id(),
              TOTCIR = sum(TOTCIR, na.rm = TRUE),
              N_LIBRARIES = n_distinct(FSCSKEY)) |>
    ungroup() |>
    mutate(
      PSU_POP_SIZE = n(),
      PSU_SAMPLING_PROB = sampling::inclusionprobabilities(a = TOTCIR, n = 100)
    ) |>
    mutate(
      FIRST_STAGE_SAMPLING_INDICATOR = sampling::UPbrewer(PSU_SAMPLING_PROB)
    ) |>
    filter(FIRST_STAGE_SAMPLING_INDICATOR == 1) |>
    select(STABR, CBSA, PSU_ID, PSU_SAMPLING_PROB, PSU_POP_SIZE)

  second_stage_sample <- library_census |>
    inner_join(y = first_stage_pps_sample,
               by = c("STABR", "CBSA")) |>
    arrange(PSU_ID, FSCSKEY)  |>
    mutate(SSU_ID = row_number()) |>
    group_by(PSU_ID) |>
    mutate(
      SSU_POP_SIZE = n(),
      SSU_SAMPLING_PROB = case_when(n() == 1 ~ 1,
                                    TRUE ~ 2/n()),
      SSU_SAMPLING_INDICATOR = sampling::srswor(n = case_when(n() == 1 ~ 1,
                                                              TRUE ~ 2),
                                                N = n())
    ) |>
    ungroup() |>
    filter(SSU_SAMPLING_INDICATOR == 1) |> select(-SSU_SAMPLING_INDICATOR) |>
    mutate(SAMPLING_PROB = PSU_SAMPLING_PROB * SSU_SAMPLING_PROB)

  library_multistage_sample <- second_stage_sample

# Rearrange columns for ease of use ----

  library_multistage_sample <- library_multistage_sample |>
    relocate(FSCSKEY, C_FSCS,
             LIBID, LIBNAME,
             SAMPLING_PROB,
             PSU_ID, SSU_ID,
             PSU_SAMPLING_PROB, PSU_POP_SIZE,
             SSU_SAMPLING_PROB, SSU_POP_SIZE,
             RESPONSE_STATUS,
             .before = 1)

  library_stsys_sample <- library_stsys_sample |>
    relocate(FSCSKEY, C_FSCS,
             LIBID, LIBNAME,
             SAMPLING_STRATUM,
             SAMPLING_SORT_ORDER,
             STRATUM_POP_SIZE,
             SAMPLING_PROB,
             RESPONSE_STATUS,
             .before = 1)

# Save the datasets to the appropriate package folder ----
  usethis::use_data(library_census, overwrite = TRUE,
                    compress = "xz")
  usethis::use_data(library_stsys_sample, overwrite = TRUE,
                    compress = "xz")
  usethis::use_data(library_multistage_sample, overwrite = TRUE,
                    compress = "xz")
