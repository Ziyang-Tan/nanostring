source("preprocess.R")
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(readr)
library(readxl)

dataDir <- "/Users/tan/OneDrive - KI.SE/ISG nanostring/data"
infoDir <- "/Users/tan/OneDrive - KI.SE/ISG nanostring/sample info"
panelInfoPath <- "/Users/tan/nanostring/doc/New ISG extended panel - for Nanostring.xlsx"

batches <- c("EXP-24-EE2602")

panel_info <- readxl::read_excel(path = panelInfoPath, sheet = 2) %>% filter(type == "Endogenous")
columnSet <- unique(unlist(panel_info %>% select(-type))) %>% na.omit()
ISG_panel <- na.omit(panel_info$`ISG score`)

dataList_nonorm <- lapply(batches, function(x) {
    # x = batches[1]
    dataPath <- paste0(dataDir, "/", x, "_RawData.csv")
    infoPath <- paste0(infoDir, "/", x, "_RNA info.xlsx")
    file <- readr::read_csv(dataPath, show_col_types = FALSE) %>%
        rename(file_name = `File Name`, gene_name = `...2`, tag_code = `...3`)
    info <- readxl::read_excel(infoPath, sheet = 2, skip = 1) %>%
        mutate(
            `Sample info` = gsub("\\(.*$", "", `Sample info`),
            `Sample ID` = as.character(`Sample ID`),
            `RNA amount (ng)` = as.character(`RNA amount (ng)`),
            `RNA vol (µl)` = as.character(`RNA vol (µl)`),
            Visit = as.character(Visit)
        )
    # normalization
    dat <- normalization(file, do_normalization = FALSE) %>%
        bind_cols(info) %>%
        add_column(batch = rep(x, dim(info)[1]))
}) # read the raw file
dataList_norm <- lapply(batches, function(x) {
    # x = batches[1]
    dataPath <- paste0(dataDir, "/", x, "_RawData.csv")
    infoPath <- paste0(infoDir, "/", x, "_RNA info.xlsx")
    file <- readr::read_csv(dataPath, show_col_types = FALSE) %>%
        rename(file_name = `File Name`, gene_name = `...2`, tag_code = `...3`)
    info <- readxl::read_excel(infoPath, sheet = 2, skip = 1) %>%
        mutate(
            `Sample info` = gsub("\\(.*$", "", `Sample info`),
            `Sample ID` = as.character(`Sample ID`),
            `RNA amount (ng)` = as.character(`RNA amount (ng)`),
            `RNA vol (µl)` = as.character(`RNA vol (µl)`),
            Visit = as.character(Visit)
        )
    # normalization
    dat <- normalization(file, do_normalization = TRUE) %>%
        bind_cols(info) %>%
        add_column(batch = rep(x, dim(info)[1]))
}) # read the raw file

dat_nonorm <- bind_rows(dataList_nonorm) %>%
    mutate(
        Group = if_else(is.na(Group), "Patients", Group),
        batch = as.factor(batch),
        `RNA amount (ng)` = as.numeric(`RNA amount (ng)`)
    ) %>%
    filter(
        Group == "Linearity test",
        QC == "Pass",
        `RNA amount (ng)` != 0,
        !is.na(`RNA amount (ng)`)
    )
ISG_nonorm <- geomean_score(dat_nonorm, ISG_panel)

dat_norm <- bind_rows(dataList_norm) %>%
    mutate(
        Group = if_else(is.na(Group), "Patients", Group),
        batch = as.factor(batch),
        `RNA amount (ng)` = as.numeric(`RNA amount (ng)`)
    ) %>%
    filter(
        Group == "Linearity test",
        QC == "Pass",
        `RNA amount (ng)` != 0,
        !is.na(`RNA amount (ng)`)
    )
ISG_norm <- geomean_score(dat_norm, ISG_panel)

ggplot(ISG_nonorm, aes(x = `RNA amount (ng)`, y = geomean)) +
    geom_point() +
    geom_smooth() +
    ylim(c(25, 160)) +
    theme_bw()
ggsave("figures/linearity test/nonorm.pdf")
ggplot(ISG_norm, aes(x = `RNA amount (ng)`, y = geomean)) +
    geom_point() +
    geom_smooth() +
    ylim(c(25, 160)) +
    theme_bw()
ggsave("figures/linearity test/norm.pdf")
dat_2norm <- dat_norm %>%
    mutate(
        across(all_of(columnSet), ~ . / `RNA amount (ng)` * 100)
    )
ISG_2norm <- geomean_score(dat_2norm, ISG_panel)
ggplot(ISG_2norm, aes(x = `RNA amount (ng)`, y = geomean)) +
    geom_point() +
    geom_smooth() +
    ylim(c(25, 160)) +
    theme_bw()
ggsave("figures/linearity test/2step_norm.pdf")
