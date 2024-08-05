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
sampleInfoPath <- "/Users/tan/Library/CloudStorage/OneDrive-KI.SE/ISG nanostring/sample info/sample and physician info.xlsx"
# metaInfoPath <- "/Users/tan/Library/CloudStorage/OneDrive-KI.SE/ISG nanostring/sample info/physician and batch info.xlsx"
sampleInfo2Path <- "/Users/tan/Library/CloudStorage/OneDrive-KI.SE/ISG nanostring/sample info/ISG IDs.xlsx"

# new panel
batches <- c(
    "EXP-21-DN5206", "EXP-21-DN5207", "EXP-21-DN5208", "EXP-21-DN5209", "EXP-21-DN5210", "EXP-21-DN5211", "EXP-21-DN5212", "EXP-21-DN5214",
    "EXP-22-DN5215", "EXP-22-DN5218", "EXP-22-DN5219", "EXP-22-DN5221", "EXP-22-DN5227", "EXP-22-run20", "EXP-22-run21", "EXP-23-DN5230",
    "EXP-23-DN5231", "EXP-23-DN5232", "EXP-23-DN5233", "EXP-23-DN5234", "EXP-23-DN5236", "EXP-24-DN5238", "EXP-24-EE2600", "EXP-24-EE2601",
    "EXP-24-EE2602", "EXP-24-EE2603", "EXP-24-EE2604", "EXP-24-EE2605"
)

panel_info <- readxl::read_excel(path = panelInfoPath, sheet = 2) %>% filter(type == "Endogenous")

sample_info_all <- left_join(
    readxl::read_excel(sampleInfoPath, sheet = 1),
    readxl::read_excel(sampleInfoPath, sheet = 2),
    by = "Referring physician"
) %>%
    mutate(`Referring physician` = if_else(is.na(`Referring physician`), "unknown", `Referring physician`)) %>%
    left_join(readxl::read_excel(sampleInfo2Path, sheet = 1) %>% select(`Patient ID`, `Personal_number`), by = "Patient ID")

dataList <- lapply(batches, function(x) {
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
    dat <- normalization(file) %>%
        bind_cols(info) %>%
        add_column(batch = rep(x, dim(info)[1]))
}) # read the raw file
dat <- bind_rows(dataList) %>%
    mutate(
        Group = if_else(is.na(Group), "Patients", Group),
        batch = as.factor(batch),
        `RNA amount (ng)` = as.numeric(`RNA amount (ng)`)
    ) %>%
    filter(
        !Group %in% c("Background", "Linearity test"),
        QC == "Pass",
        `RNA amount (ng)` != 0,
        !is.na(`RNA amount (ng)`)
    )
columnSet <- unique(unlist(panel_info %>% select(-type))) %>% na.omit()
dat <- dat %>%
    mutate(
        control = grepl("Healthy relative", Group, ignore.case = TRUE) |
            grepl("Healthy control", Group, ignore.case = TRUE),
        # across(all_of(columnSet), ~ . / `RNA amount (ng)` * 100),
        unique_ID = paste0(`Subject ID`, "_", Visit)
    )

# deal with replicates, remove for now
dat <- rbind(
    dat %>% filter(is.na(`Subject ID`)),
    dat %>% filter(!is.na(`Subject ID`)) %>% distinct(`Subject ID`, Visit, .keep_all = TRUE)
)

# with the intrinsic control, batch correction is not needed?
# write the normalized data table
write_csv(dat, "data/normalized counts latest EXP-24-EE2605.csv")


# load panel info
ISG_panel <- na.omit(panel_info$`ISG score`)
NFkb_panel <- na.omit(panel_info$`NF-kB score`)
IFNg_panel <- na.omit(panel_info$`IFN-g Score`)

# calculate scores
ISG <- geomean_score(dat, ISG_panel) %>%
    add_column(zscore = (zscore_score(dat, ISG_panel))$zscore)
datFiltered <- dat %>% filter(!batch %in% c(
    "EXP-21-DN5206", "EXP-21-DN5207", "EXP-21-DN5208", "EXP-21-DN5209",
    "EXP-21-DN5210", "EXP-21-DN5211", "EXP-21-DN5212", "EXP-21-DN5214"
)) # for those batches the panel doesn't have NFkB or IFNg genes
NFkb <- geomean_score(datFiltered, NFkb_panel) %>%
    add_column(zscore = (zscore_score(datFiltered, NFkb_panel))$zscore)
IFNg <- geomean_score(datFiltered, IFNg_panel) %>%
    add_column(zscore = (zscore_score(datFiltered, IFNg_panel))$zscore)


# batch export ISG report

exportDir <- "/Users/tan/Library/CloudStorage/OneDrive-KI.SE/ISG nanostring/ISG reports"
dir.create(exportDir, showWarnings = FALSE)
# for (folder in unique(sample_info_all$`Experiment batch`)) {
for (folder in c("EXP-24-EE2605")) {
    expath <- file.path(exportDir, folder)
    dir.create(expath, showWarnings = FALSE)
    report_families <- sample_info_all %>%
        filter(
            `Experiment batch` == folder
        ) %>%
        mutate(`Patient ID` = sub("M|F|'S'", "", `Patient ID`)) %>%
        select(`Patient ID`) %>%
        distinct() %>%
        unlist()
    lapply(report_families, function(cur_sample) {
        rmarkdown::render("ISG_report_template.Rmd",
            params = list(
                ISG = ISG,
                NFkb = NFkb,
                IFNg = IFNg,
                panel_info = panel_info,
                cur_sample = cur_sample,
                sample_info = sample_info_all %>% filter(`Patient ID` %in% c(cur_sample, paste0(cur_sample, "M"), paste0(cur_sample, "F")))
            ),
            output_file = paste0(
                expath, "/", cur_sample, "_",
                sample_info_all %>% filter(`Patient ID` == cur_sample) %>% select(Personal_number) %>% distinct() %>% unlist(),
                ".pdf"
            )
        )
    })
}

# test
cur_sample <- "ISG-6"
rmarkdown::render("ISG_report_template.Rmd",
    params = list(
        ISG = ISG,
        NFkb = NFkb,
        IFNg = IFNg,
        panel_info = panel_info,
        cur_sample = cur_sample,
        sample_info = sample_info_all %>% filter(`Patient ID` %in% c(cur_sample, paste0(cur_sample, "M"), paste0(cur_sample, "F")))
    ),
    output_file = paste0("/Users/tan/Library/CloudStorage/OneDrive-KI.SE/ISG nanostring/ISG reports test/test3.pdf")
)
