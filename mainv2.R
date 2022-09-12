source('preprocess.R')
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)

dataDir <- '/Users/tan/OneDrive - KI.SE/ISG nanostring/data'
infoDir <- '/Users/tan/OneDrive - KI.SE/ISG nanostring/sample info'
#exportDir <- '/Users/tan/OneDrive - KI.SE/ISG nanostring/ISG reports'
exportDir <- '/Users/tan/nanostring/nanostring/ISG_report'
panelInfoPath <- '/Users/tan/nanostring/doc/New ISG extended panel - for Nanostring.xlsx'
sampleInfoPath <- '/Users/tan/Library/CloudStorage/OneDrive-KI.SE/ISG nanostring/sample info/sample and physician info.xlsx'

# new panel
batches <- c('EXP-22-DN5215', 'EXP-22-DN5218', 'EXP-22-DN5219', 'EXP-22-DN5221', 'EXP-22-DN5227')

#report_samples <- list(list('cur_sample'= '20220901_30102780380723-01_Sample09_09.RCC', 'extra_label'= '20220901_30102780380723-01_Sample10_10.RCC')
#                       #list('cur_sample'= '20220901_30102780380723-01_Sample12_12.RCC', 'extra_label'= '') #ISG-22
#)

# old panel
# batches <- c('EXP-21-DN5206', 'EXP-21-DN5207',
#              'EXP-21-DN5208', 'EXP-21-DN5209',
#              'EXP-21-DN5210', 'EXP-21-DN5211',
#              'EXP-21-DN5212', 'EXP-21-DN5214',
#              'EXP-22-DN5219')
# report_samples <- list(list('cur_sample'= '20210330_30102498310821-01_Sample01_01.RCC', 'extra_label'= ''), #ISG-1
#                        list('cur_sample'= '20210617_30102498320821-01_Sample10_10.RCC', 'extra_label'= ''), #ISG-2
#                        list('cur_sample'= '20210617_30102498320821-01_Sample11_11.RCC', 'extra_label'= '20211029_30102623920622-01_Sample12_12.RCC'), #ISG-3
#                        list('cur_sample'= '20210617_30102498320821-01_Sample12_12.RCC', 'extra_label'= ''), #ISG-4
#                        list('cur_sample'= '20210721_30102467470821-01_Sample09_09.RCC', 'extra_label'= '20210721_30102467470821-01_Sample12_12.RCC'), #ISG-5
#                        list('cur_sample'= '20210723_30102467480821-01_Sample08_08.RCC', 'extra_label'= ''), #ISG-6
#                        list('cur_sample'= '20210723_30102467480821-01_Sample09_09.RCC', 'extra_label'= '20210723_30102467480821-01_Sample10_10.RCC'), #ISG-7
#                        list('cur_sample'= '20210723_30102467480821-01_Sample11_11.RCC', 'extra_label'= '20210723_30102467480821-01_Sample12_12.RCC'), #ISG-8
#                        list('cur_sample'= '20210804_30102498290821-01_Sample09_09.RCC', 'extra_label'= ''), #ISG-9
#                        list('cur_sample'= '20210804_30102498290821-01_Sample10_10.RCC', 'extra_label'= '20210804_30102498290821-01_Sample11_11.RCC'), #ISG-10
#                        list('cur_sample'= '20210804_30102498290821-01_Sample12_12.RCC', 'extra_label'= ''), #ISG-11
#                        list('cur_sample'= '20210831_30102623790622-01_Sample10_10.RCC', 'extra_label'= ''), #ISG-12
#                        list('cur_sample'= '20211029_30102623920622-01_Sample10_10.RCC', 'extra_label'= ''), #ISG-13
#                        list('cur_sample'= '20211029_30102623920622-01_Sample11_11.RCC', 'extra_label'= '20211126_30102623830622-01_Sample11_11.RCC'), #ISG-14
#                        list('cur_sample'= '20211126_30102623830622-01_Sample12_12.RCC', 'extra_label'= '') #ISG-15
#                        )

#infoPath <- '/Users/tan/nanostring/sampleInfo/2021-03-30-BRIGGS_1_C7448 info.csv'
#batchName <- '2021-03-30-BRIGGS_1_C7448'

sample_info <- left_join(
  readxl::read_excel(sampleInfoPath, sheet = 1),
  readxl::read_excel(sampleInfoPath, sheet = 2),
  by = 'Referring physician'
)


dataList <- lapply(batches, function(x){
  #x = batches[1]
  dataPath <- paste0(dataDir, '/', x, '_RawData.csv')
  infoPath <- paste0(infoDir, '/', x, '_RNA info.xlsx')
  file <- readr::read_csv(dataPath) %>% 
    rename(file_name = `File Name`, gene_name = `...2`, tag_code = `...3`)
  info <- readxl::read_excel(infoPath, sheet = 2, skip = 1) %>%
    mutate(`Sample info` = gsub('\\(.*$', '', `Sample info`)) %>%
    mutate(`Sample ID` = as.character(`Sample ID`)) %>%
    mutate(`RNA amount (ng)` = as.character(`RNA amount (ng)`)) %>%
    mutate(`RNA vol (µl)` = as.character(`RNA vol (µl)`))
  # normalization
  dat <- normalization(file) %>%
    bind_cols(info) %>%
    add_column(batch = rep(x, dim(info)[1]))
})
dat <- bind_rows(dataList) %>%
  filter(Group != 'Background') %>%
  filter(QC == 'Pass') %>%
  mutate(batch = as.factor(batch)) %>%
  select(-c('date of sampling', 'Date of analysis', 'referring physician'))
dat <- dat %>% mutate(control = grepl('Healthy relative', Group, ignore.case = T) | 
                        grepl('Healthy control', Group, ignore.case = T)) %>%
  mutate(Group = case_when(
    Group == 'Patient' ~ 'Non-interferonopathy patient',
    Group == 'Interferonopathy patient' ~ 'Interferonopathy',
    TRUE ~ Group
  )) %>%
  mutate(
    # label = case_when(
    #   !is.na(Visit) ~ paste0(`Personal info`, '_', Visit),
    #   Group == 'Interferonopathy' ~ `Personal info`,
    #   (Group == 'Non-interferonopathy patient') & 
    #     (`Personal info` != 'Non-interferonopathy patient') ~ `Personal info`
    # )
    label = case_when(
      Group == 'Interferonopathy' ~ `Personal info`
    )
  )

columnSet <- colnames(dat)[1:55]

# batch correction
#corrected <- limma::removeBatchEffect(dat %>% select(columnSet) %>% t(), 
#                                      batch = dat$batch,
#                                      design = model.matrix(~dat$control)) %>% t() %>% as_tibble()
if (length(levels(dat$batch)) > 1){
  corrected <- sva::ComBat_seq(dat %>% select(columnSet) %>% t(), 
                               batch = dat$batch,
                               covar_mod = model.matrix(~dat$control)) %>% 
    t() %>% 
    as_tibble()
  
  datCorrected <- dat %>% select(-columnSet) %>% bind_cols(corrected)
} else {
  datCorrected <- dat
}

# load panel info
panel_info <- readxl::read_excel(path = panelInfoPath, sheet = 2) %>% filter(type=='Endogenous')
ISG_panel <- na.omit(panel_info$`ISG score`)
NFkb_panel <- na.omit(panel_info$`NF-kB score`)
IFNg_panel <- na.omit(panel_info$`IFN-g Score`)

# calculate scores
ISG <- geomean_score(datCorrected, ISG_panel) %>%
  add_column(zscore = (zscore_score(datCorrected, ISG_panel))$zscore)
NFkb <- geomean_score(datCorrected, NFkb_panel) %>%
  add_column(zscore = (zscore_score(datCorrected, NFkb_panel))$zscore)
IFNg <- geomean_score(datCorrected, IFNg_panel) %>%
  add_column(zscore = (zscore_score(datCorrected, IFNg_panel))$zscore)

#ISG2 <- geomean_score(dat, columnSet) %>%
#  add_column(zscore = (zscore_score(dat, columnSet))$zscore)

#write_csv(ISG1, '/Users/tan/nanostring/nanostring/figures/correct.csv')
#write_csv(ISG2, '/Users/tan/nanostring/nanostring/figures/no correct.csv')

# dat <- dat %>% add_column(geomean_corrected = corrected[,'geomean'],
#                          zscore_corrected = corrected[,'zscore'])
# ggplot(IFNg, aes(x = geomean, y = zscore, color=Group)) +
#   geom_point() +
#   geom_text_repel(aes(label = label),
#                   force = 2,
#                   max.overlaps = 20) +
#   labs(title = 'batch correction')

#export ISG report
report_samples <- unique(sample_info$`Patient ID`)
dir.create(exportDir, showWarnings = F)
lapply(report_samples, function(cur_sample){
  # cur_sample = 'ISG-18'
  cur_data <- dat %>% filter(grepl(cur_sample, `Sample info`))
  rmarkdown::render('ISG_report_template.Rmd',
                    params = list(ISG = ISG,
                                  NFkb = NFkb,
                                  IFNg = IFNg,
                                  panel_info = panel_info, 
                                  cur_sample = cur_sample,
                                  sample_info = sample_info %>% filter(`Patient ID` == cur_sample)),
                    output_file = paste0(exportDir, '/',cur_sample,'.pdf'))
})


