source('preprocess.R')
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)

dataDir <- '/Users/tan/OneDrive - KI.SE/ISG nanostring/data'
infoDir <- '/Users/tan/OneDrive - KI.SE/ISG nanostring/sample info'
exportDir <- '/Users/tan/OneDrive - KI.SE/ISG nanostring/ISG reports'
batches <- c('EXP-21-DN5206', 'EXP-21-DN5207',
             'EXP-21-DN5208', 'EXP-21-DN5209',
             'EXP-21-DN5210', 'EXP-21-DN5211',
             'EXP-21-DN5212', 'EXP-21-DN5214')
report_samples <- list(list('cur_sample'= '20210330_30102498310821-01_Sample01_01.RCC', 'extra_label'= ''), #ISG-1
                       list('cur_sample'= '20210617_30102498320821-01_Sample10_10.RCC', 'extra_label'= ''), #ISG-2
                       list('cur_sample'= '20210617_30102498320821-01_Sample11_11.RCC', 'extra_label'= '20211029_30102623920622-01_Sample12_12.RCC'), #ISG-3
                       list('cur_sample'= '20210617_30102498320821-01_Sample12_12.RCC', 'extra_label'= ''), #ISG-4
                       list('cur_sample'= '20210721_30102467470821-01_Sample09_09.RCC', 'extra_label'= '20210721_30102467470821-01_Sample12_12.RCC'), #ISG-5
                       list('cur_sample'= '20210723_30102467480821-01_Sample08_08.RCC', 'extra_label'= ''), #ISG-6
                       list('cur_sample'= '20210723_30102467480821-01_Sample09_09.RCC', 'extra_label'= '20210723_30102467480821-01_Sample10_10.RCC'), #ISG-7
                       list('cur_sample'= '20210723_30102467480821-01_Sample11_11.RCC', 'extra_label'= '20210723_30102467480821-01_Sample12_12.RCC'), #ISG-8
                       list('cur_sample'= '20210804_30102498290821-01_Sample09_09.RCC', 'extra_label'= ''), #ISG-9
                       list('cur_sample'= '20210804_30102498290821-01_Sample10_10.RCC', 'extra_label'= '20210804_30102498290821-01_Sample11_11.RCC'), #ISG-10
                       list('cur_sample'= '20210804_30102498290821-01_Sample12_12.RCC', 'extra_label'= ''), #ISG-11
                       list('cur_sample'= '20210831_30102623790622-01_Sample10_10.RCC', 'extra_label'= ''), #ISG-12
                       list('cur_sample'= '20211029_30102623920622-01_Sample10_10.RCC', 'extra_label'= ''), #ISG-13
                       list('cur_sample'= '20211029_30102623920622-01_Sample11_11.RCC', 'extra_label'= '20211126_30102623830622-01_Sample11_11.RCC'), #ISG-14
                       list('cur_sample'= '20211126_30102623830622-01_Sample12_12.RCC', 'extra_label'= '') #ISG-15
                       )
#infoPath <- '/Users/tan/nanostring/sampleInfo/2021-03-30-BRIGGS_1_C7448 info.csv'
#batchName <- '2021-03-30-BRIGGS_1_C7448'

dataList <- lapply(batches, function(x){
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
  mutate(batch = as.factor(batch))
dat <- dat %>% mutate(control = grepl('Healthy relative', Group, ignore.case = T) | 
                        grepl('Healthy control', Group, ignore.case = T)) %>%
  mutate(Group = case_when(
    Group == 'Patient' ~ 'Non-interferonopathy patient',
    Group == 'Interferonopathy patient' ~ 'Interferonopathy',
    TRUE ~ Group
  )) %>%
  mutate(label = case_when(
    !is.na(Visit) ~ paste0(`Personal info`, '_', Visit),
    Group == 'Interferonopathy' ~ `Personal info`,
    (Group == 'Non-interferonopathy patient') & 
      (`Personal info` != 'Non-interferonopathy patient') ~ `Personal info`
  ))

columnSet <- colnames(dat)[1:30]

# batch correction
#corrected <- limma::removeBatchEffect(dat %>% select(columnSet) %>% t(), 
#                                      batch = dat$batch,
#                                      design = model.matrix(~dat$control)) %>% t() %>% as_tibble()

corrected <- sva::ComBat_seq(dat %>% select(columnSet) %>% t(), 
                             batch = dat$batch,
                             covar_mod = model.matrix(~dat$control)) %>% 
  t() %>% 
  as_tibble()

datCorrected <- dat %>% select(-columnSet) %>% bind_cols(corrected)

# calculate ISGs

ISG <- geomean_score(datCorrected, columnSet) %>%
  add_column(zscore = (zscore_score(datCorrected, columnSet))$zscore)

#ISG2 <- geomean_score(dat, columnSet) %>%
#  add_column(zscore = (zscore_score(dat, columnSet))$zscore)

#write_csv(ISG1, '/Users/tan/nanostring/nanostring/figures/correct.csv')
#write_csv(ISG2, '/Users/tan/nanostring/nanostring/figures/no correct.csv')

#dat <- dat %>% add_column(geomean_corrected = corrected[,'geomean'],
#                          zscore_corrected = corrected[,'zscore'])
ggplot(ISG, aes(x = geomean, y = zscore, color=Group)) +
  geom_point() +
  geom_text_repel(aes(label = label),
                  force = 2,
                  max.overlaps = 20) +
  labs(title = 'batch correction')

#export ISG report
dir.create(exportDir, showWarnings = F)
lapply(report_samples, function(x){
  cur_sample = x$cur_sample
  extra_label = x$extra_label
  rmarkdown::render('ISG_report_template.Rmd',
                    params = list(ISG = ISG, 
                                  extra_label = extra_label,
                                  dat = (dat %>% filter(Sample == cur_sample))[,1:30], 
                                  cur_sample = cur_sample),
                    output_file = paste0(exportDir, '/',ISG[ISG$Sample == cur_sample,]$`Subject ID`,'.pdf'))
})


