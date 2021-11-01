source('preprocess.R')
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(ggrepel)

dataDir <- '/Users/tan/nanostring/data'
infoDir <- '/Users/tan/nanostring/sampleInfo'
batches <- c('EXP-21-DN5206', 'EXP-21-DN5207',
             'EXP-21-DN5208', 'EXP-21-DN5209',
             'EXP-21-DN5210', 'EXP-21-DN5211',
             'EXP-21-DN5212')
#infoPath <- '/Users/tan/nanostring/sampleInfo/2021-03-30-BRIGGS_1_C7448 info.csv'
#batchName <- '2021-03-30-BRIGGS_1_C7448'

dataList <- lapply(batches, function(x){
  dataPath <- paste0(dataDir, '/', x, '_RawData.csv')
  infoPath <- paste0(infoDir, '/', x, '_RNA info.csv')
  file <- read_csv(dataPath) %>% 
    rename(file_name = `File Name`, gene_name = `...2`, tag_code = `...3`)
  info <- read_csv(infoPath) %>%
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
  mutate(batch = as.factor(batch))
dat <- dat %>% mutate(control = grepl('Healthy relative', Group, ignore.case = T) | 
                        grepl('Healthy control', Group, ignore.case = T)) %>%
  mutate(Group = case_when(
    Group == 'Patient' ~ 'Non-interferonopathy patient',
    Group == 'Interferonopathy patient' ~ 'Interferonopathy',
    TRUE ~ Group
  )) %>%
  mutate(label = case_when(
    Group == 'Interferonopathy'  ~ sub('\\(.*$', '', `Sample info`)
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
dir.create('ISG_report', showWarnings = F)
#for (i in 1:dim(dat)[1]){
cur_sample = '20211029_30102623920622-01_Sample12_12.RCC'
extra_label = ''
#cur_sample = ISG[i,]$Sample
rmarkdown::render('ISG_report_template.Rmd',
                  params = list(ISG = ISG, 
                                extra_label = extra_label,
                                dat = (dat %>% filter(Sample == cur_sample))[,1:30], 
                                cur_sample = cur_sample),
                  output_file = paste0('ISG_report/ISG_',cur_sample,'.pdf'))


#ggsave('/Users/tan/nanostring/nanostring/figures/correct.pdf')

# ggplot(ISG2, aes(x = geomean, y = zscore, color=Group)) +
#   geom_point() +
#   geom_text_repel(aes(label = label),
#                   force = 2,
#                   max.overlaps = 20) +
#   labs(title = 'no batch correction')
# ggsave('/Users/tan/nanostring/nanostring/figures/no correct.pdf')
#  scale_y_continuous(trans = ggallin::pseudolog10_trans)
#ggsave('/Users/tan/nanostring/nanostring/figures/no correct.pdf')

#ggplot(dat, aes(x = geomean_corrected, y = zscore_corrected, color=batch)) +
#  geom_point() +
#  scale_y_continuous(trans = ggallin::pseudolog10_trans)
#ggsave('/Users/tan/nanostring/nanostring/figures/corrected.pdf')

#write_csv(dat, '/Users/tan/nanostring/nanostring/figures/corrected.csv')
