source('preprocess.R')
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(ggrepel)

dataDir <- '/Users/tan/nanostring/data'
infoDir <- '/Users/tan/nanostring/sampleInfo'
batches <- c('2021-03-30-BRIGGS_1_C7448', 'EXP-21-DN5207','EXP-21-DN5208', 'EXP-21-DN5209')
#infoPath <- '/Users/tan/nanostring/sampleInfo/2021-03-30-BRIGGS_1_C7448 info.csv'
#batchName <- '2021-03-30-BRIGGS_1_C7448'

dataList <- lapply(batches, function(x){
  dataPath <- paste0(dataDir, '/', x, '_RawData.csv')
  infoPath <- paste0(infoDir, '/', x, '_RNA info.csv')
  file <- read_csv(dataPath) %>% 
    rename(file_name = `File Name`, gene_name = `...2`, tag_code = `...3`)
  info <- read_csv(infoPath) %>%
    mutate(`Sample ID` = as.character(`Sample ID`))
  # normalization
  dat <- normalization(file) %>%
    bind_cols(info) %>%
    add_column(batch = rep(x, dim(info)[1]))
})
dat <- bind_rows(dataList) %>%
  filter(`Sample info` != 'Water') %>%
  mutate(batch = as.factor(batch))
dat <- dat %>% mutate(control = grepl('HC', `Sample info`, ignore.case = T) | 
                        grepl('healthy control', `Sample info`, ignore.case = T)) %>%
  mutate(label = case_when(
    !control  ~ sub('\\(.*$', '', `Sample info`)
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

ISG1 <- geomean_score(datCorrected, columnSet) %>%
  add_column(zscore = (zscore_score(datCorrected, columnSet))$zscore)

ISG2 <- geomean_score(dat, columnSet) %>%
  add_column(zscore = (zscore_score(dat, columnSet))$zscore)

write_csv(ISG1, '/Users/tan/nanostring/nanostring/figures/correct.csv')
write_csv(ISG2, '/Users/tan/nanostring/nanostring/figures/no correct.csv')

#dat <- dat %>% add_column(geomean_corrected = corrected[,'geomean'],
#                          zscore_corrected = corrected[,'zscore'])
ggplot(ISG1, aes(x = geomean, y = zscore, color=batch)) +
  geom_point() +
  geom_text_repel(aes(label = label),
                  force = 2,
                  max.overlaps = 20) +
  labs(title = 'batch correction')
ggsave('/Users/tan/nanostring/nanostring/figures/correct.pdf')

ggplot(ISG2, aes(x = geomean, y = zscore, color=batch)) +
  geom_point() +
  geom_text_repel(aes(label = label),
                  force = 2,
                  max.overlaps = 20) +
  labs(title = 'no batch correction')
ggsave('/Users/tan/nanostring/nanostring/figures/no correct.pdf')
#  scale_y_continuous(trans = ggallin::pseudolog10_trans)
#ggsave('/Users/tan/nanostring/nanostring/figures/no correct.pdf')

#ggplot(dat, aes(x = geomean_corrected, y = zscore_corrected, color=batch)) +
#  geom_point() +
#  scale_y_continuous(trans = ggallin::pseudolog10_trans)
#ggsave('/Users/tan/nanostring/nanostring/figures/corrected.pdf')

#write_csv(dat, '/Users/tan/nanostring/nanostring/figures/corrected.csv')
