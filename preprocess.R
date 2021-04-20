library(dplyr)
library(readr)
library(tibble)

dataPath <- '/Users/tan/nanostring/data/2021-03-30-BRIGGS_1_C7448.csv'
file <- read_csv(dataPath) %>% 
  rename(file_name = `File Name`, gene_name = X2, tag_code = X3)

pos <- file %>% filter(file_name == 'Positive') %>% 
  mutate(across(contains('Sample'), as.numeric))
raw <- file %>% filter(file_name == 'Endogenous') %>% 
  mutate(across(contains('Sample'), as.numeric))
housekeeping <- file %>% filter(file_name == 'Housekeeping') %>%
  mutate(across(contains('Sample'), as.numeric))

norm_factor <- function(dat){
  # get the normalization factor
  # each column is a sample
  # mean of geometric mean / geometric mean
  geomean_pos <- sapply(dat, function(x){exp(mean(log(x)))})
  return(mean(geomean_pos)/geomean_pos)
}
elementwise_multiply_by_row <- function(dat, Fac){
  # dat is a tibble, Fac is a vector of multiplier, with the same number of columns as dat 
  return(t(t(dat) * Fac) %>%
           as_tibble())
}
normFac_pos <- norm_factor(dat = pos %>% select(contains('Sample')))
normFac_house <- norm_factor(dat = elementwise_multiply_by_row(dat = housekeeping %>% 
                                                                 select(contains('Sample')),
                                                               Fac = normFac_pos))

tmp <- elementwise_multiply_by_row(dat = raw %>% select(contains('Sample')),
                                   Fac = normFac_house * normFac_pos) %>% t()
colnames(tmp) <- raw$gene_name
dat <- tmp %>% as_tibble() %>% add_column(Sample = rownames(tmp))

# geomean ISG score
columnSet <- colnames(dat %>% select(-Sample))
ISG <- dat %>% 
  rowwise() %>%
  mutate(geomean = exp(mean(log(c_across(columnSet))))) %>%
  select(-columnSet)
