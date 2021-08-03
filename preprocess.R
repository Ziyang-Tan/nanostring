library(dplyr)
library(readr)
library(tibble)

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

normalization <- function(file){
  pos <- file %>% filter(file_name == 'Positive') %>% 
    mutate(across(contains('Sample'), as.numeric))
  raw <- file %>% filter(file_name == 'Endogenous') %>% 
    mutate(across(contains('Sample'), as.numeric))
  housekeeping <- file %>% filter(file_name == 'Housekeeping') %>%
    mutate(across(contains('Sample'), as.numeric))
  normFac_pos <- norm_factor(dat = pos %>% select(contains('Sample')))
  normFac_house <- norm_factor(dat = elementwise_multiply_by_row(dat = housekeeping %>% 
                                                                   select(contains('Sample')),
                                                                 Fac = normFac_pos))
  tmp <- elementwise_multiply_by_row(dat = raw %>% select(contains('Sample')),
                                     Fac = normFac_house * normFac_pos) %>% t()
  colnames(tmp) <- raw$gene_name
  dat <- tmp %>% 
    as_tibble() %>% 
    add_column(Sample = rownames(tmp))
  return(dat)
}

# geomean ISG score
geomean_score <- function(dat, columnSet){
  #columnSet <- raw$gene_name
  ISG <- dat %>%
    rowwise() %>%
    mutate(geomean = exp(mean(log(c_across(all_of(columnSet)))))) %>%
    select(-columnSet)
  return(ISG)
}


# zscore ISG score
zscore_score <- function(dat, columnSet){
  hc <- dat %>% filter(grepl('HC', `Sample info`, ignore.case = T) | 
                         grepl('healthy control', `Sample info`, ignore.case = T)) 
  hcmean <- hc %>%
    summarise(across(columnSet, mean))
  hcstd <- hc %>%
    summarise(across(columnSet, sd))
  
  zscore <-  dat %>% 
    select(columnSet) %>% 
    apply(.,1, function(x){(x-hcmean)/hcstd}) %>% 
    bind_rows() %>% 
    rowSums()
  ISG <- tibble(Sample = dat$Sample, zscore = zscore)
  return(ISG)
}

# export ISG report
#dir.create('ISG_report', showWarnings = F)
#for (i in 1:dim(dat)[1]){
#  rmarkdown::render('ISG_report_template.Rmd', 
#                    params = list(dat = dat[i,], ISG_score = ISG[i,], info = info[i,]), 
#                    output_file = paste0('ISG_report/ISG_',i,'.pdf'))
#}

