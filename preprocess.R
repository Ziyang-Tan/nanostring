library(dplyr)
library(readr)
library(tibble)

norm_factor <- function(dat) {
  # get the normalization factor
  # each column is a sample
  # mean of geometric mean / geometric mean
  geomean_pos <- sapply(dat, function(x) {
    exp(mean(log(x)))
  })
  return(mean(geomean_pos) / geomean_pos)
}
elementwise_multiply_by_row <- function(dat, Fac) {
  # dat is a tibble, Fac is a vector of multiplier, with the same number of columns as dat
  return(t(t(dat) * Fac) %>%
    as_tibble())
}

normalization <- function(file, do_normalization = TRUE) {
  pos <- file %>%
    filter(file_name == "Positive") %>%
    mutate(across(contains("Sample"), as.numeric))
  raw <- file %>%
    filter(file_name == "Endogenous") %>%
    mutate(across(contains("Sample"), as.numeric))
  housekeeping <- file %>%
    filter(file_name == "Housekeeping") %>%
    mutate(across(contains("Sample"), as.numeric))
  if (do_normalization) {
    normFac_pos <- norm_factor(dat = pos %>% select(contains("Sample")))
    normFac_house <- norm_factor(dat = elementwise_multiply_by_row(
      dat = housekeeping %>%
        select(contains("Sample")),
      Fac = normFac_pos
    ))
  } else {
    normFac_pos <- 1
    normFac_house <- 1
  }

  tmp <- elementwise_multiply_by_row(
    dat = raw %>% select(contains("Sample")),
    Fac = normFac_house * normFac_pos
  ) %>% t()
  colnames(tmp) <- raw$gene_name
  dat <- tmp %>%
    as_tibble() %>%
    add_column(Sample = rownames(tmp))
  return(dat)
}

# geomean ISG score
geomean_score <- function(dat, columnSet) {
  # columnSet <- raw$gene_name
  ISG <- dat %>%
    rowwise() %>%
    mutate(geomean = exp(mean(log(c_across(all_of(columnSet)) + 1))))
  return(ISG)
}


# zscore ISG score
zscore_score <- function(dat, columnSet) {
  hc <- dat %>% filter(control)
  hcmean <- hc %>%
    summarise(across(all_of(columnSet), mean))
  hcstd <- hc %>%
    summarise(across(all_of(columnSet), sd))

  zscore <- dat %>%
    select(all_of(columnSet)) %>%
    apply(., 1, function(x) {
      (x - hcmean) / hcstd
    }) %>%
    bind_rows() %>%
    rowSums()
  ISG <- tibble(Sample = dat$Sample, zscore = zscore)
  return(ISG)
}
