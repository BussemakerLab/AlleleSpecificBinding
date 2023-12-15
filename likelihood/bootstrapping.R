#!/usr/bin/env Rscript

library(VGAM)
library(boot)

# Read in args
arg <- as.character(commandArgs(trailingOnly = TRUE))
tf <- arg[1]
ref_ppb <- arg[2]
alt_ppb <- arg[3]
output <- arg[4]

# Read in SNPs data
snps_df <- read.csv(
  paste0("data/accB/acc_undup_", tf, ".csv"),
  stringsAsFactors = FALSE
)
snps_df$log_ref_score <- read.table(ref_ppb)$V2
snps_df$log_alt_score <- read.table(alt_ppb)$V2

# Bootstrap
logaddexp <- function(a, b) {
  pmax(a, b) + log1p(exp(-abs(a - b)))
}

fit_ll <- function(data_df, i) {
  sub_df <- data_df[i, ]

  loss <- function(rho, log_const) {
    log_ref_score <- logaddexp(sub_df$log_ref_score, log_const)
    log_alt_score <- logaddexp(sub_df$log_alt_score, log_const)
    return(
      -mean(
        dbetabinom(
          x = sub_df$ref_c,
          size = (sub_df$ref_c + sub_df$alt_c),
          prob = exp(log_ref_score - logaddexp(log_ref_score, log_alt_score)),
          rho = rho, log = TRUE
        )
      )
    )
  }

  init_log_const <- mean(c(snps_df$log_ref_score, snps_df$log_alt_score))
  return(-stats4::mle(loss, start = c(0.1, init_log_const))@min)
}


boot_ll <- boot(snps_df, fit_ll, parallel = "multicore", ncpus = 10, R = 1000)
saveRDS(boot_ll, output)
