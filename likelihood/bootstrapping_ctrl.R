#!/usr/bin/env Rscript

library(VGAM)
library(boot)

# Read in args
arg <- as.character(commandArgs(trailingOnly = TRUE))
tf <- arg[1]
output <- arg[2]


# Read in SNPs data
snps_df <- read.csv(
  paste0("data/accB/acc_undup_", tf, ".csv"),
  stringsAsFactors = FALSE
)


# Bootstrap (control model)
fit_ll_ctrl <- function(data_df, i) {
  sub_df <- data_df[i, ]

  loss <- function(rho) {
    return(
      -mean(
        dbetabinom(
          x = sub_df$ref_c,
          size = (sub_df$ref_c + sub_df$alt_c),
          prob = 0.5, rho = rho, log = TRUE
        )
      )
    )
  }

  return(-stats4::mle(loss, start = 0.1)@min)
}

boot_ll <- boot(
  snps_df, fit_ll_ctrl,
  parallel = "multicore", ncpus = 10, R = 1000
)
saveRDS(boot_ll, output)
