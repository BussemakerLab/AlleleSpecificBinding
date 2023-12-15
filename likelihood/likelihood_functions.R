library(VGAM)

## CTCF

TF <- "CTCF"
snps_df <- read.csv(paste0("data/accB/acc_undup_", TF, ".csv"), stringsAsFactors = F)


## Control Model (Beta-binomial distribution with p=0.5)

loss_ctrl <- function(rho) {
  L <- dbetabinom(x = snps_df$ref_c, size = (snps_df$ref_c + snps_df$alt_c), prob = 0.5, rho = rho, log = T)
  return(-mean(L))
}

fit_ctrl <- stats4::mle(loss_ctrl, start = 0.1)
print(paste(TF, " (Control model)"))
print(paste0("Fitted parameter rho = ", fit_ctrl@coef))
print(paste0("Maximun log-ikelihood = ", -fit_ctrl@min))


## Affinity-based model

logaddexp <- function(a, b) {
  pmax(a, b) + log1p(exp(-abs(a - b)))
}


fitLL <- function(snps_df, ref_ppb, alt_ppb) {
  snps_df$log_ref_score <- read.table(ref_ppb)$V2
  snps_df$log_alt_score <- read.table(alt_ppb)$V2

  loss <- function(rho, log_const) {
    log_ref_score <- logaddexp(snps_df$log_ref_score, log_const)
    log_alt_score <- logaddexp(snps_df$log_alt_score, log_const)
    p <- exp(log_ref_score - logaddexp(log_ref_score, log_alt_score))
    L <- dbetabinom(x = snps_df$ref_c, size = (snps_df$ref_c + snps_df$alt_c), prob = p, rho = rho, log = T)
    return(-mean(L))
  }

  init_log_const <- mean(c(snps_df$log_ref_score, snps_df$log_alt_score))
  fit_res <- stats4::mle(loss, start = c(0.1, init_log_const))
  return(fit_res)
}


# Affinities predicted using motifCentral (SELEX derived) models
fit_selex <- fitLL(snps_df,
  ref_ppb = "data/scores/CTCF_fit12715_k30_ref.tsv",
  alt_ppb = "data/scores/CTCF_fit12715_k30_alt.tsv"
)

print(paste0("Fitted parameter rho = ", fit_selex@coef[1], " const = ", exp(fit_selex@coef[2])))
print(paste0("Maximun log-ikelihood = ", -fit_selex@min))


# Affinities predicted using ChIP-seq derived binding model
fit_chip <- fitLL(snps_df,
  ref_ppb = "data/scores/CTCF_ChIP_k30_ref.tsv",
  alt_ppb = "data/scores/CTCF_ChIP_k30_alt.tsv"
)

print(paste0("Fitted parameter rho = ", fit_chip@coef[1], " const = ", exp(fit_chip@coef[2])))
print(paste0("Maximun log-ikelihood = ", -fit_chip@min))
