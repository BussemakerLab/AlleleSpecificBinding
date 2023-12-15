library(ggplot2)
TFs_info <- read.table("data/bindingModels/motifcentral_models.txt")


##### Figure 2B & Figure S1A,B #####

plotASB <- function(snps_df, cut_seq, cut_labs) {
  asb_df <- subset(snps_df, type == "ASB")
  pref_df <- data.frame(
    preferred.aff = ifelse(asb_df$ref_c >= asb_df$alt_c, exp(asb_df$ref_log_aff), exp(asb_df$alt_log_aff)),
    unpreferred.aff = ifelse(asb_df$ref_c >= asb_df$alt_c, exp(asb_df$alt_log_aff), exp(asb_df$ref_log_aff))
  )

  ggplot(pref_df, aes(unpreferred.aff, preferred.aff)) +
    geom_point(size = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "grey") +
    scale_x_continuous(trans = "log10", breaks = cut_seq, labels = cut_labs) +
    scale_y_continuous(trans = "log10", breaks = cut_seq, labels = cut_labs) +
    labs(x = "Unpreferred allele", y = "Preferred allele", title = "Predicted relative affinity") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
}


pdf("figures/Fig2B.pdf", 4, 4)
TF <- "CTCF"
snps_df <- read.csv(paste0("data/accB/acc_undup_", TF, ".csv"), stringsAsFactors = F)
snps_df$ref_log_aff <- read.table(paste0("data/scores/", TF, "_fit", TFs_info$V3[TFs_info$V2 == TF], "_k30_ref.tsv"))$V2
snps_df$alt_log_aff <- read.table(paste0("data/scores/", TF, "_fit", TFs_info$V3[TFs_info$V2 == TF], "_k30_alt.tsv"))$V2

plotASB(snps_df,
  cut_seq = c(1e-6, 1e-4, 1e-2, 1),
  cut_labs = c(expression(10^{
    -6
  }), expression(10^{
    -4
  }), expression(10^{
    -2
  }), 1)
)
dev.off()


pdf("figures/FigS1A.pdf", 4, 4)
TF <- "SPI1" # PU1
snps_df <- read.csv(paste0("data/accB/acc_undup_", TF, ".csv"), stringsAsFactors = F)
snps_df$ref_log_aff <- read.table(paste0("data/scores/", TF, "_fit", TFs_info$V3[TFs_info$V2 == TF], "_k30_ref.tsv"))$V2
snps_df$alt_log_aff <- read.table(paste0("data/scores/", TF, "_fit", TFs_info$V3[TFs_info$V2 == TF], "_k30_alt.tsv"))$V2

plotASB(snps_df,
  cut_seq = c(1e-3, 1e-2, 1e-1, 1),
  cut_labs = c(expression(10^{
    -3
  }), expression(10^{
    -2
  }), expression(10^{
    -1
  }), 1)
) +
  labs(title = "PU.1", subtitle = "Predicted relative affinity") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
dev.off()


pdf("figures/FigS1B.pdf", 4, 4)
TF <- "EBF1"
snps_df <- read.csv(paste0("data/accB/acc_undup_", TF, ".csv"), stringsAsFactors = F)
snps_df$ref_log_aff <- read.table(paste0("data/scores/", TF, "_fit", TFs_info$V3[TFs_info$V2 == TF], "_k30_ref.tsv"))$V2
snps_df$alt_log_aff <- read.table(paste0("data/scores/", TF, "_fit", TFs_info$V3[TFs_info$V2 == TF], "_k30_alt.tsv"))$V2

plotASB(snps_df,
  cut_seq = c(1e-4, 1e-3, 1e-2, 1e-1, 1),
  cut_labs = c(expression(10^{
    -4
  }), expression(10^{
    -3
  }), expression(10^{
    -2
  }), expression(10^{
    -1
  }), 1)
) +
  labs(title = TF, subtitle = "Predicted relative affinity") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
dev.off()


##### Figure 2C #####

TF <- "CTCF"
snps_df <- read.csv(paste0("data/accB/acc_undup_", TF, ".csv"), stringsAsFactors = F)
snps_df$ref_log_aff <- read.table(paste0("data/scores/", TF, "_fit", TFs_info$V3[TFs_info$V2 == TF], "_k30_ref.tsv"))$V2
snps_df$alt_log_aff <- read.table(paste0("data/scores/", TF, "_fit", TFs_info$V3[TFs_info$V2 == TF], "_k30_alt.tsv"))$V2

asb_df <- subset(snps_df, type == "ASB")
pref_df <- data.frame(
  preferred.aff = ifelse(asb_df$ref_c >= asb_df$alt_c, exp(asb_df$ref_log_aff), exp(asb_df$alt_log_aff)),
  unpreferred.aff = ifelse(asb_df$ref_c >= asb_df$alt_c, exp(asb_df$alt_log_aff), exp(asb_df$ref_log_aff))
)
# table(pref_df$preferred.aff == pref_df$unpreferred.aff)
pref_df$max_aff <- ifelse(pref_df$preferred.aff > pref_df$unpreferred.aff, pref_df$preferred.aff, pref_df$unpreferred.aff)
pref_df$concordant <- ifelse(pref_df$preferred.aff > pref_df$unpreferred.aff, TRUE, FALSE)
pref_df <- pref_df[order(pref_df$max_aff, decreasing = T), ]


binsize <- 100
nbin <- nrow(pref_df) - binsize + 1
continuous_bin_df <- data.frame(
  max_aff = pref_df$max_aff[1:nbin],
  concord = sapply(1:nbin, function(i) sum(pref_df$concordant[i:(i + binsize - 1)]) / binsize)
)

q001 <- qbinom(p = 0.01, size = binsize, prob = 0.5, lower.tail = F)

cut_seq <- c(1e-6, 1e-4, 1e-2, 1)
cut_labs <- c(expression(10^{
  -6
}), expression(10^{
  -4
}), expression(10^{
  -2
}), 1)

pdf("figures/Fig2C.pdf", 4, 4)
ggplot(continuous_bin_df, aes(max_aff, concord)) +
  geom_point(size = 0.4) +
  geom_hline(yintercept = q001 / binsize, color = "red") +
  geom_hline(yintercept = 0.5, color = "grey") +
  annotate("text", x = 0.05, y = 0.58, label = "corcordance = 62% \n(binomial p-value = 0.01)", color = "red", size = 3.5) +
  annotate("text", x = 0.05, y = 0.47, label = "corcordance = 50%", color = "grey", size = 3.5) +
  scale_x_continuous(trans = "log10", breaks = cut_seq, labels = cut_labs) +
  labs(
    x = "Predicted affinity (maximun of two alleles)",
    y = "Concordance fraction (in bin of 100 SNVs)"
  ) +
  theme_classic()
dev.off()


##### Figure S1C #####

TF <- "SPI1"
snps_df <- read.csv(paste0("data/accB/acc_undup_", TF, ".csv"), stringsAsFactors = F)
snps_df$ref_log_aff <- read.table(paste0("data/scores/", TF, "_fit", TFs_info$V3[TFs_info$V2 == TF], "_k30_ref.tsv"))$V2
snps_df$alt_log_aff <- read.table(paste0("data/scores/", TF, "_fit", TFs_info$V3[TFs_info$V2 == TF], "_k30_alt.tsv"))$V2

asb_df <- subset(snps_df, type == "ASB")
pref_df <- data.frame(
  preferred.aff = ifelse(asb_df$ref_c >= asb_df$alt_c, exp(asb_df$ref_log_aff), exp(asb_df$alt_log_aff)),
  unpreferred.aff = ifelse(asb_df$ref_c >= asb_df$alt_c, exp(asb_df$alt_log_aff), exp(asb_df$ref_log_aff))
)
# table(pref_df$preferred.aff == pref_df$unpreferred.aff)
pref_df$max_aff <- ifelse(pref_df$preferred.aff > pref_df$unpreferred.aff, pref_df$preferred.aff, pref_df$unpreferred.aff)
pref_df$concordant <- ifelse(pref_df$preferred.aff > pref_df$unpreferred.aff, TRUE, FALSE)
pref_df <- pref_df[order(pref_df$max_aff, decreasing = T), ]

binsize <- 100
nbin <- nrow(pref_df) - binsize + 1
continuous_bin_df <- data.frame(
  max_aff = pref_df$max_aff[1:nbin],
  concord = sapply(1:nbin, function(i) sum(pref_df$concordant[i:(i + binsize - 1)]) / binsize)
)

q001 <- qbinom(p = 0.01, size = binsize, prob = 0.5, lower.tail = F)

cut_seq <- c(1e-3, 1e-2, 1e-1, 1)
cut_labs <- c(expression(10^{
  -3
}), expression(10^{
  -2
}), expression(10^{
  -1
}), 1)

pdf("figures/FigS1C.pdf", 4, 4)
ggplot(continuous_bin_df, aes(max_aff, concord)) +
  geom_point(size = 0.4) +
  geom_hline(yintercept = q001 / binsize, color = "red") +
  geom_hline(yintercept = 0.5, color = "grey") +
  annotate("text", x = 0.3, y = 0.58, label = "corcordance = 62% \n(binomial p-value = 0.01)", color = "red", size = 3.5) +
  # annotate("text", x = 0.5, y = 0.47, label = "corcordance = 50%", color = "grey", size = 3.5) +
  scale_x_continuous(trans = "log10", breaks = cut_seq, labels = cut_labs) +
  labs(
    x = "Predicted affinity (maximun of two alleles)",
    y = "Concordance fraction (in bin of 100 SNVs)"
  ) +
  ggtitle("PU.1") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


##### Figure S1D #####

TF <- "EBF1"
snps_df <- read.csv(paste0("data/accB/acc_undup_", TF, ".csv"), stringsAsFactors = F)
snps_df$ref_log_aff <- read.table(paste0("data/scores/", TF, "_fit", TFs_info$V3[TFs_info$V2 == TF], "_k30_ref.tsv"))$V2
snps_df$alt_log_aff <- read.table(paste0("data/scores/", TF, "_fit", TFs_info$V3[TFs_info$V2 == TF], "_k30_alt.tsv"))$V2

asb_df <- subset(snps_df, type == "ASB")
pref_df <- data.frame(
  preferred.aff = ifelse(asb_df$ref_c >= asb_df$alt_c, exp(asb_df$ref_log_aff), exp(asb_df$alt_log_aff)),
  unpreferred.aff = ifelse(asb_df$ref_c >= asb_df$alt_c, exp(asb_df$alt_log_aff), exp(asb_df$ref_log_aff))
)
# table(pref_df$preferred.aff == pref_df$unpreferred.aff)
pref_df$max_aff <- ifelse(pref_df$preferred.aff > pref_df$unpreferred.aff, pref_df$preferred.aff, pref_df$unpreferred.aff)
pref_df$concordant <- ifelse(pref_df$preferred.aff > pref_df$unpreferred.aff, TRUE, FALSE)
pref_df <- pref_df[order(pref_df$max_aff, decreasing = T), ]

binsize <- 100
nbin <- nrow(pref_df) - binsize + 1
continuous_bin_df <- data.frame(
  max_aff = pref_df$max_aff[1:nbin],
  concord = sapply(1:nbin, function(i) sum(pref_df$concordant[i:(i + binsize - 1)]) / binsize)
)

q001 <- qbinom(p = 0.01, size = binsize, prob = 0.5, lower.tail = F)

cut_seq <- c(1e-4, 1e-3, 1e-2, 1e-1, 1)
cut_labs <- c(expression(10^{
  -4
}), expression(10^{
  -3
}), expression(10^{
  -2
}), expression(10^{
  -1
}), 1)


pdf("figures/FigS1D.pdf", 4, 4)
ggplot(continuous_bin_df, aes(max_aff, concord)) +
  geom_point(size = 0.4) +
  geom_hline(yintercept = q001 / binsize, color = "red") +
  geom_hline(yintercept = 0.5, color = "grey") +
  annotate("text", x = 0.2, y = 0.58, label = "corcordance = 62% \n(binomial p-value = 0.01)", color = "red", size = 3.5) +
  # annotate("text", x = 0.5, y = 0.47, label = "corcordance = 50%", color = "grey", size = 3.5) +
  scale_x_continuous(trans = "log10", breaks = cut_seq, labels = cut_labs) +
  labs(
    x = "Predicted affinity (maximun of two alleles)",
    y = "Concordance fraction (in bin of 100 SNVs)"
  ) +
  ggtitle(TF) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
