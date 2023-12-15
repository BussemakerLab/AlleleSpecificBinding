## Read in ASB data

library(BSgenome.Hsapiens.UCSC.hg19)

# ASB data downloaded from AlleleDB
acc <- read.table(
  "data/accB.auto.v2.1.aug16.txt",
  header = TRUE, stringsAsFactors = FALSE
)
acc$TF <- sapply(strsplit(acc$TF_indiv_accB, "_"), "[[", 1)
acc$indiv <- sapply(strsplit(acc$TF_indiv_accB, "_"), "[[", 2)

asb <- read.table(
  "data/ASB.auto.v2.1.aug16.txt",
  header = TRUE, stringsAsFactors = FALSE
)
asb$TF <- sapply(strsplit(asb$TF_indiv_ASB, "_"), "[[", 1)
asb$indiv <- sapply(strsplit(asb$TF_indiv_ASB, "_"), "[[", 2)

acc$type <- ifelse(
  paste(acc$chr, acc$start, acc$TF, acc$indiv, sep = "_")
  %in% paste(asb$chr, asb$start, asb$TF, asb$indiv, sep = "_"), "ASB", "Control"
)

## Find counts of ref and alt alleles
# SNPs with multiple alt alleles: keep the one with max count
ind <- grep(";", acc$alt)
acc$alt1 <- acc$alt
acc$alt1[ind] <- sapply(ind, function(i) {
  altc <- sapply(
    strsplit(acc$alt[[i]], split = ";")[[1]],
    function(x) acc[i, grep(x, colnames(acc)[7:10], value = TRUE)]
  )
  names(altc)[which.max(altc)]
})

acc$ref_c <- sapply(
  seq_len(nrow(acc)),
  function(i) acc[i, grep(acc$ref[i], colnames(acc)[7:10], value = TRUE)]
)
acc$alt_c <- sapply(
  seq_len(nrow(acc)),
  function(i) acc[i, grep(acc$alt1[i], colnames(acc)[7:10], value = TRUE)]
)


## Remove duplicated SNPs for each TF (keep the most significant SNP)
get_acc_tf <- function(tf) {
  df <- acc[acc$TF == tf, ]
  snps <- paste(df$chr, df$start, sep = "_")

  df_uniq <- lapply(unique(snps), function(snp) {
    df_snp <- df[snps == snp, ]
    df_snp[which.min(df_snp$p.betabinomial), ]
  })
  df_uniq <- Reduce(rbind, df_uniq)
  return(df_uniq)
}

acc_ls <- lapply(unique(acc$TF), get_acc_tf)
names(acc_ls) <- unique(acc$TF)

saveRDS(acc_ls, "data/accB/acc_ls.rds")


# Replace TF names
tf_mapping <- list(c("CTCF", "CTCF"), c("EBF", "EBF1"), c("PU1", "SPI1"))
for (tf_map in tf_mapping) {
  write.csv(
    acc_ls[[tf_map[1]]],
    paste0("data/accB/acc_undup_", tf_map[2], ".csv"),
    row.names = FALSE, quote = FALSE
  )
}

# Get sequences surrounding SNPs (k=30)
get_snp_seq <- function(df, k) {
  seq <- getSeq(Hsapiens, df$chr, df$start - k + 2, df$end + k - 1)
  seq <- as.data.frame(seq)
  names(seq) <- "ref_seq"
  seq$alt_seq <- seq$ref_seq
  substr(seq$alt_seq, k, k) <- df$alt1
  return(seq)
}

for (tf_map in tf_mapping) {
  seq_df <- get_snp_seq(acc_ls[[tf_map[1]]], k = 30)
  write.table(
    seq_df$ref_seq, paste0("data/seq/", tf_map[2], "_k30_ref.txt"),
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  write.table(
    seq_df$alt_seq, paste0("data/seq/", tf_map[2], "_k30_alt.txt"),
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}
