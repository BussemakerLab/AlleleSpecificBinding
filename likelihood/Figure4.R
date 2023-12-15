library(ggplot2)
library(reshape2)
library(scales)

##### Figure 4 #####

TF <- "CTCF"
external_df <- read.table(paste0("data/bindingModels/", TF, "_external_models.txt"), stringsAsFactors = F, header = F)


boot_pb <- readRDS(list.files("data/bootstrap/", pattern = paste0(TF, "_fit(?:.+)k30_boot"), full.names = T))
boot_ctrl <- readRDS(paste0("data/bootstrap/", TF, "_control_boot.rds"))
t_df <- data.frame(cbind(boot_ctrl$t, boot_pb$t))
colnames(t_df) <- c("Control", "motifCentral")

chip_models <- c("CUTnTag_end0", "ChIP-exo5_end1", "ChIP-seq_end5")
boot_chip_ls <- lapply(chip_models, function(m) readRDS(paste0("data/bootstrap/", TF, "_", m, "_k30_boot.rds")))
t_chip_df <- sapply(boot_chip_ls, function(bt) bt$t)
colnames(t_chip_df) <- c("CUT&Tag", "ChIP-exo5", "ChIP-seq")

boot_ex_ls <- lapply(external_df$V1, function(id) readRDS(paste0("data/bootstrap/", TF, "_external", id, "_k30_boot.rds")))
t_ex_df <- sapply(boot_ex_ls, function(bt) bt$t)
colnames(t_ex_df) <- external_df$V2

t_df <- cbind(t_df, t_chip_df, t_ex_df)
t_m <- melt(t_df, value.name = "bootstrap_ll", variable.name = "model")


res_ids <- data.frame(
  t0 = c(boot_ctrl$t0, boot_pb$t0, sapply(boot_chip_ls, function(bt) bt$t0), sapply(boot_ex_ls, function(bt) bt$t0)),
  model = colnames(t_df)
)
res_ids <- res_ids[order(res_ids$t0, decreasing = T), ]

t_m$model <- factor(t_m$model, levels = res_ids$model)

model_labels <- c(
  "CUT&Tag (this study/PyProBound)",
  "ChIP-seq (this study/PyProBound)",
  "SELEX (MotifCentral/ProBound)",
  "ChIP-exo5 (this study/PyProBound)",
  "ChIP-seq (JASPAR)",
  "ChIP-seq (HOCOMOCO)",
  "Equal-preference control"
)

col5 <- hue_pal()(7)[-c(1, 5)]
model_colors <- c(col5[1:2], hue_pal()(2)[2], col5[3:5], hue_pal()(2)[1])


pdf("figures/Fig4A.pdf", 8.5, 5)
ggplot(t_m, aes(model, bootstrap_ll, fill = model, color = model)) +
  geom_boxplot() +
  labs(y = "(bootstrap) log-likelihood", x = NULL) +
  scale_colour_manual(values = model_colors, labels = model_labels) +
  scale_fill_manual(values = model_colors, labels = model_labels) +
  theme_bw()
dev.off()
