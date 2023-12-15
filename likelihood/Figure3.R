library(ggplot2)
library(reshape2)

##### Figure 3 & Figure S2 #####

plot_boot <- function(TF) {
  boot_pb <- readRDS(list.files("data/bootstrap/", pattern = paste0(TF, "_fit(?:.+)k30_boot"), full.names = T))
  boot_ctrl <- readRDS(paste0("data/bootstrap/", TF, "_control_boot.rds"))

  t_df <- data.frame(cbind(boot_ctrl$t, boot_pb$t))
  colnames(t_df) <- c("Equal-preference control", "MotifCentral affinity-based model")
  t_m <- melt(t_df, value.name = "bootstrap_ll", variable.name = "model")
  res_ids <- data.frame(t0 = c(boot_ctrl$t0, boot_pb$t0), model = c("control model", "affinity-based model"))

  ggplot(t_m, aes(bootstrap_ll, fill = model, color = model)) +
    geom_histogram(alpha = 0.5, position = "identity") +
    geom_vline(data = res_ids, aes(xintercept = t0, color = model)) +
    labs(x = "(bootstrap) log-likelihood", title = TF) +
    theme_bw()
}


pdf("figures/Fig3.pdf", 7, 4)
plot_boot(TF = "CTCF")
dev.off()

pdf("figures/FigS2A.pdf", 7, 4)
plot_boot(TF = "EBF1")
dev.off()

pdf("figures/FigS2B.pdf", 7, 4)
plot_boot(TF = "SPI1") + ggtitle("PU.1")
dev.off()
