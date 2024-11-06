library(readr)
library(ggplot2)
library(dplyr)
library(patchwork)

mode <- "branch"
snp_keep <- 1
matvec <- TRUE

file_suffix <- mode
if (matvec) file_suffix <- paste0(file_suffix, "_matvec")
if (snp_keep < 1) file_suffix <- paste0(mode, "_prop", as.integer(snp_keep * 100))

path_to_pca <- paste0("output/chr3_pca_", file_suffix, ".csv")
pca_data <- read_csv(path_to_pca)

plot_pca <- function(data, x_axis, y_axis, color_by) {
  plot <- ggplot(data, aes_string(x = x_axis, y = y_axis, color = color_by)) +
    geom_point(size=0.8) +
    theme_classic()

  return(plot)
}

pc1 <- plot_pca(data = arrange(pca_data, PC3), x_axis = "PC1",
                y_axis = "PC2", color_by = "proband_region")
pc2 <- plot_pca(data = arrange(pca_data, PC5), x_axis = "PC3",
                y_axis = "PC4", color_by = "proband_region")
pc3 <- plot_pca(data = arrange(pca_data, PC7), x_axis = "PC5",
                y_axis = "PC6", color_by = "proband_region")
pc4 <- plot_pca(data = arrange(pca_data, PC9), x_axis = "PC7",
                y_axis = "PC8", color_by = "proband_region")


pca_plots <- (pc1 + pc2) / (pc3 + pc4) +
  plot_annotation(title = paste0(file_suffix, " PCA")) +
  plot_layout(guides = 'collect')

path_to_plot <- paste0("plots/chr3_pca_", file_suffix, ".png")
ggsave(path_to_plot, pca_plots)