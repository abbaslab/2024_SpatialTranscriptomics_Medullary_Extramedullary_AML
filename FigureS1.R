########################## Figure S1 Codes ######################################
############################ Data Import ########################################

# These are the 10X Genomics SpaceRanger output directiory
BM1 <- Load10X_Spatial("/path/to/BM1/outs/")
BM2 <- Load10X_Spatial("/path/to/BM2/outs/")
EM1 <- Load10X_Spatial("/path/to/EM1/outs/")
EM2 <- Load10X_Spatial("/path/to/EM2/outs/")
EM1_v1 <- Load10X_Spatial("/path/to/EM1_v1/outs/")
EM2_v1 <- Load10X_Spatial("/path/to/EM2_v1/outs/")

#################################################################################
########################### Figure S1G,H,I ######################################
#################################################################################



















#################################################################################
########################### Figure S1J ##########################################
#################################################################################

# Load required libraries
library(ggplot2)
library(ggpubr)  
# Create a list of data frames for each sample
sample_data <- list(
  BM1 = data.frame(x_var = BM1@meta.data[,"nCount_Spatial"], y_var = BM1@meta.data[,"nFeature_Spatial"]),
  BM2 = data.frame(x_var = BM2@meta.data[,"nCount_Spatial"], y_var = BM2@meta.data[,"nFeature_Spatial"]),
  EM2 = data.frame(x_var = EM2@meta.data[,"nCount_Spatial"], y_var = EM2@meta.data[,"nFeature_Spatial"]),
  EM1 = data.frame(x_var = EM1@meta.data[,"nCount_Spatial"], y_var = EM1@meta.data[,"nFeature_Spatial"]),
  EM2_v1 = data.frame(x_var = EM2_v1@meta.data[,"nCount_Spatial"], y_var = EM2_v1@meta.data[,"nFeature_Spatial"]),
  EM1_v1 = data.frame(x_var = EM1_v1@meta.data[,"nCount_Spatial"], y_var = EM1_v1@meta.data[,"nFeature_Spatial"])
)
# Create a list to store representative points and colors
representative_points <- list()
colors <- rainbow(length(sample_data))

# Loop through sample data and calculate representative points and colors
for (i in 1:length(sample_data)) {
  df <- sample_data[[i]]
  representative_point <- data.frame(
    Median_x = median(df$x_var, na.rm = TRUE),
    Median_y = median(df$y_var, na.rm = TRUE),
    Sample = names(sample_data)[i]
  )
  representative_points[[i]] <- representative_point
}

# Combine representative points into a single data frame
combined_data <- do.call(rbind, representative_points)

# Create the correlation plot using ggplot
gg <- ggplot(combined_data, aes(x = Median_x, y = Median_y, color = Sample)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = FALSE, color = "darkblue", size = 1) +
  xlab("Median nCount_Spatial") +
  ylab("Median nFeature_Spatial") +
  ggtitle("Correlation between Median nCount_Spatial and Median nFeature_Spatial") +
  theme_minimal() +
  theme(text = element_text(size = 14))

# Add p-values to the plot
gg <- gg + stat_cor(aes(label = paste("p =", ..p..)), method = "pearson", size = 4, color = "black")

gg
# Save the figure as a high-resolution image (e.g., PDF)
ggsave("./figs/Playground/All_correlation_plot.pdf", gg, width = 8, height = 6, dpi = 700)

# Display the figure
print(gg)
#################################################################################
########################### Figure S1K ##########################################
#################################################################################
library(SCP)
# Read CSV file that contain spot barcode information assigned with LoupeBrowser
annotations <- read.csv("/path/to/EdgeEffect/EM1_Spots.csv", stringsAsFactors = FALSE)
# Assign to Seurat object
EM1[["EdgeEffect"]] <- annotations$EdgeEffect[match(colnames(EM1), annotations$Barcode)]
# Violin Plot
FeatureStatPlot(EM1, stat.by = "percent_mito", group.by = "EdgeEffect", add_box = T, y.max = 7, palcolor = list(c("#76B7B2","#F7B6D2")))
ggsave("/path/to/EdgeEffect/EM1_VlnPlot.pdf")
# Define colors
my_cols <- c("Edge" = "#F7B6D2", "Others" = "#76B7B2")
# Spatial Map
SpatialDimPlot(S18_full, group.by = "EdgeEffect", cols = my_cols)
ggsave("/path/to/EdgeEffect/EM1_SpatialDimPlot.pdf")
################################ BM1 ############################################
# Read in the csv BM1 spot barcodes
annotations <- read.csv("/path/to/EdgeEffect/BM1_Spots.csv", stringsAsFactors = FALSE)
# Assign barcode information to Seurat object
BM1[["EdgeEffect"]] <- annotations$EdgeEffect[match(colnames(BM1), annotations$Barcode)]
# Spatial Map
SpatialDimPlot(BM1, group.by = "EdgeEffect", cols = my_cols)
ggsave("./EdgeEffect/BM1_SpatialDimPlot.pdf")
# Violin Plot
FeatureStatPlot(BM18_full, stat.by = "percent_mito", group.by = "EdgeEffect", add_box = T, y.max = 7.5, palcolor = list(c("#76B7B2","#F7B6D2")))
ggsave("./EdgeEffect/BM1_VlnPlot.pdf")
#################################################################################
########################### Figure S1L ##########################################
#################################################################################
mito_percentages <- c(
  mean(EM2@meta.data$percent_mito),
  mean(EM1@meta.data$percent_mito),
  mean(BM2@meta.data$percent_mito),
  mean(BM1@meta.data$percent_mito)
)

dv200_values <- c(24, 53, 62, 39) # DV200 values for EM2, EM1, BM2, BM1

data_for_radar <- data.frame(
  DV200 = dv200_values,
  Mito_Percentage = mito_percentages
)

rownames(data_for_radar) <- c("EM2", "EM1", "BM2", "BM1")
data_for_radar_t <- t(data_for_radar)
colnames(data_for_radar_t) <- c("EM2", "EM1", "BM2", "BM1") # This should match row names before transposition
data_for_radar_t <- as.data.frame(data_for_radar_t)
gm1 <- rbind(
  max = rep(100, ncol(data_for_radar_t)), # Max values row, since your data is in percentages
  min = rep(0, ncol(data_for_radar_t)),   # Min values row
  data_for_radar_t                        # Actual data
)
# Plot the radar chart
radarchart(gm1, axistype = 1,
           pcol = c("#4198FF", "#D02D25"), plwd = 3.5, plty = 1,
           caxislabels = seq(0, 100, by = 20), # Axis labels from 0 to 100
           cglcol = "grey", cglty = 1, cglwd = 1.5, calcex = 1.8,
           axislabcol = "grey", vlcex = 2
)

# Add a legend at the bottom
legend("topright", legend = rownames(gm1)[-c(1,2)], horiz = TRUE,
       bty = "n", pch = 20, col = c("#4198FF", "#D02D25"),
       text.col = "black", cex = 1.5, pt.cex = 2
)
#################################################################################
########################### Figure S1M ##########################################
#################################################################################
SpatialFeaturePlot(BM2, features = "nCount_Spatial", crop = F, pt.size.factor = 1.25)
ggsave("/path/to/figures/FigS1M_nCount.pdf")
SpatialFeaturePlot(BM2, features = "nFeature_Spatial", crop = F, pt.size.factor = 1.25)
ggsave("/path/to/figures/FigS1M_nFeature.pdf")
