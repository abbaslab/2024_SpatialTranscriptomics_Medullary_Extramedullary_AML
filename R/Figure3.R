########################## Figure 3 Codes ######################################
############################ Data Import #######################################
# Load Libraries
required_packages <- c("Seurat", "ggplot2", "dplyr", "SCP", "SpaCET", "ggthemes", "RColorBrewer", "ComplexHeatmap", "circlize", "grid")
lapply(required_packages, require, character.only = TRUE)

# Load RDS files
scRNA_Ref <- readRDS("/path/to/scRNA/Reference.RDS")
BM1_filtered <- readRDS("/path/to/data/BM1_filtered.RDS")
BM2_filtered <- readRDS("/path/to/data/BM2_filtered.RDS")
EM1_filtered <- readRDS("/path/to/data/EM1_filtered.RDS")
EM2_filtered <- readRDS("/path/to/data/EM2_filtered.RDS")

#################################################################################
############################# Figure 3A #########################################
#################################################################################

# Define spots that has higher score than median leukemic score
Pos_AML_Spots <- WhichCells(BM1_filtered, expression = AML > median(BM1_filtered$Seurat_AML)) # This Seurat_AML feature defined at Figure 2
# Classify spots for HLS and LLS 
BM1_filtered$AML_Spots <- ifelse(colnames(BM1_filtered) %in% Pos_AML_Spots, "HLS", "LLS")

# Spatial Map for HLS vs LLS spots
AML_Spots_col = c("HLS" ="#D8511D", "LLS" = "#212E52")
SpatialDimPlot(BM1_filtered, group.by = "AML_Spots", crop = F, pt.size.factor = 1.25, cols = AML_Spots_col)
ggsave("/path/to/Figures/HLS_LLS_Spots.pdf")

#################################################################################
############################# Figure 3B #########################################
#################################################################################

# Create the stacked bar plot
StatPlot(BM1_filtered@meta.data, stat.by = "seurat_clusters", group.by = "AML_Spots", plot_type = "trend", label = T, palcolor = seurat_cols) # Colors for unsupervised clusters are defined at Figure 2
ggsave("/path/to/Figures/HLS_LLS_vs_Clusters_StackedBarPlot.pdf")

#################################################################################
############################# Figure 3C #########################################
#################################################################################

# Define Assay as Deconvolution scores
DefaultAssay(BM1_filtered) <- "predictions"

# Run Differential Deconvolution (Gene) test with defined deconvolution assay. In this way we are getting scores from deconvolution assay for each cell type and running Wilcoxon rank sum test for those scores to get co-localization of cell types.
BM1_filtered <- RunDEtest(BM1_filtered, group_by = "AML_Spots", only.pos = F)

# Create Volcano Plot
p <- VolcanoPlot(BM1_filtered, group_by = "AML_Spots", nlabel = 5, pt.size = 6, stroke.highlight = 4, label.size = 4, sizes.highlight = 3)
p
ggsave("/path/to/Figures/BM1_CellCoLocalization.pdf", height = 8, width = 12)

#################################################################################
############################# Figure 3D #########################################
#################################################################################

# Spatial Map for Monocytes
SpatialFeaturePlot(BM1_filtered, features = "Seurat_Monocytes", crop = F, pt.size.factor = 1.25, stroke = 0)
ggsave("/path/to/Figures/BM1_Monocytes.pdf")

# Create Box Plot for cell type distribution based on clusters
FeatureStatPlot(BM1_filtered, group.by = "seurat_clusters", stat.by = "Seurat_Monocytes", plot_type = "box", palcolor = seurat_cols, y.max = 1, comparisons = list(c('2','1'), c('1', '0'), c('1', 0)))
ggsave("/path/to/Figures/Monocytes_Cluster.pdf")

# Spatial Map for GMP
SpatialFeaturePlot(BM1_filtered, features = "Seurat_GMP", crop = F, pt.size.factor = 1.25, stroke = 0)
ggsave("/path/to/Figures/BM1_SpatialMap_GMP.pdf")

# Create Box Plot for cell type distribution based on clusters
FeatureStatPlot(BM1_filtered, group.by = "seurat_clusters", stat.by = "Seurat_GMP", plot_type = "box", palcolor = seurat_cols, y.max = 1, comparisons = list(c('2','1'), c('1', '0'), c('1', 0)))
ggsave("/path/to/Figures/BM1_GMP_Clusters.pdf")

#################################################################################
############################# Figure 3E #########################################
#################################################################################

# Spatial Map for Late Erythroids
SpatialFeaturePlot(BM1_filtered, features = "Seurat_Late Erythroid", crop = F, pt.size.factor = 1.25, stroke = 0)
ggsave("/path/to/Figures/BM1_SpatialMap_LateErythroid.pdf")

# Create Box Plot for cell type distribution based on clusters
FeatureStatPlot(BM1_filtered, group.by = "seurat_clusters", stat.by = "Seurat_Late Erythroid", plot_type = "box", palcolor = seurat_cols, comparisons = list(c('2','1'), c('1', '0'), c('1', 0)))
ggsave("/path/to/Figures/BM1_LateEryth_Box.pdf")

# Spatial Map for CD8 Naive
SpatialFeaturePlot(BM1_filtered, features = "Seurat_CD8 Naive", crop = F, pt.size.factor = 1.25, stroke = 0)
ggsave("/path/to/Figures/BM1_SpatialMap_CD8_Naive.pdf")

# Create Box Plot for cell type distribution based on clusters
FeatureStatPlot(BM1_filtered, group.by = "seurat_clusters", stat.by = "Seurat_CD8 Naive", plot_type = "box", palcolor = seurat_cols, y.max = 1, comparisons = list(c('2','1'), c('1', '0'), c('1', 0)))
ggsave("/path/to/Figures/BM1_CD8_Naive_Box.pdf")

#################################################################################
############################# Figure 3F #########################################
#################################################################################

# Get the current features
current_features <- rownames(EM1_filtered[["propMatFromSpaCET"]])

# Identify the features to keep (excluding the ones want to remove)
features_to_keep <- setdiff(current_features, c("Malignant", "Malignant cell state A", "Malignant cell state B", "cDC1 CLEC9A", "Macrophage other", "Macrophage M2", "cDC3 LAMP3", "Macrophage M1", "B cell switched memory", "B cell exhausted", "B cell non-switched memory", "B cell naive", "Tfh", "Th17", "Th1", "Th2", "cDC2 CD1C", "Unidentifiable"))
features_to_keep

# Extract the data matrix and subset to include only the features want to keep
data_matrix <- GetAssayData(EM1_filtered[["propMatFromSpaCET"]], slot = "data")
data_matrix_subset <- data_matrix[features_to_keep, ]

# Update the assay with the subsetted data matrix
DefaultAssay(EM1_filtered) <- "SCT"
EM1_filtered[["SpaCET_Subset"]] <- CreateAssayObject(data = data_matrix_subset)
DefaultAssay(EM1_filtered) <- "SpaCET_Subset"

# Extract the deconvolution scores matrix
deconv_scores <- GetAssayData(EM1_filtered[["SpaCET_Subset"]], slot = "data")

# Calculate the correlation matrix
correlation_matrix <- cor(t(deconv_scores))

# Calculate the mean deconvolution score for each cell type
cell_type_abundance <- rowMeans(deconv_scores)

# Function to determine significance stars based on absolute correlation values
get_correlation_stars <- function(correlation_value) {
  if (abs(correlation_value) > 0.7) {
    return("***")
  } else if (abs(correlation_value) > 0.5) {
    return("**")
  } else {
    return("")
  }
}

# Create a matrix to hold the stars based on correlation values
correlation_stars_matrix <- matrix(nrow = nrow(correlation_matrix), ncol = ncol(correlation_matrix))
for (i in 1:nrow(correlation_matrix)) {
  for (j in 1:ncol(correlation_matrix)) {
    correlation_stars_matrix[i, j] <- get_correlation_stars(correlation_matrix[i, j])
  }
}

# Customize colors for the heatmap (Blue to White to Red)
custom_col_fun <- colorRamp2(c(-2, -1, 0, 1, 2), c("#053061", "#2166AC", "white", "#B2182B", "#67001F"))

# Create a heatmap annotation for cell type abundance
abundance_annotation <- rowAnnotation(
  Abundance = anno_barplot(
    cell_type_abundance, 
    gp = gpar(fill = "darkgray", col = "black"),
    border = TRUE, 
    axis_param = list(side = "top"),
    width = unit(2, "cm")
  )
)

# Create the heatmap with additional customizations
heatmap <- Heatmap(
  correlation_matrix,
  name = "Correlation",
  col = custom_col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Correlation",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(4, "cm")
  ),
  right_annotation = abundance_annotation,
  top_annotation = HeatmapAnnotation(
    text = anno_text(colnames(correlation_matrix), rot = 45, just = "right", gp = gpar(fontsize = 10))
  ),
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA, lwd = 0.2))
    if (correlation_stars_matrix[i, j] != "") {
      grid.text(correlation_stars_matrix[i, j], x = x, y = y, gp = gpar(fontsize = 15, col = "black"))
    }
  }
)

# Save the heatmap as a high-resolution PDF
pdf("/path/to/Figures/EM1_CoLocalization_Matrix.pdf", width = 14, height = 14)
draw(heatmap)
dev.off()

#################################################################################
############################# Figure 3G #########################################
#################################################################################

# Add deconvolution scores as features
EM1_filtered <- AddFeaturesFromAssay(EM1_filtered, assay = "propMatFromSpaCET", prefix = "Spacet_") # This function defined at Figure 2

# Spatial Map for macrophage on EM1
SpatialFeaturePlot(EM1_filtered, features = "Spacet_Macrophage", crop = F, pt.size.factor = 1.25, alpha = c(0,7), stroke = 0)
ggsave("/path/to/Figures/Macrophage_Alpha_Spatial.pdf")

# Violin Plot for Macrophages on EM1
p <- FeatureStatPlot(EM1_filtered, group.by = "seurat_clusters", stat.by = "Spacet_Macrophage", add_point = T, palcolor = seurat_cols, add_box = T, comparisons = list(c('2','1'), c('1', '0'), c('1', 0)))
panel_fix(p, height = 4, width = 3.5, save = "/path/to/Figures/EM1_Macrophages_Clusters_VlnPlot.pdf")

# Spatial Map for cDC on EM1
SpatialFeaturePlot(EM1_filtered, features = "Spacet_cDC", crop = F, pt.size.factor = 1.25, alpha = c(0,7), stroke = 0)
ggsave("/path/to/Figures/cDC_Alpha_Spatial.pdf")

# Violin Plot for cDCs on EM1
p <- FeatureStatPlot(EM1_filtered, group.by = "seurat_clusters", stat.by = "Spacet_cDC", add_point = T, palcolor = seurat_cols, add_box = T, comparisons = list(c('2','1'), c('1', '0'), c('1', 0)))
panel_fix(p, height = 4, width = 3.5, save = "/path/to/Figures/EM1_cDCs_Clusters_VlnPlot.pdf")






