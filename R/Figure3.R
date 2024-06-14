########################## Figure 3 Codes ######################################
############################ Data Import #######################################
# Load Libraries
required_packages <- c("Seurat", "ggplot2", "dplyr", "SCP", "SpaCET", "ggthemes", "RColorBrewer")
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









#################################################################################
############################# Figure 3G #########################################
#################################################################################