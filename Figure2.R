########################## Figure 2 Codes ######################################
############################ Data Import #######################################

# Load Libraries
required_packages <- c("Seurat", "ggplot2", "dplyr", "SCP", "SpaCET", "ggthemes", "RColorBrewer")
lapply(required_packages, require, character.only = TRUE)

scRNA_Ref <- readRDS("/path/to/scRNA/Reference.RDS")
BM1_filtered <- readRDS("/path/to/data/BM1_filtered.RDS")
BM2_filtered <- readRDS("/path/to/data/BM2_filtered.RDS")
EM1_filtered <- readRDS("/path/to/data/EM1_filtered.RDS")
EM2_filtered <- readRDS("/path/to/data/EM2_filtered.RDS")

#################################################################################
############################# Figure 2A #########################################
#################################################################################
# Function to assign colors for each cell type
ColAssign <- function(Var,palettes="Classic 20"){
  require(ggthemes);require(RColorBrewer)
  pal <- tableau_color_pal(palette = palettes,direction = 1,type="regular")
  if (length(Var) > 20) {
    palOut <- colorRampPalette(pal(20))(length(Var))
    names(palOut) <- Var
  } else if (length(Var) == 20) {
    palOut <- pal(20)
    names(palOut) <- Var
  } else if (length(Var) < 20) {
    palOut <- pal(20)
    palOut <- setdiff(palOut,c("#7F7F7F","#C7C7C7"))# remove grey colors
    #palOut <- sample(palOut)
    palOut <- c(palOut,c("#7F7F7F","#C7C7C7"))
    palOut <- palOut[1:length(Var)]
    names(palOut) <- Var
  }
  return(palOut)
}

# UMAP
DimPlot(scRNA_Ref,group.by="class2",cols=ColAssign(unique(scRNA_Ref$class2)),raster=F)+
  NoAxes()+ggtitle("")+
  theme(legend.title=element_text(size = 0),
        legend.text=element_text(size = 8),
        legend.key.height = unit(0.25,"cm"),
        legend.key.width = unit(0.25,"cm"))+NoAxes()
ggsave("/path/to/Figures/AllCells_DimPlot.pdf",height=12,width=18,unit="cm")
#################################################################################
############################# Figure 2B #########################################
#################################################################################
# Assign colors for spatial map
seurat_cols <- c(`0` = "#FFC312", `1` = "#C4E538", `2` = "#12CBC4")
BM1_filtered$seurat_clusters <- factor(BM1_filtered$seurat_clusters,levels=names(seurat_cols))
# Plot Spatial map for Unsupervised Clusters
SpatialDimPlot(BM1_filtered, group.by = "seurat_clusters", cols = seurat_cols, crop = F, pt.size.factor = 1.25, stroke = 0)
ggsave("/path/to/Figures/SpatialMap_UnsupervisedClusters.pdf")
# Import Pathology Annotation that defined with Loupe Browser
annotations <- read.csv("/path/to/Pathology/BM1_Pathology.csv", stringsAsFactors = FALSE)
# Add meta data of Seurat object
BM1_filtered[["Pathology"]] <- annotations$Regions[match(colnames(BM1_filtered), annotations$Barcode)]
# Define colors for Regions
region_cols <- c(`Region1` = "#ff7f0e", `Region2` = "#2ca02c", `Region3` = "#1f77b4")
# Region names changed during Adobe Illustrator editing for better representation
BM1_filtered$Regions <- factor(BM1_filtered$Regions,levels=names(region_cols))
# Spatial map for Pathology Annotation
SpatialDimPlot(BM1_filtered, group.by = "Regions", crop = F, pt.size.factor = 1.25,cols = region_cols, stroke = 0)
ggsave("/path/to/Figures/SpatialMap_PathologyAnnotation.pdf")
#################################################################################
############################# Figure 2C #########################################
#################################################################################
# Perform Deconvolution
anchors <- FindTransferAnchors(reference = scRNA_Ref, query = BM1_filtered, normalization.method = "SCT", recompute.residuals = F)
predictions.assay <- TransferData(anchorset = anchors, refdata = scRNA_Ref$class2, prediction.assay = TRUE,
                                  weight.reduction = BM1_filtered[["pca"]], dims = 1:25)
BM1_filtered[["predictions"]] <- predictions.assay

# Spatial Map for AML Cell Deconvolution
SpatialFeaturePlot(BM1_filtered, features = "AML", crop = F, pt.size.factor = 1.25, stroke = 0)
ggsave("/path/to/Figures/SpatialMap_AML_Deconvolution.pdf")

# Spatial Map for Eryhtoid Cell Deconvolution
SpatialFeaturePlot(BM1_filtered, features = "Late Erythroid", crop = F, pt.size.factor = 1.25, stroke = 0)
ggsave("/path/to/Figures/SpatialMap_Erythroid_Deconvolution.pdf")
#################################################################################
############################# Figure 2D #########################################
#################################################################################

# Create Heatmap for putative markers and compare with Unsupervised Clusters
ht <- GroupHeatmap(BM1_filtered, group.by = "Regions", cell_annotation = "seurat_clusters", 
                   features = c("HBB", "HBD", "HBA1", "HBA2", "GATA1", "GATA2",
                                "S100A12", "FCGR3A", "CD14", "CD33", "MS4A7"),
                   group_palcolor = list(region_cols), add_dot = T, add_bg = T,
                   flip = T, cell_annotation_palcolor = list(c("#FFC312", "#C4E538", "#12CBC4")),
                   show_column_names = T, cell_annotation_params = list(width = unit(12, "mm")))
ht$plot
panel_fix(ht$plot, height = 4.5, width = 7, dpi = 500, raster = F, save = "/path/to/Figures/PutativeMarkers_HeatPlot.pdf")

#################################################################################
############################# Figure 2E #########################################
#################################################################################

# Define Samples
BM1_filtered$Sample <- "BM1"
BM2_filtered$Sample <- "BM2"
BM1_filtered$Patient <- "PT1"
BM2_filtered$Patient <- "PT2"

# Store Cell Deconvolution results as feature
AddFeaturesFromAssay <- function(srt, assay, prefix){
  metadata <- srt[[assay]]@data
  transposed_data <- as.data.frame(t(metadata))
  rownames(transposed_data) <- gsub("\\.", "-", rownames(transposed_data))
  colnames(transposed_data) <- paste0(prefix, colnames(transposed_data))
  srt <- AddMetaData(srt, metadata = transposed_data)
  return(srt)
}
BM1_filtered <- AddFeaturesFromAssay(BM1_filtered, assay = "predictions", prefix = "Seurat_")
BM2_filtered <- AddFeaturesFromAssay(BM2_filtered, assay = "predictions", prefix = "Seurat_")

# Merge samples
BMs <- merge(BM1_filtered, y=BM2_filtered, add.cell.ids = c("BM1", "BM2"), project = "BoneMarrows")

# Create vector for cell type deconvolution scores stored as feature
cell_types_ordered <- c(
  "Seurat_HSC",
  "Seurat_AML",
  "Seurat_GMP",
  "Seurat_Monocytes",
  "Seurat_cDC",
  "Seurat_pDC",
  "Seurat_Platelet.Mega",
  "Seurat_CLP",
  "Seurat_B",
  "Seurat_Plasma",
  "Seurat_NK",
  "Seurat_UnconvT",
  "Seurat_CD8.Naive",
  "Seurat_CD4.Naive",
  "Seurat_CD4.Treg",
  "Seurat_CD4.Mem",
  "Seurat_CD4.Eff",
  "Seurat_CD8.Eff",
  "Seurat_CD8.Mem",
  "Seurat_Early.Erythoid",
  "Seurat_Late.Erythroid"
)

# Create DotPlot comparison
ht <- GroupHeatmap(BMs, features = cell_types_ordered, group.by = "Sample", add_dot = T, flip = T, show_column_names = T, column_names_side = "top", show_row_names = T, group_palcolor = list(c("#56CCF2", "#6FCF97")), dot_size = unit(12, "mm"))
# Print and save the plot
ht
ggsave("/path/to/Figures/BoneMarrows_DotPlot_CellTypes.pdf", height = 5, width = 10)
#################################################################################
############################# Figure 2F #########################################
#################################################################################
EM1_filtered$seurat_clusters <- factor(EM1_filtered$seurat_clusters,levels=names(seurat_cols))
# Plot Spatial map for Unsupervised Clusters
SpatialDimPlot(EM1_filtered, group.by = "seurat_clusters", cols = seurat_cols, crop = F, pt.size.factor = 1.25, stroke = 0)
ggsave("/path/to/Figures/SpatialMap_UnsupervisedClusters_EM1.pdf")
# Import Pathology Annotation that defined with Loupe Browser
annotations <- read.csv("/path/to/Pathology/EM1_Pathology.csv", stringsAsFactors = FALSE)
# Add meta data of Seurat object
EM1_filtered[["Pathology"]] <- annotations$Pathology[match(colnames(EM1_filtered), annotations$Barcode)]
# Define colors for Pathology Annotation
em1_cols <- c("Sarcoma" = "#ff793f", "Dermis" = "#218c74", "Epidermis" = "#B53471", "Gland" = "#1f77b4")
# Pathology names changed during Adobe Illustrator editing for better representation
EM1_filtered$Pathology <- factor(EM1_filtered$Pathology,levels=names(em1_cols))
# Spatial map for Pathology Annotation
SpatialDimPlot(EM1_filtered, group.by = "Pathology", crop = F, pt.size.factor = 1.25,cols = region_cols, stroke = 0)
ggsave("/path/to/Figures/SpatialMap_PathologyAnnotation_EM1.pdf")
#################################################################################
############################# Figure 2G #########################################
#################################################################################
library(SpaCET)
source("/path/to/Customization/utilities.R")

# Convert EM Seurat object to SpaCET object
EM1_spacet <- convert.Seurat(EM1_filtered)
# calculate the QC metrics
EM1_spacet <- SpaCET.quality.control(EM1_spacet)
# plot the QC metrics
SpaCET.visualize.spatialFeature(
  EM1_spacet, 
  spatialType = "QualityControl", 
  spatialFeatures=c("UMI","Gene"),imageBg = F
)
# Deconvolve EM data with Pan-Cancer dictionary from SpaCET
EM1_spacet <- SpaCET.deconvolution(EM1_spacet, cancerType="PANCAN", coreNo=4)
# show the spatial distribution of malignant cells and macrophages.
SpaCET.visualize.spatialFeature(
  EM1_spacet, 
  spatialType = "CellFraction", 
  spatialFeatures=c("Malignant","Macrophage"), imageBg = F
)
# show the spatial distribution of all cell types.
SpaCET.visualize.spatialFeature(
  EM1_spacet, 
  spatialType = "CellFraction", 
  spatialFeatures="All", 
  sameScaleForFraction = TRUE,
  pointSize = 0.1, 
  nrow=5, imageBg = T
)
# malignant cell states deconvolution
EM1_spacet <- SpaCET.deconvolution.malignant(EM1_spacet, coreNo = 4)

# Add deconvolution results to Seurat Object
EM1_filtered <- addTo.Seurat(EM1_spacet, EM1_filtered)

# Define Assay as this deconvolution 
DefaultAssay(EM1_filtered) <- "propMatFromSpaCET"

# Spatial Map
SpatialFeaturePlot(EM1_filtered, features = "Malignant", crop = F, pt.size.factor = 1.25,cols = region_cols, stroke = 0)
ggsave("/path/to/Figures/EM1_MalignantCellDeconvolution.pdf")
#################################################################################
############################# Figure 2H #########################################
#################################################################################

# Define Assay as SCTransformed data
DefaultAssay(EM1_filtered) <- "SCT"

# Create Heatmap and Pie Chart
ht <- GroupHeatmap(EM1_filtered, group.by = "Pathology", cell_annotation = "seurat_clusters",
                   features = c("CD33", "CD34", "FCGR1A", # Leukemia
                                "COL1A1", "FBN1", "VWF", #Dermis
                                "KRT14", "FLG", "LOR", # Epidermis
                                "KRT79", "KRT7", "ATP1A1"), # Gland
                   group_palcolor = list(em1_cols), add_dot = T, add_bg = T,
                   flip = T, cell_annotation_palcolor = list(c("#FFC312", "#C4E538", "#12CBC4")),
                   show_column_names = T, cell_annotation_params = list(width = unit(12, "mm")))

ht$plot
panel_fix(ht$plot, height = 4.15, width = 8, dpi = 500, raster = F, save = "/path/to/Figures/EM1_PutativeMarkersHeatMap.pdf")