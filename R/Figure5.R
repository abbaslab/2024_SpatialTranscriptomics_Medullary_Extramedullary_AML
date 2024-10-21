########################## Figure 5 Codes ######################################
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

# Subset only AML cells
aml <- subset(scRNA_Ref, subset = class2 == "AML")
# function for Seurat pipelin
analyze_seurat <- function(seurat_object, number_features=2500, res=0.8, md=0.01, k=20){
  seurat_object <- FindVariableFeatures(seurat_object,nfeatures = number_features)
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(object = seurat_object, verbose = FALSE,features = VariableFeatures(seurat_object),npcs=50)
  seurat_object <- RunHarmony(seurat_object, group.by.vars="orig.ident",assay.use ="RNA", max.iter.cluster=500)
  seurat_object <- FindNeighbors(seurat_object, dims=1:30,reduction = "harmony",k.param = k)
  seurat_object <- FindClusters(seurat_object, resolution = seq(0,1,0.1), random.seed=123)
  seurat_object <- FindClusters(object = seurat_object,resolution=res,random.seed=123,graph.name = 'SCT_snn')
  seurat_object <- RunUMAP(seurat_object,reduction = "harmony",seed.use = 123,dims=1:50,
                           a=0.5,b=1.5, verbose=FALSE)
  return(seurat_object)
}
# run the function
aml <- analyze_seurat(aml)
# Fine tune UMAP
aml <- RunUMAP(aml,reduction = "harmony",seed.use = 123,dims=1:50,
               a=3,b=1.3, verbose=FALSE)

# Subset deconvolved Visium sample based on AML score
BM1_AML <- subset(BM1_filtered, subset = AML > median(BM1_filtered$Seurat_AML))

#################################################################################
############################# Figure 5A #########################################
#################################################################################

#################################################################################
############################# Figure 5B #########################################
#################################################################################

#################################################################################
############################# Figure 5C #########################################
#################################################################################

#################################################################################
############################# Figure 5D #########################################
#################################################################################
DefaultAssay(BM18_AML) <- "SCT"

anchors <- FindTransferAnchors(reference = aml, query = BM18_AML, normalization.method = "SCT", recompute.residuals = F)
predictions.assay <- TransferData(anchorset = anchors, refdata = aml$class4, prediction.assay = TRUE,
                                  weight.reduction = BM18_AML[["pca"]], dims = 1:25)
BM18_AML[["predictions4"]] <- predictions.assay

DefaultAssay(BM18_AML) <- "predictions4"

SpatialFeaturePlot(BM18_AML, features = rownames(BM18_AML))
ggsave("./figs/Playground/Biological/PT2/BM18/Predictions4_SpatialPlot.pdf")

#################################################################################
############################# Figure 5E #########################################
#################################################################################
# Create the violin plot
p1 <- FeatureStatPlot(S18_Leukemia, group.by = "composite_score_category", stat.by = "Seurat_Committed-like", palcolor = list(c("Low" = "#2166AC", "Medium-Low" = "#67A9CF", "Medium-High" = "#F9C6AB", "High" = "#B2182B")), add_box = T, plot_type = "violin")
# Save the 
panel_fix(p1, height = 4, width = 3, dpi = 300, save = "/path/to/Figures/EM1_CommittedViolin.pdf")

#################################################################################
############################# Figure 5F #########################################
#################################################################################

#################################################################################
############################# Figure 5G #########################################
#################################################################################


