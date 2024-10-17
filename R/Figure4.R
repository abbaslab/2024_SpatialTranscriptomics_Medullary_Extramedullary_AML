########################## Figure 4 Codes ######################################
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
############################# Figure 4A #########################################
#################################################################################

#################################################################################
############################# Figure 4B #########################################
#################################################################################
# Function to add composite scores and classifications to Seurat object
addCompositeScores <- function(seurat_obj, pathways_list, image_name) {
  # Extract spatial coordinates from the specified image
  coords <- seurat_obj@images[[image_name]]@coordinates
  
  # Extract pathway scores from metadata
  pathways <- lapply(pathways_list, function(x) seurat_obj@meta.data[[x]])
  
  # Normalize and combine pathway scores
  normalized_pathways <- lapply(pathways, function(x) rescale(x, to = c(0, 1)))
  composite_score <- rowMeans(do.call(cbind, normalized_pathways))
  
  # Combine coordinates and composite scores into a data frame
  data <- data.frame(
    x = coords$imagerow,
    y = coords$imagecol,
    composite_score = composite_score
  )
  
  # Apply Jenks Natural Breaks to the composite score
  jenks_breaks <- classInt::classIntervals(data$composite_score, n = 4, style = "jenks")
  data$composite_score_category <- cut(data$composite_score, 
                                       breaks = jenks_breaks$brks, 
                                       labels = c("Low", "Medium-Low", "Medium-High", "High"))
  
  # Ensure data is in the same order as the cells in the Seurat object
  data$cell <- rownames(coords)
  rownames(data) <- data$cell
  data <- data[rownames(seurat_obj@meta.data), ]
  
  # Add composite scores and classifications to the Seurat object metadata
  seurat_obj@meta.data$composite_score <- data$composite_score
  seurat_obj@meta.data$composite_score_category <- data$composite_score_category
  
  # Plot spatial distribution of composite scores
  p1 <- SpatialFeaturePlot(seurat_obj, features = "composite_score") + 
    ggtitle("Spatial Distribution of Composite Scores")
  
  p2 <- SpatialDimPlot(seurat_obj, group.by = "composite_score_category") + 
    ggtitle("Spatial Distribution of Composite Score Categories")
  
  return(list(seurat_obj = seurat_obj, composite_plot = p1, category_plot = p2, data = data))
}
# Example usage
pathways_list <- c("HALLMARK_INFLAMMATORY_RESPONSE", 
                   "HALLMARK_IL6_JAK_STAT3_SIGNALING", 
                   "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
                   "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
                   "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
                   "HALLMARK_COMPLEMENT", 
                   "HALLMARK_IL2_STAT5_SIGNALING")
# Run function with specified pathway list
result <- addCompositeScores(BM1_filtered, pathways_list, image_name = "slice1")
#result <- addCompositeScores(EM1_filtered, pathways_list, image_name = "slice1")
# Access the updated Seurat object
BM1_filtered <- result$seurat_obj

# Spatial plot
# Rotate coordinates for 90-degree rotation
data_rotated <- data.frame(
  x = result$data$y,
  y = -result$data$x,
  composite_score = data$composite_score,
  composite_score_category = data$composite_score_category
)
# Define a bounding box for the top right part (rotated coordinates)
top_right_rotated <- data_rotated[data_rotated$x > quantile(data_rotated$x, 0.75) & data_rotated$y > quantile(data_rotated$y, 0.75), ]
# Create the plot with RdBu color scheme and transparency
spatial_plot <- ggplot(data_rotated, aes(x = x, y = y, fill = composite_score_category)) +
  geom_tile(alpha = 0.8, width = 15, height = 15) +  # Increase the size of the tiles and add transparency
  scale_fill_manual(values = c("Low" = "#2166AC", "Medium-Low" = "#67A9CF", "Medium-High" = "#D1E5F0", "High" = "#B2182B")) +  # Adjusting "Medium-High" color
  labs(title = "Composite Inflammatory Pathway Activity",
       fill = "Composite Score") +
  theme_minimal() +
  coord_fixed() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  annotate("text", x = mean(top_right_rotated$x), y = mean(top_right_rotated$y), label = "High Pathway Activity", color = "black", size = 5, fontface = "bold")
# Print
spatial_plot
# Save the plot
ggsave("/path/to/Figures/SpatialInflammationClass.pdf")
#################################################################################
############################# Figure 4C #########################################
#################################################################################
library(CellChat)

# Assign Labels function to label spots with prediction
assignLabels <- function(object, prediction = "predictions") {
  pred <- object[[prediction]]@data
  pred <- pred[1:(nrow(pred)-1), ]
  # label each spot based on the maximum prediction probability
  labels = rownames(pred)[apply(pred, 2, which.max)]
  names(labels) <- colnames(pred)
  object$labels <- factor(labels)
  Idents(object) <- "labels"
  return(object)
}

# Run the function
BM1_filtered <- assignLabels(BM1_filtered, prediction = "predictions")
# Create data input based on normalized data
data.input = Seurat::GetAssayData(BM18_filtered, layer = "data", assay = "SCT")
# Define metadata
meta = data.frame(labels = Seurat::Idents(BM18_filtered), samples = "labels", row.names = names(Seurat::Idents(BM18_filtered))) # manually create a dataframe consisting of the cell labels
meta$samples <- factor(meta$samples)
unique(meta$labels) # check the cell labels
unique(meta$samples)
# load spatial transcriptomics information
spatial.locs = Seurat::GetTissueCoordinates(BM18_filtered, scale = NULL, cols = c("imagerow", "imagecol")) 
# Spatial factors of spatial coordinates
# For 10X Visium, the conversion factor of converting spatial coordinates from Pixels to Micrometers can be computed as the ratio of the theoretical spot size (i.e., 65um) over the number of pixels that span the diameter of a theoretical spot size in the full-resolution image (i.e., 'spot_diameter_fullres' in pixels in the 'scalefactors_json.json' file). 
scalefactors = jsonlite::fromJSON(txt = file.path("./data/BM-18-002348/outs/spatial/", 'scalefactors_json.json'))
spot.size = 65 # the theoretical spot size (um) in 10X Visium
conversion.factor = spot.size/scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
d.spatial <- computeCellDistance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)
min(d.spatial[d.spatial!=0]) # this value should approximately equal 100um for 10X Visium data

# Create CellChat Object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
#cellchat
# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# set the used database in the object
cellchat@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)
# Compute the communication probability and infer cellular communication network
ptm = Sys.time()

cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE, interaction.range = 250, scale.distance = 0.01,
                              contact.dependent = TRUE, contact.range = 100)
#Users can filter out the cell-cell communication if there are only few cells in certain cell groups. By default, the minimum number of cells required in each cell group for cell-cell communication is 10.
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
pathway_analysis <- netAnalysis_computeCentrality(cellchat)
pathway_analysis <- netAnalysis_signalingRole_network(pathway_analysis)
# Visualize signaling roles of different cell types
netVisual_signalingRole_scatter(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
cellchat <- aggregateNet(cellchat)
# Creating heatmap to show differential number/strength of interaction in the cell-cell communication network
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")
# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# Sources:
## 1 = AML
## 11 = GMP
## 14 = Monocytes
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:25), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(1,14), targets.use = c(1,11), remove.isolate = FALSE) #AML-Monocyte
netVisual_bubble(cellchat, sources.use = c(1,11,14), targets.use = c(1,11,14), remove.isolate = FALSE) #AML-Monocyte-GMP
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4,5,6,7), targets.use = 1, legend.pos.x = 15)
plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = TRUE, type = "violin")
plotGeneExpression(cellchat, signaling = "CD99", enriched.only = TRUE, type = "violin")
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
# Check interaction pathways
cellchat@netP$pathways
pathways.show <- c("CXCL") 

par(mfrow=c(1,1), xpd = TRUE) # `xpd = TRUE` should be added to show the title
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#
# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 1, vertex.size.max = 0.5, alpha.image = 0.2, vertex.label.cex = 3.5, point.size = 1.5, weight.scale = T)
# Save spatial plot for CXCL
ggsave("/Path/to/Figures/BM18_CXCLPathSpatial.pdf")
# Define the cell types of interest
cell_types_of_interest <- c("AML", "Monocytes", "GMP")


# Subset the communication data
cellchat_subset@net <- cellchat@net[cells_of_interest, cells_of_interest, ]

#...

#################################################################################
############################# Figure 4D #########################################
#################################################################################

#################################################################################
############################# Figure 4E #########################################
#################################################################################
library(AUCell)
# Load genesets for read GMT file
load_genesets <- function(path, ignore_cols = 2){
  x <- scan(path, what="", sep="\n")
  y <- strsplit(x, "\t")
  #y=lapply(y,str_remove_all,'\"')
  names(y) <- sapply(y, `[[`, 1)
  for(i in 1:ignore_cols){
    y <- lapply(y, `[`, -1)
  }
  return(y)
}
# Function for running AUCell with Visium data
run_aucell<-function(sobject,path_list){
  exprMatrix=sobject@assays$Spatial$counts
  cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=2, plotStats=TRUE,splitByBlocks=TRUE)
  for (i in 1:length(path_list)){
    geneset=load_genesets(path_list[[i]],ignore_cols=2)
    cells_AUC <- AUCell_calcAUC(geneset, cells_rankings)
    score_mat <- t(SummarizedExperiment::assay(cells_AUC, 'AUC'))
    colnames(score_mat)<-make.names(colnames(score_mat))
    sobject=AddMetaData(sobject, as.data.frame(score_mat))
  }
  return(sobject)
}
# Run function with Hallmark GMT file from MsigDB
BM18_filtered <- run_aucell(BM18_filtered, "/path/to/Data/h.all.v2023.1.Hs.symbols.gmt.txt")
# Define CXCL12 and CXCR4
CXCL_Path <- c("CXCL12", "CXCR4")
# Create list to calculate AUCell score
cxcl_list <- list(
  #ChrisTest = ChrisTest,
  CXCL_Path = CXCL_Path
  #HIF1A_Path = HIF1A_Path,
  #Adipo_MSC = Adipo_MSC
)

# Run for each gene list under pathway list (Here only showed CXCL_Path but the one can increase number of gene set with this approach)
for (gene_set_name in names(cxcl_list)) {
  BM18_filtered <- calculate_AUCell(BM18_filtered, cxcl_list, gene_set_name)
  cat("Processed AUCell for", gene_set_name, "\n")
}

# Fetch data for relevant genes and pathways
cxcl12_cxcr4_scores <- FetchData(BM18_filtered, vars = "AUCell_CXCL_Path")
pi3k_data <- FetchData(BM18_filtered, vars = "HALLMARK_PI3K_AKT_MTOR_SIGNALING")
ifn_gamma_data <- FetchData(BM18_filtered, vars = "HALLMARK_INTERFERON_GAMMA_RESPONSE")
emt_data <- FetchData(BM18_filtered, vars = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
inflammation <- FetchData(BM18_filtered, vars = "composite_score")

# Combine the data into a single data frame
combined_data <- data.frame(
  CXCL12_CXCR4 = cxcl12_cxcr4_scores,
  PI3K_AKT_mTOR = pi3k_data,
  IFN_Gamma = ifn_gamma_data,
  EMT = emt_data,
  inflammation_Score <- inflammation
)
# Ensure column names are correct
colnames(combined_data) <- c("CXCL12_CXCR4", "PI3K_AKT_mTOR", "IFN_Gamma", "EMT", "inflammation_Score")

# Check if combined_data is created correctly
print(head(combined_data))

# Ensure the data frame has the expected number of rows
print(nrow(combined_data))
# Calculate correlations and p-values
cor_pi3k <- cor.test(combined_data$CXCL12_CXCR4, combined_data$PI3K_AKT_mTOR, method = "pearson")
cor_ifn <- cor.test(combined_data$CXCL12_CXCR4, combined_data$IFN_Gamma, method = "pearson")
cor_emt <- cor.test(combined_data$CXCL12_CXCR4, combined_data$EMT, method = "pearson")
cor_pi3k <- cor.test(combined_data$inflammation_Score, combined_data$PI3K_AKT_mTOR, method = "pearson")
# Print results
cat("Correlation between CXCL12-CXCR4 and PI3K/AKT/mTOR:\n")
print(cor_pi3k)
cat("\nCorrelation between CXCL12-CXCR4 and IFN Gamma:\n")
print(cor_ifn)
cat("\nCorrelation between CXCL12-CXCR4 and EMT:\n")
print(cor_emt)
# Scatter plot for CXCL12-CXCR4 vs PI3K/AKT/mTOR with inflammation class colors
p1 <- ggplot(combined_data, aes(x = CXCL12_CXCR4, y = PI3K_AKT_mTOR, color = InflammationClass)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", size = 1) +
  scale_color_manual(values = inflammation_colors, name = "Inflammation Class") +
  xlab("CXCL12-CXCR4 Co-expression Score") +
  ylab("PI3K/AKT/mTOR Signaling") +
  ggtitle("CXCL12-CXCR4 vs PI3K/AKT/mTOR Modulated by Inflammation") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )
# Print the plot
print(p1)
ggsave("/path/to/Figures/CXCL12_PI3K_Inflammation.pdf", height = 6, width = 7)


#################################################################################
############################# Figure 4F #########################################
#################################################################################
# Correlation between CXCR4 and inflammation score on EM1

# Spatial plot for CXCR4 expression on EM1
SpatialFeaturePlot(EM1_filtered, features = "CXCR4", crop = F, pt.size.factor = 1.25, stroke = 0)
#################################################################################
############################# Figure 4G #########################################
#################################################################################

# Spatial plot for PI3K/Akt/mTOR pathway
SpatialFeaturePlot(BM1_filtered, features = "HALLMARK_PI3K_AKT_MTOR_SIGNALING", crop = F, pt.size.factor = 1.25, stroke = 0)
# Spatial plot for Trans-Differentiation pathway
SpatialFeaturePlot(BM1_filtered, features = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", crop = F, pt.size.factor = 1.25, stroke = 0)

#################################################################################
############################# Figure 4H #########################################
#################################################################################

# Markers from Desai et al., 2023
cd8_dysfunction=c("RPL41", "TXNIP", "CD2", "PRF1", "PSMB9", "LIMD2", "LTB", "GIMAP1", "GIMAP4", "CX3CR1", "TRAF3IP3", "GIMAP7", "GIMAP5", "DENND2D", "CLIC1", "MYL12A", "S100A11", "RAC2", "PSMB10", "ARPC1B", "CORO1A", "PCBP1", "UCP2", "PLEK", "CD8A")
senescence <- c("CD27", "CD28", "IL2","GZMB", "B3GAT1", "KLRG1", "KLRC1", "KLRC2", "IL6", "TNF", "IFNG", "TIGIT", "HAVCR2", "CCL5", "CCL16", "CCL23")
ex_papers <- c("HNF1A", "GZMB", "CD28", "CD69", "TOX", "CD38", "CTLA4", "ENTPD1", "HAVCR2", "PDCD1", "HLA-DRA", "MKI67", "IL2RB")
ex_markers <- c("TOX", "CD38", "LAG3", "CTLA4", "TIGIT", "HAVCR2", "PDCD1", "CXCL13", "TNFRSF9", "TOX2", "ENTPD2", "LAYN","EOMES","HLA-DRA", "GZMB")
exhaustion_markers <- unique(c(ex_markers,ex_papers))
# Markers for Treg from Braun et al., 2021
Treg <- c("FOXP3", "TNFRSF4", "TNFRSF18", "TBC1D4", "IL2RA", "CTLA4", "NGFRAP1", "CORO1B", "RTKN2", "LAYN", "AC017002.1", "CISH", "SLAMF1", "NCF4", "DNPH1", "F5", "LTB", "CTSC", "BATF", "ICA1", "TIGIT", "UGP2", "PIM2", "FBLN7", "CD4", "IL32", "IKZF2", "MIR4435-1HG", "GBP2", "CARD16", "PHTF2", "GPX1", "IL1R2", "GBP5", "S100A4", "PBXIP1", "GLRX", "CLEC7A", "TBC1D8", "SPOCK2", "RPS12", "RPS27", "RPL23A", "ZFP36L2", "RARRES3", "CD99", "GIMAP1", "APMAP", "LITAF", "CLEC2B", "APOBEC3G", "GSTP1", "GIMAP5", "ANKRD32", "LYST", "MCTP2", "AC069363.1", "THEMIS", "NELL2", "CRTAM", "TC2N", "F2R", "EOMES", "PRF1", "PLEK", "GIMAP7", "AOAH", "KLRG1", "GIMAP4", "ITM2C", "CMC1", "KLRD1", "GPR56", "IFNG", "VCAM1", "LAG3", "HOPX", "ANXA1", "GZMB", "CST7", "CCL3", "GZMH", "CTSW", "CD8B", "GZMA", "GZMK", "CD8A", "CCL4", "CCL5", "NKG7", "RPS15A", "FASLG", "PTGDR", "RPS27A", "RPS4Y1", "RPS27L", "PARP8", "MGAT4A", "ACTB", "LYAR")

# Create list for each gene set
T_groups <- list(
  Poonam_exhaustion = exhaustion_markers,
  CD8_Dysfunction = cd8_dysfunction,
  T_Senescence =senescence
)
# Run for loop below for each sample
for (gene_set_name in names(T_groups)) {
  BM1_filtered <- calculate_AUCell(BM1_filtered, T_groups, gene_set_name)
  cat("Processed AUCell for", gene_set_name, "\n")
}
# Run script below for each sample to get HeatDotPlot
SCP::GroupHeatmap(S18_NONa, group.by = "composite_score_category_v2", features =  c("AUCell_Poonam_exhaustion", "AUCell_CD8_Dysfunction", "AUCell_T_Senescence", "AUCell_Treg"), exp_method = "zscore", row_names_side = "left", show_row_names = T, add_dot = T, add_bg = T, group_palcolor = list(c("Low" = "#2166AC", "Medium-Low" = "#67A9CF", "Medium-High" = "#F9C6AB", "High" = "#B2182B", "Unkown" = "#B2182B")), limits = c(-1.5,1.5))

