library(tidyverse)
library(reshape2)
library(Seurat)
library(Signac)
library(clustree)
library(ggrepel)
library(ggpubr)
library(CellChat)

# Figure 6A ---------------------------------------------------------------

DefaultAssay(full) <- 'RNA'
Idents(full) <- 'BroadCellTypes'
full.immune <- subset(full, subset = CellTypes %in% c('Immune'))

# Remove cells with low expression of CD45
mat <- GetAssayData(full.immune, slot = 'data', assay = 'RNA')
hist(mat['Ptprc',])
hist(full.immune$CD45)
full.immune[['CD45']] <- PercentageFeatureSet(full.immune, features = 'Ptprc')
full.immune <- subset(full.immune, subset = CD45 > 0.1)

# RNA analysis
DefaultAssay(full.immune) <- "RNA"
full.immune <- SCTransform(full.immune, vars.to.regress = 'percent.mt', verbose = FALSE)
full.immune <- RunPCA(full.immune, verbose = FALSE)
ElbowPlot(full.immune, ndims = 50, reduction = "pca")
full.immune <- RunUMAP(full.immune, dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
DefaultAssay(full.immune) <- "ATAC"
full.immune <- RunTFIDF(full.immune)
full.immune <- FindTopFeatures(full.immune, min.cutoff = 'q0')
full.immune <- RunSVD(full.immune)
full.immune <- RunUMAP(full.immune, reduction = 'lsi', dims = 2:30, 
                            reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# Multimodal analysis
full.immune <- FindMultiModalNeighbors(full.immune, reduction.list = list("pca", "lsi"),
                                            dims.list = list(1:30, 2:30))
full.immune <- RunUMAP(full.immune, nn.name = "weighted.nn", reduction.name = "wnn.umap",
                            reduction.key = "wnnUMAP_")
full.immune <- FindClusters(full.immune, graph.name = "wsnn", algorithm = 3, verbose = FALSE,
                                 resolution = seq(0.1,1,0.1))
clustree(full.immune, prefix = "wsnn_res.")

table(Idents(full.immune), full.immune$orig.ident)
# Remove clusters with <15 cells and recluster


DimPlot(full.immune, reduction = 'wnn.umap', group.by = 'wsnn_res.0.1',
        label = TRUE, repel = TRUE, label.box = TRUE)

# Normalize RNA assay
DefaultAssay(full.immune) <- 'RNA'
full.immune <- NormalizeData(full.immune)
full.immune <- FindVariableFeatures(full.immune, selection.method = "vst", nfeatures = 2000)
full.immune <- ScaleData(full.immune, vars.to.regress = 'percent.mt')

# Visualize dotplot of differentiation axis
axis_genes <- c('Cd247','Cd3g','Cd3e','Icos','Il7r','Cd28','Lef1','Ms4a4b',
                'Il1b','Lyz2','Csf1r','Adgre1','Mafb','Cd80','Itgam',
                'Cd79a','Btla','Ebf1','Bank1','Mef2c')
Idents(full.immune) <- 'wsnn_res.0.1'
DotPlot(object = full.immune, features = axis_genes, cols = c('navy','red'))+
  coord_flip() + RotatedAxis()

# Rename clusters
Idents(full.immune) <- "wsnn_res.0.1"
full.immune@meta.data$CellTypes <- Idents(full.immune)
Idents(full.immune) <- "CellTypes"
full.immune <- RenameIdents(full.immune,
                                 '0'='T cells', '3'='T cells',
                                 '1'='Myeloid cells', '2'='B cells')
full.immune@meta.data$CellTypes <- Idents(full.immune)

# Visualize
my_cols <- c('chartreuse3','dodgerblue','gold')
DimPlot(full.immune, reduction = 'wnn.umap', group.by = 'CellTypes', 
        cols = my_cols, label = TRUE, repel = TRUE, label.box = TRUE, label.size = 6)


# Figure 6B ---------------------------------------------------------------

# Extract and organize metadata
data <- FetchData(
  object = full.immune,
  vars = c("CellTypes", "Condition", "orig.ident")
)

data <- data %>%
  dplyr::count(orig.ident, Condition, CellTypes, name = "n") %>%
  group_by(orig.ident, Condition) %>%
  tidyr::complete(CellTypes, fill = list(n = 0)) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  group_by(Condition, CellTypes)

# Set plotting order
bar.data$Condition <- factor(
  bar.data$Condition,
  levels = c("rOVA+4NQO","4NQO","rOVA","Vehicle")   
)

bar.data$CellTypes <- factor(
  bar.data$CellTypes,
  levels = c("B cells","Myeloid cells","T cells")
)

# Plot
ggplot(bar.data,
       aes(x = Condition,
           y = mean.percent,
           fill = CellTypes)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "% of cells"
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )+
  scale_fill_manual(values = c('gold','dodgerblue','chartreuse3'),
                    guide = guide_legend(reverse = TRUE))


# Figure 6E ---------------------------------------------------------------

### CellChat analysis

# Create CellChat objects
DefaultAssay(full.immune) <- 'RNA'
Idents(full.immune) <- 'Condition'
full.immune.4NQO <- subset(full.immune, subset = Condition %in% '4NQO')
full.immune.rOVA4NQO <- subset(full.immune, subset = Condition %in% 'rOVA4NQO')

cellchat.full.immune.4NQO <- createCellChat(object = full.immune.4NQO, group.by = "CellTypes")
cellchat.full.immune.rOVA4NQO <- createCellChat(object = full.immune.rOVA4NQO, group.by = "CellTypes")

#  Extract metadeta
cellchat.full.immune.4NQO@meta$condition <- "full.immune"
cellchat.full.immune.rOVA4NQO@meta$condition <- "full.immune"

# Establish database
CellChatDB <- CellChatDB.mouse
cellchat.full.immune.4NQO@DB <- CellChatDB
cellchat.full.immune.rOVA4NQO@DB <- CellChatDB

# Subset the expression data of signaling genes
cellchat.full.immune.4NQO <- subsetData(cellchat.full.immune.4NQO)
cellchat.full.immune.rOVA4NQO <- subsetData(cellchat.full.immune.rOVA4NQO)

# Pre-process the expression data
cellchat.full.immune.4NQO <- identifyOverExpressedGenes(cellchat.full.immune.4NQO)
cellchat.full.immune.4NQO <- identifyOverExpressedInteractions(cellchat.full.immune.4NQO)
cellchat.full.immune.rOVA4NQO <- identifyOverExpressedGenes(cellchat.full.immune.rOVA4NQO)
cellchat.full.immune.rOVA4NQO <- identifyOverExpressedInteractions(cellchat.full.immune.rOVA4NQO)

# Project gene expression data onto protein-protein interaction (PPI)
cellchat.full.immune.4NQO <- smoothData(cellchat.full.immune.4NQO, adj = PPI.mouse)
cellchat.full.immune.rOVA4NQO <- smoothData(cellchat.full.immune.rOVA4NQO, adj = PPI.mouse)

# Compute the communication probability and infer cellular communication network
cellchat.full.immune.4NQO <- computeCommunProb(cellchat.full.immune.4NQO, type = "triMean")
cellchat.full.immune.4NQO <- computeCommunProbPathway(cellchat.full.immune.4NQO)
cellchat.full.immune.rOVA4NQO <- computeCommunProb(cellchat.full.immune.rOVA4NQO, type = "triMean")
cellchat.full.immune.rOVA4NQO <- computeCommunProbPathway(cellchat.full.immune.rOVA4NQO)

# Calculate the aggregated cell-cell communication network
cellchat.full.immune.4NQO <- aggregateNet(cellchat.full.immune.4NQO)
groupSize.4NQO <- as.numeric(table(cellchat.full.immune.4NQO@idents))
cellchat.full.immune.rOVA4NQO <- aggregateNet(cellchat.full.immune.rOVA4NQO)
groupSize.rOVA4NQO <- as.numeric(table(cellchat.full.immune.rOVA4NQO@idents))

# Compute the network centrality scores
cellchat.full.immune.4NQO <- netAnalysis_computeCentrality(cellchat.full.immune.4NQO, slot.name = "netP")
cellchat.full.immune.rOVA4NQO <- netAnalysis_computeCentrality(cellchat.full.immune.rOVA4NQO, slot.name = "netP")

# Merge objects
object.list <- list(NQO = cellchat.full.immune.4NQO, rOVA4NQO = cellchat.full.immune.rOVA4NQO)
cellchat.full.immune.4NQO <- mergeCellChat(object.list, add.names = names(object.list))

# Show differential number of interactions among cell types between conditions
my_cols <- c('chartreuse3','dodgerblue','gold')
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat.full.immune.4NQO, comparison = c(1, 2), measure = "count", 
                          weight.scale = T, arrow.size = 0.1, color.use = my_cols)


# Figure 6F ---------------------------------------------------------------

# Subset 4NQO(+) groups
Idents(full.immune) <- 'Condition'
immune_4NQO <- subset(full.immune, subset = Condition %in% c('4NQO','rOVA+4NQO'))

# Visualize
DotPlot(object = immune_4NQO, features = c('Il13','Il4','Il10','Maf',
                                           'Il2','Ifng','Tbx21','Stat4',
                                           'Cd4','Cd8b1','Cd8a',
                                           'Cd247','Cd3g','Cd3e','Cd3d'), 
        cols = c('navy','red'))+ coord_flip() + RotatedAxis() + 
  scale_size(range = c(1, 6)) +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))

