library(tidyverse)
library(reshape2)
library(Seurat)
library(Signac)
library(clustree)
library(ggrepel)
library(ggpubr)
library(miloR)
library(SeuratWrappers)
library(monocle3)
library(CellChat)


# Figure 2A ---------------------------------------------------------------

DefaultAssay(full) <- 'RNA'
Idents(full) <- 'BroadCellTypes'
full.epi <- subset(full, subset = CellTypes %in% c('Epithelial'))
Idents(full) <- 'Condition'
full.epi.no4NQO <- subset(full.epi, 
                          subset = Condition %in% c('Vehicle','rOVA'))

# RNA analysis
DefaultAssay(full.epi.no4NQO) <- "RNA"
full.epi.no4NQO <- SCTransform(full.epi.no4NQO, vars.to.regress = 'percent.mt', verbose = FALSE)
full.epi.no4NQO <- RunPCA(full.epi.no4NQO, verbose = FALSE)
ElbowPlot(full.epi.no4NQO, ndims = 50, reduction = "pca")
full.epi.no4NQO <- RunUMAP(full.epi.no4NQO, dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
DefaultAssay(full.epi.no4NQO) <- "ATAC"
full.epi.no4NQO <- RunTFIDF(full.epi.no4NQO)
full.epi.no4NQO <- FindTopFeatures(full.epi.no4NQO, min.cutoff = 'q0')
full.epi.no4NQO <- RunSVD(full.epi.no4NQO)
full.epi.no4NQO <- RunUMAP(full.epi.no4NQO, reduction = 'lsi', dims = 2:30, 
                           reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# Multimodal analysis
full.epi.no4NQO <- FindMultiModalNeighbors(full.epi.no4NQO, reduction.list = list("pca", "lsi"),
                                           dims.list = list(1:30, 2:30))
full.epi.no4NQO <- RunUMAP(full.epi.no4NQO, nn.name = "weighted.nn", reduction.name = "wnn.umap",
                           reduction.key = "wnnUMAP_")
full.epi.no4NQO <- FindClusters(full.epi.no4NQO, graph.name = "wsnn", algorithm = 3, verbose = FALSE,
                                resolution = seq(0.1,1,0.1))
clustree(full.epi.no4NQO, prefix = "wsnn_res.")

table(Idents(full.epi.no4NQO), full.epi.no4NQO$orig.ident)
# Remove clusters with <15 cells and recluster


my_cols <- c('mediumorchid1','chartreuse3','thistle','sandybrown','aquamarine',
             'dodgerblue','salmon','burlywood','lightgreen','slateblue1','cornflowerblue',
             'gold','firebrick','darkcyan')
DimPlot(full.epi.no4NQO, reduction = 'wnn.umap', group.by = 'wsnn_res.0.4', 
        cols = my_cols, label = TRUE, repel = TRUE, label.box = TRUE)

# Normalize RNA assay
DefaultAssay(full.epi.no4NQO) <- 'RNA'
full.epi.no4NQO <- NormalizeData(full.epi.no4NQO)
full.epi.no4NQO <- FindVariableFeatures(full.epi.no4NQO, selection.method = "vst", nfeatures = 2000)
full.epi.no4NQO <- ScaleData(full.epi.no4NQO, vars.to.regress = 'percent.mt')

# Visualize dotplot of differentiation axis
axis_genes <- c('Krtdap','Krt4','Krt13','Fabp5',  
                'Col17a1','Krt5','Krt14', 'Krt15',  
                'Mki67','Top2a')   
Idents(full.epi.no4NQO) <- 'wsnn_res.0.4'
levels(full.epi.no4NQO) <- c('13','8','4','10','15','14','9','1','0','7','2','6',
                             '12','3','5','11')
DotPlot(object = full.epi.no4NQO, features = axis_genes, cols = c('navy','red'))+
  coord_flip() + RotatedAxis()

# Rename clusters
Idents(full.epi.no4NQO) <- "wsnn_res.0.4"
full.epi.no4NQO@meta.data$CellStates <- Idents(full.epi.no4NQO)
Idents(full.epi.no4NQO) <- "CellStates"
full.epi.no4NQO <- RenameIdents(full.epi.no4NQO,
                                '13'='EP1', '8'='EP2', '4'='EP3', '10'='EP4',
                                '9'='EP5', '1'='EP6', '0'='EP7', '7'='EP8',
                                '2'='EP9', '12'='EP10', '6'='EP11',
                                '3'='EP12', '5'='EP13', '11'='EP14')
full.epi.no4NQO@meta.data$CellStates <- Idents(full.epi.no4NQO)

# Visualize
my_cols <- c('salmon','lightgreen','mediumpurple1','sandybrown','aquamarine',
             'deeppink2','chartreuse3','orchid2','khaki','powderblue',
             'dodgerblue','seagreen','thistle','chocolate2')
DimPlot(full.epi.no4NQO, reduction = 'wnn.umap', group.by = 'CellStates', 
        cols = my_cols, label = TRUE, repel = TRUE, label.box = TRUE, label.size = 6)


# Figure 2B ---------------------------------------------------------------

### RNA

DefaultAssay(full.epi.no4NQO) <- 'RNA'
Idents(full.epi.no4NQO) <- 'CellStates'
topmarkers <- FindAllMarkers(full.epi.no4NQO, assay = "RNA")
top100 <- topmarkers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC)

conditions_avex <- AggregateExpression(full.epi.no4NQO, return.seurat = T, assays = 'RNA')
conditions_avex <- ScaleData(conditions_avex, assay = 'RNA', features = top100$gene)
DoHeatmap(object = conditions_avex, features = top100$gene, assay = 'RNA', 
          group.colors = c('salmon','lightgreen','mediumpurple1','sandybrown','aquamarine',
                           'deeppink2','chartreuse3','orchid2','khaki','powderblue',
                           'dodgerblue','seagreen','thistle','chocolate2'),
          label = T, size = 3, angle = 30, draw.lines = F)+
  guides(colour=FALSE)+ 
  theme(axis.text.y = element_blank())


### ATAC

DefaultAssay(full.epi.no4NQO) <- 'ATAC'
Idents(full.epi.no4NQO) <- 'CellStates'
topmarkers <- FindAllMarkers(full.epi.no4NQO, assay = "ATAC")
top100 <- topmarkers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC)

conditions_avex <- AggregateExpression(full.epi.no4NQO, return.seurat = T, assays = 'ATAC')
conditions_avex <- ScaleData(conditions_avex, assay = 'ATAC', features = top100$gene)
DoHeatmap(object = conditions_avex, features = top100$gene, assay = 'ATAC', 
          group.colors = c('salmon','lightgreen','mediumpurple1','sandybrown','aquamarine',
                           'deeppink2','chartreuse3','orchid2','khaki','powderblue',
                           'dodgerblue','seagreen','thistle','chocolate2'),
          label = T, size = 3, angle = 30, draw.lines = F)+
  guides(colour=FALSE)+ 
  theme(axis.text.y = element_blank())



# Figure 2C ---------------------------------------------------------------

FeaturePlot(full.epi.no4NQO, reduction = "wnn.umap", features = c("Top2a","Krt14","Fabp5"))

# Figure 2D ---------------------------------------------------------------

axis_genes <- c('Krtdap','Krt4','Krt13','Fabp5',  
                'Col17a1','Krt5','Krt14', 'Krt15',  
                'Mki67','Top2a')   
Idents(full.epi.no4NQO) <- 'CellStates'
levels(full.epi.no4NQO) <- c('EP1','EP2','EP3','EP4','EP5','EP6','EP7',
                             'EP8','EP9','EP10','EP11','EP12','EP13','EP14')
DotPlot(object = full.epi.no4NQO, features = axis_genes, cols = c('navy','red'))+
  coord_flip() + RotatedAxis()


#  Figure 2E --------------------------------------------------------------

### Milo analysis

# Convert Seurat -> SingleCellExperiment
sce <- as.SingleCellExperiment(full.epi.no4NQO)


# Extract metadata & embeddings
colData(sce)$CellStates <- full.epi.no4NQO$CellStates

colData(sce)$Condition <-
  factor(full.epi.no4NQO$Condition,
         levels = c("Vehicle","rOVA"))

colData(sce)$orig.ident <- full.epi.no4NQO$orig.ident

rna.pca <- Embeddings(
  full.epi.no4NQO,
  "pca"
)[,1:50]

atac.lsi <- Embeddings(
  full.epi.no4NQO,
  "lsi"
)[,2:50]

wnn.latent <- cbind(
  rna.pca,
  atac.lsi
)

dim(wnn.latent)

reducedDim(
  sce,
  "wnn_latent"
) <- wnn.latent


# Create Milo object
milo <- Milo(sce)


# Build KNN graph & neighborhoods
milo <- buildGraph(
  milo,
  k = 30,
  d = 50,
  reduced.dim = "wnn_latent"
)

milo <- makeNhoods(
  milo,
  prop = 0.1,
  k = 30,
  d = 50,
  refined = TRUE,
  reduced_dims = "wnn_latent"
)


# Construct design matrix
milo <- countCells(
  milo,
  meta.data = as.data.frame(colData(sce)),
  sample = "orig.ident"
)

design.df <-
  full.epi.no4NQO@meta.data %>%
  select(orig.ident, Condition) %>%
  distinct()

rownames(design.df) <- design.df$orig.ident

design.df <- design.df[colnames(nhoodCounts(milo)),]

design.df$Condition <-
  factor(
    design.df$Condition,
    levels=c("Vehicle","rOVA")
  )


# Build neighborhood graph for visualization
milo <- calcNhoodDistance(
  milo,
  d = 50,
  reduced.dim = "wnn_latent"
)

milo <- buildNhoodGraph(milo)


# Differential abundance testing
reducedDim(milo, "PCA") <-
  Embeddings(full.epi.no4NQO, "pca")

da.results <- testNhoods(
  milo,
  design = ~ Condition,
  design.df = design.df
)


# Annotate neighborhoods
da.results <- annotateNhoods(
  milo,
  da.results,
  coldata_col = "CellStates"
)

sig.da <-
  da.results |>
  filter(SpatialFDR < 0.05)


# Summary by CellStates
cellstate.summary <-
  da.results |>
  group_by(CellStates) |>
  summarise(
    n_neighborhoods = n(),
    significant = sum(SpatialFDR < 0.05),
    mean_logFC = mean(logFC),
    median_logFC = median(logFC),
    min_FDR = min(SpatialFDR)
  ) |>
  arrange(min_FDR)

cellstate.summary


# Visualize
ident_order <- c("EP14", "EP13", "EP12", "EP11", "EP10", "EP9", "EP8",
                 "EP7", "EP6", "EP5", "EP4", "EP3", "EP2", "EP1")
da.results$CellStates <- factor(da.results$CellStates, levels=ident_order)

plotDAbeeswarm(
  da.results,
  group.by = "CellStates"
)+
  scale_color_gradient2(
    low = "skyblue",
    mid = "white",
    high = "darkgreen",
    midpoint = 0
  )


# Figure 2F ---------------------------------------------------------------

# Load CDS
load('frontiersclusters')

# Extract raw counts & metadata
counts <- counts(cds)

meta <- as.data.frame(colData(cds))


# Create Seurat object
cds.seurat <- CreateSeuratObject(
  counts = counts,
  meta.data = meta
)

cds.seurat <- NormalizeData(cds.seurat)
cds.seurat <- FindVariableFeatures(cds.seurat)
cds.seurat <- ScaleData(cds.seurat)
cds.seurat <- RunPCA(cds.seurat, verbose = FALSE)


# Transfer anchors
trans.anchors <- FindTransferAnchors(
  reference = cds.seurat,
  query = full.epi.no4NQO,
  dims = 1:30,
  reference.reduction = 'pca'
)

TransferID <- TransferData(
  anchorset = trans.anchors,
  refdata = cds.seurat$cluster
)

full.epi.no4NQO <- AddMetaData(full.epi.no4NQO, metadata = TransferID)


# Rename clusters
Idents(full.epi.no4NQO) <- "predicted.id"
full.epi.no4NQO <- RenameIdents(full.epi.no4NQO,
                                '4'='B1','3'='B2','1'='B3','2'='B4','9'='B5',
                                '8'='SB1','10'='SB2','7'='SB3','5'='SB4','13'='SB5','11'='SB6',
                                '6'='S1','12'='S2')
full.epi.no4NQO@meta.data$predicted.id <- Idents(full.epi.no4NQO)


# Visualize
my_cols <- c('#956D61','#C72C20','#F3D412','#21D5FE','#1442CF','#FF6900',
             '#FFBB9B','#1CC72C','#404040','#0B6B05','#18FCBB','#E264F8')
DimPlot(full.epi.no4NQO, reduction = 'wnn.umap', group.by = 'predicted.id', cols = my_cols,
        label = T, repel = T, label.box = T)


# Combined basal, suprabasal, and superficial
Idents(full.epi.no4NQO) <- "predicted.id"
full.epi.no4NQO <- RenameIdents(full.epi.no4NQO,
                                'B1'='Basal','B2'='Basal','B3'='Basal','B4'='Basal','B5'='Basal',
                                'SB1'='Suprabasal','SB2'='Suprabasal','SB3'='Suprabasal','SB4'='Suprabasal','SB5'='Suprabasal','SB6'='Suprabasal',
                                'S1'='Superficial','S2'='Superficial')
full.epi.no4NQO@meta.data$predicted.id <- Idents(full.epi.no4NQO)


# Visualize Dimplot
my_cols <- c('#E38900','#00C1AD','#DC71FA')
DimPlot(full.epi.no4NQO, reduction = 'wnn.umap', group.by = 'predicted.id', cols = my_cols,
        label = T, repel = T, label.box = T, label.size = 6)


# Visualize bar graph
data <- full.epi.no4NQO %>%
  FetchData(vars = c("predicted.id", "Condition")) %>%
  group_by(Condition, predicted.id) %>%
  tally() %>%
  group_by(Condition) %>%
  mutate(percent = n / sum(n) * 100)

ggplot(data, aes(x = Condition, y = percent, fill = predicted.id)) + 
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) + 
  scale_fill_manual(values = my_cols) +
  labs(x = "Cluster", y = "Percent of Cells (%)") +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "bold"),
        axis.title = element_text(face = "bold")) + coord_flip()


# Figure 2G ---------------------------------------------------------------

### CellChat analysis

# Create CellChat objects
DefaultAssay(full.epi.no4NQO) <- 'RNA'
Idents(full.epi.no4NQO) <- 'Condition'
full.epi.no4NQO.Vehicle <- subset(full.epi.no4NQO, subset = Condition %in% 'Vehicle')
full.epi.no4NQO.rOVA <- subset(full.epi.no4NQO, subset = Condition %in% 'rOVA')

cellchat.full.epi.Vehicle <- createCellChat(object = full.epi.no4NQO.Vehicle, group.by = "CellStates")
cellchat.full.epi.rOVA <- createCellChat(object = full.epi.no4NQO.rOVA, group.by = "CellStates")


#  Extract metadeta
cellchat.full.epi.Vehicle@meta$condition <- "full.epi.no4NQO"
cellchat.full.epi.rOVA@meta$condition <- "full.epi.no4NQO"


# Establish database
CellChatDB <- CellChatDB.mouse
cellchat.full.epi.Vehicle@DB <- CellChatDB
cellchat.full.epi.rOVA@DB <- CellChatDB

# Subset the expression data of signaling genes
cellchat.full.epi.Vehicle <- subsetData(cellchat.full.epi.Vehicle)
cellchat.full.epi.rOVA <- subsetData(cellchat.full.epi.rOVA)

# Pre-process the expression data
cellchat.full.epi.Vehicle <- identifyOverExpressedGenes(cellchat.full.epi.Vehicle)
cellchat.full.epi.Vehicle <- identifyOverExpressedInteractions(cellchat.full.epi.Vehicle)
cellchat.full.epi.rOVA <- identifyOverExpressedGenes(cellchat.full.epi.rOVA)
cellchat.full.epi.rOVA <- identifyOverExpressedInteractions(cellchat.full.epi.rOVA)

# Project gene expression data onto protein-protein interaction (PPI)
cellchat.full.epi.Vehicle <- smoothData(cellchat.full.epi.Vehicle, adj = PPI.mouse)
cellchat.full.epi.rOVA <- smoothData(cellchat.full.epi.rOVA, adj = PPI.mouse)

# Compute the communication probability and infer cellular communication network
cellchat.full.epi.Vehicle <- computeCommunProb(cellchat.full.epi.Vehicle, type = "triMean")
cellchat.full.epi.Vehicle <- computeCommunProbPathway(cellchat.full.epi.Vehicle)
cellchat.full.epi.rOVA <- computeCommunProb(cellchat.full.epi.rOVA, type = "triMean")
cellchat.full.epi.rOVA <- computeCommunProbPathway(cellchat.full.epi.rOVA)

# Calculate the aggregated cell-cell communication network
cellchat.full.epi.Vehicle <- aggregateNet(cellchat.full.epi.Vehicle)
groupSize.Vehicle <- as.numeric(table(cellchat.full.epi.Vehicle@idents))
cellchat.full.epi.rOVA <- aggregateNet(cellchat.full.epi.rOVA)
groupSize.rOVA <- as.numeric(table(cellchat.full.epi.rOVA@idents))

# Compute the network centrality scores
cellchat.full.epi.Vehicle <- netAnalysis_computeCentrality(cellchat.full.epi.Vehicle, slot.name = "netP")
cellchat.full.epi.rOVA <- netAnalysis_computeCentrality(cellchat.full.epi.rOVA, slot.name = "netP")

# Merge objects
object.list <- list(Vehicle = cellchat.full.epi.Vehicle, rOVA = cellchat.full.epi.rOVA)
cellchat.no4NQO <- mergeCellChat(object.list, add.names = names(object.list))

# Compare the overall information flow of each signaling pathway
p <- rankNet(cellchat.no4NQO, mode = "comparison", stacked = T, do.stat = TRUE)

p + scale_fill_manual(values = c("darkgreen", "skyblue")) +
  theme(
    axis.text.y = element_text(color = c("darkgreen","black","black",
                                         "skyblue","skyblue","skyblue","skyblue","skyblue",
                                         "skyblue","skyblue","skyblue","skyblue","skyblue",
                                         "skyblue","skyblue","skyblue","skyblue","skyblue",
                                         "skyblue","skyblue","skyblue","skyblue","skyblue",
                                         "skyblue","skyblue","skyblue","skyblue","skyblue"))  
  )


# Figure 2H ---------------------------------------------------------------

# visualize the contribution of each LR pairs to the communication network
netAnalysis_contribution(cellchat.full.epi.Vehicle, 
                         signaling = "NOTCH", 
                         title = "Contribution of each LR pairs")
netAnalysis_contribution(cellchat.full.epi.rOVA, 
                         signaling = "NOTCH", 
                         title = "Contribution of each LR pairs")









