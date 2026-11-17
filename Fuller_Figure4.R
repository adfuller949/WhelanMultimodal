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

# Figure 4A ---------------------------------------------------------------

DefaultAssay(full) <- 'RNA'
Idents(full) <- 'BroadCellTypes'
full.epi <- subset(full, subset = CellTypes %in% c('Epithelial'))
Idents(full) <- 'Condition'
full.epi.4NQO <- subset(full.epi, 
                          subset = Condition %in% c('4NQO','rOVA+4NQO'))

full.epi.4NQO.ds <- subset(x=full.epi.4NQO, downsample = 4700)

# RNA analysis
DefaultAssay(full.epi.4NQO.ds) <- "RNA"
full.epi.4NQO.ds <- SCTransform(full.epi.4NQO.ds, vars.to.regress = 'percent.mt', verbose = FALSE)
full.epi.4NQO.ds <- RunPCA(full.epi.4NQO.ds, verbose = FALSE)
ElbowPlot(full.epi.4NQO.ds, ndims = 50, reduction = "pca")
full.epi.4NQO.ds <- RunUMAP(full.epi.4NQO.ds, dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
DefaultAssay(full.epi.4NQO.ds) <- "ATAC"
full.epi.4NQO.ds <- RunTFIDF(full.epi.4NQO.ds)
full.epi.4NQO.ds <- FindTopFeatures(full.epi.4NQO.ds, min.cutoff = 'q0')
full.epi.4NQO.ds <- RunSVD(full.epi.4NQO.ds)
full.epi.4NQO.ds <- RunUMAP(full.epi.4NQO.ds, reduction = 'lsi', dims = 2:30, 
                           reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# Multimodal analysis
full.epi.4NQO.ds <- FindMultiModalNeighbors(full.epi.4NQO.ds, reduction.list = list("pca", "lsi"),
                                           dims.list = list(1:30, 2:30))
full.epi.4NQO.ds <- RunUMAP(full.epi.4NQO.ds, nn.name = "weighted.nn", reduction.name = "wnn.umap",
                           reduction.key = "wnnUMAP_")
full.epi.4NQO.ds <- FindClusters(full.epi.4NQO.ds, graph.name = "wsnn", algorithm = 3, verbose = FALSE,
                                resolution = seq(0.1,1,0.1))
clustree(full.epi.4NQO.ds, prefix = "wsnn_res.")

table(Idents(full.epi.4NQO.ds), full.epi.4NQO.ds$orig.ident)
# Remove clusters with <15 cells and recluster


my_cols <- c('gold','cornflowerblue','mediumpurple1','burlywood','aquamarine',
             'mediumorchid1','chartreuse3','dodgerblue','hotpink1','salmon')
DimPlot(full.epi.4NQO.ds, reduction = 'wnn.umap', group.by = 'wsnn_res.0.4', 
        cols = my_cols, label = TRUE, repel = TRUE, label.box = TRUE)

# Normalize RNA assay
DefaultAssay(full.epi.4NQO.ds) <- 'RNA'
full.epi.4NQO.ds <- NormalizeData(full.epi.4NQO.ds)
full.epi.4NQO.ds <- FindVariableFeatures(full.epi.4NQO.ds, selection.method = "vst", nfeatures = 2000)
full.epi.4NQO.ds <- ScaleData(full.epi.4NQO.ds, vars.to.regress = 'percent.mt')

# Visualize dotplot of differentiation axis
axis_genes <- c('Krtdap','Krt4','Krt13','Fabp5',  
                'Col17a1','Krt5','Krt14', 'Krt15',  
                'Mki67','Top2a')   
Idents(full.epi.4NQO.ds) <- 'wsnn_res.0.4'
levels(full.epi.4NQO.ds) <- c('7','6','9','8','0','4','2','5','1','3')
DotPlot(object = full.epi.4NQO.ds, features = axis_genes, cols = c('navy','red'))+
  coord_flip() + RotatedAxis()

# Rename clusters
Idents(full.epi.4NQO.ds) <- "wsnn_res.0.4"
full.epi.4NQO.ds@meta.data$CellStates <- Idents(full.epi.4NQO.ds)
Idents(full.epi.4NQO.ds) <- "CellStates"
full.epi.4NQO.ds <- RenameIdents(full.epi.4NQO.ds,
                                 '7'='EP1', '6'='EP2', '9'='EP3', '8'='EP4',
                                 '0'='EP5', '4'='EP6', '2'='EP7', '5'='EP8',
                                 '1'='EP9','3'='EP10')
full.epi.4NQO.ds@meta.data$CellStates <- Idents(full.epi.4NQO.ds)

# Visualize
my_cols <- c('gold','cornflowerblue','mediumpurple1','burlywood','aquamarine',
             'mediumorchid1','chartreuse3','dodgerblue','hotpink1','salmon')
DimPlot(full.epi.4NQO.ds, reduction = 'wnn.umap', group.by = 'CellStates', 
        cols = my_cols, label = TRUE, repel = TRUE, label.box = TRUE, label.size = 6)


# Figure 4B ---------------------------------------------------------------

### RNA

DefaultAssay(full.epi.4NQO.ds) <- 'RNA'
Idents(full.epi.4NQO.ds) <- 'CellStates'
topmarkers <- FindAllMarkers(full.epi.4NQO.ds, assay = "RNA")
top100 <- topmarkers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC)

conditions_avex <- AggregateExpression(full.epi.4NQO.ds, return.seurat = T, assays = 'RNA')
conditions_avex <- ScaleData(conditions_avex, assay = 'RNA', features = top100$gene)
DoHeatmap(object = conditions_avex, features = top100$gene, assay = 'RNA', 
          group.colors = c('gold','cornflowerblue','mediumpurple1','burlywood','aquamarine',
                           'mediumorchid1','chartreuse3','dodgerblue','hotpink1','salmon'),
          label = T, size = 3, angle = 30, draw.lines = F)+
  guides(colour=FALSE)+ 
  theme(axis.text.y = element_blank())


### ATAC

DefaultAssay(full.epi.4NQO.ds) <- 'ATAC'
Idents(full.epi.4NQO.ds) <- 'CellStates'
topmarkers <- FindAllMarkers(full.epi.4NQO.ds, assay = "ATAC")
top100 <- topmarkers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC)

conditions_avex <- AggregateExpression(full.epi.4NQO.ds, return.seurat = T, assays = 'ATAC')
conditions_avex <- ScaleData(conditions_avex, assay = 'ATAC', features = top100$gene)
DoHeatmap(object = conditions_avex, features = top100$gene, assay = 'ATAC', 
          group.colors = c('gold','cornflowerblue','mediumpurple1','burlywood','aquamarine',
                           'mediumorchid1','chartreuse3','dodgerblue','hotpink1','salmon'),
          label = T, size = 3, angle = 30, draw.lines = F)+
  guides(colour=FALSE)+ 
  theme(axis.text.y = element_blank())


# Figure 4C ---------------------------------------------------------------

FeaturePlot(full.epi.4NQO.ds, reduction = "wnn.umap", features = c("Top2a","Krt14","Fabp5"))


# Figure 4D ---------------------------------------------------------------

axis_genes <- c('Krtdap','Krt4','Krt13','Fabp5',  
                'Col17a1','Krt5','Krt14', 'Krt15',  
                'Mki67','Top2a')   
Idents(full.epi.4NQO.ds) <- 'CellStates'
levels(full.epi.4NQO.ds) <- c('EP1','EP2','EP3','EP4','EP5','EP6','EP7',
                             'EP8','EP9','EP10')
DotPlot(object = full.epi.4NQO.ds, features = axis_genes, cols = c('navy','red'))+
  coord_flip() + RotatedAxis()


# Figure 4E ---------------------------------------------------------------

### Milo analysis

# Convert Seurat -> SingleCellExperiment
sce <- as.SingleCellExperiment(full.epi.4NQO.ds)


# Extract metadata & embeddings
colData(sce)$CellStates <- full.epi.4NQO.ds$CellStates

colData(sce)$Condition <-
  factor(full.epi.4NQO.ds$Condition,
         levels = c("4NQO","rOVA+4NQO"))

colData(sce)$orig.ident <- full.epi.4NQO.ds$orig.ident

rna.pca <- Embeddings(
  full.epi.4NQO.ds,
  "pca"
)[,1:50]

atac.lsi <- Embeddings(
  full.epi.4NQO.ds,
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
  full.epi.4NQO.ds@meta.data %>%
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
  Embeddings(full.epi.4NQO.ds, "pca")

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
ident_order <- c("EP10", "EP9", "EP8", "EP7", "EP6", "EP5", "EP4", "EP3", 
                 "EP2", "EP1")
da.results$CellStates <- factor(da.results$CellStates, levels=ident_order)

plotDAbeeswarm(
  da.results,
  group.by = "CellStates"
)+
  scale_color_gradient2(
    low = "orangered",
    mid = "white",
    high = "darkviolet",
    midpoint = 0
  )


# Figure 4F ---------------------------------------------------------------

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
  query = full.epi.4NQO.ds,
  dims = 1:30,
  reference.reduction = 'pca'
)

TransferID <- TransferData(
  anchorset = trans.anchors,
  refdata = cds.seurat$cluster
)

full.epi.4NQO.ds <- AddMetaData(full.epi.4NQO.ds, metadata = TransferID)


# Rename clusters
Idents(full.epi.4NQO.ds) <- "predicted.id"
full.epi.4NQO.ds <- RenameIdents(full.epi.4NQO.ds,
                                '4'='B1','3'='B2','1'='B3','2'='B4','9'='B5',
                                '8'='SB1','10'='SB2','7'='SB3','5'='SB4','13'='SB5','11'='SB6',
                                '6'='S1','12'='S2')
full.epi.4NQO.ds@meta.data$predicted.id <- Idents(full.epi.4NQO.ds)


# Visualize
my_cols <- c('#956D61','#C72C20','#F3D412','#21D5FE','#1442CF','#FF6900',
             '#FFBB9B','#1CC72C','#404040','#0B6B05','#18FCBB','#E264F8')
DimPlot(full.epi.4NQO.ds, reduction = 'wnn.umap', group.by = 'predicted.id', cols = my_cols,
        label = T, repel = T, label.box = T)


# Combined basal, suprabasal, and superficial
Idents(full.epi.4NQO.ds) <- "predicted.id"
full.epi.4NQO.ds <- RenameIdents(full.epi.4NQO.ds,
                                'B1'='Basal','B2'='Basal','B3'='Basal','B4'='Basal','B5'='Basal',
                                'SB1'='Suprabasal','SB2'='Suprabasal','SB3'='Suprabasal','SB4'='Suprabasal','SB5'='Suprabasal','SB6'='Suprabasal',
                                'S1'='Superficial','S2'='Superficial')
full.epi.4NQO.ds@meta.data$predicted.id <- Idents(full.epi.4NQO.ds)


# Visualize Dimplot
my_cols <- c('#E38900','#00C1AD','#DC71FA')
DimPlot(full.epi.4NQO.ds, reduction = 'wnn.umap', group.by = 'predicted.id', cols = my_cols,
        label = T, repel = T, label.box = T, label.size = 6)


# Visualize bar graph
data <- full.epi.4NQO.ds %>%
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


