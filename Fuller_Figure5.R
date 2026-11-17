library(tidyverse)
library(reshape2)
library(Seurat)
library(Signac)
library(clustree)
library(ggrepel)
library(ggpubr)
library(org.Mm.eg.db)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(CellChat)
library(msigdbr)
library(fgsea)

# Figure 5A ---------------------------------------------------------------

# DimHeatMap
full.epi.4NQO.ds <- ProjectDim(full.epi.4NQO.ds, reduction = 'wnn.umap', assay = "RNA")
DimHeatmap(object = full.epi.4NQO.ds, reduction = 'wnn.umap', projected = T, 
           dims = 1, cells = 1000, balanced = T, fast = F, nfeatures = 100)

# To reorder
plot <- DimHeatmap(object = full.epi.4NQO.ds, reduction = 'wnn.umap', projected = T, 
                   dims = 1, cells = 1000, balanced = T, fast = F, nfeatures = 100)
current_genes <- ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range
flipped_genes <- rev(current_genes)
plot + scale_y_discrete(limits = flipped_genes)

# Generate gene programs
dim_Sox6 <- c(
  'Sox6','Antxr1','Uty','Ssbp2','Fosb','Ntm','Smoc2','Atf3','Arhgap24','Auts2',
  'Kcnq5','Vav3','Slc6a6','Hpgd','Foxp2','Car12','Ddx3y','Gpnmb','Rora','Il1r2',
  'Ccdc3','Trerf1','Eif2s3y','Cachd1','Rps6ka5','Chit1','Pappa','Krt15','Enox1','Nfia',
  'Zbtb7c','Dnajc6','Parp8','Mecom','Ttc6','Ust','Cdh13','P3h2','mt-Nd3','Plxna2',
  'Psd3','Hexb','Cnot6l','Lars2','Jun','Nr4a1','Tnfrsf21','Tgfbr3','Dgki','Lypd6b'
)
write.csv(dim_Sox6, file = "full.epi.4NQO.ds_dim_Sox6.csv")

dim_Krt17 <- c(
  'Anxa2','S100a8','Slc16a3','Rbpms','Cdc42bpa','Actb','Myo5b','S100a9','Rab31','Tmcc3',
  'Tgfbi','Grik3','Itga2','Sh3kbp1','Ero1l','Tnc','Mtss1','Col18a1','Krt6b','Ptgs2',
  'Bnip3','Prelid2','Klf7','Actn1','Tnfrsf23','Sfn','Nrg1','Aplp2','Flrt3','Phlda1',
  'Mmp13','Lamb3','Mxi1','Anxa1','Fabp5','Acsbg1','Vegfa','Gm49890','Krt16','S100a6',
  'Samd12','Krt6a','Ndrg1','Flnb','Cald1','Pcdh7','Egln3','Inhba','Lamc2','Krt17'
)
write.csv(dim_Krt17, file = "full.epi.4NQO.ds_dim_Krt17.csv")

# Add gene programs as modules
full.epi.4NQO.ds <- AddModuleScore(
  full.epi.4NQO.ds,
  features = list(dim_Sox6),
  name = "dim_Sox6"
)

full.epi.4NQO.ds <- AddModuleScore(
  full.epi.4NQO.ds,
  features = list(dim_Krt17),
  name = "dim_Krt17"
)

# Visualize
FeaturePlot(full.epi.4NQO.ds, reduction = "wnn.umap", features = c("dim_Sox61","dim_Krt171"))



# Figure 5B ---------------------------------------------------------------

VlnPlot(full.epi.4NQO.ds, features = c("dim_Sox61","dim_Krt171"),
        group.by = "predicted.id", split.by = 'Condition', pt.size = 0, 
        cols = c("orangered", "darkviolet"))

# Statistics
df <- data.frame(
  Dim = full.epi.4NQO.ds$dim_Sox61,
  Condition = full.epi.4NQO.ds$Condition,
  predicted.id = full.epi.4NQO.ds$predicted.id,
  Sample = full.epi.4NQO.ds$orig.ident
) %>%
  group_by(Sample, Condition, predicted.id) %>%
  summarise(
    MeanDim = mean(Dim),
    .groups = "drop"
  )

results <- df %>%
  filter(predicted.id %in% c("Basal", "Suprabasal", "Superficial")) %>%
  group_by(predicted.id) %>%
  summarise(
    p.value = wilcox.test(
      MeanDim ~ Condition
    )$p.value
  )

results$FDR <- p.adjust(results$p.value, method = "BH")

### Repeat for dim_Krt171


# Figure 5C ---------------------------------------------------------------

# Set genome
TxDb.Mm <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Find differentially accessible regions
DefaultAssay(full.epi.4NQO.ds) <- 'ATAC'
Idents(full.epi.4NQO.ds) <- 'Condition'
topmarkers <- FindMarkers(full.epi.4NQO.ds,
                          ident.1 = 'rOVA+4NQO',
                          ident.2 = '4NQO',
                          only.pos = T, 
                          assay = "ATAC")
signif <- topmarkers %>% filter(p_val_adj < 0.05)

# Reformat
signif_spl <- strsplit(signif$gene, split = "-" )
gene_spl <- matrix(unlist(signif_spl),ncol=3,byrow=T)
colnames(gene_spl) <- c('chr', 'start', 'end')
gene_spl <- data.frame(gene_spl)
signif_rm <- signif[ -c(7) ]
signif_spl <- cbind(signif_rm,gene_spl)
signif_granges <- GRanges(signif_spl)

# Annotate peaks
peakAnno <- annotatePeak(signif_granges, tssRegion=c(-3000, 3000), TxDb=TxDb.Mm)
write.csv(peakAnno@anno, file = 'peakanno4NQO.csv')

# Convert Entrez geneid to gene name
peakAnno <- read.csv('peakanno4NQO.csv')
myGenes <- peakAnno$geneId
myGenes <- as.character(myGenes)
gene <- select(org.Mm.eg.db, keys=myGenes, columns=c("SYMBOL"),keytype="ENTREZID")
peakAnno$gene <- gene
write.csv(peakAnno, file = 'peakanno4NQO_final.csv')

# Reformat for pie chart
data <- read.csv("peakanno4NQO_final.csv", stringsAsFactors = FALSE)
data <- data %>%
  mutate(annotation = case_when(
    str_detect(annotation, regex("intron 1 of", ignore_case = TRUE)) ~ "1st Intron",
    str_detect(annotation, regex("exon 1 of", ignore_case = TRUE))   ~ "1st Exon",
    str_detect(annotation, regex("intron", ignore_case = TRUE)) ~ "Other Intron",
    str_detect(annotation, regex("exon", ignore_case = TRUE))   ~ "Other Exon",
    TRUE ~ annotation
  ))
unique(data$annotation)

data$annotation <- factor(
  data$annotation,
  levels = c(
    "Promoter (<=1kb)",
    "Promoter (1-2kb)",
    "Promoter (2-3kb)",
    "5' UTR",
    "3' UTR",
    "1st Exon",
    "Other Exon",
    "1st Intron",
    "Other Intron",
    "Downstream (<=300bp)",
    "Distal Intergenic"
  )
)

data$cluster <- factor(
  data$cluster,
  levels = c("4NQO", "rOVA+4NQO")
)

# Summarize data
summary_df <- data %>%
  group_by(cluster, annotation) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(percent = count / sum(count))

# Visualize
ggplot(summary_df, aes(x = "", y = percent, fill = annotation)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  facet_wrap(~ cluster) +
  theme_void()


# Figure 5D ---------------------------------------------------------------

### CellChat analysis - Full dataset

# Create CellChat objects
DefaultAssay(full) <- 'RNA'
Idents(full) <- 'Condition'
full.4NQO <- subset(full, subset = Condition %in% '4NQO')
full.rOVA4NQO <- subset(full, subset = Condition %in% 'rOVA4NQO')

cellchat.full.4NQO <- createCellChat(object = full.4NQO, group.by = "BroadCellTypes")
cellchat.full.rOVA4NQO <- createCellChat(object = full.rOVA4NQO, group.by = "BroadCellTypes")

#  Extract metadeta
cellchat.full.4NQO@meta$condition <- "full"
cellchat.full.rOVA4NQO@meta$condition <- "full"

# Establish database
CellChatDB <- CellChatDB.mouse
cellchat.full.4NQO@DB <- CellChatDB
cellchat.full.rOVA4NQO@DB <- CellChatDB

# Subset the expression data of signaling genes
cellchat.full.4NQO <- subsetData(cellchat.full.4NQO)
cellchat.full.rOVA4NQO <- subsetData(cellchat.full.rOVA4NQO)

# Pre-process the expression data
cellchat.full.4NQO <- identifyOverExpressedGenes(cellchat.full.4NQO)
cellchat.full.4NQO <- identifyOverExpressedInteractions(cellchat.full.4NQO)
cellchat.full.rOVA4NQO <- identifyOverExpressedGenes(cellchat.full.rOVA4NQO)
cellchat.full.rOVA4NQO <- identifyOverExpressedInteractions(cellchat.full.rOVA4NQO)

# Project gene expression data onto protein-protein interaction (PPI)
cellchat.full.4NQO <- smoothData(cellchat.full.4NQO, adj = PPI.mouse)
cellchat.full.rOVA4NQO <- smoothData(cellchat.full.rOVA4NQO, adj = PPI.mouse)

# Compute the communication probability and infer cellular communication network
cellchat.full.4NQO <- computeCommunProb(cellchat.full.4NQO, type = "triMean")
cellchat.full.4NQO <- computeCommunProbPathway(cellchat.full.4NQO)
cellchat.full.rOVA4NQO <- computeCommunProb(cellchat.full.rOVA4NQO, type = "triMean")
cellchat.full.rOVA4NQO <- computeCommunProbPathway(cellchat.full.rOVA4NQO)

# Calculate the aggregated cell-cell communication network
cellchat.full.4NQO <- aggregateNet(cellchat.full.4NQO)
groupSize.4NQO <- as.numeric(table(cellchat.full.4NQO@idents))
cellchat.full.rOVA4NQO <- aggregateNet(cellchat.full.rOVA4NQO)
groupSize.rOVA4NQO <- as.numeric(table(cellchat.full.rOVA4NQO@idents))

# Compute the network centrality scores
cellchat.full.4NQO <- netAnalysis_computeCentrality(cellchat.full.4NQO, slot.name = "netP")
cellchat.full.rOVA4NQO <- netAnalysis_computeCentrality(cellchat.full.rOVA4NQO, slot.name = "netP")

# Merge objects
object.list <- list(NQO = cellchat.full.4NQO, rOVA4NQO = cellchat.full.rOVA4NQO)

# Compare outgoing/incoming interaction strength for all cell types
count.sum <- sapply(object.list, function(x) {rowSums(x@net$count) + 
    colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(count.sum), max(count.sum))
gg <- list()
for (i in 1:length(object.list)) {
  p <- netAnalysis_signalingRole_scatter(object.list[[i]], 
                                         title = names(object.list)[i], 
                                         weight.MinMax = weight.MinMax) 
  p <- p + scale_color_manual(values = c('lightgreen','orange','plum','deepskyblue')) +
    scale_x_continuous(limits = c(0,1.5)) +
    scale_y_continuous(limits = c(0,1.5))
  gg[[i]] <- p
}

patchwork::wrap_plots(plots = gg)


# Figure 5E ---------------------------------------------------------------

### CellChat analysis - 4NQO(+) epithelial

# Create CellChat objects
DefaultAssay(full.epi.4NQO.ds) <- 'RNA'
Idents(full.epi.4NQO.ds) <- 'Condition'
full.epi.4NQO.ds.4NQO <- subset(full.epi.4NQO.ds, subset = Condition %in% '4NQO')
full.epi.4NQO.ds.rOVA4NQO <- subset(full.epi.4NQO.ds, subset = Condition %in% 'rOVA+4NQO')

cellchat.full.epi.4NQO <- createCellChat(object = full.epi.4NQO.ds.4NQO, group.by = "CellStates")
cellchat.full.epi.rOVA4NQO <- createCellChat(object = full.epi.4NQO.ds.rOVA4NQO, group.by = "CellStates")

#  Extract metadeta
cellchat.full.epi.4NQO@meta$condition <- "full.epi.4NQO.ds"
cellchat.full.epi.rOVA4NQO@meta$condition <- "full.epi.4NQO.ds"

# Establish database
CellChatDB <- CellChatDB.mouse
cellchat.full.epi.4NQO@DB <- CellChatDB
cellchat.full.epi.rOVA4NQO@DB <- CellChatDB

# Subset the expression data of signaling genes
cellchat.full.epi.4NQO <- subsetData(cellchat.full.epi.4NQO)
cellchat.full.epi.rOVA4NQO <- subsetData(cellchat.full.epi.rOVA4NQO)

# Pre-process the expression data
cellchat.full.epi.4NQO <- identifyOverExpressedGenes(cellchat.full.epi.4NQO)
cellchat.full.epi.4NQO <- identifyOverExpressedInteractions(cellchat.full.epi.4NQO)
cellchat.full.epi.rOVA4NQO <- identifyOverExpressedGenes(cellchat.full.epi.rOVA4NQO)
cellchat.full.epi.rOVA4NQO <- identifyOverExpressedInteractions(cellchat.full.epi.rOVA4NQO)

# Project gene expression data onto protein-protein interaction (PPI)
cellchat.full.epi.4NQO <- smoothData(cellchat.full.epi.4NQO, adj = PPI.mouse)
cellchat.full.epi.rOVA4NQO <- smoothData(cellchat.full.epi.rOVA4NQO, adj = PPI.mouse)

# Compute the communication probability and infer cellular communication network
cellchat.full.epi.4NQO <- computeCommunProb(cellchat.full.epi.4NQO, type = "triMean")
cellchat.full.epi.4NQO <- computeCommunProbPathway(cellchat.full.epi.4NQO)
cellchat.full.epi.rOVA4NQO <- computeCommunProb(cellchat.full.epi.rOVA4NQO, type = "triMean")
cellchat.full.epi.rOVA4NQO <- computeCommunProbPathway(cellchat.full.epi.rOVA4NQO)

# Calculate the aggregated cell-cell communication network
cellchat.full.epi.4NQO <- aggregateNet(cellchat.full.epi.4NQO)
cellchat.full.epi.rOVA4NQO <- aggregateNet(cellchat.full.epi.rOVA4NQO)
groupSize.4NQO <- as.numeric(table(cellchat.full.epi.4NQO@idents))
groupSize.rOVA4NQO <- as.numeric(table(cellchat.full.epi.rOVA4NQO@idents))

# Compute the network centrality scores
cellchat.full.epi.4NQO <- netAnalysis_computeCentrality(cellchat.full.epi.4NQO, slot.name = "netP")
cellchat.full.epi.rOVA4NQO <- netAnalysis_computeCentrality(cellchat.full.epi.rOVA4NQO, slot.name = "netP")

# Merge objects
object.list <- list(NQO = cellchat.full.epi.4NQO, rOVA4NQO = cellchat.full.epi.rOVA4NQO)
cellchat.4NQO <- mergeCellChat(object.list, add.names = names(object.list))

# Compare the overall information flow of each signaling pathway
p <- rankNet(cellchat.no4NQO, mode = "comparison", stacked = T, do.stat = TRUE)

p + scale_fill_manual(values = c("darkviolet", "orangered")) +
  theme(
    axis.text.y = element_text(color = c("darkviolet","darkviolet","darkviolet","darkviolet",
                                         "darkviolet","darkviolet","darkviolet","darkviolet",
                                         "darkviolet","darkviolet","darkviolet","black",
                                         "darkviolet","darkviolet","darkviolet","darkviolet",
                                         "darkviolet","darkviolet","darkviolet","darkviolet",
                                         "darkviolet","black","orangered"))  
  )



# Figure 5F ---------------------------------------------------------------

netAnalysis_signalingRole_heatmap(cellchat.full.epi.4NQO, pattern = "incoming", height = 10)
netAnalysis_signalingRole_heatmap(cellchat.full.epi.rOVA4NQO, pattern = "incoming", height = 10)


# Figure 5H ---------------------------------------------------------------

### GSEA

# Merge EP2 and EP3
Idents(full.epi.4NQO.ds) <- 'CellStates'
full.epi.4NQO.ds$mergedCellStates <- as.character(Idents(full.epi.4NQO.ds))
full.epi.4NQO.ds$mergedCellStates[full.epi.4NQO.ds$mergedCellStates %in% c("Epi 2", "Epi 3")] <- "Epi 2_Epi 3"
full.epi.4NQO.ds$mergedCellStates <- factor(full.epi.4NQO.ds$mergedCellStates)
Idents(full.epi.4NQO.ds) <- "mergedCellStates"

# Find DEGs
de <- FindMarkers(
  object = full.epi.4NQO.ds,
  ident.1 = "Epi 2_Epi 3",
  assay = "RNA")
write.csv(de, file = 'full.epi.4NQO.ds_EP2EP3.csv')

# Load
deseq_results <- as.data.frame(read.csv('full.epi.4NQO.ds_EP2EP3.csv',
                                        header = T))

# Rank genes
gene_list <- deseq_results %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(X, avg_log2FC) %>%
  deframe()

# Assign gene set
all_genesets <- msigdbr(species = "Mus musculus")

specific_pathway <- all_genesets %>%
  filter(str_detect(gs_name, "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_MIGRATION"))

# Generate pathways
pathways_list <- specific_pathway %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(.$gs_name)

gobp_cs <- pathways_list
gobp_cs$GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_MIGRATION <- 
  gobp_cs$GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_MIGRATION$gene_symbol

# Run GSEA
fgsea_results <- fgsea(pathways = gobp_cs, 
                       stats = gene_list, 
                       minSize = 5, maxSize = 500)

# Visualize

plotEnrichment(
  pathway = gobp_cs[["GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_MIGRATION"]],
  stats   = gene_list
)


# Figure 5I ---------------------------------------------------------------

# Find & transfer anchors
anchors <- FindTransferAnchors(
  reference = full.epi.4NQO.ds,
  query = full.epi.no4NQO,
  dims = 1:30,
  reference.reduction = "pca"
)

pred <- TransferData(
  anchorset = anchors,
  refdata = full.epi.4NQO.ds$CellStates,
  dims = 1:30
)

# Add metadata
full.epi.no4NQO <- AddMetaData(
  object = full.epi.no4NQO,
  metadata = pred$predicted.id,
  col.name = "CancerStates"
)

# Visualize
overlap_table <- table(full.epi.no4NQO$CellStates, full.epi.no4NQO$CancerStates)
overlap_prop <- prop.table(overlap_table, margin = 1)
tab_EP1EP2EP3 <- round(100 * overlap_prop[, c("EP1", "EP2", "EP3")],1)
pheatmap(
  tab_EP1EP2EP3, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_color = "black",
  fontsize = 18,
  fontsize_number = 18,
  fontsize_row = 18,
  fontsize_col = 18,
  color = colorRampPalette(c("white", "orange2"))(100)
)