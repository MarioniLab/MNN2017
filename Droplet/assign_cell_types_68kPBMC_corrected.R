## assign cell type labels to the 68K PBMC data

# combine human preimplantation embryo data sets
# might be tricky as one data set has very few HVGs
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

pbmc <- read.table("Droplet/Results/PBMC_68k_corrected.tsv",
                   h=T, sep="\t", stringsAsFactors=FALSE)
rownames(pbmc) <- pbmc$gene_id

tcell.hvg <- read.table("Droplet/Data/Tcell/Tcell_hvg.tsv",
                        h=T, sep="\t", stringsAsFactors=FALSE)

pbmc.hvg <- read.table("Droplet/Data/PBMC/PBMC_hvg.tsv",
                       h=T, sep="\t", stringsAsFactors=FALSE)

select.hvg <- intersect(tcell.hvg$HVG, pbmc.hvg$HVG)

pbmc.tsne <- read.table("Droplet/Results/corrected_tSNE.tsv",
                        h=TRUE, sep="\t", stringsAsFactors=FALSE)

pbmc.community <- read.table("Droplet/Results/PBMC_Corrected_communities.tsv",
                             h=TRUE, sep="\t", stringsAsFactors=FALSE)

# map ensembl IDs to gene symbols to select a couple of marker genes for visualisation
all.genes <- rownames(pbmc)
ensembl <- useEnsembl(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', GRCh=37)
gene_symbol <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                     filters='ensembl_gene_id', mart=ensembl,
                     values=all.genes)

# highlight clusters with marker genes after SNN and randomwalk-based clustering
m.genes <- c("CD14", "PTPRC", "FCGR3A", "FCGR1A", "ITGAX", "ITGAM", "CD19", "HLA-DRB1", "FCGR2B", "FCGR2A",
             "CD3E", "CD4", "CD8A","CD8B", "CD28", "CD8", "TBX21", "IKAROS", "IL2RA", "CD44", "SELL",
             "CCR7", "MS4A1", "CD68", "CD163", "IL5RA", "SIGLEC8", "KLRD1", "NCR1", "CD22", "IL3RA",
             "CCR6", "IL7R", "CD27", "FOXP3", "PTCRA", "ID3", "PF4", "CCR10", "SIGLEC7", "NKG7",
             "S100A8", "CXCR3", "CCR5", "CCR3", "CCR4", "PTGDR2", "RORC", "CD1C")

m.ensembl <- gene_symbol$ensembl_gene_id[gene_symbol$external_gene_name %in% m.genes]

m.exprs <- t(pbmc[m.ensembl, 2:(dim(pbmc)[2])])
samples <- rownames(m.exprs)
m.exprs <- as.data.frame(apply(m.exprs, 2, as.numeric))
m.exprs$Sample <- samples

m.symbols <- gene_symbol$external_gene_name[gene_symbol$ensembl_gene_id %in% colnames(m.exprs)]
colnames(m.exprs) <- c(m.symbols, "Sample")

pbmc.merge <- merge(pbmc.tsne, m.exprs, by='Sample')

comm.tsne <- ggplot(pbmc.merge,
        aes(x=Dim1, y=Dim2, fill=as.factor(Community))) +
   geom_point(size=2, pch=21) +
   theme_mike() + scale_fill_Publication()
comm.tsne

ggplot(pbmc.merge,
       aes(x=Dim1, y=Dim2, colour=`MS4A1`)) +
  geom_point(size=1, alpha=0.6) + theme_mike() +
  scale_colour_distiller(palette="Reds", direction=1)

# using these sets of results I'll assign identities to the different communities
# CD3+ T cells: 1, 6, 8, 11
# CD20+ (MS4A1) B cells: 13
# CD1c+ dendritic cells: 5
# CD14+ monocytes: 3
# CD16+ macrophages: 2
# CD3- CD16+ NK cells: 4
# CD3+ CD16+ NKT cells: 10
# Unknown clusters: 7, 9, 12, 14, 15

# # write out these cell identities
# pbmc.tsne$CellType <- "Unknown"
# pbmc.tsne$CellType[pbmc.tsne$Community %in% c(1, 6, 8, 11)] <- "T_cell"
# pbmc.tsne$CellType[pbmc.tsne$Community %in% c(13)] <- "B_cell"
# pbmc.tsne$CellType[pbmc.tsne$Community %in% c(3)] <- "monocyte"
# pbmc.tsne$CellType[pbmc.tsne$Community %in% c(2)] <- "macrophage"
# pbmc.tsne$CellType[pbmc.tsne$Community %in% c(5)] <- "DCs"
# pbmc.tsne$CellType[pbmc.tsne$Community %in% c(4)] <- "NKcell"
# pbmc.tsne$CellType[pbmc.tsne$Community %in% c(10)] <- "NKT"
# 
# write.table(pbmc.tsne,
#             "~/CI_filesystem/mnt/scratcha/jmlab/morgan02/10X/MNN/data/PBMC/PBMC_meta.tsv",
#             sep="\t", quote=FALSE, row.names=FALSE)

# pbmc.merge$CellType <- "Unknown"
# pbmc.merge$CellType[pbmc.merge$Community %in% c(1, 6, 8, 11)] <- "T_cell"
# pbmc.merge$CellType[pbmc.merge$Community %in% c(13)] <- "B_cell"
# pbmc.merge$CellType[pbmc.merge$Community %in% c(3)] <- "monocyte"
# pbmc.merge$CellType[pbmc.merge$Community %in% c(2)] <- "macrophage"
# pbmc.merge$CellType[pbmc.merge$Community %in% c(5)] <- "DCs"
# pbmc.merge$CellType[pbmc.merge$Community %in% c(4)] <- "NKcell"
# pbmc.merge$CellType[pbmc.merge$Community %in% c(10)] <- "NKT"
# 
# comm.tsne <- ggplot(pbmc.merge,
#                     aes(x=Dim1, y=Dim2, fill=CellType)) +
#   geom_point(size=2, pch=21) +
#   theme_mike() + scale_fill_Publication()
# comm.tsne
# 
# ggsave(comm.tsne,
#        filename="~/Dropbox/MNNfigures/10X_data/PBMC/PBMC_68k_celltype-tSNE.png",
#        width=6.75, height=4.75, dpi=300)
