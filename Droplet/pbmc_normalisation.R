library(ggplot2)
library(ggrepel)
library(Matrix)
library(Rtsne)
library(scater)
library(scran)
library(cellrangerRkit)
library(Rtsne)
library(stringr)
library(igraph)
source("~/Dropbox/R_sessions/SingleCellFunctions/single_cell_functions.R")

pbmc_path <- "Droplet/10X_data/PBMC/"

# this has to be an absolute path, relative paths break
download_sample(sample_name="pbmc4k", sample_dir=pbmc_path,
                host="https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/")
pbmc.10x <- load_cellranger_matrix(pbmc_path)

# filter non-zero genes and normalize PBMC 10X UMI counts
pbmc.nz <- get_nonzero_genes(pbmc.10x)

pbmc.norm <- size_factor_normalize(exprs(pbmc.10x), cell.sparse=0.95,
                                   gene.sparse=0.99, cluster.size=50)

# map ensembl IDs to gene symbols to select a couple of marker genes for visualisation
all.genes <- rownames(pbmc.norm)

ensembl <- useEnsembl(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', GRCh=37)

gene_symbol <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                     filters='ensembl_gene_id', mart=ensembl,
                     values=all.genes)

# select highly variable genes (~1-2k will be more than sufficient I think)
pbmc.hvg <- find_hvg(pbmc.norm[, 1:(dim(pbmc.norm)[2]-1)], plot=FALSE)

# the mean-CV^2 plot looks OK, even with the high number of unit counts
# there are 1192 HVG at p<=0.01

write.table(pbmc.norm,
            file="~/Dropbox/MNNfigures/10X_data/PBMC/PBMC_norm.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

hvg.df <- cbind.data.frame(names(pbmc.hvg)[pbmc.hvg])
colnames(hvg.df) <- c("HVG")
write.table(hvg.df,
            file="~/Dropbox/MNNfigures/10X_data/PBMC/PBMC_hvg.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

# tSNE the shiat out of this
pbmc.select <- pbmc.norm[pbmc.hvg, ]
pbmc.select[is.na(pbmc.select)] <- 0

set.seed(42)
pbmc.map <- tsne_wrapper(pbmc.select[, 1:(dim(pbmc.norm)[2]-1)], perplexity=100)

# pull out marker genes for PBMC populations to add identities to different clusters
# also try SNNgraph for clustering!  Use just 20 dimensions from the PCA to build the graph
pbmc.snn <- buildSNNGraph(as.matrix(pbmc.norm[, 1:(dim(pbmc.norm)[2]-1)]),
                          k=5, d=20)

# use the walktrap algorithm to form communities
# try a few different values for the step parameter
pbmc.community <- cluster_walktrap(pbmc.snn, steps=4)
n.comms <- length(pbmc.community)
print(n.comms)

# vertex community membership 
pbmc.members <- pbmc.community$membership
pbmc.clusters <- cbind.data.frame(colnames(pbmc.norm[, 1:(dim(pbmc.norm)[2]-1)]), pbmc.members)
colnames(pbmc.clusters) <- c("Sample", "Community")
pbmc.clusters$Community <- as.factor(pbmc.clusters$Community)

pbmc.max <- merge(pbmc.map, pbmc.clusters, by='Sample')

comm.tsne <- ggplot(pbmc.max,
                    aes(x=Dim1, y=Dim2, colour=Community)) + 
  geom_point(size=3) + theme_classic()

ggsave(comm.tsne,
       filename="~/Dropbox/MNNfigures/10X_data/PBMC/PBMC_communities_tSNE.png",
       height=5.75, width=5.75, dpi=300)

comm.tsne

# highlight clusters with marker genes after SNN and randomwalk-based clustering
# select a bunch of marker genes for the expected cell types
m.genes <- c("CD14", "PTPRC", "FCGR3A", "ITGAX", "ITGAM", "CD19", "HLA-DRB1", "FCGR2B", "FCGR2A",
             "CD3E", "CD4", "CD8A","CD8B", "CD28", "CD8", "TBX21", "IKAROS", "IL2RA", "CD44", "SELL",
             "CCR7", "MS4A1", "CD68", "CD163", "IL5RA", "SIGLEC8", "KLRD1", "NCR1", "CD22", "IL3RA",
             "CCR6", "IL7R", "CD27", "FOXP3", "PTCRA", "ID3", "PF4", "CCR10", "SIGLEC7", "NKG7",
             "S100A8", "CXCR3", "CCR5", "CCR3", "CCR4", "PTGDR2", "RORC")

m.ensembl <- gene_symbol$ensembl_gene_id[gene_symbol$external_gene_name %in% m.genes]

m.exprs <- t(pbmc.norm[m.ensembl, 1:(dim(pbmc.norm)[2]-1)])
samples <- rownames(m.exprs)
m.exprs <- as.data.frame(apply(m.exprs, 2, as.numeric))
m.exprs$Sample <- samples

m.symbols <- gene_symbol$external_gene_name[gene_symbol$ensembl_gene_id %in% colnames(m.exprs)]
colnames(m.exprs) <- c(m.symbols, "Sample")

pbmc.merge <- merge(pbmc.max, m.exprs, by='Sample')

ggplot(pbmc.merge[pbmc.merge$Community %in% c(1:17), ],
       aes(x=Dim1, y=Dim2, colour=`CD3E`)) + 
  geom_point(size=2, alpha=0.6) + theme_classic() +
  scale_colour_distiller(palette="Reds", direction=1)

# using these sets of results I'll assign identities to the different communities
# CD3+ T cells: 4, 6, 7, 10:16
# CD20+ (MS4A1) B cells: 1, 5
# CD14+ monocytes: 2, 3
# CD11B+ neutrophils: 8
# CD16+ macrophages: 9

# unknown: 13

# write out these cell identities
pbmc.meta <- cbind.data.frame(pbmc.merge$Sample, pbmc.merge$Community)
colnames(pbmc.meta) <-c("Sample", "Community")
pbmc.meta$CellType <- "Unknown"
pbmc.meta$CellType[pbmc.meta$Community %in% c(4, 6, 7, 10:12, 14:17)] <- "T_cell"
pbmc.meta$CellType[pbmc.meta$Community %in% c(1, 5)] <- "B_cell"
pbmc.meta$CellType[pbmc.meta$Community %in% c(2, 3)] <- "monocyte"
pbmc.meta$CellType[pbmc.meta$Community %in% c(9)] <- "macrophage"
pbmc.meta$CellType[pbmc.meta$Community %in% c(8)] <- "neutrophil"

write.table(pbmc.meta,
            file="~/Dropbox/MNNfigures/10X_data/PBMC/PBMC_meta.tsv",
            row.names=F, sep="\t", quote=FALSE)
