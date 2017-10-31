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
source("SomeFuncs/convenience_functions.R")

################################################################################################
################################################################################################
# download the 10X data
tcell_path <- "Droplet/10X_data/Tcell/"

# this has to be an absolute path, relative paths break
download_sample(sample_name="t_4k", sample_dir=tcell_path,
                host="https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/")

tcell.10x <- load_cellranger_matrix(tcell_path)

# filter non-zero genes and normalize tcell 10X UMI counts
tcell.nz <- get_nonzero_genes(tcell.10x)

################################################################################################
################################################################################################
# normalize using size factors
tcell.norm <- size_factor_normalize(exprs(tcell.10x), cell.sparse=0.95,
                                   gene.sparse=0.99, cluster.size=50)

# map ensembl IDs to gene symbols to select a couple of marker genes for visualisation
all.genes <- rownames(tcell.norm)

ensembl <- useEnsembl(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', GRCh=37)

gene_symbol <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                     filters='ensembl_gene_id', mart=ensembl,
                     values=all.genes)
################################################################################################
################################################################################################
# select highly variable genes (~1-2k will be more than sufficient I think)
tcell.hvg <- find_hvg(tcell.norm[, 1:(dim(tcell.norm)[2]-1)], plot=FALSE)

# the mean-CV^2 plot looks OK, even with the high number of unit counts
# there are 1192 HVG at p<=0.01

write.table(tcell.norm,
            file="Droplet/10X_data/Tcell/Tcell_norm.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

hvg.df <- cbind.data.frame(names(tcell.hvg)[tcell.hvg])
colnames(hvg.df) <- c("HVG")
write.table(hvg.df,
            file="Droplet/10X_data/Tcell/Tcell_hvg.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

################################################################################################
################################################################################################
# tSNE embedding for visualization
tcell.select <- tcell.norm[tcell.hvg, ]
tcell.select[is.na(tcell.select)] <- 0

set.seed(42)
tcell.map <- tsne_wrapper(tcell.select[, 1:(dim(tcell.norm)[2]-1)], perplexity=100)

################################################################################################
################################################################################################
# pull out marker genes for tcell populations to add identities to different clusters
# also try SNNgraph for clustering!  Use just 20 dimensions from the PCA to build the graph
tcell.snn <- buildSNNGraph(as.matrix(tcell.norm[, 1:(dim(tcell.norm)[2]-1)]),
                          k=5, d=20)

# use the walktrap algorithm to form communities
# try a few different values for the step parameter
tcell.community <- cluster_walktrap(tcell.snn, steps=4)
n.comms <- length(tcell.community)

# vertex community membership 
tcell.members <- tcell.community$membership
tcell.clusters <- cbind.data.frame(colnames(tcell.norm[, 1:(dim(tcell.norm)[2]-1)]), tcell.members)
colnames(tcell.clusters) <- c("Sample", "Community")
tcell.clusters$Community <- as.factor(tcell.clusters$Community)

tcell.max <- merge(tcell.map, tcell.clusters, by='Sample')

comm.tsne <- ggplot(tcell.max,
                    aes(x=Dim1, y=Dim2, colour=Community)) + 
  geom_point(size=3) + theme_classic()

ggsave(comm.tsne,
       filename="Droplet/10X_data/Tcell/Tcell_communities_tSNE.png",
       height=5.75, width=5.75, dpi=300)

comm.tsne

# highlight clusters with marker genes after SNN and randomwalk-based clustering
# select a bunch of marker genes for the expected cell types
m.genes <- c("CD3E", "CD8A", "CD8B", "CD4", "CD28", "CD8", "TBX21", "IKAROS", "IL2RA", "CD44", "SELL",
             "CCR7", "CCR6", "IL7R", "CD27", "FOXP3", "PTCRA", "ID3", "PF4", "CCR10",
             "CXCR3", "CCR5", "CCR3", "CCR4", "PTGDR2", "RORC")

m.ensembl <- gene_symbol$ensembl_gene_id[gene_symbol$external_gene_name %in% m.genes]

m.exprs <- t(tcell.norm[m.ensembl, 1:(dim(tcell.norm)[2]-1)])
samples <- rownames(m.exprs)
m.exprs <- as.data.frame(apply(m.exprs, 2, as.numeric))
m.exprs$Sample <- samples

m.symbols <- gene_symbol$external_gene_name[gene_symbol$ensembl_gene_id %in% colnames(m.exprs)]
colnames(m.exprs) <- c(m.symbols, "Sample")

tcell.merge <- merge(tcell.max, m.exprs, by='Sample')

ggplot(tcell.merge[tcell.merge$Community %in% c(1:15), ],
       aes(x=Dim1, y=Dim2, colour=`CD8B`)) + 
  geom_point(size=2, alpha=0.6) + theme_classic() +
  scale_colour_distiller(palette="Reds", direction=1)

# using these sets of results I'll assign identities to the different communities
# CD4+ T cells: 1:15
# CD8+ T cells: 2, 3, 4, 9, 13
# Th1 T cells (TBX21): 6, 12
# T cells: 1:15

# write out these cell identities
tcell.meta <- cbind.data.frame(tcell.merge$Sample, tcell.merge$Community)
colnames(tcell.meta) <-c("Sample", "Community")
tcell.meta$CellType <- "T_cell"
tcell.meta$CellType[tcell.merge$CD4 >= 1] <- "CD4_Tcell"
tcell.meta$CellType[tcell.merge$Community %in% c(2:4, 9, 13) & tcell.merge$CD8B >= 1] <- "CD8_Tcell"
tcell.meta$CellType[tcell.merge$Community %in% c(6, 12) & tcell.merge$TBX21 >= 1] <- "Th1_Tcell"

write.table(tcell.meta,
            file="Droplet/10X_data/Tcell/Tcell_meta.tsv",
            row.names=F, sep="\t", quote=FALSE)
