# combine human preimplantation embryo data sets
# might be tricky as one data set has very few HVGs
library(ggplot2)
library(scran)
library(RColorBrewer)
library(Rtsne)
library(stringr)
#library(WGCNA)
library(biomaRt)
library(igraph)
source("/mnt/scratcha/jmlab/morgan02/10X/MNN/src/theme_mike.R")
source("/mnt/scratcha/jmlab/morgan02/10X/MNN/src/single_cell_functions.R")

# read in the 68K PBMC data
pbmc <- read.table("/mnt/scratcha/jmlab/morgan02/10X/MNN/data/PBMC/PBMC_norm.tsv.gz",
                   h=T, sep="\t", stringsAsFactors=FALSE)
rownames(pbmc) <- pbmc$gene_id

pbmc.hvg <- read.table("/mnt/scratcha/jmlab/morgan02/10X/MNN/data/PBMC/PBMC_hvg.tsv",
                       h=T, sep="\t", stringsAsFactors=FALSE)

# Generate tSNE on the 68K PBMCs
pbmc.select <- pbmc[pbmc.hvg$HVG, ]
pbmc.select[is.na(pbmc.select)] <- 0

print(dim(pbmc.select))
set.seed(42)
pbmc.map <- tsne_wrapper(pbmc.select[, 1:(dim(pbmc)[2]-1)])

# as.matrix errors as there are too many, columns
# fair enough, there are 68K!
pbmc.matrix <- as.matrix(pbmc.select[, 1:(dim(pbmc)[2]-1)])
pbmc.snn <- buildSNNGraph(pbmc.matrix, pc.approx=TRUE,
                          k=25, d=30)

# use the walktrap algorithm to form communities
# try a few different values for the step parameter
pbmc.community <- cluster_walktrap(pbmc.snn, steps=4)
n.comms <- length(pbmc.community)

# vertex community membership 
pbmc.members <- pbmc.community$membership
pbmc.clusters <- cbind.data.frame(colnames(pbmc[, 1:(dim(pbmc)[2]-1)]), pbmc.members)
colnames(pbmc.clusters) <- c("Sample", "Community")
pbmc.clusters$Community <- as.factor(pbmc.clusters$Community)

pbmc.max <- merge(pbmc.map, pbmc.clusters, by='Sample')

#comm.tsne <- ggplot(pbmc.max,
#                    aes(x=Dim1, y=Dim2, colour=Community)) + 
#  geom_point(size=1) + theme_classic() +
#  scale_colour_Publication()

#ggsave(comm.tsne,
#	filename="/mnt/scratcha/jmlab/morgan02/10X/MNN/plots/PBMC68k_tSNE.png",
#	height=4.75, width=4.75, dpi=300)

write.table(pbmc.max,
		file="/mnt/scratcha/jmlab/morgan02/10X/MNN/data/PBMC/PBMC_meta.tsv",
		sep="\t", quote=FALSE, row.names=FALSE)