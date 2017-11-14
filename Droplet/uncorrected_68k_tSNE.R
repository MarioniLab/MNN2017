# tSNE on combined but uncorrected 4K Tcells and 68K PBMCs
library(ggplot2)
library(scran)
library(RColorBrewer)
library(Rtsne)
library(stringr)
library(biomaRt)
library(igraph)
source("/mnt/scratcha/jmlab/morgan02/10X/MNN/src/single_cell_functions.R")

tcell <- read.table("/mnt/scratcha/jmlab/morgan02/10X/MNN/data/Tcell/Tcell_norm.tsv",
                    h=T, sep="\t", stringsAsFactors=FALSE)

rownames(tcell) <- tcell$gene_id

tcell.hvg <- read.table("/mnt/scratcha/jmlab/morgan02/10X/MNN/data/Tcell/Tcell_hvg.tsv",
                        h=T, sep="\t", stringsAsFactors=FALSE)

tcell.meta <- read.table("/mnt/scratcha/jmlab/morgan02/10X/MNN/data/Tcell/Tcell_meta.tsv",
                         h=T, sep="\t", stringsAsFactors=FALSE)
tcell.meta$Sample <- paste("tcell", tcell.meta$Sample, sep=".")
tcell.meta$Sample <- gsub(tcell.meta$Sample, pattern="-", replacement=".")

pbmc <- read.table("/mnt/scratcha/jmlab/morgan02/10X/MNN/data/PBMC/PBMC_norm.tsv.gz",
                   h=T, sep="\t", stringsAsFactors=FALSE)
rownames(pbmc) <- pbmc$gene_id

pbmc.hvg <- read.table("/mnt/scratcha/jmlab/morgan02/10X/MNN/data/PBMC/PBMC_hvg.tsv",
                       h=T, sep="\t", stringsAsFactors=FALSE)

pbmc.meta <- read.table("/mnt/scratcha/jmlab/morgan02/10X/MNN/data/PBMC/PBMC_meta.tsv",
                        h=T, sep="\t", stringsAsFactors=FALSE)
pbmc.meta$Sample <- paste("pbmc", pbmc.meta$Sample, sep=".")
pbmc.meta$Sample <- gsub(pbmc.meta$Sample, pattern="-", replacement=".")

# select the total number of overlapping genes
# and the numer of intersecting HVGs
select.genes <- intersect(rownames(tcell), rownames(pbmc))
select.hvg <- intersect(tcell.hvg$HVG, pbmc.hvg$HVG)

# there are only ~300 HVG shared between both data sets
tcell.batch <- tcell[tcell$gene_id %in% select.genes, 1:(dim(tcell)[2]-1)]
colnames(tcell.batch) <- paste("tcell", colnames(tcell.batch), sep=".")

pbmc.batch <- pbmc[pbmc$gene_id %in% select.genes, 1:(dim(pbmc)[2]-1)]
colnames(pbmc.batch) <- paste("pbmc", colnames(pbmc.batch), sep=".")

# check they have the same number of rows and ordering
dim(tcell.batch)[1] == dim(pbmc.batch)[1]

if(!all(rownames(tcell.batch) == rownames(pbmc.batch))){
  # if these don't match enforce the same ordering
  tcell.batch <- tcell.batch[rownames(pbmc.batch), ]
}

# save these for later
tcell.cells <- colnames(tcell.batch)
pbmc.cells <- colnames(pbmc.batch)

batch.meta <- do.call(rbind.data.frame, list("Tcell"=tcell.meta,
                                             "PBMC"=pbmc.meta[, c("Sample", "Community", "CellType")]))

# tSNE on the raw uncorrected data first
tcell.batch$gene_id <- rownames(tcell.batch)
pbmc.batch$gene_id <- rownames(pbmc.batch)

# merge on gene_id
all.uncorrect <- merge(pbmc.batch, tcell.batch, by='gene_id')

# assess correction using tSNE
set.seed(42)
uncorrect.tsne <- tsne_wrapper(all.uncorrect[, 2:dim(all.uncorrect)[2]])

# get study info and other meta data
uncorrect.tsne$Study <- ""
uncorrect.tsne$Study[uncorrect.tsne$Sample %in% tcell.cells] <- "Tcells"
uncorrect.tsne$Study[uncorrect.tsne$Sample %in% pbmc.cells] <- "PBMC"

tsne.merge <- merge(uncorrect.tsne, batch.meta, by='Sample')

# output this table to plot later
write.table(correct.max,
            filename="/mnt/scratcha/jmlab/morgan02/10X/MNN/data/UncorrectedtSNE_dims.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")




