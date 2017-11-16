# tSNE on combined but uncorrected 4K Tcells and 68K PBMCs
library(ggplot2)
library(scran)
library(RColorBrewer)
library(Rtsne)
library(stringr)
library(biomaRt)
library(igraph)
source("SomeFuncs/single_cell_functions.R")

tcell <- read.table("Droplet/10X_Data/Tcell/Tcell_norm.tsv",
                    h=T, sep="\t", stringsAsFactors=FALSE)

rownames(tcell) <- tcell$gene_id

tcell.hvg <- read.table("Droplet/10X_Data/Tcell/Tcell_hvg.tsv",
                        h=T, sep="\t", stringsAsFactors=FALSE)

tcell.meta <- read.table("Droplet/10X_Data/Tcell/Tcell_meta.tsv",
                         h=T, sep="\t", stringsAsFactors=FALSE)
tcell.meta$Sample <- paste("tcell", tcell.meta$Sample, sep=".")
tcell.meta$Sample <- gsub(tcell.meta$Sample, pattern="-", replacement=".")

pbmc <- read.table("Droplet/10X_Data/PBMC/PBMC_norm.tsv.gz",
                   h=T, sep="\t", stringsAsFactors=FALSE)
rownames(pbmc) <- pbmc$gene_id

pbmc.hvg <- read.table("Droplet/10X_Data/PBMC/PBMC_hvg.tsv",
                       h=T, sep="\t", stringsAsFactors=FALSE)

pbmc.meta <- read.table("Droplet/10X_Data/PBMC/PBMC_meta.tsv",
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
write.table(tsne.merge,
            file="Droplet/10X_Data/UncorrectedtSNE_dims.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

###############################################
# Plot the uncorrected, naively combined tSNE #
###############################################

tsne.merge$Study <- unlist(lapply(strsplit(as.character(tsne.merge$Sample), split='.', fixed=TRUE),
                                  FUN=function(X) paste0(X[1])))
tsne.merge$Study[tsne.merge$Study == "pbmc"] <- "PBMC"
tsne.merge$Study[tsne.merge$Study == "tcell"] <- "Tcell"

# randomly plotting points obscures the T cells because there are many fewer of them
# it is better to plot the Tcells over the PBMCs to see how they match up
study.cols <- c("#33a02c", "#cab2d6")
names(study.cols) <- c("PBMC", "Tcells")
col.vec <- labels2colors(tsne.merge$Study, colorSeq=study.cols)

png("Droplet/uncorrected_68K_PBMC-tSNE_bystudy.png",
    height=5.75, width=6.95, res=300, units="in")
par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)
plot(x=tsne.merge$Dim1, y=tsne.merge$Dim2,
     bg=col.vec, pch=21, col='black', cex=1.1, lwd=0.75,
     xlab="tSNE 1", ylab="tSNE 2",
     ylim=c(-35, 35), xlim=c(-35, 35))
legend("topright", "(x, y)", legend=names(study.cols), pch=16, col=study.cols,
       bty="n", lwd=0.75, lty=0, cex=1.5,
       inset=c(-0.35, 0))
dev.off()

## now with the cell type labels

cell.hex <- c("#386cb0", "#fdb462",
              "#7fc97f","#ef3b2c","#662506",
              "#a6cee3","#fb9a99","#984ea3",
              "#ffff33", "#b2df8a", "#BDBABD")
cell.cols <- labels2colors(tsne.merge$CellType, colorSeq=cell.hex)
names(cell.cols) <- tsne.merge$CellType

png("Droplet/uncorrected_68K_PBMC-tSNE_bycell.png",
    height=5.75, width=6.95, res=300, units="in")
par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)
plot(x=tsne.merge$Dim1, y=tsne.merge$Dim2,
     bg=alpha(cell.cols, 0.75), pch=21, col='black', cex=1.1, lwd=0.75,
     xlab="tSNE 1", ylab="tSNE 2",
     ylim=c(-35, 35), xlim=c(-35, 35))
legend("topright", legend=unique(names(cell.cols)),
       pch=16, col=unique(cell.cols), ncol=1,
       bty="n", lwd=0.75, lty=0, cex=1,
       inset=c(-0.35, 0))
dev.off()




