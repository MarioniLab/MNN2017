# mnncorrect on 68K PBMC and 4k T cells
library(ggplot2)
library(WGCNA)
library(RColorBrewer)
library(scran)
library(igraph)

pbmc.tsne <- read.table("Droplet/Results/corrected_tSNE.tsv",
                        h=TRUE, sep="\t", stringsAsFactors=FALSE)

pbmc.tsne$Study <- unlist(lapply(strsplit(pbmc.tsne$Sample, split=".", fixed=T),
                                 FUN=function(X) paste0(X[1])))
pbmc.tsne$Study[pbmc.tsne$Study == "pbmc"] <- "PBMC"
pbmc.tsne$Study[pbmc.tsne$Study == "tcell"] <- "Tcells"

pbmc.meta <- read.table("Droplet/10X_Data/PBMC/PBMC_meta.tsv",
                        h=TRUE, sep="\t", stringsAsFactors=FALSE)
pbmc.meta$Sample <- paste("pbmc", pbmc.meta$Sample, sep=".")

tcell.meta <- read.table("Droplet/10X_Data/Tcell/Tcell_meta.tsv",
                         sep="\t", h=TRUE, stringsAsFactors=FALSE)
tcell.meta$Sample <- paste("tcell", tcell.meta$Sample, sep=".")
tcell.meta$Sample <- gsub(tcell.meta$Sample, pattern="-", replacement=".")

all.meta <- do.call(rbind.data.frame, list("Tcell"=tcell.meta,
                                           "PBMC"=pbmc.meta[, c("Sample", "Community", "CellType")]))
pbmc.merge <- merge(pbmc.tsne, all.meta, by='Sample')

# randomly plotting points obscures the T cells because there are many fewer of them
# it is better to plot the Tcells over the PBMCs to see how they match up
study.cols <- c("#33a02c", "#cab2d6")
names(study.cols) <- c("PBMC", "Tcells")
col.vec <- labels2colors(pbmc.merge$Study, colorSeq=study.cols)

png("Droplet/mnnCorrect_68K_PBMC-tSNE_bystudy.png",
    height=5.75, width=6.95, res=300, units="in")
par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)
plot(x=pbmc.merge$Dim1, y=pbmc.merge$Dim2,
     bg=col.vec, pch=21, col='black', cex=1.1, lwd=0.75,
     xlab="tSNE 1", ylab="tSNE 2",
     ylim=c(-35, 35), xlim=c(-35, 35))
legend("topright", "(x, y)", legend=names(study.cols), pch=16, col=study.cols,
       bty="n", lwd=0.75, lty=0, cex=1.5,
       inset=c(-0.35, 0))
dev.off()

cell.hex <- c("#386cb0", "#fdb462",
               "#7fc97f","#ef3b2c","#662506",
               "#a6cee3","#fb9a99","#984ea3",
               "#ffff33", "#b2df8a", "#BDBABD")
cell.cols <- labels2colors(pbmc.merge$CellType, colorSeq=cell.hex)
names(cell.cols) <- pbmc.merge$CellType

png("Droplet/mnnCorrect_68K_PBMC-tSNE_bycell.png",
     height=5.75, width=6.95, res=300, units="in")
par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)
plot(x=pbmc.merge$Dim1, y=pbmc.merge$Dim2,
     bg=alpha(cell.cols, 0.75), pch=21, col='black', cex=1.1, lwd=0.75,
     xlab="tSNE 1", ylab="tSNE 2",
     ylim=c(-35, 35), xlim=c(-35, 35))
legend("topright", legend=unique(names(cell.cols)),
       pch=16, col=unique(cell.cols), ncol=1,
       bty="n", lwd=0.75, lty=0, cex=1,
       inset=c(-0.35, 0))
dev.off()


# construct an SNN graph on the tSNE dimensions
tsne.snn <- buildSNNGraph(t(pbmc.merge[, c("Dim1", "Dim2")]), k=100)
tsne.community <- cluster_walktrap(tsne.snn, steps=8)

tsne.cluster <- cbind.data.frame(pbmc.merge$Sample, tsne.community$membership)
colnames(tsne.cluster) <- c("Sample", "NewCommunity")

pbmc.uber <- merge(pbmc.merge, tsne.cluster, by='Sample')

# # get marker expression to identify cell types 
# # map ensembl IDs to gene symbols to select a couple of marker genes for visualisation
# all.genes <- correct.df$gene_id
# 
# ensembl <- useEnsembl(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', GRCh=37)
# 
# gene_symbol <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
#                      filters='ensembl_gene_id', mart=ensembl,
#                      values=all.genes)
# 
# m.genes <- c("CD14", "PTPRC", "FCGR3A", "ITGAX", "ITGAM", "CD19", "HLA-DRB1", "FCGR2B", "FCGR2A",
#              "CD3E", "CD4", "CD8A","CD8B", "CD28", "CD8", "TBX21", "IKAROS", "IL2RA", "CD44", "SELL",
#              "CCR7", "MS4A1", "CD68", "CD163", "IL5RA", "SIGLEC8", "KLRD1", "NCR1", "CD22", "IL3RA",
#              "CCR6", "IL7R", "CD27", "FOXP3", "PTCRA", "ID3", "PF4", "CCR10", "SIGLEC7", "NKG7",
#              "S100A8", "CXCR3", "CCR5", "CCR3", "CCR4", "PTGDR2", "RORC")
# 
# m.ensembl <- gene_symbol$ensembl_gene_id[gene_symbol$external_gene_name %in% m.genes]

# ggplot(pbmc.uber, aes(x=Dim1, y=Dim2, fill=as.factor(NewCommunity))) +
#   geom_point(shape=21, size=2)

# plot the marker genes used to identify cell types against the communities to assign cell types

pbmc.merge$NewCell <- "Unknown"
pbmc.merge$NewCell[pbmc.merge$Community %in% c(1:4, 6:9, 11)] <- "T_cell"
pbmc.merge$NewCell[pbmc.merge$Community == 13] <- "B_cell"
pbmc.merge$NewCell[pbmc.merge$Community == 5] <- "DCs"
pbmc.merge$NewCell[pbmc.merge$Community == 2] <- "macrophage"
pbmc.merge$NewCell[pbmc.merge$Community == 3] <- "monocyte"
pbmc.merge$NewCell[pbmc.merge$Community == 4] <- "NKcell"
pbmc.merge$NewCell[pbmc.merge$Community == 10] <- "NKT"

