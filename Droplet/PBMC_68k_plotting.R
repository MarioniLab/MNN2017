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