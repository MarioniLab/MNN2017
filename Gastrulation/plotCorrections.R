library(scran)
require(Rtsne)
load("mesoandwolf.Rdata") 
dir.create("results", showWarnings=FALSE)

#########################################################

# Organizing cell type and cell stage labels in the mesoderm data. 
meta.file <- file.path("meso", "metadata.txt")
if (!file.exists(meta.file)) { 
    download.file("http://gastrulation.stemcells.cam.ac.uk/data/metadata.txt", meta.file)
}
meta.data <- read.table(meta.file, header=TRUE, stringsAsFactors=FALSE, row.names=1) 
stopifnot(identical(rownames(meta.data), colnames(data.meso)))

meso.stage <- meta.data$embryoStage
meso.stage[meso.stage=="HF"]<-"E7.75"
meso.stage[meso.stage=="NP"]<-"E7.5"
meso.stage[meso.stage=="PS"]<-"E7.0"

meso.cluster <- meta.data$cluster

# Keeping only cells in the first two stages.
select.meso <- which(meso.stage=="E6.5" | meso.stage=="E7.0")
meso.stage <- meso.stage[select.meso]
meso.cluster <- meso.cluster[select.meso]
data.meso.2 <- data.meso[,select.meso]

# Removing some clusters for clarity.
excl.meso <- which(meso.cluster=="pink" | meso.cluster=="magenta" | meso.cluster=="red")
meso.stage <- meso.stage[-excl.meso]
meso.cluster <- meso.cluster[-excl.meso]
data.meso.2 <- data.meso.2[,-excl.meso]

# Cleaning up the names.
meso.cluster[meso.cluster=="turquoise"]<-"Epiblast-PS"
meso.cluster[meso.cluster=="yellow"]<-"Blood prog."
meso.cluster[meso.cluster=="blue"]<-"Mesoderm prog."
meso.cluster[meso.cluster=="red"]<-"Endothelial"
meso.cluster[meso.cluster=="darkorange"]<-"Post-mesoderm"
new.labs.meso <- paste0(meso.cluster, " (", meso.stage, ")")

#########################################################

# Organizing cell type and cell stage labels in Wolf's data. 
meta.file <- file.path("wolf", "Cell_clusters.txt")
datainfo <- read.table(meta.file, sep="\t", stringsAsFactors = FALSE, header=TRUE)
data.wolf.2 <- data.wolf
wolf.stage <- datainfo$Stage
wolf.lineage <- datainfo$Lineage

# Excluding cells in VE.
keep <- wolf.stage %in% c("E5.5", "E6.5", "E6.75") & !grepl("^VE", wolf.lineage)
data.wolf.2 <- data.wolf.2[,keep.wolf]
wolf.lineage <- wolf.lineage[keep.wolf]
wolf.stage <- wolf.stage[keep.wolf]

# Cleaning up the names.
wolf.lineage[grep("^epiblast", wolf.lineage)] <- "Epiblast"
wolf.lineage[wolf.lineage=="PS"] <- "PS"
new.labs.wolf <- paste0(wolf.lineage, " (", wolf.stage, ")")

#########################################################
# Setting up a central plotting function.

raw.all <- cbind(data.wolf.2,data.meso.2)
first.batch <- rep(c(TRUE, FALSE), c(ncol(data.wolf.2), ncol(data.meso.2)))

typing <- c(new.labs.wolf, new.labs.meso)
cols <- rep("grey80", length(typing))
cols[typing=="Epiblast-PS (E6.5)"] <- "orange"
cols[typing=="Epiblast (E4.5)"] <- "lightblue"
cols[typing=="Epiblast (E5.5)"] <- "dodgerblue"
cols[typing=="Epiblast (E6.5)"] <- "blue"
cols[typing=="Epiblast (E6.75)"] <- "darkblue"

plotFUN <- function(fname, Y, col=NULL, xlab="PC1",ylab="PC2", ...) { 
    if (is.null(col)) { 
         col <- cols
    }
    pdf(file=fname,width=9,height=7)
    par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=1.2,cex.main=1.4,cex.lab=1.4)
    shuffler <- sample(nrow(Y)) 
    plot(Y[shuffler,1], Y[shuffler,2], pch=21, cex=2, col="black", bg=col[shuffler], 
         type=ifelse(!is.null(subset), "n", "p"), ..., xlab=xlab, ylab=ylab)

    dev.off()
}

#########################################################

# Raw
pc.unc <- prcomp(t(raw.all[any.hvg,]), rank=2)
plotFUN("results/uncGast.pdf", pc.unc$x)

# MNN
mnn.out <-mnnCorrect(data.wolf.2, data.meso.2, subset.row=any.hvg, k=20, sigma=0.1, cos.norm=TRUE, svd.dim=NA)
X.mnn <-cbind(mnn.out$corrected[[1]], mnn.out$corrected[[2]])
pc.mnn <-prcomp(t(X.mnn), rank=2)
plotFUN("results/mnnGast.pdf", pc.mnn$x)

# limma
library(limma)
X.lm <- removeBatchEffect(raw.all[any.hvg,], factor(first.batch))
pc.lm <-prcomp(t(X.lm), rank=2)
plotFUN("results/limmaGast.pdf", pc.lm$x)

# ComBat (needs addition of some small value, for some reason).
library(sva)
X.com <- ComBat(raw.all[any.hvg,] + 0.01, first.batch, mod=NULL, prior.plots = FALSE)
pc.com <-prcomp(t(X.com), rank=2)
plotFUN("results/combatGast.pdf", pc.com$x)

## Legend.
#pdf(file="leg_detailed.pdf",width=900,height=700)
#par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
##plot(tsne.c$Y[1,1],tsne.c$Y[1,2], pch=3,cex=4,col=alpha(allcolors[1],0.6),main="legend",xlim=c(-20,30),ylim=c(-13,13),xlab="tSNE 1",ylab="tSNE 2")
#plot(0,0,type="n")
#forleg<-table(celltypes,allcolors)
#legend("bottomright", "(x,y)", legend = leg.txt, col =unique(allcolors) , pch = 17,cex = 2.5,bty = "n",lwd = 3,lty=0)   #, trace = TRUE)
#legend("bottomleft", "(x,y)", legend = leg.txt ,pt.bg=forcoloringleg, pch = 21,cex = 2.5,bty = "n",lwd = 3,lty=0)   #, trace = TRUE)
##legend(x=7,y=15, legend = c("CEL-Seq 1","CEL-Seq 2","SMART-Seq 1","SMART-Seq 2"), col ="black" , pch = c(4,1,18,3),cex = 2.5,bty = "n")   #, trace = TRUE)
#dev.off()

#########################################################

# Colouring by expression of certain genes.
library(RColorBrewer)
rbPal <- colorRampPalette(c('grey95', 'blue'))

# Mesp1
col.mesp <- rbPal(100)[cut(raw.all["ENSMUSG00000030544",], breaks = 100)]
plotFUN("results/mnn_Mesp1.pdf", pc.mnn$x, col=col.mesp)

# T (Brachyury)
col.t <- rbPal(100)[cut(raw.all["ENSMUSG00000062327",], breaks = 100)]
plotFUN("results/mnn_T.pdf", pc.mnn$x, col=col.t)

# Mettl7a3
col.mettl <- rbPal(100)[cut(raw.all["ENSMUSG00000024406",], breaks = 100)]
plotFUN("results/mnn_Mettl7a3.pdf", pc.mnn$x, col=col.mettl)

# Legend.
pdf(file="results/leg_exprs.pdf")
plot(0, 0, type="n",xlab="", ylab="", bty="n", axes=FALSE)
y <- 1:100/100 - 0.5
rect(-0.1, y, 0.1, y+0.015, col=rbPal(100), border=NA)
text(0, min(y), pos=1, "Low")
text(0, max(y), pos=3, "High")
dev.off()

#########################################################
# End.

