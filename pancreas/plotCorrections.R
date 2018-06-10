# This script performs batch correction with each method and creates t-SNE plots
# for the pancreas data sets; this corresponds to Figure 4 of the manuscript.

library(scran)

dir.create("results", showWarning=FALSE)
sceA <- readRDS("sce.gse81076.rds")
sceB <- readRDS("sce.gse85241.rds")
sceC <- readRDS("sce.gse86473.rds")
sceD <- readRDS("sce.emtab5601.rds")

decA <- readRDS("dec.gse81076.rds")
decB <- readRDS("dec.gse85241.rds")
decC <- readRDS("dec.gse86473.rds")
decD <- readRDS("dec.emtab5601.rds")

# Keeping the top 5000 HVGs with the largest average biological component.
rownames(decC) <- rownames(sceC) <- scater::uniquifyFeatureNames(rownames(sceC), rowData(sceC)$Symbol)
universe <- Reduce(intersect, list(rownames(decA), rownames(decB), rownames(decC), rownames(decD)))
combined.bio <- decA[universe,"bio"] + decB[universe,"bio"] + decC[universe,"bio"] + decD[universe,"bio"]
chosen <- universe[order(combined.bio, decreasing=TRUE)[1:5000]]

# Adjusting the scale to minimize differences in variance.
nout <- multiBatchNorm(sceA[chosen,], sceB[chosen,], sceC[chosen,], sceD[chosen,])
sceA <- nout[[1]]
sceB <- nout[[2]]
sceC <- nout[[3]]
sceD <- nout[[4]]

# Organizing cell type colors.
cell.types <- c(rep("unknown", ncol(sceA)), rep("unknown", ncol(sceB)), sceC$CellType, sceD$CellType)
type.colors <- c(
    Acinar="#ffff00",
    Alpha="#ff0000",
    Beta="#c400ff",
    Delta="#ff7800",
    Ductal="#00f5ff",
    Gamma="#0000ff",
    Other="#000000",
    unknown="grey80"
)
cell.types[!cell.types %in% names(type.colors)] <- "Other"
allcolors <- type.colors[cell.types]

# Organizing batch colors.
batch.id <- rep(1:4, c(ncol(sceA), ncol(sceB), ncol(sceC), ncol(sceD)))
batchcolor <- c("lavender", "lightcoral", "goldenrod1", "lightblue")[batch.id]

# Making a plotting function.
plotFUN <- function(fname, Y, subset=NULL, ..., xlab="tSNE 1",ylab="tSNE 2",main="") {
    if (is.null(subset)) {
        subset <- seq_len(nrow(Y))
    }
    png(fname,width=900,height=700)
    par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
    plot(Y[,1], Y[,2], cex=2, pch=21, bg=allcolors[subset], xlab=xlab, ylab=ylab, main=main) 
    dev.off()
}

plotFUNb <- function(fname, Y, subset=NULL, ...) {
    if (is.null(subset)) {
        subset <- seq_len(nrow(Y))
    }
    png(fname,width=900,height=700)
    par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
    plot(Y[,1], Y[,2], cex=2, pch=21, bg=batchcolor[subset], ...)
    dev.off()
}

######################
# No correction.

combined <- cbind(logcounts(sceA), logcounts(sceB), logcounts(sceC), logcounts(sceD))
t.unc <- as.matrix(t(combined))

## Generating a t-SNE plot.
library(Rtsne)
tsne.unc <- Rtsne(t.unc, perplexity = 30)
plotFUN("results/tsne_unc_type.png", tsne.unc$Y, main="Uncorrected", xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/tsne_unc_batch.png", tsne.unc$Y, main="Uncorrected", xlab="tSNE 1",ylab="tSNE 2")

rm(t.unc)
gc()

######################## 
# Performing the correction with MNN (turned down the sigma to improve mixing).

mnn.out <- mnnCorrect(logcounts(sceA), logcounts(sceB), logcounts(sceC), logcounts(sceD),
    k=20, sigma=0.1, cos.norm.in=TRUE, cos.norm.out=TRUE, var.adj=TRUE)
t.mnn <- as.matrix(t(do.call(cbind, mnn.out$corrected)))

#png(file="results/angles.png",width=900,height=700)
#par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
#hist(mnn.out$angles[[2]],xlab="Angle",ylab="Frequency",main="")
#dev.off()

# Generating a t-SNE plot.
set.seed(0)
tsne.mnn <- Rtsne(t.mnn, perplexity = 30)
plotFUN("results/tsne_mnn_type.png", tsne.mnn$Y, main="MNN", xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/tsne_mnn_batch.png", tsne.mnn$Y, main="MNN", xlab="tSNE 1",ylab="tSNE 2")
gc()

rm(t.mnn)
gc()

######################## 
# Performing the correction with faster MNN.

set.seed(1000)
mnn.out2 <- fastMNN(logcounts(sceA), logcounts(sceB), logcounts(sceC), logcounts(sceD), k=20, approximate=TRUE, cos.norm=TRUE)
t.mnn <- mnn.out2$corrected

# Generating a t-SNE plot.
set.seed(0)
tsne.mnn <- Rtsne(t.mnn, perplexity = 30)
plotFUN("results/tsne_mnn2_type.png", tsne.mnn$Y, main="Fast MNN", xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/tsne_mnn2_batch.png", tsne.mnn$Y, main="Fast MNN", xlab="tSNE 1",ylab="tSNE 2")

rm(t.mnn)
gc()

######################## 
# Performing the correction with limma.

library(limma)
X.lm <- removeBatchEffect(as.matrix(combined), factor(batch.id))
t.lm <- t(X.lm)

## Generating a t-SNE plot.
set.seed(0)
tsne.lm <- Rtsne(t.lm, perplexity = 30)
plotFUN("results/tsne_limma_type.png", tsne.lm$Y, main="limma", xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/tsne_limma_batch.png", tsne.lm$Y, main="limma", xlab="tSNE 1",ylab="tSNE 2")

rm(t.lm)
gc()

######################## 
# Performing the correction with ComBat.

library(sva)
Z <- factor(batch.id)
X.combat <- ComBat(as.matrix(combined), Z, mod=NULL,prior.plots = FALSE)
t.combat <- t(X.combat)

## Generating a t-SNE plot.
set.seed(0)
tsne.combat <- Rtsne(t.combat, perplexity = 30)
plotFUN("results/tsne_combat_type.png", tsne.combat$Y, main="ComBat", xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/tsne_combat_batch.png", tsne.combat$Y, main="ComBat", xlab="tSNE 1",ylab="tSNE 2")

rm(t.combat)
gc()

######################## 
# Performing the correction with CCA. 
# Only using the two batches with the cell type labels,
# as I don't know how to do multiple corrections.

library(Seurat)
SeuC <- CreateSeuratObject(logcounts(sceC))
SeuD <- CreateSeuratObject(logcounts(sceD))
SeuC@meta.data$group <- "group1"
SeuD@meta.data$group <- "group2"

SeuC <- ScaleData(SeuC)
SeuD <- ScaleData(SeuD)
Y <- RunCCA(SeuC, SeuD, genes.use=rownames(sceA), do.normalize=FALSE)
suppressWarnings(Y <- AlignSubspace(Y, grouping.var="group", dims.align=1:20))

t.cca <- Y@dr$cca.aligned@cell.embeddings

## Generating a t-SNE plot.
set.seed(0)
tsne.cca <- Rtsne(t.cca, perplexity = 30)
plotFUN("results/tsne_cca_type.png", tsne.cca$Y, subset=(batch.id%in%c(3,4)), main="CCA", xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/tsne_cca_batch.png", tsne.cca$Y, subset=(batch.id %in% c(3,4)), main="CCA", xlab="tSNE 1",ylab="tSNE 2")

rm(t.cca)
gc()

######################## 
# Performing the correction with CCA, starting from the raw counts and normalizing it their way.

library(Seurat)
SeuC <- CreateSeuratObject(counts(sceC))
SeuD <- CreateSeuratObject(counts(sceD))
SeuC@meta.data$group <- "group1"
SeuD@meta.data$group <- "group2"

SeuC <- NormalizeData(SeuC)
SeuD <- NormalizeData(SeuD)
SeuC <- ScaleData(SeuC)
SeuD <- ScaleData(SeuD)
Y <- RunCCA(SeuC, SeuD, genes.use=rownames(sceA))
suppressWarnings(Y <- AlignSubspace(Y, grouping.var="group", dims.align=1:20))

t.cca <- Y@dr$cca.aligned@cell.embeddings

## Generating a t-SNE plot.
set.seed(0)
tsne.cca2 <- Rtsne(t.cca, perplexity = 30)
plotFUN("results/tsne_cca2_type.png", tsne.cca2$Y, subset=(batch.id%in%c(3,4)), main="CCA native", xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/tsne_cca2_batch.png", tsne.cca2$Y, subset=(batch.id%in%c(3,4)), main="CCA native", xlab="tSNE 1",ylab="tSNE 2")

rm(t.cca)
gc()

######################## 
# END
