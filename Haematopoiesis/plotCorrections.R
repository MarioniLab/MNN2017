# This script performs batch correction with each method and creates t-SNE plots
# for the haematopoietic data sets; this corresponds to Figure 3 of the manuscript.

library(scran)
dir.create("results", showWarning=FALSE)
sceA <- readRDS("haem_data_A.rds")
sceF <- readRDS("haem_data_F.rds")

# Load data.
sceF$CellType[sceF$CellType=="other"] <- "Unsorted"
sceA$CellType[sceA$CellType=="ERY"] <- "MEP"
sce <- cbind(sceF, sceA)
sce$Batch <- rep(c(TRUE, FALSE), c(ncol(sceF), ncol(sceA)))

# Adding colours.
color.legendF <- c(MEP="orange", GMP="chartreuse4", CMP="magenta", 
                   HSPC="cyan", LTHSC="dodgerblue", MPP="blue", LMPP="light blue", Unsorted="grey")
colmatF <- col2rgb(color.legendF) 
colmatA <- colmatF + 100 # A lighter shade.
colmatA[colmatA > 255] <- 255

color.legendA <- setNames(rgb(colmatA[1,], colmatA[2,], colmatA[3,], maxColorValue=255), names(color.legendF))
allcolors <- c(color.legendF[sceF$CellType], color.legendA[sceA$CellType])
first.batch <- rep(c(TRUE, FALSE), c(ncol(sceF), ncol(sceA)))

# Only keeping common cell types for PCA.
pca.retain <- sce$CellType %in% c("MEP", "GMP", "CMP") 

# Making a plotting function.
plotFUN <- function(fname, Y, subset=NULL, ..., xlab="tSNE 1",ylab="tSNE 2",main="") {
    if (is.null(subset)) {
        subset <- seq_len(nrow(Y))
    }
    png(fname,width=900,height=700)
    par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
    plot(Y[,1], Y[,2], cex=2,
         pch=ifelse(first.batch, 21, 1)[subset], 
         col=ifelse(first.batch, "black", allcolors)[subset],
         bg=allcolors[subset], xlab=xlab, ylab=ylab, main=main) 
    dev.off()
}

batchcolor <- c("lavender","lightcoral")[first.batch + 1]
plotFUNb <- function(fname, Y, subset=NULL, ...) {
    if (is.null(subset)) {
        subset <- seq_len(nrow(Y))
    }
    png(fname,width=900,height=700)
    par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
    plot(Y[,1], Y[,2], cex=2,
         pch=ifelse(first.batch, 21, 1)[subset], 
         col=ifelse(first.batch, "black", batchcolor)[subset],
         bg=batchcolor[subset], ...)
    dev.off()
}

######################
# No correction.

t.unc <- as.matrix(t(logcounts(sce)))

## Generating a t-SNE plot.
library(Rtsne)
tsne.unc <- Rtsne(t.unc, perplexity = 90)
plotFUN("results/tsne_unc_type.png", tsne.unc$Y, main="Uncorrected", xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/tsne_unc_batch.png", tsne.unc$Y, main="Uncorrected", xlab="tSNE 1",ylab="tSNE 2")

gc()

# Generating a PCA plot.
pca.unc <- prcomp(t.unc[pca.retain,], rank=2)
plotFUN("results/pca_unc_type.png", pca.unc$x, subset=pca.retain, main="Uncorrected", xlab="PC 1",ylab="PC 2")
plotFUNb("results/pca_unc_batch.png", pca.unc$x, subset=pca.retain, main="Uncorrected", xlab="PC 1",ylab="PC 2")

rm(t.unc)
gc()

######################## 
# Performing the correction with MNN (turned down the sigma to improve mixing).

mnn.out <- mnnCorrect(logcounts(sceF), logcounts(sceA), k=20, sigma=0.1,cos.norm.in=TRUE, cos.norm.out=TRUE, var.adj=TRUE,compute.angle=TRUE)
t.mnn <- as.matrix(t(do.call(cbind, mnn.out$corrected)))

#png(file="results/angles.png",width=900,height=700)
#par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
#hist(mnn.out$angles[[2]],xlab="Angle",ylab="Frequency",main="")
#dev.off()

# Generating a t-SNE plot.
set.seed(0)
tsne.mnn <- Rtsne(t.mnn, perplexity = 90)
plotFUN("results/tsne_mnn_type.png", tsne.mnn$Y, main="MNN", xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/tsne_mnn_batch.png", tsne.mnn$Y, main="MNN", xlab="tSNE 1",ylab="tSNE 2")
gc()

# Generating a PCA plot.
pca.mnn <- prcomp(t.mnn[pca.retain,], rank=2)
plotFUN("results/pca_mnn_type.png", pca.mnn$x, subset=pca.retain, main="MNN", xlab="PC 1",ylab="PC 2")
plotFUNb("results/pca_mnn_batch.png", pca.mnn$x, subset=pca.retain, main="MNN", xlab="PC 1",ylab="PC 2")

## Generating diffusion map plots.
#library(destiny)
#dm<-DiffusionMap(t.mnn,n_local = 150)
#plotFUN("results/mnnFAdm12.png", dm@eigenvectors[,1:2], main="MNN",  xlab="DC 1",ylab="DC 2")
#plotFUN("results/mnnFAdm23.png", cbind(dm@eigenvectors[,2],dm@eigenvectors[,3]), main="MNN",  xlab="DC 2",ylab="DC 4")
#plotFUN("results/mnnFAdm13.png", cbind(dm@eigenvectors[,1],dm@eigenvectors[,3]), main="MNN",  xlab="DC 2",ylab="DC 4")

rm(t.mnn)
gc()

######################## 
# Performing the correction with faster MNN.

set.seed(1000)
mnn.out2 <- fastMNN(logcounts(sceF), logcounts(sceA), k=20, approximate=TRUE, cos.norm=TRUE)
t.mnn <- mnn.out2$corrected

# Generating a t-SNE plot.
set.seed(0)
tsne.mnn <- Rtsne(t.mnn, perplexity = 90)
plotFUN("results/tsne_mnn2_type.png", tsne.mnn$Y, main="Fast MNN", xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/tsne_mnn2_batch.png", tsne.mnn$Y, main="Fast MNN", xlab="tSNE 1",ylab="tSNE 2")
gc()

# Generating a PCA plot.
pca.mnn <- prcomp(t.mnn[pca.retain,], rank=2)
plotFUN("results/pca_mnn2_type.png", pca.mnn$x, subset=pca.retain, main="Fast MNN", xlab="PC 1",ylab="PC 2")
plotFUNb("results/pca_mnn2_batch.png", pca.mnn$x, subset=pca.retain, main="Fast MNN", xlab="PC 1",ylab="PC 2")

rm(t.mnn)
gc()

######################## 
# Performing the correction with limma.

library(limma)
X.lm <- removeBatchEffect(as.matrix(logcounts(sce)), factor(first.batch))
t.lm <- t(X.lm)

## Generating a t-SNE plot.
set.seed(0)
tsne.lm <- Rtsne(t.lm, perplexity = 90)
plotFUN("results/tsne_limma_type.png", tsne.lm$Y, main="limma", xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/tsne_limma_batch.png", tsne.lm$Y, main="limma", xlab="tSNE 1",ylab="tSNE 2")

gc()

# Generating a PCA plot.
pca.lm <- prcomp(t.lm[pca.retain,], rank=2)
plotFUN("results/pca_limma_type.png", pca.lm$x, subset=pca.retain, main="limma", xlab="PC 1",ylab="PC 2")
plotFUNb("results/pca_limma_batch.png", pca.lm$x, subset=pca.retain, main="limma", xlab="PC 1",ylab="PC 2")

rm(t.lm)
gc()

######################## 
# Performing the correction with ComBat.

library(sva)
Z <- factor(first.batch)
X.combat <- ComBat(as.matrix(logcounts(sce)), Z, mod=NULL,prior.plots = FALSE)
t.combat <- t(X.combat)

## Generating a t-SNE plot.
set.seed(0)
tsne.combat <- Rtsne(t.combat, perplexity = 90)
plotFUN("results/tsne_combat_type.png", tsne.combat$Y, main="ComBat", xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/tsne_combat_batch.png", tsne.combat$Y, main="ComBat", xlab="tSNE 1",ylab="tSNE 2")

gc()

# Generating a PCA plot.
pca.combat <- prcomp(t.combat[pca.retain,], rank=2)
plotFUN("results/pca_combat_type.png", pca.combat$x, subset=pca.retain, main="ComBat", xlab="PC 1",ylab="PC 2")
plotFUNb("results/pca_combat_batch.png", pca.combat$x, subset=pca.retain, main="ComBat", xlab="PC 1",ylab="PC 2")

rm(t.combat)
gc()

######################## 
# Performing the correction with CCA.

library(Seurat)
SeuA <- CreateSeuratObject(logcounts(sceA))
SeuF <- CreateSeuratObject(logcounts(sceF))
SeuA@meta.data$group <- "group1"
SeuF@meta.data$group <- "group2"

SeuA <- ScaleData(SeuA)
SeuF <- ScaleData(SeuF)
Y <- RunCCA(SeuA, SeuF, genes.use=rownames(sceA), do.normalize=FALSE)
suppressWarnings(Y <- AlignSubspace(Y, grouping.var="group", dims.align=1:20))

t.cca <- Y@dr$cca.aligned@cell.embeddings

## Generating a t-SNE plot.
set.seed(0)
tsne.cca <- Rtsne(t.cca, perplexity = 90)
plotFUN("results/tsne_cca_type.png", tsne.cca$Y, main="CCA", xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/tsne_cca_batch.png", tsne.cca$Y, main="CCA", xlab="tSNE 1",ylab="tSNE 2")

gc()

# Generating a PCA plot.
pca.cca <- prcomp(t.cca[pca.retain,], rank=2)
plotFUN("results/pca_cca_type.png", pca.cca$x, subset=pca.retain, main="CCA", xlab="PC 1",ylab="PC 2")
plotFUNb("results/pca_cca_batch.png", pca.cca$x, subset=pca.retain, main="CCA", xlab="PC 1",ylab="PC 2")

rm(t.cca)
gc()

######################## 
# Performing the correction with CCA, starting from the raw counts and normalizing it their way.

library(Seurat)
SeuA <- CreateSeuratObject(counts(sceA))
SeuF <- CreateSeuratObject(counts(sceF))
SeuA@meta.data$group <- "group1"
SeuF@meta.data$group <- "group2"

SeuA <- NormalizeData(SeuA)
SeuF <- NormalizeData(SeuF)
SeuA <- ScaleData(SeuA)
SeuF <- ScaleData(SeuF)
Y <- RunCCA(SeuA, SeuF, genes.use=rownames(sceA))
suppressWarnings(Y <- AlignSubspace(Y, grouping.var="group", dims.align=1:20))

t.cca <- Y@dr$cca.aligned@cell.embeddings

## Generating a t-SNE plot.
set.seed(0)
tsne.cca2 <- Rtsne(t.cca, perplexity = 90)
plotFUN("results/tsne_cca2_type.png", tsne.cca2$Y, main="CCA native", xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/tsne_cca2_batch.png", tsne.cca2$Y, main="CCA native", xlab="tSNE 1",ylab="tSNE 2")
gc()

# Generating a PCA plot.
pca.cca2 <- prcomp(t.cca[pca.retain,], rank=2)
plotFUN("results/pca_cca2_type.png", pca.cca2$x, subset=pca.retain, main="CCA native", xlab="PC 1",ylab="PC 2")
plotFUNb("results/pca_cca2_batch.png", pca.cca2$x, subset=pca.retain, main="CCA native", xlab="PC 1",ylab="PC 2")

rm(t.cca)
gc()


######################## 
# Making the legend (using PDF for better resolution).

pdf(file="results/legend.pdf", width=7, height=7)
plot(0,0,type="n", bty="n", axes=FALSE, xlab="", ylab="")
legend(x=-1, y=1, legend=names(color.legendF), pch=21, cex=2.5, col="black", pt.bg=color.legendF, bty="n")
legend(x=0, y=1, legend=names(color.legendA)[1:3], pch=21, cex=2.5, col="black", pt.bg=color.legendA[1:3], bty="n")
dev.off()

######################## 
# END
