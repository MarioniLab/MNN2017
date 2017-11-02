# Batch correction and t-SNE plots of haematopoietic data sets, as in Figure 3.
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(scran)
dir.create("results", showWarning=FALSE)

# Load the output of "preparedata.R".
load("logdataFandA_all.RData") 
colnames(logDataF3)[colnames(logDataF3)=="other"] <- "Unsorted"
colnames(logDataA3)[colnames(logDataA3)=="ERY"] <- "MEP"
raw.all <- cbind(logDataF3, logDataA3)
first.batch <- rep(c(TRUE, FALSE), c(ncol(logDataF3), ncol(logDataA3)))

# Adding colours.
#base.color <- "grey"
color.legendF <- c(MEP="gold", GMP="chartreuse4", CMP="magenta", 
                   HSPC="cyan", LTHSC="dodgerblue", MPP="blue", LMPP="light blue", Unsorted="grey")
colmatF <- col2rgb(color.legendF) 

colmatA <- colmatF + 100 # A lighter shade.
colmatA[colmatA > 255] <- 255

#colmatF<-colmatF+200
#colmatF[colmatF > 255] <- 255
color.legendA <- setNames(rgb(colmatA[1,], colmatA[2,], colmatA[3,], maxColorValue=255), names(color.legendF))
allcolors <- c(color.legendF[colnames(logDataF3)], color.legendA[colnames(logDataA3)])
batch<-c( rep(1,ncol(logDataF3)),rep(2,ncol(logDataA3)) )

# Only keeping common cell types for PCA.
celltypes <- c(colnames(logDataF3), colnames(logDataA3))
pca.retain <- celltypes %in% c("MEP", "GMP", "CMP") 

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

batchcolor=c("lavender","lightcoral")
plotFUNb <- function(fname, Y, subset=NULL, ...) {
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  png(fname,width=900,height=700)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=ifelse(first.batch, 21, 1)[subset], 
       col=ifelse(first.batch, "black", batchcolor[batch[subset]]),
       bg=batchcolor[batch[subset]], ...)#,  xlab="tSNE 1",ylab="tSNE 2")
  dev.off()
}

######################
# No correction.

X.unc <- raw.all
t.unc <- t(X.unc)

## Generating a t-SNE plot.
require(Rtsne)
set.seed(0)
all.dists.unc <- as.matrix(dist(t.unc))
tsne.unc <- Rtsne(all.dists.unc, is_distance=TRUE, perplexity = 90)
plotFUN("results/uncFA.png", tsne.unc$Y, main="Uncorrected",  xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/uncFAb.png", tsne.unc$Y, main="Uncorrected",  xlab="tSNE 1",ylab="tSNE 2")

rm(all.dists.unc)
gc()

# Generating a PCA plot.
pca.unc <- prcomp(t.unc[pca.retain,], rank=2)
pca.unc$x[ (pca.unc$x<(-0.08))]<- (-0.08)
plotFUN("results/pca_raw.png", pca.unc$x, subset=pca.retain, main="Uncorrected", ylim=c(-0.1, max(pca.unc$x[,2])),  xlab="PC 1",ylab="PC 2")

rm(t.unc)
gc()
######################## 
# Performing the correction with MNN (turned down the sigma to improve mixing).

mnn.out <- mnnCorrect(logDataF3, logDataA3, k=20, sigma=0.1, cos.norm=TRUE, svd.dim=0)
X.mnn <- cbind(mnn.out$corrected[[1]], mnn.out$corrected[[2]])
t.mnn <- t(X.mnn)

# plot histogram of angles between batch vectors and 2 svds of the reference batch.
png(file="results/angles.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
hist(unlist(mnn.out$ang.with.ref),xlab="Angle",ylab="Frequency",main="")
dev.off()
# Generating a t-SNE plot.
set.seed(0)
all.dists.mnn <- as.matrix(dist(t.mnn))
tsne.mnn <- Rtsne(all.dists.mnn, is_distance=TRUE, perplexity = 90)
plotFUN("results/mnnFA.png", tsne.mnn$Y, main="MNN",  xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/mnnFAb3.png", tsne.mnn$Y, main="MNN",  xlab="tSNE 1",ylab="tSNE 2")

set.seed(0)
tsne.mnn2 <- Rtsne(t.mnn, perplexity = 90)
plotFUN("results/mnnFA_conventsne.png", tsne.mnn2$Y, main="MNN",  xlab="tSNE 1",ylab="tSNE 2")

rm(all.dists.mnn)
gc()

# Generating a PCA plot.
pca.mnn <- prcomp(t.mnn[pca.retain,], rank=2)
pca.mnn$x[ (pca.mnn$x<(-0.08))]<- (-0.08)
plotFUN("results/pca_mnn.png", pca.mnn$x, subset=pca.retain, main="MNN", ylim=c(-0.08,0.05),  xlab="PC 1",ylab="PC 2")

# Generating diffusion map plots.
library(destiny)
dm<-DiffusionMap(t.mnn,n_local = 150)
plotFUN("results/mnnFAdm12.png", dm@eigenvectors[,1:2], main="MNN",  xlab="DC 1",ylab="DC 2")
plotFUN("results/mnnFAdm24.png", cbind(dm@eigenvectors[,2],dm@eigenvectors[,4]), main="MNN",  xlab="DC 2",ylab="DC 4")

rm(t.mnn)
gc()

######################## 
# Performing the correction with limma.

library(limma)
X.lm <- removeBatchEffect(raw.all, factor(first.batch))
t.lm <- t(X.lm)

## Generating a t-SNE plot.
set.seed(0)
all.dists.lm <- as.matrix(dist(t.lm))
tsne.lm <- Rtsne(all.dists.lm, is_distance=TRUE, perplexity = 90)
plotFUN("results/lmfitFA.png", tsne.lm$Y, main="limma",  xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/lmfitFAb.png", tsne.lm$Y, main="limma",  xlab="tSNE 1",ylab="tSNE 2")

rm(all.dists.lm)
gc()

# Generating a PCA plot.
pca.lm <- prcomp(t.lm[pca.retain,], rank=2)
pca.lm$x[ (pca.lm$x[,2] > (0.08)),2]<- (0.08)

plotFUN("results/pca_lm.png", pca.lm$x, subset=pca.retain, main="limma", ylim=c(min(pca.lm$x[,2]), 0.05),  xlab="PC 1",ylab="PC 2")

rm(t.lm)
gc()

######################## 
# Performing the correction with ComBat.

library(sva)
Z <- factor(first.batch)
X.combat <- ComBat(raw.all,Z,mod=NULL,prior.plots = FALSE)
t.combat <- t(X.combat)

## Generating a t-SNE plot.
set.seed(0)
all.dists.combat <- as.matrix(dist(t.combat))
tsne.combat<-Rtsne(all.dists.combat, is_distance=TRUE,perplexity = 90)
plotFUN("results/combatFA.png", tsne.combat$Y, main="ComBat",  xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/combatFAb.png", tsne.combat$Y, main="ComBat",  xlab="tSNE 1",ylab="tSNE 2")

rm(all.dists.combat)
gc()

# Generating a PCA plot.
pca.combat <- prcomp(t.combat[pca.retain,], rank=2)
plotFUN("results/pca_com.png", pca.combat$x, subset=pca.retain, main="ComBat", ylim=c(min(pca.combat$x[,2]), 0.05),  xlab="PC 1",ylab="PC 2")

rm(t.combat)
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
