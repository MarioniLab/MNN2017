# Batch correction and t-SNE plots of haematopoietic data sets, as in Figure 3.

library(scran)
dir.create("results", showWarning=FALSE)

# Load the output of "preparedata.R".
load("logdataFandA_all.RData") 
raw.all <- cbind(logDataF3, logDataA3)
first.batch <- rep(c(TRUE, FALSE), c(ncol(logDataF3), ncol(logDataA3)))

# Organize cell type labels.
celltypes <- c(colnames(logDataF3), colnames(logDataA3))
celltypes[celltypes=="other"]<- "Unsorted"
celltypes[celltypes=="ERY"]<- "MEP"

# Adding colours.
base.color <- "grey"
color.legend <- c(MEP="red",
                  GMP="blue",
                  CMP="goldenrod", 
                  HSPC=base.color, LTHSC=base.color, MPP=base.color, LMPP=base.color, Unsorted=base.color)
colmat <- col2rgb(color.legend[celltypes])
colmat[,!first.batch] <- colmat[,!first.batch]+200 # A lighter shade.
colmat[colmat > 255] <- 255
allcolors <- rgb(colmat[1,], colmat[2,], colmat[3,], maxColorValue=255)

# Making a plotting function.
plotFUN <- function(fname, Y, subset=NULL, ..., xlab="tSNE 1",ylab="tSNE 2") {
    # Setting up the plotting device.
    png(fname,width=900,height=700)
    par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)

    # Plot the first, smaller batch ON TOP OF the second, larger batch.
    reorder <- rev(seq_len(nrow(Y)))
    plot(Y[reorder,1], Y[reorder,2], cex=2, pch=21, col="black", 
         bg=allcolors[reorder], ..., xlab=xlab, ylab=ylab) 
    dev.off()
}

######################
# No correction.

X.unc <- raw.all
t.unc <- t(X.unc)

## Generating a t-SNE plot.
#require(Rtsne)
#set.seed(0)
#all.dist.unc <- dist(t.unc)
#tsne.unc <- Rtsne(all.dist.unc, is_distance=TRUE, perplexity = 30)
#plotFUN("results/uncFA.png", tsne.unc$Y, main="Uncorrected")
#plotFUN("results/uncFAb.png", tsne.unc$Y, main="Uncorrected", by.batch=TRUE)
#
#gc()

# Generating a PCA plot.
pca.unc <- prcomp(t.unc, rank=2)
plotFUN("results/pca_raw.png", pca.unc$x, main="Uncorrected", xlab="PC1", ylab="PC2")

rm(t.unc)
gc()

######################## 
# Performing the correction with MNN (turned down the sigma to improve mixing).

mnn.out <- mnnCorrect(logDataF3, logDataA3, k=20, sigma=0.1, cos.norm=TRUE, svd.dim=NA)
X.mnn <- cbind(mnn.out$corrected[[1]], mnn.out$corrected[[2]])
t.mnn <- t(X.mnn)

## Generating a t-SNE plot.
#set.seed(0)
#all.dist.mnn <- as.matrix(dist(t.mnn))
#tsne.mnn <- Rtsne(all.dist.mnn, is_distance=TRUE, perplexity = 30)
#plotFUN("results/mnnFA.png", tsne.mnn$Y, main="MNN")
#plotFUN("results/mnnFAb.png", tsne.mnn$Y, main="MNN", by.batch=TRUE)
#
#gc()

# Generating a PCA plot.
pca.mnn <- prcomp(t.mnn, rank=2)
plotFUN("results/pca_mnn.png", pca.mnn$x, main="MNN", xlab="PC1", ylab="PC2")

rm(t.mnn)
gc()

######################## 
# Performing the correction with limma.

library(limma)
X.lm <- removeBatchEffect(raw.all, factor(first.batch))
t.lm <- t(X.lm)

## Generating a t-SNE plot.
#set.seed(0)
#all.dist.lm <- as.matrix(dist(t.lm))
#tsne.lm <- Rtsne(all.dist.lm, is_distance=TRUE, perplexity = 30)
#plotFUN("results/lmfitFA.png", tsne.lm$Y, main="limma")
#plotFUN("results/lmfitFAb.png", tsne.lm$Y, main="limma", by.batch=TRUE)
#
#gc()

# Generating a PCA plot.
pca.lm <- prcomp(t.lm, rank=2)
plotFUN("results/pca_lm.png", pca.lm$x, main="limma", xlab="PC1", ylab="PC2")

rm(t.lm)
gc()

######################## 
# Performing the correction with ComBat.

library(sva)
Z <- factor(first.batch)
X.combat <- ComBat(raw.all,Z,mod=NULL,prior.plots = FALSE)
t.combat <- t(X.combat)

## Generating a t-SNE plot.
#set.seed(0)
#all.dists.combat <- as.matrix(dist(t.combat))
#tsne.combat <- Rtsne(all.dists.combat, is_distance=TRUE, perplexity=30)
#plotFUN("results/combatFA.png", tsne.combat$Y, main="ComBat")
#plotFUN("results/combatFAb.png", tsne.combat$Y, main="ComBat", by.batch=TRUE)
#
#gc()

# Generating a PCA plot.
pca.combat <- prcomp(t.combat, rank=2)
plotFUN("results/pca_com.png", pca.combat$x, main="ComBat", xlab="PC1", ylab="PC2")

rm(t.combat)
gc()

######################## 
# Making the legend (using PDF for better resolution).

pdf(file="results/legFA.pdf", width=7, height=7)
plot(0,0,type="n", bty="n", axes=FALSE, xlab="", ylab="")
legend(x=-1, y=1, legend=names(color.legend), pch=21, cex=2.5, col="black", pt.bg=color.legend, bty="n")
legend(x=0, y=1, legend=names(color.legend)[1:3], col = color.legend[1:3], pch = 1, cex = 2.5, bty = "n", lwd = 3, lty=0) 
dev.off()

######################## 
# END
