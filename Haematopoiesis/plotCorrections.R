# Batch correction and t-SNE plots of haematopoietic data sets, as in Figure 3.

library(scran)
dir.create("results", showWarning=FALSE)

# Load the output of "preparedata.R".
load("logdataFandA_all.RData") 
colnames(logDataF3)[colnames(logDataF3)=="other"] <- "Unsorted"
colnames(logDataA3)[colnames(logDataA3)=="ERY"] <- "MEP"
raw.all <- cbind(logDataF3, logDataA3)
first.batch <- rep(c(TRUE, FALSE), c(ncol(logDataF3), ncol(logDataA3)))

# Adding colours.
base.color <- "grey"
color.legendF <- c(MEP="red", GMP="blue", CMP="orange", 
                   HSPC=base.color, LTHSC=base.color, MPP=base.color, LMPP=base.color, Unsorted=base.color)
colmatF <- col2rgb(color.legendF)
colmatA <- colmatF + 200 # A lighter shade.
colmatA[colmatA > 255] <- 255
color.legendA <- setNames(rgb(colmatA[1,], colmatA[2,], colmatA[3,], maxColorValue=255), names(color.legendF))
allcolors <- c(color.legendF[colnames(logDataF3)], color.legendA[colnames(logDataA3)])

# Making a plotting function.
plotFUN <- function(fname, Y, subset=NULL, ..., xlab="tSNE 1",ylab="tSNE 2") {
    # Setting up the plotting device.
    pdf(fname,width=9,height=7)
    par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=1.2,cex.main=1.4,cex.lab=1.4)

    # Plot the first, smaller batch ON TOP OF the second, larger batch.
    set.seed(0)
    reorder <- sample(nrow(Y))
    plot(Y[reorder,1], Y[reorder,2], cex=2, pch=21, col="black", lwd=0.1,
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
#plotFUN("results/uncFA.pdf", tsne.unc$Y, main="Uncorrected")
#plotFUN("results/uncFAb.pdf", tsne.unc$Y, main="Uncorrected", by.batch=TRUE)
#
#gc()

# Generating a PCA plot.
pca.unc <- prcomp(t.unc, rank=2)
plotFUN("results/pca_raw.pdf", pca.unc$x, main="Uncorrected", xlab="PC1", ylab="PC2")

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
#plotFUN("results/mnnFA.pdf", tsne.mnn$Y, main="MNN")
#plotFUN("results/mnnFAb.pdf", tsne.mnn$Y, main="MNN", by.batch=TRUE)
#
#gc()

# Generating a PCA plot.
pca.mnn <- prcomp(t.mnn, rank=2)
plotFUN("results/pca_mnn.pdf", pca.mnn$x, main="MNN", xlab="PC1", ylab="PC2")

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
#plotFUN("results/lmfitFA.pdf", tsne.lm$Y, main="limma")
#plotFUN("results/lmfitFAb.pdf", tsne.lm$Y, main="limma", by.batch=TRUE)
#
#gc()

# Generating a PCA plot.
pca.lm <- prcomp(t.lm, rank=2)
plotFUN("results/pca_lm.pdf", pca.lm$x, main="limma", xlab="PC1", ylab="PC2")

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
#plotFUN("results/combatFA.pdf", tsne.combat$Y, main="ComBat")
#plotFUN("results/combatFAb.pdf", tsne.combat$Y, main="ComBat", by.batch=TRUE)
#
#gc()

# Generating a PCA plot.
pca.combat <- prcomp(t.combat, rank=2)
plotFUN("results/pca_com.pdf", pca.combat$x, main="ComBat", xlab="PC1", ylab="PC2")

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
