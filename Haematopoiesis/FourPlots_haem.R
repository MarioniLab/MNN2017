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
color.legend <- c(MEP="red", GMP="dark green", CMP="magenta", HSPC="cyan", LTHSC="dodgerblue", MPP="blue", LMP="light blue", Unsorted="grey")
allcolors <- color.legend[celltypes]

plotFUN <- function(fname, Y, subset=NULL, ...) {
    if (is.null(subset)) {
        subset <- seq_len(nrow(Y))
    }
    png(fname,width=900,height=700)
    par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
    plot(Y[,1], Y[,2], cex=2,
         pch=ifelse(first.batch, 21, 1)[subset], 
         col=ifelse(first.batch, "black", allcolors)[subset],
         bg=allcolors[subset], ...,  xlab="tSNE 1",ylab="tSNE 2")
    dev.off()
}

# Only keeping common cell types for PCA.
pca.retain <- celltypes %in% c("MEP", "GMP", "CMP") 

######################
# No correction.

X.unc <- raw.all
t.unc <- t(X.unc)

# Generating a t-SNE plot.
require(Rtsne)
set.seed(0)
all.dists.unc <- as.matrix(dist(t.unc))
tsne.unc <- Rtsne(all.dists.unc, is_distance=TRUE)#, perplexity = 90)
plotFUN("results/uncFA.png", tsne.unc$Y, main="Uncorrected")

rm(all.dists.unc)
gc()

# Generating a PCA plot.
pca.unc <- prcomp(t.unc[pca.retain,], rank=2)
plotFUN("results/pca_raw.png", pca.unc$x, subset=pca.retain, main="Uncorrected", ylim=c(-0.1, max(pca.unc$x[,2])))

rm(t.unc)
gc()

######################## 
# Performing the correction with MNN (turned down the sigma to improve mixing).

mnn.out <- mnnCorrect(logDataF3, logDataA3, k=20, sigma=0.01, cos.norm=TRUE, svd.dim=2)
X.mnn <- cbind(mnn.out$corrected[[1]], mnn.out$corrected[[2]])
t.mnn <- t(X.mnn)

# Generating a t-SNE plot.
set.seed(0)
all.dists.mnn <- as.matrix(dist(t.mnn))
tsne.mnn <- Rtsne(all.dists.mnn, is_distance=TRUE)#, perplexity = 90)
plotFUN("results/mnnFA.png", tsne.mnn$Y, main="MNN")

rm(all.dists.mnn)
gc()

# Generating a PCA plot.
pca.mnn <- prcomp(t.mnn[pca.retain,], rank=2)
plotFUN("results/pca_mnn.png", pca.mnn$x, subset=pca.retain, main="MNN")

rm(t.mnn)
gc()

######################## 
# Performing the correction with limma.

library(limma)
X.lm <- removeBatchEffect(raw.all, factor(first.batch))
t.lm <- t(X.lm)

# Generating a t-SNE plot.
set.seed(0)
all.dists.lm <- as.matrix(dist(t.lm))
tsne.lm <- Rtsne(all.dists.lm, is_distance=TRUE)#, perplexity = 90)
plotFUN("results/lmfitFA.png", tsne.lm$Y, main="limma")

rm(all.dists.lm)
gc()

# Generating a PCA plot.
pca.lm <- prcomp(t.lm[pca.retain,], rank=2)
plotFUN("results/pca_lm.png", pca.lm$x, subset=pca.retain, main="limma", ylim=c(min(pca.lm$x[,2]), 0.05))

rm(t.lm)
gc()

######################## 
# Performing the correction with ComBat.

library(sva)
Z <- factor(first.batch)
X.combat <- ComBat(raw.all,Z,mod=NULL,prior.plots = FALSE)
t.combat <- t(X.combat)

# Generating a t-SNE plot.
set.seed(0)
all.dists.combat <- as.matrix(dist(t.combat))
tsne.combat<-Rtsne(all.dists.combat, is_distance=TRUE)#,perplexity = 90)
plotFUN("results/combatFA.png", tsne.combat$Y, main="ComBat")

rm(all.dists.combat)
gc()

# Generating a PCA plot.
pca.combat <- prcomp(t.combat[pca.retain,], rank=2)
plotFUN("results/pca_com.png", pca.combat$x, subset=pca.retain, main="ComBat", ylim=c(min(pca.combat$x[,2]), 0.05))

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
