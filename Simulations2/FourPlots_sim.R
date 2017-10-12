# This script compares the performance of various batch correction methods.

library(scran)
require(Rtsne)
library(limma)
library(sva)
library(scales)

load("Sim.RData")
dir.create("figs", showWarnings=FALSE)

plotFUN <- function(fname, Y, batch.id, cols, xlab="tSNE 1", ylab="tSNE 2", ...) {
    png(fname, width=900, height=700)
    par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
    plot(Y[,1], Y[,2], 
         pch=c(16, 2)[batch.id],
         cex=c(2.5, 3.5)[batch.id],
         col=alpha(cols, 0.6),
         xlab=xlab, ylab=ylab, ...)
    dev.off()
}

##########################################
# Running through all of the methods.

for (easy in c(FALSE, TRUE)) {
    if (easy) {
        B2 <- B2ii
        clust2 <- clust2ii
        prefix <- "easy_"
    } else {
        B2 <- B2i
        clust2 <- clust2i
        prefix <- ""
    }

    # Uncorrected.
    raw.all <- cbind(B1, B2)
    clust.cols <- c(clust1, clust2)
    batch.id <- rep(1:2, c(ncol(B1), ncol(B2)))

    all.dists2.unc <- as.matrix(dist(t(raw.all)))
    set.seed(0)
    tsne.unc <- Rtsne(all.dists2.unc, is_distance=TRUE)#, perplexity = 0.9)
    plotFUN(paste0("figs/", prefix, "unc.png"), Y=tsne.unc$Y, batch.id=batch.id, cols=clust.cols, main="Uncorrected",xlim=c(-35,30),ylim=c(-40,45))

    # MNN corrected.
    Xmnn <- mnnCorrect(B1, B2, k=20, sigma=1, cos.norm=FALSE, svd.dim=2)
    corre <- cbind(Xmnn$corrected[[1]],Xmnn$corrected[[2]])
    all.dists2.c <- as.matrix(dist(t(corre)))

    set.seed(0)
    tsne.c <- Rtsne(all.dists2.c, is_distance=TRUE)#, perplexity = 0.9)
    plotFUN(paste0("figs/", prefix, "mnn.png"), Y=tsne.c$Y, batch.id=batch.id, col=clust.cols, main="MNN corrected",xlim=c(-30,40),ylim=c(-45,45))

    # limma.
    Xlm <- removeBatchEffect(raw.all, factor(batch.id))
    all.dists2.lm <- as.matrix(dist(t(Xlm)))

    set.seed(0)
    tsne.lm <- Rtsne(all.dists2.lm, is_distance=TRUE)#, perplexity = 0.9)
    plotFUN(paste0("figs/", prefix, "lmfit.png"), Y=tsne.lm$Y, batch.id=batch.id, col=clust.cols, main="limma corrected",xlim=c(-30,30),ylim=c(-35,45))

    # ComBat.
    cleandat <- ComBat(raw.all, factor(batch.id), mod=NULL, prior.plots = FALSE)
    all.dists.combat <- as.matrix(dist(t(cleandat)))
    
    set.seed(0)
    tsne.combat<-Rtsne(all.dists.combat, is_distance=TRUE)
    plotFUN(paste0("figs/", prefix, "combat.png"), Y=tsne.combat$Y, batch.id=batch.id, col=clust.cols, main="ComBat corrected",xlim=c(-35,30),ylim=c(-35,30))
}

##########################################
# Printing out the legend 

png(file="leg.png",width=900,height=700)
plot.new()
legend("topleft", legend = c("Cell type 1", "Cell type 2", "Cell type 3", "Batch 1", "Batch 2"),
       col = c("brown1", "dark green", "blue", "black", "black"), 
       pch = c(15, 15, 15, 16, 2),
       cex = 2.7,bty = "n")   
dev.off()
