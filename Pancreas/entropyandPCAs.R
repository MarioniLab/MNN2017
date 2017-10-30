# Code for PCA plots and entropy of batch mixing for pancreas data sets (Suppl. Fig. 3)
# Please set the working directory to the location of this repository, e.g. ~/MNN2017/

load("Pancreas/Data/ObjectsForPlotting.RDS")  #load the output of "FourPlots_panc.R"
library(scales)
library(RANN)

batch <- colnames(corrected.df)[1:(dim(corrected.df)[2]-1)]

# set plotting colors
library(ggplot2)
library(cowplot)
library(RColorBrewer)
colors4 <- brewer.pal(4, "Set3")
names(colors4) <- unique(all.meta$Study)

# set plotting symbols
forpch <- numeric(length(batch))
forpch[batch %in% all.meta$Sample[all.meta$Study == "GSE81076"]] <- 3
forpch[batch %in% all.meta$Sample[all.meta$Study == "GSE85241"]] <- 18
forpch[batch %in% all.meta$Sample[all.meta$Study == "GSE86473"]] <- 1
forpch[batch %in% all.meta$Sample[all.meta$Study == "E-MTAB-5061"]] <- 4

## PCA of uncorrected data
data <- raw.all[common.hvgs, ]
pca.raw <- prcomp(t(data), center=TRUE)
pca.raw$x <- scran:::cosine.norm(pca.raw$x)  #for orthonormal pca basis. Basis normalization is not done in prcomp
pcs.raw <- as.data.frame(pca.raw$x)
pcs.raw$Sample <- colnames(data)
raw.merge <- merge(pcs.raw, all.meta, by='Sample')

# randomize order of cells for plotting
set.seed(2)
ix <- sample(1:nrow(pcs.raw), nrow(pcs.raw)) 

# plot of points by cell type and shape by study
ggplot(raw.merge, aes(x=V1, y=V2, 
                      fill=Study,
                      group=Study)) +
  geom_point(size=1.5, shape=21) + theme_classic() +
  scale_fill_manual(values=colors4) +
  scale_y_continuous(limits=c(-0.04, 0.04)) +
  scale_x_continuous(limits=c(-0.04, 0.04)) +
  labs(x="PC 1", y="PC 2", title="Uncorrected")

# png(file="pca_raw.png", width=900, height=700)
# par(mfrow=c(1,1), mar=c(6,6,4,2), cex.axis=2, cex.main=3, cex.lab=2.5, pch=16)
# plot(pca.raw$x[ix,1], pca.raw$x[ix,2],# col=alpha(colors4[batch[ix]], 0.6),
#      pch=16, main="Uncorrected", xlab="PC 1", ylab="PC 2", cex=1.5)
# dev.off()

###### PCA of MNN corrected
pca.mnn <- prcomp(t(corrected.df[, 1:(dim(corrected.df)[2]-1)]), center=TRUE)
pca.mnn$x <- scran:::cosine.norm(pca.mnn$x) 
pcs.mnn <- as.data.frame(pca.mnn$x)
pcs.mnn$Sample <- colnames(corrected.df[, 1:(dim(corrected.df)[2]-1)])
mnn.merge <- merge(pcs.mnn, all.meta, by='Sample')

# plot of points by cell type and shape by study
ggplot(mnn.merge, aes(x=V1, y=V2, 
                      fill=Study,
                      group=Study)) +
  geom_point(size=1.5, shape=21) + theme_classic() +
  scale_fill_manual(values=colors4) +
  scale_y_continuous(limits=c(-0.04, 0.04)) +
  scale_x_continuous(limits=c(-0.04, 0.04)) +
  labs(x="PC 1", y="PC 2", title="MNN corrected")

# png(file="pca_mnn.png", width=900, height=700)
# par(mfrow=c(1,1), mar=c(6,6,4,2), cex.axis=2, cex.main=3, cex.lab=2.5, pch=16)
# plot(pca.mnn$x[ix,1], pca.mnn$x[ix,2], col=alpha(colors4[batch[ix]], 0.6),
#      main="MNN corrected", xlab="PC 1", ylab="PC 2", cex=1.5)
# dev.off()

###### PCA of limma corrected
pca.lm <- prcomp(lm.mat, center=TRUE)
pca.lm$x <- scran:::cosine.norm(pca.lm$x) 
pcs.lm <- as.data.frame(pca.lm$x)
pcs.lm$Sample <- rownames(lm.mat)
lm.merge <- merge(pcs.lm, all.meta, by='Sample')

# plot of points by batch
ggplot(lm.merge, aes(x=V1, y=V2, 
                      fill=Study,
                      group=Study)) +
  geom_point(size=1.5, shape=21) + theme_classic() +
  scale_fill_manual(values=colors4) +
  scale_y_continuous(limits=c(-0.04, 0.04)) +
  scale_x_continuous(limits=c(-0.04, 0.04)) +
  labs(x="PC 1", y="PC 2", title="limma corrected")


# png(file="pca_lm.png", width=900, height=700)
# par(mfrow=c(1,1), mar=c(6,6,4,2), cex.axis=2, cex.main=3, cex.lab=2.5, pch=16)
# plot(pca.lm$x[ix,1], pca.lm$x[ix,2], col=alpha(colors4[batch[ix]], 0.6),
#      main="limma corrected", xlab="PC 1", ylab="PC 2", cex=1.5)
# dev.off()

######PCA of ComBat corrected
pca.com <- prcomp(combat.mat, center=TRUE)
pca.com$x <- scran:::cosine.norm(pca.com$x)
pcs.com <- as.data.frame(pca.com$x)
pcs.com$Sample <- rownames(combat.mat)
com.merge <- merge(pcs.com, all.meta, by='Sample')

# plot of points by batch
ggplot(com.merge, aes(x=V1, y=V2, 
                     fill=Study,
                     group=Study)) +
  geom_point(size=1.5, shape=21) + theme_classic() +
  scale_fill_manual(values=colors4) +
  scale_y_continuous(limits=c(-0.04, 0.04)) +
  scale_x_continuous(limits=c(-0.04, 0.04)) +
  labs(x="PC 1", y="PC 2", title="ComBat corrected")


png(file="pca_com.png", width=900, height=700)
par(mfrow=c(1,1), mar=c(6,6,4,2), cex.axis=2, cex.main=3, cex.lab=2.5, pch=16)
plot(pca.com$x[ix,1], pca.com$x[ix,2], col=alpha(colors4[batch[ix]], 0.6),
     main="ComBat corrected", xlab="PC 1", ylab="PC 2", cex=1.5)
dev.off()

# ########### the legend
# png(file="pancpca_leg.png", width=900, height=700)
# par(mfrow=c(1,1), mar=c(6,6,4,2), cex.axis=2, cex.main=3, cex.lab=2.5, pch=16)
# plot(1 , 2, col=alpha(colors4[batch], 0.6),
#      main="limma corrected", xlab="PC 1", ylab="PC 2", cex=1.5,
#      xlim=c(-0.025,0.04))
# legend("topleft", legend=c("CEL-Seq", "CEL-Seq2", "SMART-Seq2 (I)", "SMART-Seq (II)"),
#        col=unique(batch), pch=16, cex=2.5, bty="n")
# dev.off()
###################Entropy of batch mixing for uncorrected and batch corrected data using different methods

source("SomeFuncs/BatchMixingEntropy.R")
entrop.raw <- BatchEntropy(pca.raw$x[,1:2], batch)
entrop.mnn <- BatchEntropy(pca.mnn$x[,1:2], batch)
entrop.lm <- BatchEntropy(pca.lm$x[,1:2], batch)
entrop.com <- BatchEntropy(pca.com$x[,1:2], batch)

En <- cbind(entrop.raw, entrop.mnn, entrop.lm, entrop.com)

# box plot of entropy of batch mixing 
png(file="pcabatch_entropy.png", width=900, height=700)
par(mfrow=c(1,1), mar=c(8,8,5,3), cex.axis=3, cex.main=2, cex.lab=3)
boxplot(En,main="", names=c("Raw","MNN","limma","ComBat"),
        lwd=4, ylab="Entropy of batch mixing")
dev.off()

