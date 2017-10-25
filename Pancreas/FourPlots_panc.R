# Code for the t-SNE plots of pancreas data sets and the Silhouette coefficients before and after batch correction by different methods (main text Figure 4).
# PLEASE SET THE WORKING DIRECTORY TO ./MNN2017/

# setwd("results")
require(WGCNA)
require(Rtsne)
library(scales)
library(scran)
library(ggplot2)
library(cowplot)
############
# GSE81076 #
############
# do a bit of tidying to make sure expression matrices are all conformable and
# only contain data with cell type assignments

# explicitly read in the study data
datah1 <- read.table("Pancreas/Data/GSE81076_SFnorm.tsv",
                     h=TRUE, sep="\t", stringsAsFactors=FALSE)

# gene IDs are in the last column
genes1 <- datah1$gene_id
rownames(datah1) <- genes1
datah1 <- datah1[, 1:(dim(datah1)[2]-1)]

hvg1 <- read.table("Pancreas/Data/GSE81076-HVG.tsv",
                   h=TRUE, sep="\t", stringsAsFactors=F)
HVG1 <- hvg1$gene_id

meta1 <- read.table("Pancreas/Data/GSE81076_marker_metadata.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)

# select meta data for which there is expression data
meta1 <- meta1[meta1$Sample %in% colnames(datah1), ]

# standardize the cell labels
## NB there are 112 samples with unassigned cell types, remove these
celltypes1 <- meta1$CellType
no.label1 <- meta1$Sample[meta1$CellType == ""]
datah1 <- datah1[, !colnames(datah1) %in% no.label1]
meta1 <- meta1[!meta1$Sample %in% no.label1, ]
celltypes1 <- celltypes1[celltypes1 != ""]
samples1 <- colnames(datah1)

# check all dimensions match up
if(dim(datah1)[2] == dim(meta1)[1]) {dim(datah1)[2] == length(celltypes1)}

############
# GSE85241 #
############
# do a bit of tidying to make sure expression matrices are all conformable and
# only contain data with cell type assignments

# explicitly read in the study data
datah2 <- read.table("Pancreas/Data/GSE85241_SFnorm.tsv",
                     h=TRUE, sep="\t", stringsAsFactors=FALSE)

# gene IDs are in the last column
genes2 <- datah2$gene_id
rownames(datah2) <- genes2
datah2 <- datah2[, 1:(dim(datah2)[2]-1)]

hvg2 <- read.table("Pancreas/Data/GSE85241-HVG.tsv",
                   h=TRUE, sep="\t", stringsAsFactors=F)
HVG2 <- hvg2$gene_id

meta2 <- read.table("Pancreas/Data/GSE85241_marker_metadata.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)

# select meta data for which there is expression data
meta2 <- meta2[meta2$Sample %in% colnames(datah2), ]

# standardize the cell labels
## NB there are not any samples with unassigned cell types
celltypes2 <- meta2$CellType
no.label2 <- meta2$Sample[meta2$CellType == ""]
datah2 <- datah2[, !colnames(datah2) %in% no.label2]
meta2 <- meta2[!meta2$Sample %in% no.label2, ]
celltypes2 <- celltypes2[celltypes2 != ""]
samples2 <- colnames(datah2)

# check all dimensions match up
if(dim(datah2)[2] == dim(meta2)[1]) {dim(datah2)[2] == length(celltypes2)}


############
# GSE86473 #
############
# do a bit of tidying to make sure expression matrices are all conformable and
# only contain data with cell type assignments

# explicitly read in the study data
datah3 <- read.table("Pancreas/Data/GSE86473_SFnorm.tsv",
                     h=TRUE, sep="\t", stringsAsFactors=FALSE)

# gene IDs are in the last column
genes3 <- datah3$gene_id
rownames(datah3) <- genes3
datah3 <- datah3[, 1:(dim(datah3)[2]-1)]

hvg3 <- read.table("Pancreas/Data/GSE86473-HVG.tsv",
                   h=TRUE, sep="\t", stringsAsFactors=F)
HVG3 <- hvg3$gene_id

meta3 <- read.table("Pancreas/Data/GSE86473_metadata.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)

# select meta data for which there is expression data
meta3 <- meta3[meta3$Sample %in% colnames(datah3), ]

# standardize the cell labels
## NB there are not any samples with unassigned cell types
celltypes3 <- meta3$CellType
no.label3 <- meta3$Sample[meta3$CellType == ""]
datah3 <- datah3[, !colnames(datah3) %in% no.label3]
meta3 <- meta3[!meta3$Sample %in% no.label3, ]
celltypes3 <- celltypes3[celltypes3 != ""]
samples3 <- colnames(datah3)

# check all dimensions match up
if(dim(datah3)[2] == dim(meta3)[1]) {dim(datah3)[2] == length(celltypes3)}


###############
# E-MTAB-5061 #
###############
# do a bit of tidying to make sure expression matrices are all conformable and
# only contain data with cell type assignments

# explicitly read in the study data
datah4 <- read.table("Pancreas/Data/E-MTAB-5061_SFnorm.tsv",
                     h=TRUE, sep="\t", stringsAsFactors=FALSE)

# gene IDs are in the last column
genes4 <- datah4$gene_id
rownames(datah4) <- genes4
datah4 <- datah4[, 1:(dim(datah4)[2]-1)]

hvg4 <- read.table("Pancreas/Data/E-MTAB-5061-HVG.tsv",
                   h=TRUE, sep="\t", stringsAsFactors=F)
HVG4 <- hvg4$gene_id

meta4 <- read.table("Pancreas/Data/E-MTAB-5061_metadata.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)

# select meta data for which there is expression data
meta4 <- meta4[meta4$Sample %in% colnames(datah4), ]

# standardize the cell labels
## NB there are not any samples with unassigned cell types
celltypes4 <- meta4$CellType
no.label4 <- meta4$Sample[meta4$CellType == ""]
datah4 <- datah4[, !colnames(datah4) %in% no.label4]
meta4 <- meta4[!meta4$Sample %in% no.label4, ]
celltypes4 <- celltypes4[celltypes4 != ""]
samples4 <- colnames(datah4)

# check all dimensions match up
if(dim(datah4)[2] == dim(meta4)[1]) {dim(datah4)[2] == length(celltypes4)}

# create one big meta data frame
all.meta <- do.call(rbind.data.frame, list("b1"=meta1[, c("Sample", "CellType", "Protocol", "Study")],
                                           "b2"=meta2[, c("Sample", "CellType", "Protocol", "Study")],
                                           "b3"=meta3[, c("Sample", "CellType", "Protocol", "Study")],
                                           "b4"=meta4[, c("Sample", "CellType", "Protocol", "Study")]))

##########################################################################
##########################################################################
# Merge all four data matrices together based on a common set of gene IDs
common.genes <- intersect(genes1, intersect(genes2, intersect(genes3, genes4)))
r.datah1 <- datah1[common.genes, ]
r.datah2 <- datah2[common.genes, ]
r.datah3 <- datah3[common.genes, ]
r.datah4 <- datah4[common.genes, ]

# check the orders of names are the same
all(rownames(r.datah1) == rownames(r.datah2) &
      rownames(r.datah3) == rownames(r.datah4) &
      rownames(r.datah1) == rownames(r.datah4))

# merge uncorrected data based on gene IDs
r.datah1$gene_id <- rownames(r.datah1)
r.datah2$gene_id <- rownames(r.datah2)
r.datah3$gene_id <- rownames(r.datah3)
r.datah4$gene_id <- rownames(r.datah4)

raw.all <- Reduce(x=list("b1"=r.datah1, "b2"=r.datah2,
                         "b3"=r.datah3, "b4"=r.datah4),
                  f=function(x, y) merge(x, y, by='gene_id'))
rownames(raw.all) <- raw.all$gene_id
raw.all <- raw.all[, 2:dim(raw.all)[2]]

# take intersection of all hvgs
common.hvgs <- intersect(HVG1, intersect(HVG2, intersect(HVG3, HVG4)))

# assign small weird cell types from GSE85241 and E=MTAB-5061 to 'other'
all.meta$CellType[all.meta$CellType == "PP"] <- "Gamma"
all.meta$CellType[grepl(all.meta$CellType, pattern="Mesenchyme")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Co-ex")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Endo")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Epsi")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Mast")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="MHC")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Uncl")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Not")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="PSC")] <- "other"

# setup the matrix of highly variable gene expression across all cells
# this needs to drop any all-zero rows/columns
raw.hvg <- raw.all[rownames(raw.all) %in% common.hvgs, ]

## set cell type colorings
allcolors <- c("#CD0303", "#F59500", 
               "#C400FE", "#17C4CC",
               "#352FBF", "#716D6E",
               "#D4D407")
# allcolors <- labels2colors(all.meta$CellType)
# allcolors[allcolors == "red"] <- "deeppink"
# allcolors[allcolors == "yellow"] <- "orange1" # "darkgoldenrod1"

##########################################################################
##########################################################################
## Uncorrected data 
# perform tSNE on uncorrected data using HVGs
uncorrected <- as.matrix(t(raw.hvg[rownames(raw.hvg) %in% common.hvgs, ]))
set.seed(0)
tsne.unc <- Rtsne(uncorrected, distance=FALSE, perplexity=100)

# input and output order are the same
unc.tsne <- data.frame(tsne.unc$Y)
colnames(unc.tsne) <- c("Dim1", "Dim2")
unc.tsne$Sample <- colnames(raw.hvg)
unc.merge <- merge(unc.tsne, all.meta, by='Sample')
cell.colors <- unique(allcolors)
names(cell.colors) <- unique(unc.merge$CellType)

# plot of points by study
unc.bybatch <- ggplot(unc.merge, aes(x=Dim1, y=Dim2, colour=Study)) +
  geom_point(size=2) + theme_classic() +
  scale_y_continuous(limits=c(-40, 40)) +
  scale_x_continuous(limits=c(-40, 40)) +
  guides(colour=FALSE)

# plot of points by cell type and shape by study
unc.bycell <- ggplot(unc.merge, aes(x=Dim1, y=Dim2, 
                                    colour=CellType,
                                    shape=Study,
                                    group=CellType)) +
  geom_point(size=2) + theme_classic() +
  scale_colour_manual(values=cell.colors) +
  scale_y_continuous(limits=c(-40, 40)) +
  scale_x_continuous(limits=c(-40, 40)) +
  guides(colour=FALSE, shape=FALSE) +
  labs(x="tSNE 1", y="tSNE 2", main="Uncorrected")
  
## MNN batch correction
hvg.common <- common.hvgs[common.hvgs %in% common.genes]
Xmnn <- mnnCorrect(as.matrix(r.datah1[, 1:(dim(r.datah1)[2]-1)]),
                   as.matrix(r.datah2[, 1:(dim(r.datah2)[2]-1)]),
                   as.matrix(r.datah3[, 1:(dim(r.datah3)[2]-1)]),
                   as.matrix(r.datah4[, 1:(dim(r.datah4)[2]-1)]),
                   subset.row=hvg.common, 
                   svd.dim=NA,
                   k=20, sigma=0.1,
                   cos.norm=TRUE)

# combine corrected matrices together
corrected.df <- do.call(cbind.data.frame, Xmnn$corrected)
corrected.mat <- as.matrix(t(corrected.df))

set.seed(0)
tsne.c <- Rtsne(corrected.mat, is_distance=FALSE, perplexity=100)

# input and output order are the same
mnn.tsne <- data.frame(tsne.c$Y)
colnames(mnn.tsne) <- c("Dim1", "Dim2")
mnn.tsne$Sample <- colnames(raw.hvg)
mnn.merge <- merge(mnn.tsne, all.meta, by='Sample')
cell.colors <- unique(allcolors)
names(cell.colors) <- unique(mnn.merge$CellType)

# plot by batch ID
correct.bybatch <- ggplot(mnn.merge, aes(x=Dim1, y=Dim2, colour=Study, shape=CellType)) +
  geom_point(size=2) + theme_classic() +
  scale_y_continuous(limits=c(-40, 40)) +
  scale_x_continuous(limits=c(-40, 40)) +
  guides(colour=FALSE, shape=FALSE)

# plot by cell type, and shape by batch
correct.bycell <- ggplot(mnn.merge, aes(x=Dim1, y=Dim2, colour=CellType, shape=Study)) +
  geom_point(size=2) + theme_classic() +
  scale_colour_manual(values=cell.colors) +
  scale_y_continuous(limits=c(-40, 40)) +
  scale_x_continuous(limits=c(-40, 40)) +
  #guides(colour=FALSE, shape=FALSE) +
  labs(x="tSNE 1", y="tSNE 2", main="MNN corrected")

##########################################################################
##########################################################################
## limma batch correction
library(limma)

# construct a factor containing the batch IDs
batch.vector <- factor(all.meta$Study)
length(batch.vector) == dim(raw.hvg)[2]

# remove batch effect based on batch id
Xlm <- removeBatchEffect(raw.hvg, batch.vector)
lm.mat <- as.matrix(t(Xlm))

set.seed(0)
tsne.lm <- Rtsne(lm.mat, is_distance=FALSE, perplexity = 100)

# input and output order are the same
lm.tsne <- data.frame(tsne.lm$Y)
colnames(lm.tsne) <- c("Dim1", "Dim2")
lm.tsne$Sample <- colnames(raw.hvg)
lm.merge <- merge(lm.tsne, all.meta, by='Sample')
cell.colors <- unique(allcolors)
names(cell.colors) <- unique(lm.merge$CellType)

# plot by batch ID
lm.bybatch <- ggplot(lm.merge, aes(x=Dim1, y=Dim2, colour=Study, shape=CellType)) +
  geom_point(size=2) + theme_classic() +
  scale_y_continuous(limits=c(-40, 40)) +
  scale_x_continuous(limits=c(-40, 40)) +
  guides(colour=FALSE, shape=FALSE)

# plot by cell type, and shape by batch ID
lm.bycell<- ggplot(lm.merge, aes(x=Dim1, y=Dim2, colour=CellType, shape=Study)) +
  geom_point(size=2) + theme_classic() +
  scale_colour_manual(values=cell.colors) +
  scale_y_continuous(limits=c(-40, 40)) +
  scale_x_continuous(limits=c(-40, 40)) +
  guides(colour=FALSE, shape=FALSE) +
  labs(x="tSNE 1", y="tSNE 2", main="limma corrected")

##########################################################################
##########################################################################
## ComBat correction
library(sva)
cleandat.combat <- ComBat(raw.hvg, all.meta$Study,
                          mod=NULL, prior.plots=FALSE)
combat.mat <- as.matrix(t(cleandat.combat))

set.seed(0)
tsne.combat <- Rtsne(combat.mat, is_distance=FALSE, perplexity=100)

# input and output order are the same
combat.tsne <- data.frame(tsne.combat$Y)
colnames(combat.tsne) <- c("Dim1", "Dim2")
combat.tsne$Sample <- colnames(raw.hvg)
combat.merge <- merge(combat.tsne, all.meta, by='Sample')
cell.colors <- unique(allcolors)
names(cell.colors) <- unique(combat.merge$CellType)

# plot by batch ID
combat.bybatch <- ggplot(combat.merge, aes(x=Dim1, y=Dim2, colour=Study, shape=CellType)) +
  geom_point(size=2) + theme_classic() +
  scale_y_continuous(limits=c(-40, 40)) +
  scale_x_continuous(limits=c(-40, 40)) +
  guides(colour=FALSE, shape=FALSE)

# plot by cell type and shape by batch ID
combat.bycell <- ggplot(combat.merge, aes(x=Dim1, y=Dim2, colour=CellType, shape=Study)) +
  geom_point(size=2) + theme_classic() +
  scale_colour_manual(values=cell.colors) +
  scale_y_continuous(limits=c(-40, 40)) +
  scale_x_continuous(limits=c(-40, 40)) +
  guides(colour=FALSE, shape=FALSE) +
  labs(x="tSNE 1", y="tSNE 2", main="ComBat corrected")

##########################################################################
##########################################################################
# write out a single file for all of the corrected data and a 
# file for the combined meta data
corrected.df$gene_id <- common.hvgs
write.table(corrected.df, file="Pancreas/Data/mnnCorrected.tsv", row.names=FALSE, sep="\t", quote=FALSE)

## ouput a single plot of all 4 tSNEs coloured by cell type
plot_grid(unc.bycell, correct.bycell,
          lm.bycell, combat.bycell,
          ncol=2)
##########################################################################
##########################################################################
## compute Silhouette coefficients on t-SNE coordinates
require(cluster)
ct.fac <- factor(all.meta$CellType)
# only calculate silhouette on alpha cells.

dd.unc <- as.matrix(dist(tsne.unc$Y)) 
score_sil <- silhouette(as.numeric(ct.fac), dd.unc)
sil_unc <- score_sil[, 3]  #for uncorrected data

dd.c <- as.matrix(dist(tsne.c$Y))
score_sil <- silhouette(as.numeric(ct.fac), dd.c)
sil_c<-score_sil[,3] #for MNN corrected data

dd.lm <- as.matrix(dist(tsne.lm$Y))
score_sil <- silhouette(as.numeric(ct.fac), dd.lm)
sil_lm<-score_sil[,3] #for limma corrected data

dd.com <- as.matrix(dist(tsne.combat$Y))
score_sil <- silhouette(as.numeric(ct.fac), dd.com)
sil_com<-score_sil[,3] #for ComBat corrected data

### boxplot of Silhouette coefficients
sils <- cbind(sil_unc, sil_c, sil_lm, sil_com)

png(file="sils_alltypes_tsnespace.png", width=900, height=700)
par(mfrow=c(1, 1), mar=c(4, 6, 2, 2), cex.axis=1.5, cex.main=2, cex.lab=2)
boxplot(sils, main="", names=c("Raw", "MNN", "limma", "ComBat"),
        ylim=c(-0.4, 0.4),
	      lwd=2, ylab="Silhouette coefficient")#, col="Yellow", ylab="Alpha dists")
dev.off()

##################