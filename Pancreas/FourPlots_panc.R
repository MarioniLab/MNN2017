# Code for the t-SNE plots of pancreas data sets and the Silhouette coefficients before and after batch correction by different methods (main text Figure 4).
# PLEASE SET THE WORKING DIRECTORY TO ./MNN2017/

require(WGCNA)
require(Rtsne)
library(scales)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
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
HVG1 <- hvg1$gene_id[hvg1$HVG == 1]

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
HVG2 <- hvg2$gene_id[hvg2$HVG == 1]

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
HVG3 <- hvg3$gene_id[hvg3$HVG == 1]

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
HVG4 <- hvg4$gene_id[hvg4$HVG == 1]

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

# to get a set of informative highly variable genes we will meta-analyse the
# p-values from the indivdual chi-squared tests on from each batch
# using the mean p-value, because it has a breakpoint of 0, we will
# then take the set of HVGs across batches with a combined pval <= 0.01
# if there are many large pvalues, i.e. 1 the geometric mean will better handle
# this situation so no pvalue > 1
# we can only meta-analyse the genes which are in common to all data sets
# define a geomtric mean function, no inbuilt function in R
# credit to https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in?answertab=votes#tab-top
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# geometric mean on v.small fractional numbers might become quite
# unstable.  Could just use take geometric mean on a log scale
# this still uses the intersection of all genes, if it doesn't
# appear in the list of genes at all, set the pval - 1.
all.genes <- unique(c(genes1, genes2, genes3, genes4))
common.hvgs <- c()
for(i in seq_along(common.genes)){
  g <- common.genes[i]
  pvals <- c(hvg1$pval[hvg1$gene_id == g],
             hvg2$pval[hvg2$gene_id == g],
             hvg3$pval[hvg3$gene_id == g],
             hvg4$pval[hvg4$gene_id == g])
  meta.p <- gm_mean(-log10(pvals))
  if(10**(-meta.p) <= 0.001){
    common.hvgs <- c(common.hvgs, g)
  }
}

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

# construct a variable that captures interaction between cell types and batch
# colour by cell type, shade by batch
# note not all cell types are present in GSE86473
acinar <- c("#ffff00")#, "#d4d432", "#d4d496")
alpha <- c("#ff0000")#, "#9b3232", "#cc6464", "#cc9696")
beta <- c("#c400ff")#, "#c496ff", "#7e00b9", "#60329b")
delta <- c("#ff7800")#, "#ffaa64", "#c84b23", "#c87d64")
ductal <- c("#00f5ff")#, "#00c8c8", "#006464")
gamma <- c("#0000ff")#, "#6464e1", "#afafe1", "#000087")
other <- c("#000000")#, "#323232", "#646464")

interact.cols <- c(acinar, alpha, beta, delta, ductal, gamma, other)
names(interact.cols) <- levels(all.meta$CellType)

batch.cols <- brewer.pal(4, "Set3")
names(batch.cols) <- unique(all.meta$Study)

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

# plot of points by cell type and shape by study
unc.bycell <- ggplot(unc.merge, aes(x=Dim1, y=Dim2, 
                                    fill=CellType,
                                    group=CellType)) +
  geom_point(size=2, shape=21) + theme_classic() +
  scale_fill_manual(values=interact.cols) +
  scale_y_continuous(limits=c(-50, 50)) +
  scale_x_continuous(limits=c(-50, 50)) +
  labs(x="tSNE 1", y="tSNE 2", title="Uncorrected") +
  guides(fill=FALSE) +
  theme(axis.title=element_text(size=24), axis.text=element_text(size=16))

ggsave(unc.bycell,
       filename="Pancreas/Uncorrected_pancreas-tSNE.png",
       height=5.75, width=9.3035, dpi=300)

# plot of points by cell type and shape by study
unc.bybatch <- ggplot(unc.merge, aes(x=Dim1, y=Dim2, 
                                    fill=Study,
                                    group=Study)) +
  geom_point(size=2, shape=21) + theme_classic() +
  scale_fill_manual(values=batch.cols) +
  scale_y_continuous(limits=c(-50, 50)) +
  scale_x_continuous(limits=c(-50, 50)) +
  labs(x="tSNE 1", y="tSNE 2", title="Uncorrected") +
  guides(fill=FALSE) +
  theme(axis.title=element_text(size=24), axis.text=element_text(size=16))

ggsave(unc.bybatch,
       filename="Pancreas/Uncorrected_by_batch_pancreas-tSNE.png",
       height=5.75, width=9.3035, dpi=300)

##########################################################################
##########################################################################
## MNN batch correction
# source the devel version of mnnCorrect for the purposes of the ms
# source("SomeFuncs/devel_mnnCorrect.R")
hvg.common <- common.hvgs[common.hvgs %in% common.genes]
Xmnn <- mnnCorrect(as.matrix(r.datah1[, 1:(dim(r.datah1)[2]-1)]),
                   as.matrix(r.datah2[, 1:(dim(r.datah2)[2]-1)]),
                   as.matrix(r.datah3[, 1:(dim(r.datah3)[2]-1)]),
                   as.matrix(r.datah4[, 1:(dim(r.datah4)[2]-1)]),
                   hvg.genes=hvg.common, 
                   svd.dim=FALSE,
                   k=20, sigma=0.1,
                   cos.norm="inside_and_output")

# combine corrected matrices together
corrected.df <- do.call(cbind.data.frame, Xmnn$corrected)
corrected.mat <- as.matrix(t(corrected.df))
# attach original column names, output order is the same as input order
colnames(corrected.df) <- c(colnames(r.datah1[, 1:(dim(r.datah1)[2]-1)]), colnames(r.datah2[, 1:(dim(r.datah2)[2]-1)]),
                            colnames(r.datah3[, 1:(dim(r.datah3)[2]-1)]), colnames(r.datah4[, 1:(dim(r.datah4)[2]-1)]))

set.seed(0)
tsne.c <- Rtsne(corrected.mat, is_distance=FALSE, perplexity=30)

# input and output order are the same
mnn.tsne <- data.frame(tsne.c$Y)
colnames(mnn.tsne) <- c("Dim1", "Dim2")
mnn.tsne$Sample <- colnames(raw.hvg)
mnn.merge <- merge(mnn.tsne, all.meta, by='Sample')

# plot of points by cell type and shape by study
mnn.bycell <- ggplot(mnn.merge, aes(x=Dim1, y=Dim2, 
                                    fill=CellType,
                                    group=CellType)) +
  geom_point(size=2, shape=21, alpha=1) + theme_classic() +
  scale_fill_manual(values=interact.cols) +
  scale_y_continuous(limits=c(-50, 50)) +
  scale_x_continuous(limits=c(-50, 50)) +
  labs(x="tSNE 1", y="tSNE 2", title="MNN corrected") +
  guides(fill=FALSE) +
  theme(axis.title=element_text(size=24), axis.text=element_text(size=16))

ggsave(mnn.bycell,
       filename="Pancreas/MNNcorrect_pancreas-tSNE.png",
       height=5.75, width=9.3035, dpi=300)

# plot of points by cell type and shape by study
mnn.bybatch <- ggplot(mnn.merge, aes(x=Dim1, y=Dim2, 
                                     fill=Study,
                                     group=Study)) +
  geom_point(size=2, shape=21) + theme_classic() +
  scale_fill_manual(values=batch.cols) +
  scale_y_continuous(limits=c(-50, 50)) +
  scale_x_continuous(limits=c(-50, 50)) +
  labs(x="tSNE 1", y="tSNE 2", title="MNN corrected") +
  guides(fill=FALSE) +
  theme(axis.title=element_text(size=24), axis.text=element_text(size=16))

ggsave(mnn.bybatch,
       filename="Pancreas/MNNcorrect_by_batch_pancreas-tSNE.png",
       height=5.75, width=9.3035, dpi=300)
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
tsne.lm <- Rtsne(lm.mat, is_distance=FALSE, perplexity = 30)

# input and output order are the same
lm.tsne <- data.frame(tsne.lm$Y)
colnames(lm.tsne) <- c("Dim1", "Dim2")
lm.tsne$Sample <- colnames(raw.hvg)
lm.merge <- merge(lm.tsne, all.meta, by='Sample')

# plot by cell type, and shape by batch ID
lm.bycell <- ggplot(lm.merge, aes(x=Dim1, y=Dim2, 
                                    fill=CellType,
                                    group=CellType)) +
  geom_point(size=2, shape=21, alpha=1) + theme_classic() +
  scale_fill_manual(values=interact.cols) +
  scale_y_continuous(limits=c(-50, 50)) +
  scale_x_continuous(limits=c(-50, 50)) +
  labs(x="tSNE 1", y="tSNE 2", title="limma corrected") +
  guides(fill=FALSE) +
  theme(axis.title=element_text(size=24), axis.text=element_text(size=16))

ggsave(lm.bycell,
       filename="Pancreas/limma_pancreas-tSNE.png",
       height=5.75, width=9.3035, dpi=300)

# plot of points by cell type and shape by study
lm.bybatch <- ggplot(lm.merge, aes(x=Dim1, y=Dim2, 
                                     fill=Study,
                                     group=Study)) +
  geom_point(size=2, shape=21) + theme_classic() +
  scale_fill_manual(values=batch.cols) +
  scale_y_continuous(limits=c(-50, 50)) +
  scale_x_continuous(limits=c(-50, 50)) +
  labs(x="tSNE 1", y="tSNE 2", title="limma corrected") +
  guides(fill=FALSE) +
  theme(axis.title=element_text(size=24), axis.text=element_text(size=16))

ggsave(lm.bybatch,
       filename="Pancreas/limma_by_batch_pancreas-tSNE.png",
       height=5.75, width=9.3035, dpi=300)
##########################################################################
##########################################################################
## ComBat correction
library(sva)
cleandat.combat <- ComBat(raw.hvg, all.meta$Study,
                          mod=NULL, prior.plots=FALSE)
combat.mat <- as.matrix(t(cleandat.combat))

set.seed(0)
tsne.combat <- Rtsne(combat.mat, is_distance=FALSE, perplexity=30)

# input and output order are the same
combat.tsne <- data.frame(tsne.combat$Y)
colnames(combat.tsne) <- c("Dim1", "Dim2")
combat.tsne$Sample <- colnames(raw.hvg)
combat.merge <- merge(combat.tsne, all.meta, by='Sample')

# plot by cell type, and shape by batch ID
combat.bycell <- ggplot(combat.merge, aes(x=Dim1, y=Dim2, 
                                  fill=CellType,
                                  group=CellType)) +
  geom_point(size=2, shape=21, alpha=1) + theme_classic() +
  scale_fill_manual(values=interact.cols) +
  scale_y_continuous(limits=c(-50, 50)) +
  scale_x_continuous(limits=c(-50, 50)) +
  labs(x="tSNE 1", y="tSNE 2", title="ComBat corrected") +
  guides(fill=FALSE) +
  theme(axis.title=element_text(size=24), axis.text=element_text(size=16))

ggsave(combat.bycell,
       filename="Pancreas/ComBat_pancreas-tSNE.png",
       height=5.75, width=9.3035, dpi=300)

# plot of points by cell type and shape by study
combat.bybatch <- ggplot(combat.merge, aes(x=Dim1, y=Dim2,
                                        fill=Study,
                                        group=Study)) +
  geom_point(size=2, shape=21) + theme_classic() +
  scale_fill_manual(values=batch.cols) +
  scale_y_continuous(limits=c(-50, 50)) +
  scale_x_continuous(limits=c(-50, 50)) +
  labs(x="tSNE 1", y="tSNE 2", title="ComBat corrected") +
  guides(fill=FALSE) +
  theme(axis.title=element_text(size=24), axis.text=element_text(size=16))

ggsave(combat.bybatch,
       filename="Pancreas/ComBat_by_batch_pancreas-tSNE.png",
       height=5.75, width=9.3035, dpi=300)

##########################################################################
##########################################################################
# write out a single file for all of the corrected data and a 
# file for the combined meta data
corrected.df$gene_id <- common.genes
write.table(corrected.df, file="Pancreas/Data/mnnCorrected.tsv",
            row.names=FALSE, sep="\t", quote=FALSE)

##################
# output a few objects for use in other scripts for supplementary figures
celltype1 <- all.meta$CellType[all.meta$Study == "GSE81076"]
celltype2 <- all.meta$CellType[all.meta$Study == "GSE85241"]
celltype3 <- all.meta$CellType[all.meta$Study == "GSE86473"]
celltype4 <- all.meta$CellType[all.meta$Study == "EMTAB5061"]

# serialize objects to an RDS file
save(raw.all, corrected.df, lm.mat, combat.mat,
     datah1, datah2, datah3, datah4,
     celltype1, celltype2, celltype3, celltype4,
     HVG1, HVG2, HVG3, HVG4,
     common.hvgs,
     all.meta, interact.cols,
     file="Pancreas/Data/ObjectsForPlotting.RDS")