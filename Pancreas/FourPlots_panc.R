#Code for the t-SNE plots of pancreas data sets and the Silhouette coefficients before and after batch correction by different methods (main text Figure 4).
#this.dir <- dirname(parent.frame(2)$ofile)
#setwd(this.dir)

#load("Pancreas/raw_complete4DataSets.RData") #load the output of "preparedata.R"

# need to create 'results' directory
#dir.create(paste0(this.dir, "/results"))
# setwd("results")
require(WGCNA)
require(Rtsne)
library(scales)
library(scran)
library(ggplot2)

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

# sample 1000 cells from each batch for t-SNE plotting for computational reasons
set.seed(2)
sub1 <- colnames(r.datah1)[sample(1:dim(r.datah1)[2], 1000)]
sub2 <- colnames(r.datah2)[sample(1:dim(r.datah2)[2], 1000)]
sub3 <- colnames(r.datah3)[sample(1:dim(r.datah3)[2], 1000)]
sub4 <- colnames(r.datah4)[sample(1:dim(r.datah4)[2], 1000)]
subsamples <- c(sub1, sub2, sub3, sub4)

sub.datah1 <- r.datah1[common.genes, sub1]
sub.datah1$gene_id <- rownames(sub.datah1)

sub.datah2 <- r.datah2[common.genes, sub2]
sub.datah2$gene_id <- rownames(sub.datah2)

sub.datah3 <- r.datah3[common.genes, sub3]
sub.datah3$gene_id <- rownames(sub.datah3)

sub.datah4 <- r.datah4[common.genes, sub4]
sub.datah4$gene_id <- rownames(sub.datah4)

all.sub <- Reduce(x=list("b1"=sub.datah1, "b2"=sub.datah2,
                         "b3"=sub.datah3, "b4"=sub.datah4),
                  f=function(x, y) merge(x, y, by='gene_id'))
rownames(all.sub) <- all.sub$gene_id
all.sub <- all.sub[, 2:dim(all.sub)[2]]

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

# take union of all hvgs together
all.hvgs <- unique(c(HVG1, HVG2, HVG3, HVG4))
common.hvgs <- intersect(HVG1, intersect(HVG2, intersect(HVG3, HVG4)))
allsamples <- c(samples1, samples2, samples3, samples4)

# tidy cell type to just keep the major lineages
celltypes <- c(celltypes1, celltypes2, celltypes3, celltypes4)
celltypes[celltypes == "PP"] <- "Gamma"
celltypes[grepl(celltypes, pattern="Mesenchyme")] <- "other"
celltypes[grepl(celltypes, pattern="Co-ex")] <- "other"
celltypes[grepl(celltypes, pattern="Endo")] <- "other"
celltypes[grepl(celltypes, pattern="Epsi")] <- "other"
celltypes[grepl(celltypes, pattern="Mast")] <- "other"
celltypes[grepl(celltypes, pattern="MHC")] <- "other"
celltypes[grepl(celltypes, pattern="Uncl")] <- "other"
celltypes[grepl(celltypes, pattern="Not")] <- "other"
celltypes[grepl(celltypes, pattern="PSC")] <- "other"

# celltypes[celltypes == "co-e"] <- "other"
# celltypes[celltypes == "endo"] <- "other"
# celltypes[celltypes == "epsi"] <- "other"
# celltypes[celltypes == "mast"] <- "other"
# celltypes[celltypes == "mhc "] <- "other"
# celltypes[celltypes == "uncl"] <- "other"
# celltypes[celltypes == "psc "] <- "other"
# 
# celltypes[celltypes == "alph"] <- "alpha"
# celltypes[celltypes == "delt"] <- "delta"
# celltypes[celltypes == "gamm"] <- "gamma"
# celltypes[celltypes == "psc "] <- "other"

# sub sample the cell types
sub.type1 <- meta1$CellType[meta1$Sample %in% sub1]
sub.type2 <- meta2$CellType[meta2$Sample %in% sub2]
sub.type3 <- meta3$CellType[meta3$Sample %in% sub3]
sub.type4 <- meta4$CellType[meta4$Sample %in% sub4]

subtypes <- c(sub.type1, sub.type2, sub.type3, sub.type4)
subtypes[subtypes == "PP"] <- "Gamma"
subtypes[grepl(subtypes, pattern="Mesenchyme")] <- "other"
subtypes[grepl(subtypes, pattern="Co-ex")] <- "other"
subtypes[grepl(subtypes, pattern="Endo")] <- "other"
subtypes[grepl(subtypes, pattern="Epsi")] <- "other"
subtypes[grepl(subtypes, pattern="Mast")] <- "other"
subtypes[grepl(subtypes, pattern="MHC")] <- "other"
subtypes[grepl(subtypes, pattern="Uncl")] <- "other"
subtypes[grepl(subtypes, pattern="Not")] <- "other"
subtypes[grepl(subtypes, pattern="PSC")] <- "other"

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

# check cell types and samples are the same length
length(allsamples) == length(celltypes)

# check subsamples cell types are the same length as the subsampled IDs
length(subtypes) == length(subsamples)

# setup the matrix of highly variable gene expression across all cells
# this needs to drop any all-zero rows/columns
raw.hvg <- raw.all[rownames(raw.all) %in% all.hvgs, ]
sub.hvg <- all.sub[rownames(all.sub) %in% all.hvgs, ]

# get a vector of study IDs
batch.id <- c(meta1$Study, meta2$Study, meta3$Study, meta4$Study)
sub.batch <- c(meta1$Study[meta1$Sample %in% sub1],
               meta2$Study[meta2$Sample %in% sub2],
               meta3$Study[meta3$Sample %in% sub3],
               meta4$Study[meta4$Sample %in% sub4])

##### set cell type colorings
allcolors <- labels2colors(celltypes)
allcolors[allcolors == "red"] <- "deeppink"
allcolors[allcolors == "yellow"] <- "orange1" # "darkgoldenrod1"

subcolours <- labels2colors(subtypes)
subcolours[subcolours == "red"] <- "deeppink"
subcolours[subcolours == "yellow"] <- "orange1"

N <- c(1000, 2000, 3000, 4000)

# N<-c(length(celltype2),length(celltype2)+length(celltype3),
#      length(celltype2)+length(celltype3)+length(celltype4),
#      length(celltype2)+length(celltype3)+length(celltype4)+length(celltype1))

#### Uncorrected data 
# this doesn't need to be a distance matrix!
#all.dists2.unc <- as.matrix(dist(t(raw.all[all.hvgs, allsamples])))
#all.dists2.unc <- as.matrix(dist(t(sub.hvg)))
all.dists2.unc <- as.matrix(t(raw.hvg))

par(mfrow=c(1, 1))

set.seed(0)
tsne.unc <- Rtsne(all.dists2.unc, distance=FALSE, perplexity = 100)

# input and output order are the same
unc.tsne <- data.frame(tsne.unc$Y)
colnames(unc.tsne) <- c("Dim1", "Dim2")
unc.tsne$Sample <- colnames(raw.hvg)
unc.merge <- merge(unc.tsne, all.meta, by='Sample')
cell.colors <- unique(allcolors)
names(cell.colors) <- unique(unc.merge$CellType)

ggplot(unc.merge, aes(x=Dim1, y=Dim2, colour=Study, shape=CellType)) +
  geom_point(size=2) + theme_classic()

ggplot(unc.merge, aes(x=Dim1, y=Dim2, colour=CellType, shape=Study)) +
  geom_point(size=1) + theme_classic() +
  scale_colour_manual(values=cell.colors) +
  scale_y_continuous(limits=c(-50, 50)) +
  scale_x_continuous(limits=c(-50, 50))
  
# png(file="unc4321.png", width=900, height=700)
# par(mfrow=c(1,1), mar=c(6,6,4,2), cex.axis=2, cex.main=3, cex.lab=2.5)
# 
# # create a plot with colour as cell type and shape as batch
# plot(unc.tsne$Dim1,
#      unc.tsne$Dim2,
#      pch=c(rep(3, dim(unc.tsne[unc.tsne$Study == "GSE81076",])[2]),
#            rep(18, dim(unc.tsne[unc.tsne$Study == "GSE85241",])[2]),
#            rep(1, dim(unc.tsne[unc.tsne$Study == "GSE86473",])[2]),
#            rep(4, dim(unc.tsne[unc.tsne$Study == "E-MTAB-5061",])[2])),
#      cex=1, 
#      col=alpha(unc.tsne$CellColour, 0.6),
#      main="Uncorrected", xlab="tSNE 1", ylab="tSNE 2")


# points(tsne.unc$Y[(N[1]+1):N[2], 1], tsne.unc$Y[(N[1]+1):N[2], 2], pch=18,cex=4, col=alpha(allcolors[(N[1]+1):N[2]], 0.6))
# points(tsne.unc$Y[(N[2]+1):N[3], 1], tsne.unc$Y[(N[2]+1):N[3], 2], pch=1,cex=4, col=alpha(allcolors[(N[2]+1):N[3]], 0.6))
# points(tsne.unc$Y[(N[3]+1):N[4], 1], tsne.unc$Y[(N[3]+1):N[4], 2], pch=4,cex=4, col=alpha(allcolors[(N[3]+1):N[4]], 0.6))
# dev.off()

### MNN batch correction
# HVGs need to be a subset of the genes commonly expressed across all data sets
hvg.common <- all.hvgs[all.hvgs %in% common.genes]
Xmnn <- mnnCorrect(as.matrix(r.datah1[, 1:(dim(r.datah1)[2]-1)]),
                   as.matrix(r.datah2[, 1:(dim(r.datah2)[2]-1)]),
                   as.matrix(r.datah3[, 1:(dim(r.datah3)[2]-1)]),
                   as.matrix(r.datah4[, 1:(dim(r.datah4)[2]-1)]),
                   subset.row=hvg.common,
                   k=30, sigma=0.5, svdim=0,
                   cos.norm=TRUE)

# corre <- cbind(Xmnn$corrected[[1]], Xmnn$corrected[[2]], Xmnn$corrected[[3]], Xmnn$corrected[[4]])
corre <- do.call(cbind.data.frame, Xmnn$corrected)
#all.dists2.c <- as.matrix(dist(t(corre[HVG, allsamples])))
all.dists2.c <- as.matrix(t(corre))

set.seed(0)
tsne.c <- Rtsne(all.dists2.c, is_distance=FALSE, perplexity = 100)

# input and output order are the same
mnn.tsne <- data.frame(tsne.c$Y)
colnames(mnn.tsne) <- c("Dim1", "Dim2")
mnn.tsne$Sample <- colnames(raw.hvg)
mnn.merge <- merge(mnn.tsne, all.meta, by='Sample')
cell.colors <- unique(allcolors)
names(cell.colors) <- unique(mnn.merge$CellType)

ggplot(mnn.merge, aes(x=Dim1, y=Dim2, colour=Study, shape=CellType)) +
  geom_point(size=1) + theme_classic() +
  scale_y_continuous(limits=c(-50, 50)) +
  scale_x_continuous(limits=c(-50, 50))

ggplot(mnn.merge, aes(x=Dim1, y=Dim2, colour=CellType, shape=Study)) +
  geom_point(size=1) + theme_classic() +
  scale_colour_manual(values=cell.colors) +
  scale_y_continuous(limits=c(-50, 50)) +
  scale_x_continuous(limits=c(-50, 50))


# png(file="mnn4321_s.png", width=900, height=700)
# par(mfrow=c(1,1), mar=c(6,6,4,2), cex.axis=2, cex.main=3, cex.lab=2.5)
# plot(tsne.c$Y[1:N[1],1], tsne.c$Y[1:N[1],2], pch=3, cex=4, col=alpha(allcolors[1:N[1]], 0.6),
# 			 main="MNN corrected", xlim=c(-25,25), ylim=c(-25,20), xlab="tSNE 1", ylab="tSNE 2")
# points(tsne.c$Y[(N[1]+1):N[2],1], tsne.c$Y[(N[1]+1):N[2],2], pch=18,cex=4, col=alpha(allcolors[(N[1]+1):N[2]], 0.6))
# points(tsne.c$Y[(N[2]+1):N[3],1], tsne.c$Y[(N[2]+1):N[3],2], pch=1,cex=4, col=alpha(allcolors[(N[2]+1):N[3]], 0.6))
# points(tsne.c$Y[(N[3]+1):N[4],1], tsne.c$Y[(N[3]+1):N[4],2], pch=4,cex=4, col=alpha(allcolors[(N[3]+1):N[4]], 0.6))
# dev.off()

### limma batch correction
library(limma)
# construct a factor containing the batch IDs
batch.vector <- factor(all.meta$Study)

length(batch.vector) == dim(raw.hvg)[2]

Xlm <- removeBatchEffect(raw.hvg, batch.vector)
all.dists2.lm <- as.matrix(t(Xlm))

set.seed(0)
tsne.lm <- Rtsne(all.dists2.lm, is_distance=FALSE, perplexity = 100)

# input and output order are the same
lm.tsne <- data.frame(tsne.lm$Y)
colnames(lm.tsne) <- c("Dim1", "Dim2")
lm.tsne$Sample <- colnames(raw.hvg)
lm.merge <- merge(lm.tsne, all.meta, by='Sample')
cell.colors <- unique(allcolors)
names(cell.colors) <- unique(lm.merge$CellType)

ggplot(lm.merge, aes(x=Dim1, y=Dim2, colour=Study, shape=CellType)) +
  geom_point(size=2) + theme_classic()

ggplot(lm.merge, aes(x=Dim1, y=Dim2, colour=CellType, shape=Study)) +
  geom_point(size=1) + theme_classic() +
  scale_colour_manual(values=cell.colors) +
  scale_y_continuous(limits=c(-50, 50)) +
  scale_x_continuous(limits=c(-50, 50))


# png(file="lmfit4321.png", width=900, height=700)
# par(mfrow=c(1,1), mar=c(6,6,4,2), cex.axis=2, cex.main=3, cex.lab=2.5)
# plot(tsne.lm$Y[1:N[1],1], tsne.lm$Y[1:N[1],2], pch=3,cex=4, col=alpha(allcolors[1:N[1]], 0.6),
# 			  main="limma corrected", xlim=c(-15,20), ylim=c(-25,20), xlab="tSNE 1" ,ylab="tSNE 2")
# points(tsne.lm$Y[(N[1]+1):N[2], 1], tsne.lm$Y[(N[1]+1):N[2], 2], pch=18,cex=4, col=alpha(allcolors[(N[1]+1):N[2]], 0.6))
# points(tsne.lm$Y[(N[2]+1):N[3], 1], tsne.lm$Y[(N[2]+1):N[3], 2], pch=1,cex=4, col=alpha(allcolors[(N[2]+1):N[3]], 0.6))
# points(tsne.lm$Y[(N[3]+1):N[4], 1], tsne.lm$Y[(N[3]+1):N[4], 2], pch=4,cex=4, col=alpha(allcolors[(N[3]+1):N[4]], 0.6))
# dev.off()

#### ComBat correction
library(sva)

cleandat.combat <- ComBat(raw.hvg, all.meta$Study,
                          mod=NULL, prior.plots=FALSE)

all.dists.combat <- as.matrix(t(cleandat.combat))
set.seed(0)
tsne.combat <- Rtsne(all.dists.combat, is_distance=FALSE, perplexity=100)

# input and output order are the same
combat.tsne <- data.frame(tsne.combat$Y)
colnames(combat.tsne) <- c("Dim1", "Dim2")
combat.tsne$Sample <- colnames(raw.hvg)
combat.merge <- merge(combat.tsne, all.meta, by='Sample')
cell.colors <- unique(allcolors)
names(cell.colors) <- unique(combat.merge$CellType)

ggplot(combat.merge, aes(x=Dim1, y=Dim2, colour=Study, shape=CellType)) +
  geom_point(size=2) + theme_classic()

ggplot(combat.merge, aes(x=Dim1, y=Dim2, colour=CellType, shape=Study)) +
  geom_point(size=1) + theme_classic() +
  scale_colour_manual(values=cell.colors) +
  scale_y_continuous(limits=c(-50, 50)) +
  scale_x_continuous(limits=c(-50, 50))


# png(file="combat4321.png", width=900, height=700)
# par(mfrow=c(1,1), mar=c(6,6,4,2), cex.axis=2, cex.main=3, cex.lab=2.5)
# plot(tsne.combat$Y[1:N[1], 1], tsne.combat$Y[1:N[1],2], pch=3, cex=4, 
# 			   col=alpha(allcolors[1:N[1]],0.6), 
# 			   main="ComBat corrected", xlim=c(-12,15), ylim=c(-10,10),
# 			   xlab="tSNE 1", ylab="tSNE 2")
# points(tsne.combat$Y[(N[1]+1):N[2],1], tsne.combat$Y[(N[1]+1):N[2],2],
# 				       pch=18, cex=4, col=alpha(allcolors[(N[1]+1):N[2]], 0.6))
# points(tsne.combat$Y[(N[2]+1):N[3],1], tsne.combat$Y[(N[2]+1):N[3],2],
# 				       pch=1, cex=4, col=alpha(allcolors[(N[2]+1):N[3]], 0.6))
# points(tsne.combat$Y[(N[3]+1):N[4],1], tsne.combat$Y[(N[3]+1):N[4],2],
# 				       pch=4, cex=4, col=alpha(allcolors[(N[3]+1):N[4]], 0.6))
# #plot(tsne.combat$Y[,1],tsne.combat$Y[,2], pch=16, col=labels2colors(celltypes), main="ComBat corrected")
# dev.off()

### save results
save(file="completedata_correcteds.RData", raw.all, Xmnn,
					   Xlm, cleandat.combat, allsamples, celltypes)
# ###### the legend
# png(file="leg_detailed4321.png", width=900, height=700)
# par(mfrow=c(1,1), mar=c(6,6,4,2), cex.axis=2, cex.main=3, cex.lab=2.5)
# plot(1, 2, pch=3, cex=4, col=alpha(allcolors[1],0.6),
# 	main="legend", xlim=c(-20,30), ylim=c(-13,13), xlab="tSNE 1", ylab="tSNE 2")
# forleg <- table(celltypes,allcolors)
# leg.txt <- unique(celltypes[allsamples])
# legend("bottomright", "(x,y)", legend = leg.txt, col =unique(allcolors) , pch = 4,cex = 2.5,bty = "n",lwd = 3,lty=0)   #, trace = TRUE)
# legend("bottomleft", "(x,y)", legend = leg.txt, col =unique(allcolors) , pch = 1,cex = 2.5,bty = "n",lwd = 3,lty=0)   #, trace = TRUE)
# legend("topright", "(x,y)", legend = leg.txt, col =unique(allcolors) , pch = 18,cex = 2.5,bty = "n",lwd = 3,lty=0)   #, trace = TRUE)
# legend("topleft", "(x,y)", legend = leg.txt, col =unique(allcolors) , pch = 3,cex = 2.5,bty = "n",lwd = 3,lty=0)   #, trace = TRUE)
# dev.off()

# ##################### write.table
# Ns <- c(ncol(datah4), ncol(datah3), ncol(datah2), ncol(datah1))
# correh4 <- corre[, 1:Ns[1]]
# correh3 <- corre[, (Ns[1] + 1):(Ns[1] + Ns[2])]
# correh2 <- corre[, (Ns[1] + Ns[2] + 1):(Ns[1] + Ns[2] + Ns[3])]
# correh1 <- corre[, (Ns[1] + Ns[2] + Ns[3] + 1):(Ns[1] + Ns[2] + Ns[3] + Ns[4])]
# 
# colnames(correh4) <- celltype4
# colnames(correh3) <- celltype3
# colnames(correh2) <- celltype2
# colnames(correh1) <- celltype1
# 
# colnames(correh4) <- colnames(datah4)
# colnames(correh3) <- colnames(datah3)
# colnames(correh2) <- colnames(datah2)
# colnames(correh1) <- colnames(datah1)
 
write.table(file="C_Smartseq_GSE86473.txt", correh4, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="C_Smartseq_EMATB5061.txt", correh3, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="C_CELseq_SSE85241.txt", correh2, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="C_CELseq_GSE81076.txt", correh1, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
##########

########## compute Silhouette coefficients on t-SNE coordinates
# ct.fac <- factor(celltypes[allsamples])
ct.fac <- factor(all.meta$CellType)

dd.unc <- as.matrix(dist(tsne.unc$Y)) 
score_sil <- (cluster::silhouette(as.numeric(ct.fac), dd.unc))
sil_unc <- score_sil[, 3]  #for uncorrected data

dd.c <- as.matrix(dist(tsne.c$Y))
score_sil <- (cluster::silhouette(as.numeric(ct.fac), dd.c))
sil_c<-score_sil[,3] #for MNN corrected data

dd.lm <- as.matrix(dist(tsne.lm$Y))
score_sil <- (cluster::silhouette(as.numeric(ct.fac), dd.lm))
sil_lm<-score_sil[,3] #for limma corrected data

dd.com <- as.matrix(dist(tsne.combat$Y))
score_sil <- (cluster::silhouette(as.numeric(ct.fac), dd.com))
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