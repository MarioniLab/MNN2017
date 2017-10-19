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

# merge uncorrected data
raw.all <- cbind(r.datah1, r.datah2, r.datah3, r.datah4)

# sample 1000 cells from each batch for t-SNE plotting for computational reasons
set.seed(2)
#samples1 <- sample(1:dim(datah1)[2], 1000) 
#samples2 <- sample(1:dim(datah2)[2], 1000)
#samples3 <- sample(1:dim(datah3)[2], 1000)
#samples4 <- sample(1:dim(datah4)[2], 1000)
allsamples <- c(samples1, samples2, samples3, samples4)

# tidy cell type to just keep the major lineages
celltypes <- c(celltype1, celltype2, celltype3, celltype4)
celltypes[celltypes == "pp"] <- "gamma"
celltypes[celltypes == "mese"] <- "other"
celltypes[celltypes == "co-e"] <- "other"
celltypes[celltypes == "endo"] <- "other"
celltypes[celltypes == "epsi"] <- "other"
celltypes[celltypes == "mast"] <- "other"
celltypes[celltypes == "mhc "] <- "other"
celltypes[celltypes == "uncl"] <- "other"
celltypes[celltypes == "psc "] <- "other"

celltypes[celltypes == "alph"] <- "alpha"
celltypes[celltypes == "delt"] <- "delta"
celltypes[celltypes == "gamm"] <- "gamma"
celltypes[celltypes == "psc "] <- "other"

##### set cell type colorings
allcolors <- labels2colors(celltypes[allsamples])
allcolors[allcolors == "red"] <- "deeppink"
allcolors[allcolors == "yellow"] <- "orange1" # "darkgoldenrod1"

N <- c(1000, 2000, 3000, 4000)

# N<-c(length(celltype2),length(celltype2)+length(celltype3),
#      length(celltype2)+length(celltype3)+length(celltype4),
#      length(celltype2)+length(celltype3)+length(celltype4)+length(celltype1))

#### Uncorrected data 
all.dists2.unc <- as.matrix(dist(t(raw.all[HVG, allsamples])))
par(mfrow=c(1, 1))

set.seed(0)
tsne.unc <- Rtsne(all.dists2.unc, is_distance=TRUE)#, perplexity = 5)

png(file="unc4321.png", width=900, height=700)
par(mfrow=c(1,1), mar=c(6,6,4,2), cex.axis=2, cex.main=3, cex.lab=2.5)
plot(tsne.unc$Y[1:N[1], 1], tsne.unc$Y[1:N[1], 2], pch=3, cex=4, col=alpha(allcolors[1:N[1]],0.6),
			main="Uncorrected", xlim=c(-30,25), ylim=c(-20,30), xlab="tSNE 1", ylab="tSNE 2")

points(tsne.unc$Y[(N[1]+1):N[2], 1], tsne.unc$Y[(N[1]+1):N[2], 2], pch=18,cex=4, col=alpha(allcolors[(N[1]+1):N[2]], 0.6))
points(tsne.unc$Y[(N[2]+1):N[3], 1], tsne.unc$Y[(N[2]+1):N[3], 2], pch=1,cex=4, col=alpha(allcolors[(N[2]+1):N[3]], 0.6))
points(tsne.unc$Y[(N[3]+1):N[4], 1], tsne.unc$Y[(N[3]+1):N[4], 2], pch=4,cex=4, col=alpha(allcolors[(N[3]+1):N[4]], 0.6))
dev.off()

### MNN batch correction
inquiry_genes <- row.names(datah4)
Xmnn <- mnnCorrect(datah4, datah3, datah2, datah1,
     	#inquiry.genes=inquiry_genes,
        hvg.genes=HVG, k=20, sigma=0.1,
	cos.norm=TRUE, svd.dim=0) # batch correction is throwing an error because of non-numeric values?

corre <- cbind(Xmnn$corrected[[1]], Xmnn$corrected[[2]], Xmnn$corrected[[3]], Xmnn$corrected[[4]])
all.dists2.c <- as.matrix(dist(t(corre[HVG, allsamples])))

set.seed(0)
tsne.c <- Rtsne(all.dists2.c, is_distance=TRUE)#, perplexity = 5)

png(file="mnn4321_s.png", width=900, height=700)
par(mfrow=c(1,1), mar=c(6,6,4,2), cex.axis=2, cex.main=3, cex.lab=2.5)
plot(tsne.c$Y[1:N[1],1], tsne.c$Y[1:N[1],2], pch=3, cex=4, col=alpha(allcolors[1:N[1]], 0.6),
			 main="MNN corrected", xlim=c(-25,25), ylim=c(-25,20), xlab="tSNE 1", ylab="tSNE 2")
points(tsne.c$Y[(N[1]+1):N[2],1], tsne.c$Y[(N[1]+1):N[2],2], pch=18,cex=4, col=alpha(allcolors[(N[1]+1):N[2]], 0.6))
points(tsne.c$Y[(N[2]+1):N[3],1], tsne.c$Y[(N[2]+1):N[3],2], pch=1,cex=4, col=alpha(allcolors[(N[2]+1):N[3]], 0.6))
points(tsne.c$Y[(N[3]+1):N[4],1], tsne.c$Y[(N[3]+1):N[4],2], pch=4,cex=4, col=alpha(allcolors[(N[3]+1):N[4]], 0.6))
dev.off()

### limma batch correction
library(limma)
Xlm <- removeBatchEffect(raw.all, factor(colnames(raw.all)))
all.dists2.lm <- as.matrix(dist(t(Xlm[HVG, allsamples])))

set.seed(0)
tsne.lm<-Rtsne(all.dists2.lm, is_distance=TRUE)#, perplexity = 0.9)
png(file="lmfit4321.png", width=900, height=700)
par(mfrow=c(1,1), mar=c(6,6,4,2), cex.axis=2, cex.main=3, cex.lab=2.5)
plot(tsne.lm$Y[1:N[1],1], tsne.lm$Y[1:N[1],2], pch=3,cex=4, col=alpha(allcolors[1:N[1]], 0.6),
			  main="limma corrected", xlim=c(-15,20), ylim=c(-25,20), xlab="tSNE 1" ,ylab="tSNE 2")
points(tsne.lm$Y[(N[1]+1):N[2], 1], tsne.lm$Y[(N[1]+1):N[2], 2], pch=18,cex=4, col=alpha(allcolors[(N[1]+1):N[2]], 0.6))
points(tsne.lm$Y[(N[2]+1):N[3], 1], tsne.lm$Y[(N[2]+1):N[3], 2], pch=1,cex=4, col=alpha(allcolors[(N[2]+1):N[3]], 0.6))
points(tsne.lm$Y[(N[3]+1):N[4], 1], tsne.lm$Y[(N[3]+1):N[4], 2], pch=4,cex=4, col=alpha(allcolors[(N[3]+1):N[4]], 0.6))
dev.off()

#### ComBat correction
library(sva)

Z <- colnames(raw.all)
cleandat.combat <- ComBat(raw.all, Z, mod=NULL, prior.plots=FALSE)

all.dists.combat <- as.matrix(dist(t(cleandat.combat[HVG, allsamples])))
set.seed(0)
tsne.combat <- Rtsne(all.dists.combat, is_distance=TRUE)

png(file="combat4321.png", width=900, height=700)
par(mfrow=c(1,1), mar=c(6,6,4,2), cex.axis=2, cex.main=3, cex.lab=2.5)
plot(tsne.combat$Y[1:N[1], 1], tsne.combat$Y[1:N[1],2], pch=3, cex=4, 
			   col=alpha(allcolors[1:N[1]],0.6), 
			   main="ComBat corrected", xlim=c(-12,15), ylim=c(-10,10),
			   xlab="tSNE 1", ylab="tSNE 2")
points(tsne.combat$Y[(N[1]+1):N[2],1], tsne.combat$Y[(N[1]+1):N[2],2],
				       pch=18, cex=4, col=alpha(allcolors[(N[1]+1):N[2]], 0.6))
points(tsne.combat$Y[(N[2]+1):N[3],1], tsne.combat$Y[(N[2]+1):N[3],2],
				       pch=1, cex=4, col=alpha(allcolors[(N[2]+1):N[3]], 0.6))
points(tsne.combat$Y[(N[3]+1):N[4],1], tsne.combat$Y[(N[3]+1):N[4],2],
				       pch=4, cex=4, col=alpha(allcolors[(N[3]+1):N[4]], 0.6))
#plot(tsne.combat$Y[,1],tsne.combat$Y[,2], pch=16, col=labels2colors(celltypes), main="ComBat corrected")
dev.off()

### save results
save(file="completedata_correcteds.RData", raw.all, Xmnn,
					   Xlm, cleandat.combat, allsamples, celltypes)
###### the legend
png(file="leg_detailed4321.png", width=900, height=700)
par(mfrow=c(1,1), mar=c(6,6,4,2), cex.axis=2, cex.main=3, cex.lab=2.5)
plot(1, 2, pch=3, cex=4, col=alpha(allcolors[1],0.6),
	main="legend", xlim=c(-20,30), ylim=c(-13,13), xlab="tSNE 1", ylab="tSNE 2")
forleg <- table(celltypes,allcolors)
leg.txt <- unique(celltypes[allsamples])
legend("bottomright", "(x,y)", legend = leg.txt, col =unique(allcolors) , pch = 4,cex = 2.5,bty = "n",lwd = 3,lty=0)   #, trace = TRUE)
legend("bottomleft", "(x,y)", legend = leg.txt, col =unique(allcolors) , pch = 1,cex = 2.5,bty = "n",lwd = 3,lty=0)   #, trace = TRUE)
legend("topright", "(x,y)", legend = leg.txt, col =unique(allcolors) , pch = 18,cex = 2.5,bty = "n",lwd = 3,lty=0)   #, trace = TRUE)
legend("topleft", "(x,y)", legend = leg.txt, col =unique(allcolors) , pch = 3,cex = 2.5,bty = "n",lwd = 3,lty=0)   #, trace = TRUE)
dev.off()

##################### write.table
Ns <- c(ncol(datah4), ncol(datah3), ncol(datah2), ncol(datah1))
correh4 <- corre[, 1:Ns[1]]
correh3 <- corre[, (Ns[1] + 1):(Ns[1] + Ns[2])]
correh2 <- corre[, (Ns[1] + Ns[2] + 1):(Ns[1] + Ns[2] + Ns[3])]
correh1 <- corre[, (Ns[1] + Ns[2] + Ns[3] + 1):(Ns[1] + Ns[2] + Ns[3] + Ns[4])]

colnames(correh4) <- celltype4
colnames(correh3) <- celltype3
colnames(correh2) <- celltype2
colnames(correh1) <- celltype1

colnames(correh4) <- colnames(datah4)
colnames(correh3) <- colnames(datah3)
colnames(correh2) <- colnames(datah2)
colnames(correh1) <- colnames(datah1)
 
write.table(file="C_Smartseq_GSE86473.txt", correh4, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="C_Smartseq_EMATB5061.txt", correh3, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="C_CELseq_SSE85241.txt", correh2, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="C_CELseq_GSE81076.txt", correh1, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
##########

########## compute Silhouette coefficients on t-SNE coordinates
ct.fac <- factor(celltypes[allsamples])

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
par(mfrow=c(1,1), mar=c(8,8,5,3), cex.axis=3, cex.main=2, cex.lab=3)
boxplot(sils, main="", names=c("Raw", "MNN", "limma", "ComBat"),
	      lwd=4, ylab="Silhouette coefficient")#, col="Yellow", ylab="Alpha dists")
dev.off()

##################