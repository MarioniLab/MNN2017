#Code for the t-SNE plots of pancreas data sets and the Silhouette coefficients before and after batch correction by different methods (main text Figure 4).
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

load("raw_complete4DataSets.RData") #load the output of "preparedata.R"

# need to create 'results' directory
#dir.create(paste0(this.dir, "/results"))
setwd("results")
require(WGCNA)
library(scales)
library(scran)


######uncorrected data
raw.all <- cbind(datah4, datah3, datah2, datah1)
colnames(raw.all) <- c(rep(4, dim(datah4)[2]), rep(3, dim(datah3)[2]), rep(2, dim(datah2)[2]), rep(1, dim(datah1)[2]))

###### sample 1000 cells from each batch for t-SNE plotting 
set.seed(2)
samples1 <- sample(1:dim(datah1)[2], 1000) 
samples2 <- sample(1:dim(datah2)[2], 1000)
samples3 <- sample(1:dim(datah3)[2], 1000)
samples4 <- sample(1:dim(datah4)[2], 1000)
allsamples <- c(samples4, dim(datah4)[2]+samples3, dim(datah4)[2] + dim(datah3)[2] + samples2,
	   dim(datah4)[2] + dim(datah3)[2] + dim(datah2)[2] + samples1)

######## tidy cell type
celltypes <- c(celltype4, celltype3, celltype2, celltype1)
celltypes[celltypes == "pp"] <- "other"
celltypes[celltypes == "mese"] <- "delt"
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
all.dists2.unc <- as.matrix(dist(t(raw.all[hvg_genes, allsamples])))
require(Rtsne)
par(mfrow=c(1, 1))

set.seed(0)
tsne.unc <- Rtsne(all.dists2.unc, is_distance=TRUE)#, perplexity = 5)

png(file="unc4321.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.unc$Y[1:N[1],1],tsne.unc$Y[1:N[1],2], pch=3,cex=4,col=alpha(allcolors[1:N[1]],0.6),main="Uncorrected",xlim=c(-30,25),ylim=c(-20,30),xlab="tSNE 1",ylab="tSNE 2")
points(tsne.unc$Y[(N[1]+1):N[2],1],tsne.unc$Y[(N[1]+1):N[2],2], pch=18,cex=4,col=alpha(allcolors[(N[1]+1):N[2]],0.6))
points(tsne.unc$Y[(N[2]+1):N[3],1],tsne.unc$Y[(N[2]+1):N[3],2], pch=1,cex=4,col=alpha(allcolors[(N[2]+1):N[3]],0.6))
points(tsne.unc$Y[(N[3]+1):N[4],1],tsne.unc$Y[(N[3]+1):N[4],2], pch=4,cex=4,col=alpha(allcolors[(N[3]+1):N[4]],0.6))
dev.off()

### MNN batch correction
inquiry_genes <- row.names(datah4)
Xmnn <- mnnCorrect(datah4, datah3, datah2, datah1,
     	inquiry.genes=inquiry_genes,
        hvg.genes=hvg_genes, k=20, sigma=0.1,
	cos.norm=TRUE,svd.dim=0) # batch correction is throwing an error because of non-numeric values?

corre <- cbind(Xmnn$corrected[[1]],Xmnn$corrected[[2]],Xmnn$corrected[[3]],Xmnn$corrected[[4]])
all.dists2.c <- as.matrix(dist(t(corre[hvg_genes,allsamples])))

set.seed(0)
tsne.c <- Rtsne(all.dists2.c, is_distance=TRUE)#, perplexity = 5)

png(file="mnn4321_s.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.c$Y[1:N[1],1],tsne.c$Y[1:N[1],2], pch=3,cex=4,col=alpha(allcolors[1:N[1]],0.6),main="MNN corrected",xlim=c(-25,25),ylim=c(-25,20),xlab="tSNE 1",ylab="tSNE 2")
points(tsne.c$Y[(N[1]+1):N[2],1],tsne.c$Y[(N[1]+1):N[2],2], pch=18,cex=4,col=alpha(allcolors[(N[1]+1):N[2]],0.6))
points(tsne.c$Y[(N[2]+1):N[3],1],tsne.c$Y[(N[2]+1):N[3],2], pch=1,cex=4,col=alpha(allcolors[(N[2]+1):N[3]],0.6))
points(tsne.c$Y[(N[3]+1):N[4],1],tsne.c$Y[(N[3]+1):N[4],2], pch=4,cex=4,col=alpha(allcolors[(N[3]+1):N[4]],0.6))
dev.off()

### limma batch correction
library(limma)
Xlm <- removeBatchEffect(raw.all, factor(colnames(raw.all)))
all.dists2.lm <- as.matrix(dist(t(Xlm[hvg_genes,allsamples])))

set.seed(0)
tsne.lm<-Rtsne(all.dists2.lm, is_distance=TRUE)#, perplexity = 0.9)
png(file="lmfit4321.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.lm$Y[1:N[1],1],tsne.lm$Y[1:N[1],2], pch=3,cex=4,col=alpha(allcolors[1:N[1]],0.6),main="limma corrected",xlim=c(-15,20),ylim=c(-25,20),xlab="tSNE 1",ylab="tSNE 2")
points(tsne.lm$Y[(N[1]+1):N[2],1],tsne.lm$Y[(N[1]+1):N[2],2], pch=18,cex=4,col=alpha(allcolors[(N[1]+1):N[2]],0.6))
points(tsne.lm$Y[(N[2]+1):N[3],1],tsne.lm$Y[(N[2]+1):N[3],2], pch=1,cex=4,col=alpha(allcolors[(N[2]+1):N[3]],0.6))
points(tsne.lm$Y[(N[3]+1):N[4],1],tsne.lm$Y[(N[3]+1):N[4],2], pch=4,cex=4,col=alpha(allcolors[(N[3]+1):N[4]],0.6))
dev.off()

#### ComBat correction
library(sva)

Z <- colnames(raw.all)
cleandat.combat <- ComBat(raw.all,Z,mod=NULL,prior.plots = FALSE)

all.dists.combat <- as.matrix(dist(t(cleandat.combat[hvg_genes, allsamples])))
set.seed(0)
tsne.combat <- Rtsne(all.dists.combat, is_distance=TRUE)

png(file="combat4321.png", width=900,height=700)
par(mfrow=c(1,1), mar=c(6,6,4,2), cex.axis=2, cex.main=3, cex.lab=2.5)
plot(tsne.combat$Y[1:N[1], 1], 
			   tsne.combat$Y[1:N[1],2], pch=3, cex=4, 
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
write.table(file="C_CELLseq_SSE85241.txt", correh2, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="C_CELLseq_GSE81076.txt", correh1, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
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