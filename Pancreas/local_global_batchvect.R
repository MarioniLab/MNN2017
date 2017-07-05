#Code for the comparison of locally variable batch vector correction by MNN and global batch vector correction by MNN (Suppl. Fig. 6).

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

load("raw_complete4DataSets.RData") #load the output of "preparedata.R"
setwd(paste0(this.dir,"/results/loc-glob-tests/"))

require(WGCNA)
library(scales)
library(scran)

raw.all<-cbind(datah4,datah3,datah2,datah1)
colnames(raw.all)<-c(rep(2,dim(datah4)[2]),rep(3,dim(datah3)[2]),rep(4,dim(datah2)[2]),rep(1,dim(datah1)[2]))


set.seed(2)
samples1<-sample(1:dim(datah1)[2],1000) 
samples2<-sample(1:dim(datah2)[2],1000)
samples3<-sample(1:dim(datah3)[2],1000)
samples4<-sample(1:dim(datah4)[2],1000)
allsamples<-c(samples4, dim(datah4)[2]+samples3, dim(datah4)[2]+dim(datah3)[2]+samples2, dim(datah4)[2]+dim(datah3)[2]+dim(datah2)[2]+samples1)


celltypes<-c(celltype4,celltype3,celltype2,celltype1)
celltypes[celltypes=="pp"]="other"
celltypes[celltypes=="mese"]="delt"
celltypes[celltypes=="co-e"]="other"
celltypes[celltypes=="endo"]="other"
celltypes[celltypes=="epsi"]="other"
celltypes[celltypes=="mast"]="other"
celltypes[celltypes=="mhc "]="other"
celltypes[celltypes=="uncl"]="other"
celltypes[celltypes=="psc "]="other"

celltypes[celltypes=="alph"]="alpha"
celltypes[celltypes=="delt"]="delta"
celltypes[celltypes=="gamm"]="gamma"
celltypes[celltypes=="psc "]="other"

allcolors<-labels2colors(celltypes[allsamples])
allcolors[allcolors=="red"]<-"deeppink"
allcolors[allcolors=="yellow"]<-"orange1"#"darkgoldenrod1"

N<-c(1000,2000,3000,4000)

####### MNN correction; allowing local batch vector (sigma >0)

inquiry_genes=row.names(datah4)
Xmnn1<-mnnCorrect(datah4,datah3,datah2,datah1,inquiry.genes=inquiry_genes, hvg.genes=hvg_genes, k=20, sigma=0.1, cos.norm=TRUE,svd.dim=0)
corre1<-cbind(Xmnn1$corrected[[1]],Xmnn1$corrected[[2]],Xmnn1$corrected[[3]],Xmnn1$corrected[[4]])
all.dists2.c1 <- as.matrix(dist(t(corre1[hvg_genes,allsamples])))

require(Rtsne)
set.seed(0)
tsne.c1<-Rtsne(all.dists2.c1, is_distance=TRUE)#, perplexity = 5)
#tsne.c<-Rtsne(t(Xmnn$corrected))
png(file="mnn4321_local.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.c1$Y[1:N[1],1],tsne.c1$Y[1:N[1],2], pch=3,cex=4,col=alpha(allcolors[1:N[1]],0.6),main="MNN corrected",xlim=c(-20,25),ylim=c(-20,17),xlab="tSNE 1",ylab="tSNE 2")
points(tsne.c1$Y[(N[1]+1):N[2],1],tsne.c1$Y[(N[1]+1):N[2],2], pch=18,cex=4,col=alpha(allcolors[(N[1]+1):N[2]],0.6))
points(tsne.c1$Y[(N[2]+1):N[3],1],tsne.c1$Y[(N[2]+1):N[3],2], pch=1,cex=4,col=alpha(allcolors[(N[2]+1):N[3]],0.6))
points(tsne.c1$Y[(N[3]+1):N[4],1],tsne.c1$Y[(N[3]+1):N[4],2], pch=4,cex=4,col=alpha(allcolors[(N[3]+1):N[4]],0.6))
dev.off()

####### MNN correction; single global batch vector (sigma=0)
Xmnn2<-mnnCorrect(datah4,datah3,datah2,datah1,inquiry.genes=inquiry_genes, hvg.genes=hvg_genes, k=20, sigma=0, cos.norm=TRUE,svd.dim=0)
corre2<-cbind(Xmnn2$corrected[[1]],Xmnn2$corrected[[2]],Xmnn2$corrected[[3]],Xmnn2$corrected[[4]])
all.dists2.c2 <- as.matrix(dist(t(corre2[hvg_genes,allsamples])))

set.seed(0)
tsne.c2<-Rtsne(all.dists2.c2, is_distance=TRUE)#, perplexity = 5)
png(file="mnn4321_global.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.c2$Y[1:N[1],1],tsne.c2$Y[1:N[1],2], pch=3,cex=4,col=alpha(allcolors[1:N[1]],0.6),main="MNN corrected",xlim=c(-15,20),ylim=c(-25,25),xlab="tSNE 1",ylab="tSNE 2")
points(tsne.c2$Y[(N[1]+1):N[2],1],tsne.c2$Y[(N[1]+1):N[2],2], pch=18,cex=4,col=alpha(allcolors[(N[1]+1):N[2]],0.6))
points(tsne.c2$Y[(N[2]+1):N[3],1],tsne.c2$Y[(N[2]+1):N[3],2], pch=1,cex=4,col=alpha(allcolors[(N[2]+1):N[3]],0.6))
points(tsne.c2$Y[(N[3]+1):N[4],1],tsne.c2$Y[(N[3]+1):N[4],2], pch=4,cex=4,col=alpha(allcolors[(N[3]+1):N[4]],0.6))
dev.off()

############# compute Silhouette coefficients on t-SNE coordinates
library(kBET)

ct.fac <- factor(celltypes[allsamples])
dd<-all.dists2.c1
score_sil <- (cluster::silhouette(as.numeric(ct.fac), dd))
sil_c1<-score_sil[,3]

dd<-all.dists2.c2
score_sil <- (cluster::silhouette(as.numeric(ct.fac), dd))
sil_c2<-score_sil[,3]


###### box plot of Silhouette coefficients
sils<-cbind(sil_c1,sil_c2)

png(file="sils_local_global.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(8,8,5,3),cex.axis=3,cex.main=2,cex.lab=3)
boxplot(sils,main="",names=c("Local","Public"),lwd=4,ylab="Silhouette coefficient")#,col="Yellow",ylab="Alpha dists")
dev.off()

################## Compute PCAs for local and global
source("../../../SomeFuncs/cosine-norm.R") #need cosine-norm to orthonormalize PCA basis

batch0<-c(rep(4,dim(datah4)[2]),rep(3,dim(datah3)[2]),rep(2,dim(datah2)[2]),rep(1,dim(datah1)[2]))
library(RColorBrewer)
colors4<-brewer.pal(4,"Set3")

#colors4<-c("#2ca25f","#8856a7","#43a2ca","#e34a33")

data<-corre1[hvg_genes,] #local

pca.mnn1 <- prcomp(t(data), center=TRUE)#,scale. = TRUE)
pca.mnn1$x<-cosine.norm(pca.mnn1$x)

data<-corre2[hvg_genes,] #global

pca.mnn2 <- prcomp(t(data), center=TRUE)#,scale. = TRUE)
pca.mnn2$x<-cosine.norm(pca.mnn2$x)

#randomize order of cells for plotting
set.seed(2)
ix<-sample(1:nrow(pca.mnn1$x),nrow(pca.mnn1$x))
  
png(file="pca_mnn_local.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
plot(pca.mnn1$x[ix,1],pca.mnn1$x[ix,2],col=alpha(colors4[batch0[ix]],0.6),main="MNN corrected",xlab="PC 1",ylab="PC 2",cex=1.5)#,xlim=c(-0.015,0.025))
dev.off()

png(file="pca_mnn_global.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
plot(pca.mnn2$x[ix,1],pca.mnn2$x[ix,2],col=alpha(colors4[batch0[ix]],0.6),main="MNN corrected",xlab="PC 1",ylab="PC 2",cex=1.5)#,xlim=c(-0.015,0.025))
dev.off()
#### calculate entropy of batch mixings
source("../../../SomeFuncs/BatchMixingEntropy.R")
entrop.loc<-BatchEntropy(pca.mnn1$x[,1:2],batch0)
entrop.glob<-BatchEntropy(pca.mnn2$x[,1:2],batch0)

##box plots
En<-cbind(entrop.loc,entrop.glob)
boxplot(En, names=c("Local","Global"),ylab="Entorpy of batch mixing")

png(file="entropy_loc_global.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(8,8,5,3),cex.axis=3,cex.main=2,cex.lab=3)
boxplot(En,main="",names=c("Local","Global"),lwd=4,ylab="Entorpy of batch mixing")
dev.off()
