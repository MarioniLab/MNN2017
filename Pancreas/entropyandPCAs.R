#Code for PCA plots and entropy of batch mixing for pancreas data sets (Suppl. Fig. 3)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

load("raw_complete4DataSets.RData")  #load the output of "preparedata.R"
setwd(paste0(this.dir,"/results"))
load("completedata_correcteds.RData") #load the output of "FourPlots_panc.R"

library(scales)
library(RANN)

source("../../SomeFuncs/cosine-norm.R")

batch<-c(rep(4,dim(datah4)[2]),rep(3,dim(datah3)[2]),rep(2,dim(datah2)[2]),rep(1,dim(datah1)[2]))

#set plotting colors
library(RColorBrewer)
colors4<-brewer.pal(4,"Set3")
#set plotting symbols
forpch<-batch
forpch[batch==4]<-3
forpch[batch==3]<-18
forpch[batch==2]<-1
forpch[batch==1]<-4

##PCA of uncorrected data
data<-raw.all[hvg_genes,]#
pca.raw <- prcomp(t(data), center=TRUE)#,scale. = TRUE)
pca.raw$x<-cosine.norm(pca.raw$x)  #for orthonormal pca basis. Basis normalization is not done in prcomp

#randomize order of cells for plotting
set.seed(2)
ix<-sample(1:nrow(pca.raw$x),nrow(pca.raw$x)) 

png(file="pca_raw.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
plot(pca.raw$x[ix,1],pca.raw$x[ix,2],col=alpha(colors4[batch[ix]],0.6),pch=16,main="Uncorrected",xlab="PC 1",ylab="PC 2",cex=1.5)#,xlim=c(-0.03,0.02))
dev.off()

######PCA of MNN corrected
corre<-cbind(Xmnn$corrected[[1]],Xmnn$corrected[[2]],Xmnn$corrected[[3]],Xmnn$corrected[[4]])
data<-corre[hvg_genes,]  #PCA on log-scale output
pca.mnn <- prcomp(t(data), center=TRUE)#,scale. = TRUE)
pca.mnn$x<-cosine.norm(pca.mnn$x) 

png(file="pca_mnn.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
plot(pca.mnn$x[ix,1],pca.mnn$x[ix,2],col=alpha(colors4[batch[ix]],0.6),main="MNN corrected",xlab="PC 1",ylab="PC 2",cex=1.5)#,xlim=c(-0.015,0.025))
dev.off()

######PCA of limma corrected
data<-Xlm[hvg_genes,]
pca.lm <- prcomp(t(data), center=TRUE)#,scale. = TRUE)
pca.lm$x<-cosine.norm(pca.lm$x) 

png(file="pca_lm.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
plot(pca.lm$x[ix,1],pca.lm$x[ix,2],col=alpha(colors4[batch[ix]],0.6),main="limma corrected",xlab="PC 1",ylab="PC 2",cex=1.5)#,xlim=c(-0.025,0.04))
dev.off()

######PCA of ComBat corrected
data<-cleandat.combat[hvg_genes,]
pca.com <- prcomp(t(data), center=TRUE)#,scale. = TRUE)
pca.com$x<-cosine.norm(pca.com$x) 

png(file="pca_com.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
plot(pca.com$x[ix,1],pca.com$x[ix,2],col=alpha(colors4[batch[ix]],0.6),main="ComBat corrected",xlab="PC 1",ylab="PC 2",cex=1.5)#,xlim=c(-0.02,0.035))
dev.off()

##save(file="panc4321_all.RData",list=ls())
########### the legend
png(file="pancpca_leg.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
plot(1,2,col=alpha(colors4[batch],0.6),main="limma corrected",xlab="PC 1",ylab="PC 2",cex=1.5,xlim=c(-0.025,0.04))
legend("topleft", legend =c("CEL-Seq","CEL-Seq2","SMART-Seq2 (I)","SMART-Seq (II)"), col = unique(batch), pch =16,cex = 2.5,bty = "n")
dev.off()
###################Entropy of batch mixing for uncorrected and batch corrected data using different methods

source("../../SomeFuncs/BatchMixingEntropy.R")
#batch<-c(rep(4,dim(datah4)[2]),rep(3,dim(datah3)[2]),rep(2,dim(datah2)[2]),rep(1,dim(datah1)[2]))
entrop.raw<-BatchEntropy(pca.raw$x[,1:2],batch)
entrop.mnn<-BatchEntropy(pca.mnn$x[,1:2],batch)
entrop.lm<-BatchEntropy(pca.lm$x[,1:2],batch)
entrop.com<-BatchEntropy(pca.com$x[,1:2],batch)

En<-cbind(entrop.raw,entrop.mnn,entrop.lm,entrop.com)

#box plot of entropy of batch mixing 
png(file="pcabatch_entropy.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(8,8,5,3),cex.axis=3,cex.main=2,cex.lab=3)
boxplot(En,main="",names=c("Raw","MNN","limma","ComBat"),lwd=4,ylab="Entorpy of batch mixing")
dev.off()

