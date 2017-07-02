#checked
library(scran)
require(WGCNA)
library(scales)

setwd("/Users/laleh/projects/deposit2/Simulations2/")
load("Sim.RData")

raw.all<-cbind(B1,B2)
clustlabels<-cbind(clust1,clust2)
colnames(raw.all)<-c(rep(1,dim(B1)[2]),rep(2,dim(B2)[2]))
require(Rtsne)
par(mfrow=c(1,1))

all.dists2.unc <- as.matrix(dist(t(raw.all)))
set.seed(0)
tsne.unc<-Rtsne(all.dists2.unc, is_distance=TRUE)#, perplexity = 0.9)
########################
N=1000
setwd("/Users/laleh/projects/deposit2/Simulations2/figs")
########################

png(file="unc.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.unc$Y[1:N,1],tsne.unc$Y[1:N,2], pch=16,cex=2.5,col=alpha((clustlabels[1:N]),0.6),main="Uncorrected",xlim=c(-35,30),ylim=c(-40,45),xlab="tSNE 1",ylab="tSNE 2")
points(tsne.unc$Y[(N+1):(2*N),1],tsne.unc$Y[(N+1):(2*N),2], pch=2,cex=3.5,col=(clustlabels[(N+1):(2*N)]))
dev.off()

#############MNN
Xmnn<-mnnCorrect(B1,B2,inquiry.genes=NULL, hvg.genes=NULL, k=20, sigma=1, cos.norm=FALSE,svd.dim=2)

corre<-cbind(Xmnn$corrected[[1]],Xmnn$corrected[[2]])
all.dists2.c <- as.matrix(dist(t(corre)))

set.seed(0)
tsne.c<-Rtsne(all.dists2.c, is_distance=TRUE)#, perplexity = 0.9)

png(file="mnn.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.c$Y[1:N,1],tsne.c$Y[1:N,2], pch=16,cex=2.5,col=alpha((clustlabels[1:N]),0.6),main="MNN corrected",xlim=c(-30,40),ylim=c(-45,45),xlab="tSNE 1",ylab="tSNE 2")
points(tsne.c$Y[(N+1):(2*N),1],tsne.c$Y[(N+1):(2*N),2], pch=2,cex=3.5,col=(clustlabels[(N+1):(2*N)]))
dev.off()

#############limma
library(limma)
Xlm <- removeBatchEffect(raw.all, factor(colnames(raw.all)))
all.dists2.lm <- as.matrix(dist(t(Xlm)))

set.seed(0)
tsne.lm<-Rtsne(all.dists2.lm, is_distance=TRUE)#, perplexity = 0.9)
png(file="lmfit.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.lm$Y[1:N,1],tsne.lm$Y[1:N,2], pch=16,cex=2.5,col=alpha((clustlabels[1:N]),0.6),main="limma corrected",xlim=c(-30,30),ylim=c(-35,45),xlab="tSNE 1",ylab="tSNE 2")
points(tsne.lm$Y[(N+1):(2*N),1],tsne.lm$Y[(N+1):(2*N),2], pch=2,cex=3.5,col=(clustlabels[(N+1):(2*N)]))
dev.off()

####ComBat correction
library(sva)
Z <- colnames(raw.all)
#mod<-model.matrix(~X)
cleandat <- ComBat(raw.all,Z,mod=NULL,prior.plots = FALSE)
all.dists.combat <- as.matrix(dist(t(cleandat)))
set.seed(0)
tsne.combat<-Rtsne(all.dists.combat, is_distance=TRUE)

png(file="combat.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.combat$Y[1:N,1],tsne.combat$Y[1:N,2], pch=16,cex=2.5,col=alpha((clustlabels[1:N]),0.6),main="ComBat corrected",xlim=c(-35,30),ylim=c(-35,30),xlab="tSNE 1",ylab="tSNE 2")
points(tsne.combat$Y[(N+1):(2*N),1],tsne.combat$Y[(N+1):(2*N),2], pch=2,cex=3.5,col=(clustlabels[(N+1):(2*N)]))
dev.off()
#############legend
png(file="leg.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.c$Y[1,1],tsne.c$Y[1,2], pch=16,cex=2.5,col=alpha((clustlabels[1]),0.6),main="MNN corrected",xlim=c(-40,40),ylim=c(-45,45),xlab="tSNE 1",ylab="tSNE 2")
#points(tsne.c$Y[(N+1):(2*N),1],tsne.c$Y[(N+1):(2*N),2], pch=2,cex=3.5,col=(clustlabels[(N+1):(2*N)]))
leg.txt1=c("Cell type 1","Celltype 2","Celltype 3")
leg.txt2=c("Batch1","Batch2")
legend("topleft", legend = leg.txt1, col = c("brown1","dark green","blue") , pch = 15,
       cex = 2.7,bty = "n")   #, trace = TRUE)
legend(x=-42,y=25, legend = leg.txt2, col = "black" , pch = c(16,2),
       cex = 2.7,bty = "n")   #, trace = TRUE)
dev.off()
