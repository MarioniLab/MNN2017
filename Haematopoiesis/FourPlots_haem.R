#Batch correction and t-SNA plots of haematopoietic data sets as in Fig.3
library(scran)
require(WGCNA)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

load("logdataFandA_all.RData") #load the output of "preparedata.R"
setwd(paste0(this.dir,"/results"))

raw.all<-cbind(logdataF3,logdataA3)
colnames(raw.all)<-c(rep(1,dim(logdataF3)[2]),rep(2,dim(logdataA3)[2]))
celltypeF<-colnames(logdataF3)
celltypeA<-substr( (colnames(logdataA3)),start = 1,stop = 3)  

##organize cell type labels
celltypes<-c(celltypeF,celltypeA)
celltypes[celltypes=="other"]<- "Unsorted"
celltypes[celltypes=="ERY"]<- "MEP"
allcolors<-celltypes
allcolors[celltypes=="MPP"]<-"blue"
allcolors[celltypes=="LTH"]<-"dodgerblue"
allcolors[celltypes=="LMP"]<-"light blue"
allcolors[celltypes=="HSP"]<-"cyan"
allcolors[celltypes=="Unsorted"]<-"grey"
allcolors[celltypes=="MEP"]<-"red" 
#allcolors[celltypes=="ERY"]<-"red"
allcolors[celltypes=="GMP"]<-"dark green"
allcolors[celltypes=="CMP"]<-"magenta"


N<-c(length(celltypeF),length(celltypeF)+length(celltypeA))

######################Raw (uncorrected data)
all.dists2.unc <- as.matrix(dist(t(raw.all)))

require(Rtsne)
par(mfrow=c(1,1))

set.seed(0)
tsne.unc<-Rtsne(all.dists2.unc, is_distance=TRUE)#, perplexity = 90)
#tsne.unc<-Rtsne(t(Xmnn$raw))#, perplexity = 0.9)

png(file="uncFA.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.unc$Y[1:N[1],1],tsne.unc$Y[1:N[1],2], pch=21,col="black",cex=2,bg=allcolors[1:N[1]],main="Uncorrected",xlim=c(-35,35),ylim=c(-35,40),xlab="tSNE 1",ylab="tSNE 2")
points(tsne.unc$Y[(N[1]+1):N[2],1],tsne.unc$Y[(N[1]+1):N[2],2], pch=1,cex=2,col=allcolors[(N[1]+1):N[2]])
dev.off()

######################## MNN corrected
row.names(logdataF3)<-row.names(logdataA3)
Xmnn<-mnnCorrect(logdataF3,logdataA3,inquiry.genes=row.names(logdataF3), hvg.genes=row.names(logdataF3), k=20, sigma=0.01, cos.norm=TRUE,svd.dim=2)

corre<-cbind(Xmnn$corrected[[1]],Xmnn$corrected[[2]])
all.dists2.c <- as.matrix(dist(t(corre)))

set.seed(0)
tsne.c<-Rtsne(all.dists2.c, is_distance=TRUE)#, perplexity = 90)
png(file="mnnFA.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.c$Y[1:N[1],1],tsne.c$Y[1:N[1],2],pch=21,col="black",cex=2.5,bg=allcolors[1:N[1]],main="MNN corrected",xlim=c(-30,35),ylim=c(-45,40),xlab="tSNE 1",ylab="tSNE 2")
points(tsne.c$Y[(N[1]+1):N[2],1],tsne.c$Y[(N[1]+1):N[2],2], pch=1,cex=2,col=allcolors[(N[1]+1):N[2]])
dev.off()

##### legend
leg.txt2<-c("MEP", "GMP" ,"CMP", "HSP","LTH" ,   "MPP" , "LMP",  "Unsorted")
forcoloringleg<-c("red", "dark green","magenta", "cyan", "dodgerblue","blue","light blue","grey" )

png(file="legFA.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.c$Y[1,1],tsne.c$Y[1,2], pch=1,cex=2,col=allcolors[1],xlim=c(-20,-12),ylim=c(0,3))
leg.txt<-unique(celltypes)
legend(x=-20, y=3, legend =leg.txt2, pch=21, cex = 2.5, pt.bg=forcoloringleg ,bty = "n")   #, trace = TRUE)
legend(x=-18,y=3, legend =leg.txt2[1:3], col = forcoloringleg[1:3], pch = 1,cex = 2.5,bty = "n",lwd = 3,lty=0) 
dev.off()

#################limma
library(limma)
Xlm <- removeBatchEffect(raw.all, factor(colnames(raw.all)))
all.dists2.lm <- as.matrix(dist(t(Xlm)))

set.seed(0)
tsne.lm<-Rtsne(all.dists2.lm, is_distance=TRUE)#, perplexity = 90)
png(file="lmfitFA.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.lm$Y[1:N[1],1],tsne.lm$Y[1:N[1],2], pch=21,col="black",cex=2,bg=allcolors[1:N[1]],main="limma corrected",xlim=c(-40,35),ylim=c(-40,45),xlab="tSNE 1",ylab="tSNE 2")
points(tsne.lm$Y[(N[1]+1):N[2],1],tsne.lm$Y[(N[1]+1):N[2],2], pch=1,cex=2,col=allcolors[(N[1]+1):N[2]])
dev.off()

####ComBat 

library(sva)
Z <- colnames(raw.all)
cleandat.combat <- ComBat(raw.all,Z,mod=NULL,prior.plots = FALSE)
all.dists.combat <- as.matrix(dist(t(cleandat.combat)))
set.seed(0)
tsne.combat<-Rtsne(all.dists.combat, is_distance=TRUE)#,perplexity = 90)


png(file="combatFA.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.combat$Y[1:N[1],1],tsne.combat$Y[1:N[1],2], pch=21,col="black",cex=2,bg=allcolors[1:N[1]],main="ComBat corrected",xlim=c(-30,40),ylim=c(-45,35),xlab="tSNE 1",ylab="tSNE 2")
points(tsne.combat$Y[(N[1]+1):N[2],1],tsne.combat$Y[(N[1]+1):N[2],2], pch=1,cex=2,col=allcolors[(N[1]+1):N[2]])
dev.off()

#save corrected data
save(list = ls(), file = "FAcorrected.RData")
