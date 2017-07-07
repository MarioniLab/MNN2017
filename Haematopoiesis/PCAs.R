#Generates PCA plots of haematopoietic data sets as in Suppl. Fig.2
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
load("FAcorrected.RData") #load the output of "FourPlots_haem.R"

batch0<-c(rep(1,length(logdataF3)),rep(2,length(logdataA3)))

##select the shared cell types between the two data sets
selectid<-which(celltypes=="GMP" | celltypes=="CMP" | celltypes=="MEP")
batch<-batch0[selectid]
selecttype<-celltypes[selectid]

#set colors and symbols for plotting
forcol<-selecttype
forcol[selecttype=="GMP"]<-"dark green"
forcol[selecttype=="CMP"]<-"magenta"
forcol[selecttype=="MEP"]<-"red"

forpch<-batch
forpch[batch==1]<-16
forpch[batch==2]<-1

########Raw
  
data<-raw.all[,selectid]
N<-c(sum(batch==1),dim(data)[2])

pca.data <- prcomp(t(data), center=TRUE)
pca.data$x<-cosine.norm(pca.data$x)  #for orthonormal PC basis

pca.data$x[(pca.data$x[,2]< (-0.1)),2]<-(-0.1)
png(file="pca_raw.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(pca.data$x[1:N[1],1],pca.data$x[1:N[1],2], pch=21,col="black",cex=2,bg=forcol[1:N[1]],main="Uncorrected",xlab="PC1",ylab="PC2"
     ,xlim=c(-0.07,0.023),ylim=c(-0.06,0.06))
points(pca.data$x[(N[1]+1):N[2],1],pca.data$x[(N[1]+1):N[2],2], pch=1,cex=2,col=forcol[(N[1]+1):N[2]])
dev.off()
##########################MNN
data<-corre[,selectid]
pca.data <- prcomp(t(data), center=TRUE)
pca.data$x<-cosine.norm(pca.data$x) #for orthonormal PC basis

png(file="pca_mnn.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(pca.data$x[1:N[1],1],pca.data$x[1:N[1],2], pch=21,col="black",cex=2,bg=forcol[1:N[1]],main="MNN corrected",xlab="PC1",ylab="PC2"
     ,xlim=c(-0.06,0.02),ylim=c(-0.04,0.08))
points(pca.data$x[(N[1]+1):N[2],1],pca.data$x[(N[1]+1):N[2],2], pch=1,cex=2,col=forcol[(N[1]+1):N[2]])
dev.off()

##########################Xlm
data<-Xlm[,selectid]
pca.data <- prcomp(t(data), center=TRUE)
pca.data$x<-cosine.norm(pca.data$x) #for orthonormal PC basis
pca.data$x[(pca.data$x[,2]< -(0.1)),2]<--(0.1)

png(file="pca_lm.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(pca.data$x[1:N[1],1],pca.data$x[1:N[1],2], pch=21,col="black",cex=2,bg=forcol[1:N[1]],main="limma corrected",xlab="PC1",ylab="PC2"
     ,xlim=c(-0.025,0.065),ylim=c(-0.1,0.05))
points(pca.data$x[(N[1]+1):N[2],1],pca.data$x[(N[1]+1):N[2],2], pch=1,cex=2,col=forcol[(N[1]+1):N[2]])
dev.off()

##########################combat
data<-cleandat.combat[,selectid]
pca.data <- prcomp(t(data), center=TRUE)
pca.data$x<-cosine.norm(pca.data$x) #for orthonormal PC basis
pca.data$x[(pca.data$x[,2]> (0.1)),2]<-(0.1)

png(file="pca_com.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(pca.data$x[1:N[1],1],pca.data$x[1:N[1],2], pch=21,col="black",cex=2,bg=forcol[1:N[1]],main="ComBat corrected",xlab="PC1",ylab="PC2"
     ,xlim=c(-0.02,0.1),ylim=c(-0.1,0.07))
points(pca.data$x[(N[1]+1):N[2],1],pca.data$x[(N[1]+1):N[2],2], pch=1,cex=2,col=forcol[(N[1]+1):N[2]])
dev.off()

