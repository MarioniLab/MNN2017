
BatchEntropy <- function(dataset, batch0, L=100, M=100, k=500) {
#entropy of batch mixing
# L is the number bootstrapping times
# M is the number of randomly picked cells    
# k is the number of nearest neighbours of cell (from all batches) to check   
  
require(RANN)  
nbatches<-length(unique(batch0))

entropy<-matrix(0,L,1)
set.seed(0) 
for (boot in 1:L) {
  bootsamples<-sample(1:nrow(dataset),M)
  W21<-nn2(dataset,query=dataset[bootsamples,],k)
  
  for (i in 1:length(bootsamples)){
    
    for (j in 1:nbatches) {
    xi<-max(1,sum(batch0[W21$nn.idx[i,]]==j))
    entropy[boot]<-entropy[boot]+xi*log(xi)
    }
  }
}

return( (-1)*entropy/length(bootsamples) )
}


###usage:
#batch0<-c(rep(4,dim(datah4)[2]),rep(3,dim(datah3)[2]),rep(2,dim(datah2)[2]),rep(1,dim(datah1)[2]))
#entrop.raw<-BatchEntropy(pca.raw$x[,1:2],batch0)

##repeat function with each pca.x and the corresponding entrop.x before the next line and making the boxplots
#En<-cbind(entrop.raw,entrop.mnn,entrop.lm,entrop.com)
#boxplot(En, names=c("RAW","MNN","limma","ComBat"),ylab="Entorpy of batch mixing")

#png(file="pcabatch_entropy0p1.png",width=900,height=700)
#par(mfrow=c(1,1),mar=c(8,8,5,3),cex.axis=3,cex.main=2,cex.lab=3)
#boxplot(En,main="",names=c("Raw","MNN","limma","ComBat"),lwd=4,ylab="Entorpy of batch mixing")#,col="Yellow",ylab="Alpha dists")
#dev.off()
