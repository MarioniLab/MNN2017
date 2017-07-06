

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

N <- 1000  #no. cells

###############Batch1
### Set proportions of each of the three cell types in batch 1, and center and variance of each of the Gaussians in 2D 
components <- sample(1:3,prob=c(0.3,0.5,0.2),size=N,replace=TRUE)
xmus <- c(0,3,10)
xsds <- c(1,0.1,1)

ymus <- c(2,5,1)
ysds <- c(1,0.1,1)

########Sample N cells from Gaussian dist. according to the components (batch1)
samples1<-matrix(0,nrow=N,ncol=2)
set.seed(0)
samples1[,1] <- rnorm(n=N,mean=xmus[components],sd=xsds[components])
samples1[,2] <- rnorm(n=N,mean=ymus[components],sd=ysds[components])

########### Get the true cluter labels for batch 1
clust1<-components
clust1[components==1]<-"blue"
clust1[components==2]<-"brown1"
clust1[components==3]<-"dark green"
par(mfrow=c(1,1))
plot(samples1, pch=16,cex=1.5,col=(clust1))

############ Projection to D dimensional space
D<-100
set.seed(0)
R1D<-matrix(rnorm(D*N),  nrow=D, ncol=2)

A1<-as.matrix(samples1) %*% t(as.matrix(R1D))

#######Add normal noise
rawd1 <- matrix(rnorm(D*N), ncol=D, nrow=N)
#A1<-as.data.frame(A1 %*% T1 + rawd1)
A1<-as.data.frame(A1 + rawd1)

row.names(A1)<-c(1:nrow(A1))
colnames(A1)<-c(1:D)
###############Batch2
### Set proportions of each of the three cell types in batch 2
components <- sample(1:3,prob=c(0.3,0.5,0.2),size=N,replace=TRUE)
########Sample N cells from Gaussian dist. according to the components (batch2)

samples2<-matrix(0,nrow=N,ncol=2)
set.seed(0)
samples2[,1] <- rnorm(n=N,mean=xmus[components],sd=xsds[components])
samples2[,2] <- rnorm(n=N,mean=ymus[components],sd=ysds[components])

########################### Get the true cluter labels for batch 2
clust2<-components
clust2[components==1]<-"blue"
clust2[components==2]<-"brown1"
clust2[components==3]<-"dark green"

par(mfrow=c(1,1))
plot(samples2, pch=16,cex=1.5,col=(clust2))
############ Projection to D dimensional space

A2<-as.matrix(samples2) %*% t(as.matrix(R1D))
T2<-matrix(0,D,D)
diag(T2)<-1
####Add normal noise
rawd2 <- matrix(rep(rnorm(1*D),each=N), nrow=N, ncol=D ) + matrix(rnorm(N*D), nrow=N, ncol=D) 
A2<-as.data.frame(A2 %*% T2+  rawd2 )
row.names(A2)<-c(1:nrow(A2))
colnames(A2)<-c(1:D)
##########save batches and celltype lables 

B2<-as.data.frame(t(A2))
B1<-as.data.frame(t(A1))

save(file="easySim.RData",B1,B2,clust1,clust2)

