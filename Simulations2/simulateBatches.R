
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
source('batchpartitions.R')
set.seed(10)
parts<-  batchpartitions(t(samples1),nclusters=4)  
clust1<-vector()
clust1[parts$part1]<-"dark green"#2
clust1[parts$part2]<-"blue"#3
clust1[parts$part3]<-"blue"#3
clust1[parts$part4]<-"brown1"#1
par(mfrow=c(1,1))
#plot(samples1, pch=16,cex=1.5,col=labels2colors(clust1))
plot(samples1, pch=16,cex=1.5,col=(clust1))
#clustlabels<-cbind(clust1,clust2)
############ Projection to D dimensional space
D<-100
set.seed(0)
R1D<-matrix(rnorm(D*N),  nrow=D, ncol=2)

A1<-as.matrix(samples1) %*% t(as.matrix(R1D))
#######Add normal noise
rawd1 <- matrix(rnorm(D*N), ncol=D, nrow=N)
A1<-as.data.frame(A1 + rawd1)

row.names(A1)<-c(1:nrow(A1))
colnames(A1)<-c(1:D)
###############Batch2
### Set proportions of each of the three cell types in batch 2
components <- sample(1:3,prob=c(0.65,0.3,0.05),size=N,replace=TRUE) 

########Sample N cells from Gaussian dist. according to the components (batch2)
samples2<-matrix(0,nrow=N,ncol=2)
set.seed(0)
samples2[,1] <- rnorm(n=N,mean=xmus[components],sd=xsds[components])
samples2[,2] <- rnorm(n=N,mean=ymus[components],sd=ysds[components])

########################### Get the true cluter labels for batch 2
set.seed(10)
parts<-  batchpartitions(t(samples2),nclusters=5)  
clust2<-vector()
clust2[parts$part1]<-"brown1"#1
clust2[parts$part2]<-"dark green"#2
clust2[parts$part3]<-"blue"#3
clust2[parts$part4]<-"dark green"#2
clust2[parts$part5]<-"dark green"#2
#plot(samples2, pch=16,cex=1.5,col=labels2colors(clust2))
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
#save(file="Sim.RData",B1,B2,clust1,clust2)

