

setwd("/Users/laleh/projects/deposit2/Simulations2")


N <- 1000

components <- sample(1:3,prob=c(0.3,0.5,0.2),size=N,replace=TRUE)
xmus <- c(0,3,10)
xsds <- c(1,0.1,1)

ymus <- c(2,5,1)
ysds <- c(1,0.1,1)

samples1<-matrix(0,nrow=N,ncol=2)
set.seed(0)
samples1[,1] <- rnorm(n=N,mean=xmus[components],sd=xsds[components])
samples1[,2] <- rnorm(n=N,mean=ymus[components],sd=ysds[components])

###########True cluster1
source('batchpartitions.R')
set.seed(10)
parts<-  batchpartitions(t(samples1),"kmeans",nclusters=4)  
clust1<-vector()
clust1[parts$part1]<-"dark green"#2
clust1[parts$part2]<-"blue"#1
clust1[parts$part3]<-"blue"#3
clust1[parts$part4]<-"brown1"#3
par(mfrow=c(1,1))
plot(samples1, pch=16,cex=1.5,col=(clust1))

############
D<-100
set.seed(0)
R1D<-matrix(rnorm(D*N),  nrow=D, ncol=2)

A1<-as.matrix(samples1) %*% t(as.matrix(R1D))

#######batch effects
rawd1 <- matrix(rnorm(D*N), ncol=D, nrow=N)
#A1<-as.data.frame(A1 %*% T1 + rawd1)
A1<-as.data.frame(A1 + rawd1)

row.names(A1)<-c(1:nrow(A1))
colnames(A1)<-c(1:D)
###
components <- sample(1:3,prob=c(0.3,0.5,0.2),size=N,replace=TRUE)

samples2<-matrix(0,nrow=N,ncol=2)
set.seed(0)
samples2[,1] <- rnorm(n=N,mean=xmus[components],sd=xsds[components])
samples2[,2] <- rnorm(n=N,mean=ymus[components],sd=ysds[components])

################True cluster2
set.seed(10)
parts<-  batchpartitions(t(samples2),"kmeans",nclusters=5)  
clust2<-vector()
clust2[parts$part1]<-"blue"#3
clust2[parts$part2]<-"dark green"#2
clust2[parts$part3]<-"brown1"#1
clust2[parts$part4]<-"blue"#3
clust2[parts$part5]<-"dark green"#2
plot(samples2, pch=16,cex=1.5,col=(clust2))
#######
A2<-as.matrix(samples2) %*% t(as.matrix(R1D))

T2<-matrix(0,D,D)
diag(T2)<-1

rawd2 <- matrix(rep(rnorm(1*D),each=N), nrow=N, ncol=D ) + matrix(rnorm(N*D), nrow=N, ncol=D) 

A2<-as.data.frame(A2 %*% T2+  rawd2 )
row.names(A2)<-c(1:nrow(A2))
colnames(A2)<-c(1:D)
##########

B2<-as.data.frame(t(A2))
B1<-as.data.frame(t(A1))

#######
#save(file="easySim.RData",B1,B2,clust1,clust2)

