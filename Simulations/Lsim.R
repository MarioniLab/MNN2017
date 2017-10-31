# L shaped data simulations (supplement figures)
#of non-orthogonal batch effect + small angle rotation

cx=0
cy=0
n=200

set.seed(0)
hor.x=runif(n, min = 0, max = 2)
hor.y=cy+rnorm(n, mean = 0, sd = 0.1)
hor<-cbind(hor.x,hor.y)


ver.x=cx+rnorm(n, mean = 0, sd = 0.1)
ver.y=runif(n, min = 0, max = 2)
ver<-cbind(ver.x,ver.y)
L1 <- rbind(hor,ver)

plot(L1,xlim=c(-1,3),ylim=c(-1,2.5),col="red")

hor.x=runif(n, min = 0, max = 2)
hor.y=cy+rnorm(n, mean = 0, sd = 0.1)
hor<-cbind(hor.x,hor.y)

ver.x=cx+rnorm(n, mean = 0, sd = 0.1)
ver.y=runif(n, min = 0, max = 2)
ver<-cbind(ver.x,ver.y)

#simulation of non-orth. batch vector
shift=matrix(rep(c(0.5,-0.5) ,each=2*n),nrow=2*n,ncol=2)
L2 <- rbind(hor,ver)+shift
#L2 <- L1+shift

points(L2)

Xmnn<-mnnCorrect2(t(L1),t(L2),k=10, sigma=0.1, cos.norm=FALSE,svd.dim=0,k.clara = 0,withQC = FALSE,varCare=TRUE)
corre<-t(do.call(cbind,Xmnn$corrected))

par(mfrow=c(1,2),mar=c(6,6,4,2),cex.axis=1,cex.main=1,cex.lab=1)
plot(L1,xlim=c(-1.5,3),ylim=c(-1,2.5),col="red", main="Raw", xlab="x",ylab="y")
points(L2)
plot(corre[1:400,],col="red",main="MNN corrected",xlim=c(-1.5,3),ylim=c(-1,2.5), xlab="x",ylab="y")
points(corre[401:800,])
########Simulation of rotation
theta=20*pi/180
rotM=matrix(0,nrow=2,ncol=2)
rotM[1,]=c(cos(theta),-sin(theta))
rotM[2,]=c(sin(theta),cos(theta))

L2r=t(rotM%*% t(L1))
colnames(L2r)<-colnames(L2)
Xmnn<-mnnCorrect2(t(L1),t(L2r),k=10, sigma=0.1, cos.norm=FALSE,svd.dim=0,k.clara = 0,withQC = FALSE,varCare=TRUE)
corre<-t(do.call(cbind,Xmnn$corrected))

par(mfrow=c(1,2),mar=c(6,6,4,2),cex.axis=1,cex.main=1,cex.lab=1)
plot(L1,xlim=c(-1.5,3),ylim=c(-1,2.5),col="red", main="Raw", xlab="x",ylab="y")
points(L2r)
plot(corre[1:400,],col="red",main="MNN corrected",xlim=c(-1.5,3),ylim=c(-1,2.5), xlab="x",ylab="y")
points(corre[401:800,])