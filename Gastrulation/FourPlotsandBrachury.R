#checked
library(scran)
library(scales)
require(WGCNA)
require(Rtsne)

#SET WORKING DIRECTORY ######
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

load("mesoandwolf.Rdata") #load the output of preparedata.R

setwd(paste0(this.dir,"/meso"))

##organize cell type * cell stage labels
meta.data<-read.table("metadata.txt", head=T)
meso.cluster<-as.character(meta.data$cluster)
meso.stage<-as.character(meta.data$embryoStage)

meso.stage[meso.stage=="HF"]<-"E7.75"
meso.stage[meso.stage=="NP"]<-"E7.5"
meso.stage[meso.stage=="PS"]<-"E7.0"

select.meso<-which(meso.stage=="E6.5" | meso.stage=="E7.0")
meso.stage<-meso.stage[select.meso]
meso.cluster<-meso.cluster[select.meso]
data.meso.2<-data.meso.2[,select.meso]
excl.meso<-which(meso.cluster=="pink" | meso.cluster=="magenta" | meso.cluster=="red")
data.meso.2<-data.meso.2[,-excl.meso]
meso.cluster<-meso.cluster[-excl.meso]
meso.stage<-meso.stage[-excl.meso]
new.labs.meso<-paste0(substr(meso.cluster,start=1,stop=3),meso.stage)

setwd("../wolf")
datainfo<-read.table("Cell_clusters.txt",sep="\t",stringsAsFactors = F, head=T)

batch<-datainfo$Plate

select.wolf<-which(datainfo$Stage=="E5.5" | datainfo$Stage=="E6.5"
                  |datainfo$Stage=="E6.75")

wolf.stage<-datainfo$Stage[select.wolf]
wolf.lineage<-datainfo$Lineage[select.wolf]
data.wolf.2<-data.wolf.2[,select.wolf]

excl.wolf<-which(substr(wolf.lineage,start=1,stop=2)=="VE")
data.wolf.2<-data.wolf.2[,-excl.wolf]
wolf.lineage<-wolf.lineage[-excl.wolf]
wolf.stage<-wolf.stage[-excl.wolf]

new.labs.wolf<-paste0(substr(wolf.lineage,start=1,stop=3),wolf.stage)
new.labs.all<-c(new.labs.wolf,new.labs.meso)

###################
datawolf<-data.wolf.2
datameso<-data.meso.2

raw.all<-cbind(datawolf,datameso)
colnames(raw.all)<-c(rep(2,dim(datawolf)[2]),rep(1,dim(datameso)[2]))

batch<-colnames(raw.all)

new.labs.all<-c(new.labs.wolf,new.labs.meso)
celltypes<-new.labs.all

#celltypes[substr(celltypes,start=1,stop=2)=="VE"]<-"VE"
#celltypes[substr(celltypes,start=1,stop=3)=="mag"]<-"VE"
#celltypes[substr(celltypes,start=1,stop=3)=="pin"]<-"Ext.end"
celltypes[substr(celltypes,start=1,stop=3)=="tur"]<-"epi-psE6.5"
celltypes[substr(celltypes,start=1,stop=3)=="yel"]<-"blood.prE7.0"
celltypes[substr(celltypes,start=1,stop=3)=="blu"]<-"meso.prE7.0"
celltypes[substr(celltypes,start=1,stop=3)=="red"]<-"Endoth"
celltypes[substr(celltypes,start=1,stop=3)=="dar"]<-"post.mesE7.0"

allcolors<-labels2colors(celltypes)

allcolors[allcolors=="yellow"]<-"darkorange"
allcolors[allcolors=="greenyellow"]<-"darkgreen"

N<-c(dim(datawolf)[2],dim(datawolf)[2]+dim(datameso)[2])
setwd("/Users/laleh/projects/deposit2/Gast/results")
#################Raw
all.dists2.unc <- as.matrix(dist(t(raw.all)))
set.seed(0)
tsne.unc<-Rtsne(all.dists2.unc, is_distance=TRUE, perplexity = 90)
#tsne.unc<-Rtsne(t(Xmnn$raw))#, perplexity = 0.9)
png(file="uncGast.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.unc$Y[1:N[1],1],tsne.unc$Y[1:N[1],2], pch=17,cex=2.5,col=alpha(allcolors[1:N[1]],0.6),main="Uncorrected"
     ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-7,7),ylim=c(-12,12))
points(tsne.unc$Y[(N[1]+1):N[2],1],tsne.unc$Y[(N[1]+1):N[2],2], pch=21,cex=2.5,col="black",bg=(allcolors[(N[1]+1):N[2]]))
dev.off()
##############MNN
#Xmnn<-mnnCorrect(datawolf,datameso, k=20, sigma=0.1, cos.norm=TRUE,svd.dim=1)
Xmnn<-mnnCorrect(datawolf,datameso,inquiry.genes=row.names(datawolf), hvg.genes=hvg_genes, k=20, sigma=0.1, cos.norm=TRUE,svd.dim=2)
corre<-cbind(Xmnn$corrected[[1]],Xmnn$corrected[[2]])
all.dists2.c <- as.matrix(dist(t(corre[hvg_genes,])))

set.seed(0)
tsne.c<-Rtsne(all.dists2.c, is_distance=TRUE, perplexity = 90)
#tsne.c<-Rtsne(t(Xmnn$corrected))
# png(file="MNNGast.png",width=900,height=700)
# par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
# plot(tsne.c$Y[1:N[1],1],tsne.c$Y[1:N[1],2],pch=17,cex=2.5,col=alpha(allcolors[1:N[1]],0.6),main="MNN corrected"
#      ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-10,10),ylim=c(-10,10))
# points(tsne.c$Y[(N[1]+1):N[2],1],tsne.c$Y[(N[1]+1):N[2],2], pch=21,cex=2.5,col="black",bg=allcolors[(N[1]+1):N[2]])
# dev.off()

png(file="MNNGast.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.c$Y[(N[1]+1):N[2],1],tsne.c$Y[(N[1]+1):N[2],2], pch=21,cex=2.5,col="black",bg=allcolors[(N[1]+1):N[2]],main="MNN corrected"
     ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-10,10),ylim=c(-10,10))
points(tsne.c$Y[1:N[1],1],tsne.c$Y[1:N[1],2],pch=17,cex=2.5,col=alpha(allcolors[1:N[1]],0.6) )
dev.off()


###############limma
library(limma)
Xlm <- removeBatchEffect(raw.all, factor(colnames(raw.all)))
all.dists2.lm <- as.matrix(dist(t(Xlm)))

set.seed(0)
tsne.lm<-Rtsne(all.dists2.lm, is_distance=TRUE, perplexity = 90)
#tsne.unc<-Rtsne(t(Xmnn$raw))#, perplexity = 0.9)
png(file="limmaGast.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.lm$Y[1:N[1],1],tsne.lm$Y[1:N[1],2], pch=17,cex=2.5,col=alpha(allcolors[1:N[1]],0.6),main="limma corrected"
     ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-6,5),ylim=c(-7,9))
points(tsne.lm$Y[(N[1]+1):N[2],1],tsne.lm$Y[(N[1]+1):N[2],2], pch=21,cex=2.5,col="black",bg=(allcolors[(N[1]+1):N[2]]))
dev.off()

###########ComBat
library(sva)
Z <- colnames(raw.all)
cleandat.combat <- ComBat(raw.all+0.01,Z,mod=NULL,prior.plots = FALSE)

all.dists.combat <- as.matrix(dist(t(cleandat.combat)))
set.seed(0)
tsne.combat<-Rtsne(all.dists.combat, is_distance=TRUE)

png(file="CombatGast.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.combat$Y[1:N[1],1],tsne.combat$Y[1:N[1],2], pch=17,cex=2.5,col=alpha(allcolors[1:N[1]],0.6),main="ComBat corrected"
     ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-9,9),ylim=c(-20,12))
points(tsne.combat$Y[(N[1]+1):N[2],1],tsne.combat$Y[(N[1]+1):N[2],2], pch=21,cex=2.5,col="black",bg=(allcolors[(N[1]+1):N[2]]))
dev.off()

#########
celltypes[celltypes=="PSE6.5"]="psE6.5"
leg.txt<-unique(celltypes)
#leg.txt<-unique(celltypes)
forcoloringleg<-unique(allcolors)


png(file="leg_detailed.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
#plot(tsne.c$Y[1,1],tsne.c$Y[1,2], pch=3,cex=4,col=alpha(allcolors[1],0.6),main="legend",xlim=c(-20,30),ylim=c(-13,13),xlab="tSNE 1",ylab="tSNE 2")
plot(1,2, pch=3,cex=4,col=alpha(allcolors[1],0.6),main="legend",xlim=c(-20,30),ylim=c(-13,13),xlab="tSNE 1",ylab="tSNE 2")
forleg<-table(celltypes,allcolors)
legend("bottomright", "(x,y)", legend = leg.txt, col =unique(allcolors) , pch = 17,cex = 2.5,bty = "n",lwd = 3,lty=0)   #, trace = TRUE)
legend("bottomleft", "(x,y)", legend = leg.txt ,pt.bg=forcoloringleg, pch = 21,cex = 2.5,bty = "n",lwd = 3,lty=0)   #, trace = TRUE)
#legend(x=7,y=15, legend = c("CEL-Seq 1","CEL-Seq 2","SMART-Seq 1","SMART-Seq 2"), col ="black" , pch = c(4,1,18,3),cex = 2.5,bty = "n")   #, trace = TRUE)
dev.off()

##############################
library(RColorBrewer)
rbPal <- colorRampPalette(c('blue','yellow','red'))

datCol1 <- rbPal(100)[as.numeric(cut(raw.all[which(row.names(raw.all)=="Mesp1"),],breaks = 100))]
datCol2 <- rbPal(100)[as.numeric(cut(raw.all[which(row.names(raw.all)=="T"),],breaks = 100))]
datCol3 <- rbPal(100)[as.numeric(cut(raw.all[which(row.names(raw.all)=="Mettl7a3"),],breaks = 100))]

png(file="Raw_Mesp1.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.unc$Y[1:N[1],1],tsne.unc$Y[1:N[1],2], pch=17,cex=2.5,col=alpha(datCol1[1:N[1]],0.6),main=" Mesp1 uncorrected"
     ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-7,7),ylim=c(-12,12))
points(tsne.unc$Y[(N[1]+1):N[2],1],tsne.unc$Y[(N[1]+1):N[2],2], pch=20,cex=3,col=alpha(datCol1[(N[1]+1):N[2]],0.6))
dev.off()

png(file="Raw_T.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.unc$Y[1:N[1],1],tsne.unc$Y[1:N[1],2], pch=17,cex=2.5,col=alpha(datCol2[1:N[1]],0.6),main="T uncorrected"
     ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-7,7),ylim=c(-12,12))
points(tsne.unc$Y[(N[1]+1):N[2],1],tsne.unc$Y[(N[1]+1):N[2],2], pch=20,cex=3,col=alpha(datCol2[(N[1]+1):N[2]],0.6))
dev.off()

png(file="Raw_Mettl7a3.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.unc$Y[1:N[1],1],tsne.unc$Y[1:N[1],2], pch=17,cex=2.5,col=alpha(datCol3[1:N[1]],0.6),main="Mettl7a3 uncorrected"
     ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-7,7),ylim=c(-12,12))
points(tsne.unc$Y[(N[1]+1):N[2],1],tsne.unc$Y[(N[1]+1):N[2],2], pch=20,cex=3,col=alpha(datCol3[(N[1]+1):N[2]],0.6))
dev.off()

png(file="MNN_Mesp1.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.c$Y[1:N[1],1],tsne.c$Y[1:N[1],2], pch=17,cex=2.5,col=alpha(datCol1[1:N[1]],0.6),main=" Mesp1 MNN corrected"
     ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-10,10),ylim=c(-10,10))
points(tsne.c$Y[(N[1]+1):N[2],1],tsne.c$Y[(N[1]+1):N[2],2], pch=20,cex=3,col=alpha(datCol1[(N[1]+1):N[2]],0.6))
dev.off()

png(file="MNN_T.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.c$Y[1:N[1],1],tsne.c$Y[1:N[1],2], pch=17,cex=2.5,col=alpha(datCol2[1:N[1]],0.6),main="T MNN corrected"
     ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-10,10),ylim=c(-10,10))
points(tsne.c$Y[(N[1]+1):N[2],1],tsne.c$Y[(N[1]+1):N[2],2], pch=20,cex=3,col=alpha(datCol2[(N[1]+1):N[2]],0.6))
dev.off()

png(file="MNN_Mettl7a3.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.c$Y[1:N[1],1],tsne.c$Y[1:N[1],2], pch=17,cex=2.5,col=alpha(datCol3[1:N[1]],0.6),main="Mettl7a3 MNN corrected"
     ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-10,10),ylim=c(-10,10))
points(tsne.c$Y[(N[1]+1):N[2],1],tsne.c$Y[(N[1]+1):N[2],2], pch=20,cex=3,col=alpha(datCol3[(N[1]+1):N[2]],0.6))
dev.off()

png(file="limma_Mesp1.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.lm$Y[1:N[1],1],tsne.lm$Y[1:N[1],2], pch=17,cex=2.5,col=alpha(datCol1[1:N[1]],0.6),main=" Mesp1 limma corrected"
     ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-6,5),ylim=c(-7,9))
points(tsne.lm$Y[(N[1]+1):N[2],1],tsne.lm$Y[(N[1]+1):N[2],2], pch=20,cex=3,col=alpha(datCol1[(N[1]+1):N[2]],0.6))
dev.off()

png(file="limma_T.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.lm$Y[,1],tsne.lm$Y[,2], pch=17,cex=2.5,col=alpha(datCol2[1:N[1]],0.6),main="T limma corrected"
     ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-6,5),ylim=c(-7,9))
points(tsne.lm$Y[(N[1]+1):N[2],1],tsne.lm$Y[(N[1]+1):N[2],2], pch=20,cex=3,col=alpha(datCol2[(N[1]+1):N[2]],0.6))
dev.off()

png(file="limma_Mettl7a3.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.lm$Y[1:N[1],1],tsne.lm$Y[1:N[1],2], pch=17,cex=2.5,col=alpha(datCol3[1:N[1]],0.6),main="Mettl7a3 limma corrected"
     ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-6,5),ylim=c(-7,9))
points(tsne.lm$Y[(N[1]+1):N[2],1],tsne.lm$Y[(N[1]+1):N[2],2], pch=20,cex=3,col=alpha(datCol3[(N[1]+1):N[2]],0.6))
dev.off()

png(file="Combat_Mesp1.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.combat$Y[1:N[1],1],tsne.combat$Y[1:N[1],2], pch=17,cex=2.5,col=alpha(datCol1[1:N[1]],0.6),main="Mesp1 ComBat corrected"
     ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-9,9),ylim=c(-20,12))
points(tsne.combat$Y[(N[1]+1):N[2],1],tsne.combat$Y[(N[1]+1):N[2],2], pch=20,cex=3,col=alpha(datCol1[(N[1]+1):N[2]],0.6))
dev.off()

png(file="Combat_T.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.combat$Y[1:N[1],1],tsne.combat$Y[1:N[1],2], pch=17,cex=2.5,col=alpha(datCol2[1:N[1]],0.6),main="T ComBat corrected"
     ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-9,9),ylim=c(-20,12))
points(tsne.combat$Y[(N[1]+1):N[2],1],tsne.combat$Y[(N[1]+1):N[2],2], pch=20,cex=3,col=alpha(datCol2[(N[1]+1):N[2]],0.6))
dev.off()

png(file="Combat_Mettl7a3.png",width=900,height=700)
par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
plot(tsne.combat$Y[1:N[1],1],tsne.combat$Y[1:N[1],2], pch=17,cex=2.5,col=alpha(datCol3[1:N[1]],0.6),main="Mettl7a3 ComBat corrected"
     ,xlab="tSNE 1",ylab="tSNE 2",xlim=c(-9,9),ylim=c(-20,12))
points(tsne.combat$Y[(N[1]+1):N[2],1],tsne.combat$Y[(N[1]+1):N[2],2], pch=20,cex=3,col=alpha(datCol3[(N[1]+1):N[2]],0.6))
dev.off()

#############
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

png(file="Combat_leg.png",width=900,height=700)
color.bar(colorRampPalette(c("blue", "yellow", "red"))(100), -1)
dev.off()

#########
