#checked
#SET WORKING DIRECTORY ######

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

############naive size factor corection
c_by_sizefactor<-function(X){
  G<-dim(X)[1]
  n<-dim(X)[2]
  cellnorm<-colSums(X)
  cellnorm.nonzero<-cellnorm
  cellnorm.nonzero[which(cellnorm==0)]<-1
  Y<-X/matrix(rep(cellnorm.nonzero,each=G),G,n)
}
###########
library(utils)
#download and read the counts and meta data of Nestorowa et al. 2016
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81682&format=file&file=GSE81682%5FHTSeq%5Fcounts%2Etxt%2Egz", "GSE81682_HTSeq_counts.txt")
dataF<-read.table("GSE81682_HTSeq_counts.txt", head=T)
featuresF<-dataF[,1]
dataF<-dataF[,-1]
row.names(dataF)<-featuresF
dim(dataF)

download.file("http://blood.stemcells.cam.ac.uk/data/all_cell_types.txt", "metaF.txt")
metaF<-read.table("metaF.txt",stringsAsFactors = F, head=T)
metainds<-match(row.names(metaF),colnames(dataF))
narow<-which(is.na(metainds))
metainds<-metainds[-narow]
metaF<-metaF[-narow,]

metatypeF<-vector(mode="character",length=dim(metaF)[1])
for (i in 1:dim(metaF)[1]) {
whichcol<-which(metaF[i,]==1)
  if (length(whichcol) >0) {
  metatypeF[i]<-substr(colnames(metaF)[whichcol[1]],start=1,stop=3) 
  }
  else{metatypeF[i]<-"other"
  }
}

metatypeF[metatypeF=="ESL"]<-"HSP"
leftouts<- setdiff(c(1:ncol(dataF)),metainds )
metainds<-c(metainds,leftouts)
metatypeF<-c(metatypeF,substr(colnames(dataF)[leftouts], start=1,stop=3))
dataF<-dataF[,metainds]
metatypeF[metatypeF=="LT."]<-"LTH"
metatypeF[metatypeF=="Pro"]<-"other"
colnames(dataF)<-metatypeF

#download and read the counts and meta data of Paul et al. 2015
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE72857&format=file&file=GSE72857%5Fumitab%2Etxt%2Egz","umitab_Amit.txt")
dataA<-read.table("umitab_Amit.txt", head=T)
metaA<-read.csv2("MAP.csv",sep=",",stringsAsFactors = F, head=T)
row.names(metaA)<-metaA[,1]
dim(dataA)

metainds<-match(row.names(metaA),colnames(dataA))
dataA<-dataA[,metainds]

#organize cell type labels
metatypeA<-vector(mode="character",length=dim(metaA)[1])
metatypeA[metaA[,2]<7]<-"ERY"
metatypeA[metaA[,2]>6 &metaA[,2]<12]<-"CMP"
metatypeA[metaA[,2]>11]<-"GMP"
Amit.celltypes<-metatypeA

colnames(dataA)<-unlist(Amit.celltypes)

##download list of highly variable genes identified by Nestrowa et al. 2016
download.file("http://blood.stemcells.cam.ac.uk/data/coordinates_gene_counts_flow_cytometry.txt.gz","coordinates_gene_counts_flow_cytometry.txt")
TFs <- read.table("coordinates_gene_counts_flow_cytometry.txt",
                  nrows =  1)
###### organize and match gene names between the two data sets
source("../SomeFuncs/match_gene_names.R")
features<-as.character(unlist(TFs))
features<-features[28:4800]
out<-matchgenenamestoENS(features,idknown = TRUE)
dataF2<-dataF[out$ensembl_gene_id,]
dataA2<-dataA[out$mgi_symbol,]

colnames(dataA2)<-colnames(dataA)
##tidy up 
 narow<-which(is.na(dataF2[,2])) #there was a narow!
if (length(narow) >0) {
 dataF2<-dataF2[-narow,]
 dataA2<-dataA2[-narow,]
}
 narow<-which(is.na(dataA2[,2]))  #there were several narow!
 if (length(narow) >0) {
 dataF2<-dataF2[-narow,]
 dataA2<-dataA2[-narow,]
 }
####size factor corrections + log transform
dataF2<-(c_by_sizefactor(dataF2)) 
dataA2<-(c_by_sizefactor(dataA2))  # 
logdataF3<-(log(1+dataF2))
logdataA3<-(log(1+dataA2))

save(list=c("logdataA3","logdataF3"),file="logdataFandA_all.RData")
