#Read pancreas data and meta data and set of highly variable genes + preprocessing of data befor batch correction; match gene names accross data sets.  
require(WGCNA)

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

##read data files
datah1<-read.table("CELseq/GSE81076-norm.tsv",sep="\t",stringsAsFactors = F, head=T)
datah2<-read.table("CELseq/GSE85241-norm.tsv",sep="\t",stringsAsFactors = F, head=T)
datah3<-read.table("Smartseq2/E-MTAB-5061-norm.tsv",sep="\t",stringsAsFactors = F, head=T)
datah4<-read.table("Smartseq2/GSE86473-norm.tsv",sep="\t",stringsAsFactors = F, head=T)

HVG1<-read.table("CELseq/GSE81076-HVG.tsv")
HVG2<-read.table("CELseq/GSE85241-HVG.tsv")
HVG3<-read.table("Smartseq2/E-MTAB-5061-HVG.tsv")
HVG4<-read.table("Smartseq2/GSE86473-HVG.tsv")


meta1<-read.table("CELseq/GSE81076_marker_metadata.tsv",sep="\t",stringsAsFactors = F, head=T)
meta2<-read.table("CELseq/GSE85241_marker_metadata.tsv",sep="\t",stringsAsFactors = F, head=T)
meta3<-read.table("Smartseq2/E-MTAB-5061_marker_metadata.tsv",sep="\t",stringsAsFactors = F, head=T)
meta4<-read.table("Smartseq2/GSE86473_marker_metadata.tsv",sep="\t",stringsAsFactors = F, head=T)

#gene names
genes1<-as.character(datah1[,1190]) #last columns is the genes name
duplrows<-which(duplicated(genes1))
datah1<-datah1[,-1190]
row.names(datah1)<-genes1

genes2<-as.character(datah2[,2406]) #last columns is the genes name
datah2<-datah2[,-2406]
datah2<-datah2[,-1]  # column1 does not have celltype label
row.names(datah2)<-genes2

genes3<-as.character(datah3[,2189]) #last columns is the genes name
datah3<-datah3[,-2189]
row.names(datah3)<-genes3

genes4<-as.character(datah4[,1451]) #last columns is the genes name
datah4<-datah4[,-1451]
row.names(datah4)<-genes4

############process cell type labels

celltype1<- vector("character", length = dim(datah1)[2])
for (i in 1:dim(datah1)[2]) {
  ci<-which(meta1$Sample==colnames(datah1)[i])
  celltype1[i] <-substr( tolower(meta1$CellType[ci]),start = 1,stop = 4)
}
celltype2<- vector("character", length = dim(datah2)[2])
for (i in 1:dim(datah2)[2]) {
  ci<-which(meta2$Sample==colnames(datah2)[i])
  celltype2[i] <-substr( tolower(meta2$CellType[ci]),start = 1,stop = 4)
}

celltype3<- vector("character", length = dim(datah3)[2])
for (i in 1:dim(datah3)[2]) {
  ci<-which(meta3$Sample==colnames(datah3)[i])
  celltype3[i] <-substr( tolower(meta3$CellType[ci]),start = 1,stop = 4)
}

celltype4<- vector("character", length = dim(datah4)[2])
for (i in 1:dim(datah4)[2]) {
  ci<-which(meta4$Sample==colnames(datah4)[i])
  celltype4[i] <-substr( tolower(meta4$CellType[ci]),start = 1,stop = 4)
}
###prepare batches with identical row names (matched gene names) 
inquiry_genes<- intersect(genes1,intersect(genes2,intersect(genes3,genes4)))
datah1<-datah1[inquiry_genes,]
datah2<-datah2[inquiry_genes,]
datah3<-datah3[inquiry_genes,]
datah4<-datah4[inquiry_genes,]
###find of set of highly variable gene names which are present in all data sets 
HVG<-unique(c(as.character(unlist(HVG1)),as.character(unlist(HVG2)),as.character(unlist(HVG3)),as.character(unlist(HVG4)))) #union of highly variable genes
common_genes<- intersect(genes1,intersect(genes2,intersect(genes3,intersect(HVG,genes4))))
hvg_genes<-common_genes

##further cleaning
badcol<-which(celltype2=="")
if (length(badcol)>0){
datah2<-datah2[,-badcol] #datah2 has a nonlabled column
celltype2<-celltype2[-badcol]
}

narow<-which(is.na(datah1[,2])) #there was a narow!
if (length(narow) >0) {
  datah1<-datah1[-narow,]
}

narow<-which(is.na(datah2[,2])) #there was a narow!
if (length(narow) >0) {
  datah2<-datah2[-narow,]
}

narow<-which(is.na(datah3[,2])) #there was a narow!
if (length(narow) >0) {
  datah3<-datah3[-narow,]
}

narow<-which(is.na(datah4[,2])) #there was a narow!
if (length(narow) >0) {
  datah4<-datah4[-narow,]
}
########
save(datah1, datah2,datah3,datah4,celltype1,celltype2,celltype3,celltype4,inquiry_genes,hvg_genes,file="raw_complete4DataSets.RData")

