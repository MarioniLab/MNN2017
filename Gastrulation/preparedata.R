#SET WORKING DIRECTORY ######
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(scales)
require(WGCNA)
require(scran)

#######
source("../SomeFuncs/match_gene_names.R")
############

setwd(paste0(this.dir,"/meso"))

#load and read data&metadata scialdon3 2016
#download.file("http://gastrulation.stemcells.cam.ac.uk/data/counts.gz","counts.txt")
#download.file("http://gastrulation.stemcells.cam.ac.uk/data/metadataTal1.txt","metadata.txt")

data.meso<-read.table("counts.txt", head=T)
meta.data<-read.table("metadata.txt", head=T)

meso.cluster<-as.character(meta.data$cluster)
meso.stage<-as.character(meta.data$embryoStage)

row.names(meta.data)<-meta.data[,"cellName"]
#load highly variable genes (HVG)
load("high_var_genes.RData")# HVG across all cells, got from the authors
hvg.meso<-high.var.genes.all

##match the genes names to ENS
out<-matchgenenamestoENS(hvg.meso,idknown = TRUE)
hvg.meso<-unique(out$mgi_symbol)
#################scran size factor normalisation for mesodermal data
sce = newSCESet(countData = data.meso)
#filter low abundance genes
sce = sce[calcAverage(sce)>0.1,]
clusts <- quickCluster(sce, min.size=120)
#number of cells in each cluster should be at least twice that of the largest 'sizes'
min.clust = min(table(clusts))/2
new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce = computeSumFactors(sce, clusters = clusts, sizes = new_sizes)
#extract size factors
sfs = sizeFactors(sce)
#normalise data
sce = normalise(sce)
#dataF2 <- data.frame(exprs(sce))
data.meso <- t(t(counts(sce))/sizeFactors(sce))
data.meso<-log(1+data.meso)
#######download and read Wolf's data
setwd(paste0(this.dir,"/wolf"))
#download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100597&format=file&file=GSE100597%5Fcount%5Ftable%5FQC%5Ffiltered%2Etxt%2Egz","Total_table_NoTrof.bed")
data.wolf<-read.table("Total_table_NoTrof.bed",sep="\t",stringsAsFactors = F, head=T)
datainfo<-read.table("Cell_clusters.txt",sep="\t",stringsAsFactors = F, head=T)  #got from the authors
batch<-datainfo$Plate
######################
source('../../SomeFuncs/highly_var_genes.R')
genes<-row.names(data.wolf)

hvg.wolf<-find.high.var.biol.genes(data.wolf, genes, red.line=TRUE, 0.01,minMean=1, plot=TRUE)
hvg.wolf<-hvg.wolf$genes.high.var
hvg.wolf<-hvg.wolf[-which(is.na(hvg.wolf))]

#clean gene names for matching
genes<-gsub( "_.*$", "", genes)
hvg.wolf<-gsub( "_.*$", "", hvg.wolf)
#######scran size factor normalisation for wolf's data
sce = newSCESet(countData = data.wolf)
#filter low abundance genes
sce = sce[calcAverage(sce)>0.1,]
clusts <- quickCluster(sce, min.size=120)
#number of cells in each cluster should be at least twice that of the largest 'sizes'
min.clust = min(table(clusts))/2
new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce = computeSumFactors(sce, clusters = clusts, sizes = new_sizes)
#extract size factors
sfs = sizeFactors(sce)
#normalise data
sce = normalise(sce)
#dataF2 <- data.frame(exprs(sce))
data.wolf <- t(t(counts(sce))/sizeFactors(sce))
data.wolf<-log(1+data.wolf)

Stage<-datainfo$Stage
Lineage<-datainfo$Lineage

genes.wolf<-row.names(data.wolf)

## match gene names between the two data sets and pick the common genes
out<-matchgenenamestoENS(row.names(data.meso),idknown = TRUE)
dupls<-which(duplicated(out$ensembl_gene_id))
if (length(dupls)>0) {
out<-out[-dupls,]
data.meso<-data.meso[-dupls,]}
dupls<-which(duplicated(out$mgi_symbol))
if (length(dupls)>0) {
out<-out[-dupls,]
data.meso<-data.meso[-dupls,]
}
genes.wolf<-gsub( "_.*$", "", genes.wolf)
duplwolf<-which(duplicated(genes.wolf))
if (length(duplwolf)>0) {
genes.wolf<-genes.wolf[-duplwolf]
data.wolf<-data.wolf[-duplwolf,]
}
row.names(data.wolf)<-genes.wolf
commongenes<-intersect(out$mgi_symbol,row.names(data.wolf))

out2<-matchgenenamestoENS(commongenes,idknown=FALSE)
dupls<-which(duplicated(out2$mgi_symbol))
if (length(dupls)>0) {
out2<-out2[-dupls,]}
dupls<-which(duplicated(out2$ensembl_gene_id))

data.wolf.2<-data.wolf[out2$mgi_symbol,]
data.meso.2<-data.meso[match(out2$ensembl_gene_id,row.names(data.meso)),]

row.names(data.meso.2)<-row.names(data.wolf.2)
nas<-which(is.na(rowSums(data.meso.2)))
data.meso.2<-data.meso.2[-nas,]
data.wolf.2<-data.wolf.2[-nas,]

hvg.meso<-intersect(hvg.meso,row.names(data.meso.2)) #
hvg_genes<-intersect(row.names(data.meso.2),union(hvg.meso,hvg.wolf)) #common hvg_genes
save(file="../mesoandwolf.Rdata",data.meso.2,data.wolf.2,hvg_genes)
