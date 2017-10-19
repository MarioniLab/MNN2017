# Read pancreas data and meta data and set of highly variable genes + preprocessing of data befor batch correction; match gene names accross data sets.  
# this script needs to start from the raw counts matrix (or even download it directly if possible)

#### This script is superfluous as we can work directly from the flat files

#this.dir <- dirname(parent.frame(2)$ofile)
#setwd(this.dir)

# read data files
datah1 <- read.table("Pancreas/Data/GSE81076_SFnorm.tsv", sep="\t", stringsAsFactors=FALSE, head=TRUE)
datah2 <- read.table("Pancreas/Data/GSE85241_SFnorm.tsv", sep="\t", stringsAsFactors=FALSE, head=TRUE)
datah3 <- read.table("Pancreas/Data/GSE86473_SFnorm.tsv", sep="\t", stringsAsFactors=FALSE, head=TRUE)
datah4 <- read.table("Pancreas/Data/E-MTAB-5061_SFnorm.tsv", sep="\t", stringsAsFactors=FALSE, head=TRUE)

# read in highly variable gene files
HVG1 <- read.table("Pancreas/Data/GSE81076-HVG.tsv", sep="\t", h=TRUE, stringsAsFactors=FALSE)
HVG2 <- read.table("Pancreas/Data/GSE85241-HVG.tsv", sep="\t", h=TRUE, stringsAsFactors=FALSE)
HVG3 <- read.table("Pancreas/Data/GSE86473-HVG.tsv", sep="\t", h=TRUE, stringsAsFactors=FALSE)
HVG4 <- read.table("Pancreas/Data/E-MTAB-5061-HVG.tsv", sep="\t", h=TRUE, stringsAsFactors=FALSE)

# read in meta data with cell type labels
# only cells that have passed previous QC steps are included in these meta data files
meta1 <- read.table("Pancreas/Data/GSE81076_marker_metadata.tsv", sep="\t", stringsAsFactors = FALSE, head=TRUE)
meta2 <- read.table("Pancreas/Data/GSE85241_marker_metadata.tsv", sep="\t", stringsAsFactors = FALSE, head=TRUE)
meta3 <- read.table("Pancreas/Data/GSE86473_metadata.tsv", sep="\t", stringsAsFactors = FALSE, head=TRUE)
meta4 <- read.table("Pancreas/Data/E-MTAB-5061_metadata.tsv", sep="\t", stringsAsFactors = FALSE, head=TRUE)

# subset metadata to include cells for which data are available, i.e. passed QC, keep gene ID column
datah1 <- datah1[, c(intersect(colnames(datah1), meta1$Sample), "gene_id")]
datah2 <- datah2[, c(intersect(colnames(datah2), meta2$Sample), "gene_id")]
datah3 <- datah3[, c(intersect(colnames(datah3), meta3$Sample), "gene_id")]
datah4 <- datah4[, c(intersect(colnames(datah4), meta4$Sample), "gene_id")]

# last columns is the genes name
# remove any duplicate names and set the rownames to geneIDs
genes1 <- as.character(datah1[, dim(datah1)[2]])
duplrows <- which(duplicated(genes1))
# trim off the last column that contains the gene IDs
datah1 <- datah1[, -dim(datah1)[2]]
row.names(datah1) <- genes1

# last columns is the genes name
genes2 <- as.character(datah2[, dim(datah2)[2]])
datah2 <- datah2[, -dim(datah2)[2]]
datah2 <- datah2[, -1]  # column1 does not have celltype label
row.names(datah2) <- genes2

# last columns is the genes name
genes3 <- as.character(datah3[, dim(datah3)[2]])
datah3 <- datah3[, -dim(datah3)[2]]
row.names(datah3) <- genes3

# last columns is the genes name
genes4 <- as.character(datah4[, dim(datah4)[2]])
datah4 <- datah4[, -dim(datah4)[2]]
row.names(datah4) <- genes4

# subset the first 4 characters of the cell type labels for plotting
celltype1 <- vector("character", length=dim(datah1)[2])
for (i in 1:dim(datah1)[2]) {
  ci <- which(meta1$Sample==colnames(datah1)[i])
  celltype1[i] <- substr(tolower(meta1$CellType[ci]), start=1, stop=4)
}

celltype2 <- vector("character", length=dim(datah2)[2])
for (i in 1:dim(datah2)[2]) {
  ci <- which(meta2$Sample == colnames(datah2)[i])
  celltype2[i] <- substr(tolower(meta2$CellType[ci]), start=1, stop=4)
}

celltype3 <- vector("character", length=dim(datah3)[2])
for (i in 1:dim(datah3)[2]) {
  ci <- which(meta3$Sample == colnames(datah3)[i])
  celltype3[i] <- substr(tolower(meta3$CellType[ci]), start=1, stop=4)
}

celltype4 <- vector("character", length=dim(datah4)[2])
for (i in 1:dim(datah4)[2]) {
  ci <- which(meta4$Sample == colnames(datah4)[i])
  celltype4[i] <- substr(tolower(meta4$CellType[ci]), start=1, stop=4)
}

# test without these inquiry genes to make sure it all still works
# prepare batches with identical row names (matched gene names) 
#inquiry_genes <- intersect(genes1, intersect(genes2, intersect(genes3, genes4)))
#datah1 <- datah1[inquiry_genes, ]
#datah2 <- datah2[inquiry_genes, ]
#datah3 <- datah3[inquiry_genes, ]
#datah4 <- datah4[inquiry_genes, ]

# find of set of highly variable gene names which are present in all data sets 
# this should be the union of highly variable genes, not the intersection?
HVG <- unique(c(HVG1$gene_id, HVG2$gene_id, HVG3$gene_id, HVG4$gene_id))

common_genes <- intersect(genes1, intersect(genes2, intersect(genes3, genes4)))
# hvg_genes <- common_genes

## further cleaning
# remove cells with no cell label
badcol <- which(celltype2 == "")
if (length(badcol) > 0){
datah2 <- datah2[, -badcol] #datah2 has a nonlabled column
celltype2 <- celltype2[-badcol]
}

# remove NA rows or just set values to 0?
narow <- which(is.na(datah1[, 2])) #there was a narow!
if (length(narow) > 0) {
  datah1 <- datah1[-narow, ]
}

narow <- which(is.na(datah2[, 2])) #there was a narow!
if (length(narow) > 0) {
  datah2 <- datah2[-narow, ]
}

narow <- which(is.na(datah3[, 2])) #there was a narow!
if (length(narow) > 0) {
  datah3 <- datah3[-narow, ]
}

narow <- which(is.na(datah4[, 2])) #there was a narow!
if (length(narow) > 0) {
  datah4 <- datah4[-narow,]
}

# only write out the matrices that conform based on gene IDs
datah1 <- datah1[common_genes, ]
datah2 <- datah2[common_genes, ]
datah3 <- datah3[common_genes, ]
datah4 <- datah4[common_genes, ] 

# save the R data objects of the normalized expression matrices, cell labels and highly variable gene
save(datah1, datah2, datah3, datah4,
	     celltype1, celltype2, celltype3, celltype4,
	     HVG,
	     file="Pancreas/raw_complete4DataSets.RData")
