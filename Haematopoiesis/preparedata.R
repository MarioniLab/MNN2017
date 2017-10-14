# This script prepares data for the haematopoiesis analysis.
# It involves two batches of publicly available data.

##########################################
##########################################

# Download and read the counts, metadata of Nestorowa et al. 2016
fname <- "GSE81682_HTSeq_counts.txt.gz"
if (!file.exists(fname)) { download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81682&format=file&file=GSE81682%5FHTSeq%5Fcounts%2Etxt%2Egz", fname) }
dataF <- read.table(fname, header=TRUE, row.names=1, check.names=FALSE)
dim(dataF)

fname <- "metaF.txt"
if (!file.exists(fname)) { download.file("http://blood.stemcells.cam.ac.uk/data/all_cell_types.txt", fname) }
metaF <- read.table(fname, stringsAsFactors = FALSE, header=TRUE, check.names=FALSE)
metainds <- match(colnames(dataF), rownames(metaF))
missing.meta <- is.na(metainds)
metaF <- metaF[metainds,] # This will contain NA's... which is okay, at this point, to preserve length.

# Defining the cell type based on the metadata.
metatypeF <- rep("other", nrow(metaF))
for (col in rev(colnames(metaF))) { # reverse, so earlier columns end up overwriting later ones.
    chosen <- metaF[,col]==1
    metatypeF[chosen] <- sub("[0-9]?_.*", "", col)
}
metatypeF[metatypeF=="ESLAM"] <- "HSPC"

# Filling in metadata from the cell sorting label, if metadata was missing.
metatypeF[missing.meta] <- sub("_.*", "", colnames(dataF)[missing.meta])
metatypeF[metatypeF=="LT-HSC"] <- "LTHSC"
metatypeF[metatypeF=="Prog"] <- "other"
colnames(dataF)<-metatypeF

##########################################
##########################################

# Download and read the counts and meta data of Paul et al. 2015
fname <- "umitab_Amit.txt.gz"
if (!file.exists(fname)) { download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE72857&format=file&file=GSE72857%5Fumitab%2Etxt%2Egz", fname) }
dataA <- read.table(fname, header=TRUE, row.names=1)
metaA <- read.csv2("MAP.csv",sep=",",stringsAsFactors = FALSE, head=TRUE, row.names=1)
dim(dataA)

# Only selecting cells that are in the metadata.
metainds <- match(rownames(metaA), colnames(dataA))
dataA <- dataA[,metainds]

# Organizing cell type labels.
metatypeA <- character(nrow(metaA))
metatypeA[metaA[,1]<7] <- "ERY"
metatypeA[metaA[,1]>6 & metaA[,1]<12] <- "CMP"
metatypeA[metaA[,1]>11] <- "GMP"
colnames(dataA) <- metatypeA

##########################################
##########################################

# Download list of highly variable genes identified by Nestrowa et al. 2016
fname <- "coordinates_gene_counts_flow_cytometry.txt.gz"
if (!file.exists(fname)) { download.file("http://blood.stemcells.cam.ac.uk/data/coordinates_gene_counts_flow_cytometry.txt.gz", fname) }
TFs <- read.table(fname, nrows=1, stringsAsFactors=FALSE)
features <- as.character(unlist(TFs))
features <- features[grep("ENSMUS", features)]

# Pull down IDs from BioMaRt.
library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org" )
out <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), values = features, mart = mart,filters = "ensembl_gene_id")

# Select features that are HVGs _and_ present in both data sets.
mF <- match(out$ensembl_gene_id, rownames(dataF))
mA <- pmatch(out$mgi_symbol, rownames(dataA)) # partial, due to use of concatenated gene symbols.
keep <- !is.na(mF) & !is.na(mA)

dataA2 <- dataA[mA[keep],]
dataF2 <- dataF[mF[keep],]
rownames(dataA2) <- rownames(dataF2)

# Apply library size normalization and perform log-transformation.
libF <- colSums(dataF2)
dataF2 <- t(t(dataF2)/libF)
logDataF3 <- as.matrix(log(1 + dataF2))

libA <- colSums(dataA2)
dataA2 <- t(t(dataA2)/libA)
logDataA3 <- as.matrix(log(1 + dataA2))

save(logDataA3, logDataF3, file="logdataFandA_all.RData")

###########
# END
