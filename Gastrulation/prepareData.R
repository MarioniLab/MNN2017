require(scran)
library(scater)

####################################
# Create a directory for the mesoderm data.

dir.create("meso", showWarnings=FALSE)
count.file <- file.path("meso", "counts.txt")

if (!file.exists(count.file)) { 
    download.file("http://gastrulation.stemcells.cam.ac.uk/data/counts.gz","meso/counts.gz")
    system("cd meso; tar -xvf counts.gz") # god knows why they named it like this.
}
data.meso <- read.table(count.file, header=TRUE, row.names=1)

# Compute size factors.
sce.meso <- SingleCellExperiment(list(counts=as.matrix(data.meso)))

high.ab <- calcAverage(sce.meso)>1
clusts <- quickCluster(sce.meso, method="igraph", subset.row=high.ab, min.size=100)
sce.meso <- computeSumFactors(sce.meso, clusters = clusts, subset.row=high.ab)

pdf(file.path("meso", "norm.pdf"))
plot(colSums(counts(sce.meso)), sizeFactors(sce.meso))
dev.off()

sce.meso <- normalise(sce.meso)

# Discovering highly variable genes using the Brennecke method.
meso.var <- technicalCV2(sce.meso, spike.type=NA, min.bio.disp=0)
hvg.meso <- meso.var$FDR <= 0.05 & !is.na(meso.var$FDR)

pdf(file.path("meso", "hvg.pdf"))
plot(meso.var$mean, meso.var$cv2, log="xy", pch=16, cex=0.5)
o <- order(meso.var$mean)
lines(meso.var$mean[o], meso.var$trend[o], col="red")
dev.off()

####################################
# Create a directory for Wolf's data.

dir.create("wolf", showWarnings=FALSE)
count.file <- "wolf/GSE100597_count_table_QC_filtered.txt.gz"
if (!file.exists(count.file)) { 
    download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100597&format=file&file=GSE100597%5Fcount%5Ftable%5FQC%5Ffiltered%2Etxt%2Egz", count.file)
}

data.wolf <- read.table(count.file, sep="\t", stringsAsFactors = FALSE, header = TRUE, row.names=1)

# Compute size factors.
sce.wolf <- SingleCellExperiment(list(counts=as.matrix(data.wolf)))

high.ab <- calcAverage(sce.wolf)>1
clusts <- quickCluster(sce.wolf, method="igraph", subset.row=high.ab)
sce.wolf <- computeSumFactors(sce.wolf, clusters = clusts, subset.row=high.ab)

pdf(file.path("wolf", "norm.pdf"))
plot(colSums(counts(sce.wolf)), sizeFactors(sce.wolf))
dev.off()

sce.wolf <- normalise(sce.wolf)

# Discovering highly variable genes using the Brennecke method.
wolf.var <- technicalCV2(sce.wolf, spike.type=NA, min.bio.disp=0)
hvg.wolf <- wolf.var$FDR <= 0.05 & !is.na(wolf.var$FDR)

pdf(file.path("wolf", "hvg.pdf"))
plot(wolf.var$mean, wolf.var$cv2, log="xy", pch=16, cex=0.5)
o <- order(wolf.var$mean)
lines(wolf.var$mean[o], wolf.var$trend[o], col="red")
dev.off()

####################################
# Match gene names to ensembl IDs in the Wolf data set.

data.meso <- exprs(sce.meso)
data.wolf <- exprs(sce.wolf)

library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org")
cleaned.names <- gsub( "_.*$", "", rownames(data.wolf))
wolf.anno <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), values = cleaned.names, mart = mart, filters = "mgi_symbol")
wolf.ens <- wolf.anno$ensembl_gene_id[match(cleaned.names, wolf.anno$mgi_symbol)] # need to get back to input order!
wolf.ens <- ifelse(is.na(wolf.ens), cleaned.names, wolf.ens)
rownames(data.wolf) <- wolf.ens

# Picking genes that are present in both data sets, _and_ are also HVGs.
in.both <- intersect(wolf.ens, rownames(data.meso)) 
any.hvg <- union(rownames(data.wolf)[hvg.wolf], rownames(data.meso)[hvg.meso])
any.hvg <- intersect(in.both, any.hvg) 

data.meso <- data.meso[in.both,]
data.wolf <- data.wolf[in.both,]
save(file="mesoandwolf.Rdata", data.meso, data.wolf, any.hvg)
