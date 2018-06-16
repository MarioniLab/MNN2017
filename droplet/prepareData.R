library(scran)
library(scater)
library(DropletUtils)

##########################################
##########################################

sce.68 <- read10xCounts("raw_data/filtered_matrices_mex/hg19/") 

library(org.Hs.eg.db)
symb <- mapIds(org.Hs.eg.db, keys=rownames(sce.68), keytype="ENSEMBL", column="SYMBOL")
rowData(sce.68)$Symbol <- symb

# Adding locational annotation (using a slightly off-version ensembl, but chromosome assignment shouldn't change).
library(EnsDb.Hsapiens.v86)
loc <- mapIds(EnsDb.Hsapiens.v86, keys=rownames(sce.68), keytype="GENEID", column="SEQNAME")
rowData(sce.68)$Chr <- loc

# Brief quality control.
sce.68 <- calculateQCMetrics(sce.68, compact=TRUE, feature_controls=list(Mt=which(loc=="MT")))
lowlib <- isOutlier(sce.68$scater_qc$all$log10_total_counts, type="lower", nmads=3)
lowfeat <- isOutlier(sce.68$scater_qc$all$log10_total_features_by_counts, type="lower", nmads=3)
highmito <- isOutlier(sce.68$scater_qc$feature_control_Mt$pct_counts, type="higher", nmads=3)
discard <- lowlib | lowfeat | highmito

sce.68 <- sce.68[,!discard]

# Performing normalization, breaking the problem up into smaller blocks and subclustering within them.
blocks <- rep(seq_len(10), length.out=ncol(sce.68))
clusters <- quickCluster(sce.68, min.mean=0.1, block=blocks, method="igraph", BPPARAM=MulticoreParam(2))
##table(clusters)
##table(clusters, blocks)
sce.68 <- computeSumFactors(sce.68, clusters=clusters, min.mean=0.1, BPPARAM=MulticoreParam(2))
##plot(sce.68$scater_qc$all$total_counts, sizeFactors(sce.68), log="xy")
sce.68 <- normalize(sce.68)

# Identifying (relatively) highly variable genes.
fit <- trendVar(sce.68, use.spikes=FALSE, loess.args=list(span=0.1))
##plot(fit$mean, fit$var)
##curve(fit$trend(x), col="red",add=TRUE)
dec.68 <- decomposeVar(fit=fit)

saveRDS(file="dec68.rds", dec.68)
saveRDS(file="sce68.rds", sce.68)
