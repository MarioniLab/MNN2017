library(scran)
library(scater)
library(DropletUtils)

##########################################
##########################################

# Pre-processing the 68K PBMC dataset.
sce.68 <- read10xCounts("raw_data/pbmc68k/hg19/") 

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
##summary(discard)
sce.68 <- sce.68[,!discard]

# Performing normalization, breaking the problem up into smaller blocks and subclustering within them.
blocks <- rep(seq_len(10), length.out=ncol(sce.68))
clusters <- quickCluster(sce.68, min.mean=0.1, block=blocks, method="igraph", block.BPPARAM=MulticoreParam(2))
##table(clusters)
##table(clusters, blocks)

sce.68 <- computeSumFactors(sce.68, clusters=clusters, min.mean=0.1, BPPARAM=MulticoreParam(2))
##plot(sce.68$scater_qc$all$total_counts, sizeFactors(sce.68), log="xy")

saveRDS(file="sce.pbmc68k.rds", sce.68)

# Cleaning out the memory.
rm(list=ls())
gc()

##########################################
##########################################

# Pre-processing the 4K T-cell dataset.
sce.4 <- read10xCounts("raw_data/t4k/GRCh38") 
symb <- mapIds(org.Hs.eg.db, keys=rownames(sce.4), keytype="ENSEMBL", column="SYMBOL")
rowData(sce.4)$Symbol <- symb

# Adding locational annotation.
library(EnsDb.Hsapiens.v86)
loc <- mapIds(EnsDb.Hsapiens.v86, keys=rownames(sce.4), keytype="GENEID", column="SEQNAME")
rowData(sce.4)$Chr <- loc

# Brief quality control.
sce.4 <- calculateQCMetrics(sce.4, compact=TRUE, feature_controls=list(Mt=which(loc=="MT")))
lowlib <- isOutlier(sce.4$scater_qc$all$log10_total_counts, type="lower", nmads=3)
lowfeat <- isOutlier(sce.4$scater_qc$all$log10_total_features_by_counts, type="lower", nmads=3)
highmito <- isOutlier(sce.4$scater_qc$feature_control_Mt$pct_counts, type="higher", nmads=3)
discard <- lowlib | lowfeat | highmito
##summary(discard)

sce.4 <- sce.4[,!discard]

# Performing normalization.
clusters <- quickCluster(sce.4, min.mean=0.1, method="igraph")
##table(clusters)

sce.4 <- computeSumFactors(sce.4, clusters=clusters, min.mean=0.1, BPPARAM=MulticoreParam(2))
##plot(sce.4$scater_qc$all$total_counts, sizeFactors(sce.4), log="xy")

saveRDS(file="sce.t4k.rds", sce.4)
