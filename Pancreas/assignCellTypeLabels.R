# this script assigns cell-type labels to cells based on the original implementation
# described in the manuscript for each study
library(Rtsne)
library(ggplot2)
library(cluster)

############
# GSE81076 #
############
# read in the normalized expression data
gse81076.norm <- read.table("Pancreas/Data/GSE81076_SFnorm.tsv",
                            h=TRUE, sep="\t", stringsAsFactors=FALSE)

# remove gene id column for hvg discovery, last column
rownames(gse81076.norm) <- gse81076.norm$gene_id

gse81076.hvg_df <- read.table("Pancreas/Data/GSE81076-HVG.tsv", 
                              sep="\t", h=TRUE, stringsAsFactors=FALSE)

gse81076.hvg <- gse81076.norm[gse81076.norm$gene_id %in% gse81076.hvg_df$gene_id,
                              1:(dim(gse81076.norm)[2]-1)]

gse81076.meta <- read.table("Pancreas/Data/GSE81076_metadata.tsv",
                            sep="\t", h=TRUE, stringsAsFactors=FALSE)

# get low dimensional embedding of genes by tSNE
gse81076.tsne <- Rtsne(t(gse81076.hvg), perplexity=30)
gse81076.map <- data.frame(gse81076.tsne$Y)
colnames(gse81076.map) <- c("Dim1", "Dim2")
gse81076.map$Sample <- colnames(gse81076.hvg)

# merge metadata and clustering info
gse81076.uber <- merge(gse81076.map, gse81076.meta,
                       by='Sample')

marker_genes <- rownames(gse81076.norm)[grepl(rownames(gse81076.norm), pattern="(^GCG)|(^KRT19)|(^INS)|
                                            |(^SST)|(^PPY)|(^PRSS1)|(^COL1A1)|(^GHRL)|(^ESAM)") ]

marker_exprs <- gse81076.norm[marker_genes, 1:(dim(gse81076.norm)[2]-1)]

set.seed(42)
kmed <- pam(t(gse81076.hvg), 9)
panc.kmed <- data.frame(cbind(kmed$clustering))
colnames(panc.kmed) <- "Kmediods"
panc.kmed$Sample <- rownames(panc.kmed)
panc.kmed$Kmediods <- as.factor(panc.kmed$Kmediods)
panc.meta <- merge(gse81076.uber, panc.kmed, by='Sample')

marker_df <- data.frame(t(marker_exprs))
marker_df$Sample <- rownames(marker_df)

# can we cluster the cells on just the marker genes?
# hierarchical clustering doesn't work so well...

mark.dim <- dim(marker_df)
marker.kmed <- pam(marker_df[, 1:mark.dim[2]-1], 9)
marker.cluster <- data.frame(cbind(marker.kmed$clustering))
colnames(marker.cluster) <- c("markClust")
marker.cluster$Sample <- rownames(marker.cluster)
marker.cluster$markClust <- as.factor(marker.cluster$markClust)

marker_merge <- merge(panc.meta, marker_df, by='Sample')
marker_uber <- merge(marker_merge, marker.cluster, by='Sample')

# uncomment the following plots to see the overlay of each marker gene
# on the tSNE.  This illustrates the cell type label inference.
# # mediods 8 = alpha cells
# ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=GCG,
#                          shape=Kmediods)) + 
#   geom_point(size=1) + theme_classic() +
#   scale_colour_continuous(low='white', high='red') + 
#   scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))
# 
# # mediod 6 = beta cells, INS highest
# ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=INS,
#                          shape=Kmediods)) + 
#   geom_point(size=1) + theme_classic() +
#   scale_colour_continuous(low='white', high='red') + 
#   scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))
# 
# # mediod 7 = delta cells, SST highest
# ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=SST,
#                          shape=Kmediods)) + 
#   geom_point(size=1) + theme_classic() +
#   scale_colour_continuous(low='white', high='red') +
#   scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))
# 
# # mediod 7 & 8 = PP cells, PPY highest
# ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=PPY,
#                          shape=Kmediods)) + 
#   geom_point(size=1) + theme_classic() +
#   scale_colour_continuous(low='white', high='red') +
#   scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))
# 
# # mediod 1 & 2 = acinar cells, PRSS1 highest
# ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=PRSS1,
#                          shape=Kmediods)) + 
#   geom_point(size=1) + theme_classic() +
#   scale_colour_continuous(low='white', high='red') +
#   scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))
# 
# # mediod 4 = ductal cells, KRT19 highest
# ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=KRT19,
#                          shape=Kmediods)) + 
#   geom_point(size=1) + theme_classic() +
#   scale_colour_continuous(low='white', high='red') +
#   scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))
# 
# # small pop of mediod 9 = mesenchyme cells, COL1A1 highest, same as PP cells
# ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=COL1A1,
#                          shape=Kmediods)) + 
#   geom_point(size=1) + theme_classic() +
#   scale_colour_continuous(low='white', high='red') +
#   scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))

table(marker_uber$Kmediods, marker_uber$markClust)
marker_uber$CellType <- ""
marker_uber$CellType[marker_uber$Kmediods %in% c(7, 8) & marker_uber$PPY >= 5] <- "PP"
marker_uber$CellType[marker_uber$Kmediods == 6] <- "Beta"
marker_uber$CellType[marker_uber$Kmediods %in% c(8)] <- "Alpha"
marker_uber$CellType[marker_uber$Kmediods == 7] <- "Delta"
marker_uber$CellType[marker_uber$Kmediods %in% c(1, 2)] <- "Acinar"
marker_uber$CellType[marker_uber$Kmediods %in% c(3, 4) & marker_uber$KRT19 >= 1.5] <- "Ductal"
marker_uber$CellType[marker_uber$Kmediods %in% c(9) & marker_uber$COL1A1 >= 2] <- "Mesenchyme"
table(marker_uber$CellType)

write.table(marker_uber,
            file="Pancreas/Data/GSE81076_marker_metadata.tsv",
            sep="\t", quote=F, row.names=F, col.names=T)

############
# GSE85241 #
############
rm(list=ls())
gc()

# read in the normalized expression data
gse85241.norm <- read.table("Pancreas/Data/GSE85241_SFnorm.tsv",
                            h=TRUE, sep="\t", stringsAsFactors=FALSE)

# remove gene id column for hvg discovery, last column
rownames(gse85241.norm) <- gse85241.norm$gene_id

gse85241.hvg_df <- read.table("Pancreas/Data/GSE85241-HVG.tsv", 
                              sep="\t", h=TRUE, stringsAsFactors=FALSE)

gse85241.hvg <- gse85241.norm[gse85241.norm$gene_id %in% gse85241.hvg_df$gene_id,
                              1:(dim(gse85241.norm)[2]-1)]

gse85241.meta <- read.table("Pancreas/Data/GSE85241_metadata.tsv",
                            sep="\t", h=TRUE, stringsAsFactors=FALSE)

# get low dimensional embedding of genes by tSNE
gse85241.tsne <- Rtsne(t(gse85241.hvg), perplexity=30)
gse85241.map <- data.frame(gse85241.tsne$Y)
colnames(gse85241.map) <- c("Dim1", "Dim2")
gse85241.map$Sample <- colnames(gse85241.hvg)

# merge metadata and clustering info
gse85241.uber <- merge(gse85241.map, gse85241.meta,
                       by='Sample')

marker_genes <- rownames(gse85241.norm)[grepl(rownames(gse85241.norm), pattern="(^GCG)|(^KRT19)|(^INS)|
                                              |(^SST)|(^PPY)|(^PRSS1)|(^COL1A1)|(^GHRL)|(^ESAM)") ]

marker_exprs <- gse85241.norm[marker_genes, 1:(dim(gse85241.norm)[2]-1)]

set.seed(42)
kmed <- pam(t(gse85241.hvg), 9)
panc.kmed <- data.frame(cbind(kmed$clustering))
colnames(panc.kmed) <- "Kmediods"
panc.kmed$Sample <- rownames(panc.kmed)
panc.kmed$Kmediods <- as.factor(panc.kmed$Kmediods)
panc.meta <- merge(gse85241.uber, panc.kmed, by='Sample')

marker_df <- data.frame(t(marker_exprs))
marker_df$Sample <- rownames(marker_df)

mark.dim <- dim(marker_df)
marker.kmed <- pam(marker_df[, 1:mark.dim[2]-1], 9)
marker.cluster <- data.frame(cbind(marker.kmed$clustering))
colnames(marker.cluster) <- c("markClust")
marker.cluster$Sample <- rownames(marker.cluster)
marker.cluster$markClust <- as.factor(marker.cluster$markClust)

marker_merge <- merge(panc.meta, marker_df, by='Sample')
marker_uber <- merge(marker_merge, marker.cluster, by='Sample')

# uncomment the following plots to see the overlay of each marker gene
# on the tSNE.  This illustrates the cell type label inference.
# # mediods 1 & 6 = alpha cells
# ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=GCG,
#                          shape=Kmediods)) +
#   geom_point(size=1) + theme_classic() +
#   scale_colour_continuous(low='white', high='red') +
#   scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))
# 
# # mediod 9 = beta cells, INS highest
# ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=INS,
#                          shape=Kmediods)) +
#   geom_point(size=1) + theme_classic() +
#   scale_colour_continuous(low='white', high='red') +
#   scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))
# 
# # mediod 3 = delta cells, SST highest
# ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=SST,
#                          shape=Kmediods)) +
#   geom_point(size=1) + theme_classic() +
#   scale_colour_continuous(low='white', high='red') +
#   scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))
# 
# # mediod 4 = PP cells, PPY highest
# ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=PPY,
#                          shape=Kmediods)) +
#   geom_point(size=1) + theme_classic() +
#   scale_colour_continuous(low='white', high='red') +
#   scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))
# 
# # mediod 5 = acinar cells, PRSS1 highest
# ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=PRSS1,
#                          shape=Kmediods)) +
#   geom_point(size=1) + theme_classic() +
#   scale_colour_continuous(low='white', high='red') +
#   scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))
# 
# # mediod 4 = ductal cells, KRT19 highest
# ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=KRT19,
#                          shape=Kmediods)) +
#   geom_point(size=1) + theme_classic() +
#   scale_colour_continuous(low='white', high='red') +
#   scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))
# 
# # small pop of mediod 9 = mesenchyme cells, COL1A1 highest, same as PP cells
# ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=COL1A1,
#                          shape=Kmediods)) +
#   geom_point(size=1) + theme_classic() +
#   scale_colour_continuous(low='white', high='red') +
#   scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))


table(marker_uber$Kmediods, marker_uber$markClust)
marker_uber$CellType <- ""
marker_uber$CellType[marker_uber$Kmediods %in% c(4)] <- "PP"
marker_uber$CellType[marker_uber$Kmediods %in%  c(9)] <- "Beta"
marker_uber$CellType[marker_uber$Kmediods %in% c(1, 7)] <- "Alpha"
marker_uber$CellType[marker_uber$Kmediods %in% c(3)] <- "Delta"
marker_uber$CellType[marker_uber$Kmediods %in% c(5, 8)] <- "Acinar"
marker_uber$CellType[marker_uber$Kmediods %in% c(6) ] <- "Ductal"
marker_uber$CellType[marker_uber$Kmediods %in% c(2)] <- "Mesenchyme"
table(marker_uber$CellType)

write.table(marker_uber,
            file="Pancreas/Data/GSE85241_marker_metadata.tsv",
            sep="\t", quote=F, row.names=F, col.names=T)

rm(list=ls())
gc()

