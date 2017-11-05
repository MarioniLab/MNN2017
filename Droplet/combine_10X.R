# combine human preimplantation embryo data sets
# might be tricky as one data set has very few HVGs

library(ggplot2)
library(scran)
library(gplots)
library(RColorBrewer)
library(Rtsne)
library(stringr)
library(WGCNA)
library(biomaRt)
library(igraph)
source("SomeFuncs/convenience_functions.R")

tcell <- read.table("Droplet/10X_data/Tcell/Tcell_norm.tsv",
                    h=T, sep="\t", stringsAsFactors=FALSE)

rownames(tcell) <- tcell$gene_id

tcell.hvg <- read.table("Droplet/10X_data/Tcell/Tcell_hvg.tsv",
                        h=T, sep="\t", stringsAsFactors=FALSE)

tcell.meta <- read.table("Droplet/10X_data/Tcell/Tcell_meta.tsv",
                        h=T, sep="\t", stringsAsFactors=FALSE)
tcell.meta$Sample <- paste("tcell", tcell.meta$Sample, sep=".")
tcell.meta$Sample <- gsub(tcell.meta$Sample, pattern="-", replacement=".")

pbmc <- read.table("Droplet/10X_data/PBMC/PBMC_norm.tsv",
                   h=T, sep="\t", stringsAsFactors=FALSE)
rownames(pbmc) <- pbmc$gene_id

pbmc.hvg <- read.table("Droplet/10X_data/PBMC/PBMC_hvg.tsv",
                       h=T, sep="\t", stringsAsFactors=FALSE)

pbmc.meta <- read.table("Droplet/10X_data/PBMC/PBMC_meta.tsv",
                        h=T, sep="\t", stringsAsFactors=FALSE)
pbmc.meta$Sample <- paste("pbmc", pbmc.meta$Sample, sep=".")
pbmc.meta$Sample <- gsub(pbmc.meta$Sample, pattern="-", replacement=".")

# select the total number of overlapping genes
# and the numer of intersecting HVGs
select.genes <- intersect(rownames(tcell), rownames(pbmc))
select.hvg <- intersect(tcell.hvg$HVG, pbmc.hvg$HVG)

# there are only ~300 HVG shared between both data sets
tcell.batch <- tcell[tcell$gene_id %in% select.genes, 1:(dim(tcell)[2]-1)]
colnames(tcell.batch) <- paste("tcell", colnames(tcell.batch), sep=".")

pbmc.batch <- pbmc[pbmc$gene_id %in% select.genes, 1:(dim(pbmc)[2]-1)]
colnames(pbmc.batch) <- paste("pbmc", colnames(pbmc.batch), sep=".")

# check they have the same number of rows and ordering
dim(tcell.batch)[1] == dim(pbmc.batch)[1]

all(rownames(tcell.batch) == rownames(pbmc.batch))

# save these for later
tcell.cells <- colnames(tcell.batch)
pbmc.cells <- colnames(pbmc.batch)

batch.meta <- do.call(rbind.data.frame, list("Tcell"=tcell.meta,
                                             "PBMC"=pbmc.meta))

# tSNE on the raw uncorrected data first
tcell.batch$gene_id <- rownames(tcell.batch)
pbmc.batch$gene_id <- rownames(pbmc.batch)

merge.batch <- merge(tcell.batch, pbmc.batch, by='gene_id')
uncorrect.tsne <- tsne_wrapper(merge.batch[, 2:dim(merge.batch)[2]])
uncorrect.tsne$Study <- ""
uncorrect.tsne$Study[uncorrect.tsne$Sample %in% tcell.cells] <- "Tcells"
uncorrect.tsne$Study[uncorrect.tsne$Sample %in% pbmc.cells] <- "PBMC"

uncorrect.tsne.merge <- merge(uncorrect.tsne, batch.meta, by='Sample')

label_cols <- c("#386cb0", "#fdb462",
                "#7fc97f","#ef3b2c","#662506",
                "#a6cee3","#fb9a99","#984ea3",
                "#ffff33", "#b2df8a")
names(label_cols) <- unique(uncorrect.tsne.merge$CellType)

batch.cols <- c("#33a02c", "#cab2d6")
names(batch.cols) <- unique(uncorrect.tsne.merge$Study)

unc.by_study <- ggplot(uncorrect.tsne.merge,
       aes(x=Dim1, y=Dim2, fill=Study)) +
  geom_point(size=3, shape=21) + theme_classic() +
  labs(x="tSNE Dimension 1", y="tSNE Dimension 2") +
  guides(colour=guide_legend(title='Batch', nrow=3)) +
  scale_fill_manual(values=batch.cols) +
  theme(axis.title=element_text(size=20, colour='black'),
        axis.text=element_text(size=16, colour='black'))

ggsave(unc.by_study,
       filename="Droplet/Uncorrected_by_study.png",
       height=4.75, width=6.75, dpi=300)

unc.by_cell <- ggplot(uncorrect.tsne.merge,
                       aes(x=Dim1, y=Dim2, fill=CellType)) +
  geom_point(size=3, shape=21) + theme_classic() +
  labs(x="tSNE Dimension 1", y="tSNE Dimension 2") +
  guides(colour=guide_legend(title='Batch', nrow=3)) +
  scale_fill_manual(values=label_cols) +
  theme(axis.title=element_text(size=20, colour='black'),
        axis.text=element_text(size=16, colour='black'))

ggsave(unc.by_cell,
       filename="Droplet/Uncorrected_by_celltype.png",
       height=4.75, width=6.75, dpi=300)

# remove the gene_id column before correction, etc
tcell.batch <- tcell.batch[, 1:(dim(tcell.batch)[2]-1)]
pbmc.batch <- pbmc.batch[, 1:(dim(pbmc.batch)[2]-1)]

correct.data <- mnnCorrect(as.matrix(pbmc.batch),
                           as.matrix(tcell.batch),
                           k=20, subset.row=select.hvg)

pbmc.correct <- as.data.frame(correct.data$corrected[[1]])
colnames(pbmc.correct) <- colnames(pbmc.batch)
pbmc.correct$gene_id <- rownames(pbmc.correct)

tcell.correct <- as.data.frame(correct.data$corrected[[2]])
colnames(tcell.correct) <- colnames(tcell.batch)
tcell.correct$gene_id <- rownames(tcell.correct)

correct.df <- merge(pbmc.correct, tcell.correct, by='gene_id')
rownames(correct.df) <- correct.df$gene_id

# assess correction using tSNE
set.seed(42)
correct.tsne <- tsne_wrapper(correct.df[, 2:dim(correct.df)[2]])

# get study info and other meta data
correct.tsne$Study <- ""
correct.tsne$Study[correct.tsne$Sample %in% tcell.cells] <- "Tcells"
correct.tsne$Study[correct.tsne$Sample %in% pbmc.cells] <- "PBMC"

tsne.merge <- merge(correct.tsne, batch.meta, by='Sample')

# looks pretty good!  Not bad for only 306 genes!!!
by.label <- ggplot(tsne.merge, aes(x=Dim1, y=Dim2, fill=CellType)) +
  geom_point(size=3, shape=21) + theme_classic() +
  labs(x="tSNE Dimension 1", y="tSNE Dimension 2") +
  guides(colour=guide_legend(title='Cell Label', nrow=3)) +
  scale_fill_manual(values=label_cols) +
  theme(axis.title=element_text(size=20, colour='black'),
        axis.text=element_text(size=16, colour='black'))

ggsave(by.label, 
       filename="Droplet/10X_data/PBMC_Tcell-4K_bylabel-tSNE.png",
       height=5.5, width=7.5, dpi=300)

# let's construct an SNN graph out of this and see how that looks
correct.snn <- buildSNNGraph(as.matrix(correct.df[, 2:dim(correct.df)[2]]),
                             k=5, d=20)

correct.community <- cluster_walktrap(correct.snn, steps=4)
n.comms <- length(correct.community)
print(n.comms)

# vertex community membership
correct.members <- correct.community$membership
correct.clusters <- cbind.data.frame(colnames(correct.df[, 2:dim(correct.df)[2]]), correct.members)
colnames(correct.clusters) <- c("Sample", "Community")
correct.clusters$Community <- as.factor(correct.clusters$Community)
correct.max <- merge(tsne.merge, correct.clusters, by='Sample')


comm.tsne <- ggplot(correct.max,
                    aes(x=Dim1, y=Dim2, fill=Study)) + 
  geom_point(size=3, shape=21) + theme_classic() +
  guides(colour=guide_legend(title='Data Set')) +
  labs(x="tSNE Dimension 1", y="tSNE Dimension 2") +
  scale_fill_manual(values=batch.cols) +
  theme(axis.title=element_text(size=20, colour='black'),
        axis.text=element_text(size=16, colour='black'))

ggsave(comm.tsne,
       filename="Droplet/10X_data/PBMC_Tcell-4K-tSNE.png",
       height=5.5, width=7.5, dpi=300)

# now it would be a good idea to check the corresponding cell labels
table(correct.max$Community.y, correct.max$Study)

# how many T cells from the 4K fall into the correct cluster?
table(correct.max$CellType, correct.max$Study)

# how many T cells from the 4K fall into the correct cluster?
# map the community clusters to actual cell types,
# e.g. B_cell = 5
# simplify all T cell subsets in to 'T cells'
table(correct.max$CellType, correct.max$Community.y)

correct.max$Community.CellType <- "Unassigned"
correct.max$Community.CellType[correct.max$Community.y %in% c(5)] <- "B_cell"
correct.max$Community.CellType[correct.max$Community.y %in% c(13)] <- "macrophage"
correct.max$Community.CellType[correct.max$Community.y %in% c(9)] <- "neutrophil"
correct.max$Community.CellType[correct.max$Community.y %in% c(1, 2)] <- "monocyte"
correct.max$Community.CellType[correct.max$Community.y %in% c(3, 4, 6, 7, 8, 11, 12, 14:16)] <- "T_cell"
table(correct.max$Study, correct.max$Community.CellType)
table(correct.max$CellType, correct.max$Community.CellType)

new.cols <- c("#386cb0", "#518A87",
              "#7fc97f","#ef3b2c","#662506",
              "#a6cee3")
names(new.cols) <- unique(correct.max$Community.CellType)

by.newlabel <- ggplot(correct.max,
                       aes(x=Dim1, y=Dim2, fill=Community.CellType)) + 
  geom_point(size=3, shape=21) + theme_classic() +
  guides(colour=guide_legend(title='Data Set')) +
  labs(x="tSNE Dimension 1", y="tSNE Dimension 2") +
  scale_fill_manual(values=new.cols) +
  theme(axis.title=element_text(size=20, colour='black'),
        axis.text=element_text(size=16, colour='black'))

