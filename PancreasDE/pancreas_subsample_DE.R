### subsample pancreas data after combination to demonstrate the extra DE genes are not just an
### artifact of the batch effect correction
library(ggplot2)
library(reshape2)
library(limma)
library(RColorBrewer)
library(goseq)
library(org.Hs.eg.db)
library(flashClust)
library(VennDiagram)
library(biomaRt)
library(stringr)
library(WGCNA)
library(ggrepel)
library(cowplot)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")
source("~/Dropbox/R_sessions/SingleCellFunctions/single_cell_functions.R")
source("~/Dropbox/R_sessions/Clustering/CosineKNN.R")
ensembl <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl', GRCh=37)

## this script handles the subsampling of cells down to the smallest batch after batch correction

gse81076.correct <- read.table("~/Dropbox/pancreas/CorrectedData/C_CELLseq_GSE81076.txt",
                               h=T, stringsAsFactors=F)

gse85241.correct <- read.table("~/Dropbox/pancreas/CorrectedData/C_CELLseq_SSE85241.txt",
                               h=T, stringsAsFactors=F)

emtab5061.correct <- read.table("~/Dropbox/pancreas/CorrectedData/C_Smartseq_EMATB5061.txt",
                                h=T, stringsAsFactors=F)

gse86473.correct <- read.table("~/Dropbox/pancreas/CorrectedData/C_Smartseq_GSE86473.txt",
                               h=T, stringsAsFactors=F)

gse81076.raw <- read.table("~/Dropbox/pancreas/GSE81076-norm.tsv.gz", sep="\t",
                           h=T, stringsAsFactors=F)

gse85241.raw <- read.table("~/Dropbox/pancreas/GSE85241-norm.tsv.gz", sep="\t",
                           h=T, stringsAsFactors=F)

emtab5061.raw <- read.table("~/Dropbox/pancreas/E-MTAB-5061-norm.tsv.gz", sep="\t",
                            h=T, stringsAsFactors=F)

gse86473.raw <- read.table("~/Dropbox/pancreas/GSE86473-norm.tsv.gz", sep="\t",
                           h=T, stringsAsFactors=F)

all.genes <- unique(c(rownames(gse81076.correct), rownames(gse85241.correct), rownames(emtab5061.correct), rownames(gse86473.correct)))

gene_symbol <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                     filters='external_gene_name', mart=ensembl,
                     values=all.genes)

# meta data
gse81076.meta <- read.table("~/Dropbox/pancreas/GSE81076_marker_metadata.tsv", sep="\t",
                            stringsAsFactors=F, h=T)

gse85241.meta <- read.table("~/Dropbox/pancreas/GSE85241_marker_metadata.tsv", sep="\t",
                            stringsAsFactors=F, h=T)

emtab5061.meta <- read.table("~/Dropbox/pancreas/E-MTAB-5061_marker_metadata.tsv", sep="\t",
                             stringsAsFactors=F, h=T)

gse86473.meta <- read.table("~/Dropbox/pancreas/GSE86473_marker_metadata.tsv", sep="\t",
                            stringsAsFactors=F, h=T)
gse86473.meta$Study <- "GSE86473"

gse81076.correct$gene_id <- rownames(gse81076.correct)

gse85241.correct$gene_id <- rownames(gse85241.correct)

gse86473.correct$gene_id <- rownames(gse86473.correct)

emtab5061.correct$gene_id <- rownames(emtab5061.correct)

set.seed(42)
merge.list <- list("GSE81076"=gse81076.correct,
                   "GSE85241"=gse85241.correct,
                   "EMTAB5061"=emtab5061.correct,
                   "GSE86473"=gse86473.correct)

merged.correct <- Reduce(f=function(x, y) merge(x, y, by='gene_id'),
                         x=merge.list)
rownames(merged.correct) <- merged.correct$gene_id
merged.correct <- merged.correct[, 2:dim(merged.correct)[2]]

merged.meta <- data.frame(do.call(rbind, list("GSE81076"=gse81076.meta[, c("Sample", "Study", "CellType")],
                                              "GSE85241"=gse85241.meta[, c("Sample", "Study", "CellType")],
                                              "EMTAB5061"=emtab5061.meta[, c("Sample", "Study", "CellType")],
                                              "GSE86473"=gse86473.meta[, c("Sample", "Study", "CellType")])))
merged.meta$CellType <- gsub(tolower(merged.meta$CellType), pattern=" cell", replacement="")

table(merged.meta$CellType)

# use the original set of HVG
gse86473.hvg <- read.table("~/Dropbox/pancreas/GSE86473-HVG.tsv", h=F, sep="\t",
                           stringsAsFactors=F)

gse81076.hvg <- read.table("~/Dropbox/pancreas/GSE81076-HVG.tsv", h=F, sep="\t",
                           stringsAsFactors=F)

emtab5061.hvg <- read.table("~/Dropbox/pancreas/E-MTAB-5061-HVG.tsv", h=F, sep="\t",
                            stringsAsFactors=F)

gse85241.hvg <- read.table("~/Dropbox/pancreas/GSE85241-HVG.tsv", h=F, sep="\t",
                           stringsAsFactors=F)

union.hvg <- unique(c(unlist(gse86473.hvg), unlist(gse81076.hvg), unlist(emtab5061.hvg), unlist(gse85241.hvg)))
all.hvg <- rownames(merged.correct) %in% union.hvg
names(all.hvg) <- rownames(merged.correct)

table(all.hvg)

# remove the small scrappy cell types that are poorly defined
inc.types <- c("pp", "alpha", "beta", "acinar", "delta", "ductal", "mesenchyme", "psc", "gamma")
inc.cells <- merged.meta$Sample[merged.meta$CellType %in% inc.types]

# cluster cells on the hvg
all.cells <- merged.correct[all.hvg, intersect(inc.cells, colnames(merged.correct))]
all.cells[is.na(all.cells)] <- 0

all.cells <- all.cells[, intersect(colnames(all.cells), merged.meta$Sample)]
merged.meta <- merged.meta[merged.meta$Sample %in% intersect(colnames(all.cells), merged.meta$Sample), ]
rownames(merged.meta) <- merged.meta$Sample

# get the cell IDs for each study
gse81076.cells <- merged.meta$Sample[merged.meta$Study == "GSE81076"]
gse86473.cells <- merged.meta$Sample[merged.meta$Study == "GSE86473"]
emtab5061.cells <- merged.meta$Sample[merged.meta$Study == "EMTAB5061"]
gse85241.cells <- merged.meta$Sample[merged.meta$Study == "GSE85241"]

# do a clustering on Spearman correlation distance
all.cor <- cor(all.cells, method='spearman')
all.cor[is.na(all.cor)] <- 0
all.spearman.dist <- as.dist(1 - all.cor)
all.spearman.tree <- flashClust(all.spearman.dist, method='average')
all.spearman.cut <- cutreeDynamicTree(dendro=all.spearman.tree, minModuleSize = 50, 
                                      deepSplit=F)
all.spearman.cols <- labels2colors(all.spearman.cut)
all.cluster <- data.frame(cbind(colnames(all.cells), all.spearman.cols))
colnames(all.cluster) <- c("Sample", "Cluster")

set.seed(42)
all.dist <- dist(t(all.cells), method='euclidean')
all.tsne <- tsne_wrapper(all.dist, perplexity=30, is.dist=TRUE)

all.map <- merge(all.tsne, all.cluster, by='Sample')
all.merge <- merge(all.map, merged.meta, by='Sample')

# try a k-means with k=8
# try with PCs first, then tSNE
set.seed(42)
#all.kmeans <- kmeans(all.dist, centers=12)
all.kmeans <- kmeans(all.tsne[, 1:2], centers=8)
kmean.df <- data.frame(cbind(colnames(all.cells), all.kmeans$cluster))
colnames(kmean.df) <- c("Sample", "Kmeans")

kmean.merge <- merge(all.merge, kmean.df, by='Sample')

marker.genes <- c("SST", "PPY", "PRSS1", "KRT19", "INS", "GHRL", "GCG", "ESAM", "COL1A1")
marker.exprs <- data.frame(t(all.cells[marker.genes, ]))
marker.exprs$Sample <- rownames(marker.exprs)

meta.marker <- merge(kmean.merge, marker.exprs, by='Sample')

meta.marker$new.labels <- "Unk"
meta.marker$new.labels[meta.marker$Kmeans == 1 & meta.marker$INS >= 4] <- "Beta"
meta.marker$new.labels[meta.marker$Kmeans == 2 & meta.marker$INS >= 4] <- "Beta"
meta.marker$new.labels[meta.marker$Kmeans == 3 & meta.marker$INS >= 4] <- "Beta"
meta.marker$new.labels[meta.marker$Kmeans == 4 & meta.marker$GCG >= 4] <- "Alpha"
meta.marker$new.labels[meta.marker$Kmeans == 5 & meta.marker$PRSS1 >= 4] <- "Acinar"
meta.marker$new.labels[meta.marker$Kmeans == 6 & meta.marker$PPY > 5] <- "Gamma"
meta.marker$new.labels[meta.marker$Kmeans == 7 & meta.marker$KRT19 >= 4] <- "Ductal"
meta.marker$new.labels[meta.marker$Kmeans == 7 & meta.marker$COL1A1 >= 4] <- "Mesenchyme"
meta.marker$new.labels[meta.marker$Kmeans == 8 & meta.marker$SST > 5] <- "Delta"

table(meta.marker$CellType, meta.marker$new.labels)
print(randIndex(table(meta.marker$CellType, meta.marker$new.labels), adjust=TRUE))


table(meta.marker$CellType, meta.marker$Study)

# subsample pp to 15 cells and delta to 51 cells after batch correction
# demonstrate good mixing from each study


# select the same raw genes that are in the corrected data
# need to use the raw values, and block on batch, something very wrong with the corrected
# values on the log space
correct.genes <- rownames(merged.correct)
gse86473.raw$gene_id <- gse86473.raw$hgnc_symbol
emtab5061.raw$gene_id <- emtab5061.raw$hgnc_symbol

# explicitly select the same cells that are subsampled from the corrected data
rownames(gse86473.raw) <- gse86473.raw$gene_id
gse86473.tomerge <- gse86473.raw[, colnames(gse86473.raw) %in% meta.marker$Sample]
gse86473.tomerge$gene_id <- rownames(gse86473.tomerge)
gse86473.tomerge[is.na(gse86473.tomerge)] <- 0

rownames(gse85241.raw) <- gse85241.raw$gene_id
gse85241.tomerge <- gse85241.raw[, colnames(gse85241.raw) %in% meta.marker$Sample]
gse85241.tomerge$gene_id <- rownames(gse85241.tomerge)
gse85241.tomerge[is.na(gse85241.tomerge)] <- 0

rownames(emtab5061.raw) <- emtab5061.raw$gene_id
emtab5061.tomerge <- emtab5061.raw[, colnames(emtab5061.raw) %in% meta.marker$Sample]
emtab5061.tomerge$gene_id <- rownames(emtab5061.tomerge)
emtab5061.tomerge[is.na(emtab5061.tomerge)] <- 0

rownames(gse81076.raw) <- gse81076.raw$gene_id
gse81076.tomerge <- gse81076.raw[, colnames(gse81076.raw) %in% meta.marker$Sample]
gse81076.tomerge$gene_id <- rownames(gse81076.tomerge)
gse81076.tomerge[is.na(gse81076.tomerge)] <- 0

raw.list <- list("GSE86473"=gse86473.tomerge,
                 "GSE85241"=gse85241.tomerge,
                 "EMTAB5061"=emtab5061.tomerge,
                 "GSE81076"=gse81076.tomerge)

merged.raw <- Reduce(x=raw.list,
                     f=function(x, y) merge(x, y, by='gene_id'))

rownames(merged.raw) <- merged.raw$gene_id
correct.exprs <- merged.raw
correct.exprs <- correct.exprs[, 2:(dim(correct.exprs)[2])]

# try with original cell labels and newly assigned cell labels
# on the down sampled cells
n.gamma <- 15
n.delta <- 51

# randomly sample 15 PP/Gamma cells and 51 Delta cells
all.gamma <- meta.marker$Sample[meta.marker$new.labels == "Gamma"]
all.delta <- meta.marker$Sample[meta.marker$new.labels == "Delta"]

# perform 100 iterations of this sampling
# report the number and identity of DE genes
# check the intersection with the original combined DE dataset
correct.de <- read.table("~/Dropbox/MNNfigures/Pancreas_Corrected_DE_genes.tsv",
                         h=T, sep="\t", stringsAsFactors=F)
de.symbols <- unique(gene_symbol$external_gene_name[gene_symbol$ensembl_gene_id %in% correct.de$x])

n.iters <- 101
de.list <- list()
i <- 1
while(i < n.iters){
  tryCatch({
    sample.gamma <- sample(all.gamma, n.gamma)
    sample.delta <- sample(all.delta, n.delta)
    sample.exprs <- correct.exprs[, colnames(correct.exprs) %in% c(sample.gamma, sample.delta)]
    sample.meta <- meta.marker[meta.marker$Sample %in% c(sample.gamma, sample.delta), ]
    
    table(sample.meta$new.labels, sample.meta$Study)
    
    sample.meta$interact.group <- as.factor(apply(sample.meta,
                                                  1, FUN=function(X) paste(X["Study"], X["new.labels"], sep=".")))
    rownames(sample.meta) <- sample.meta$Sample
    
    sample.design <- model.matrix(~ 0 + interact.group, data=sample.meta)
    colnames(sample.design) <- levels(sample.meta$interact.group)
    
    # fit a linear model
    sample.fit <- limma::lmFit(sample.exprs[, rownames(sample.meta)],
                               sample.design)
  
    contrast.mat <- limma::makeContrasts(GSE81076.Delta-GSE81076.Gamma, GSE85241.Delta-GSE85241.Gamma,
                                         EMTAB5061.Delta-EMTAB5061.Gamma, GSE86473.Delta-GSE86473.Gamma,
                                         levels=sample.design)
    # get the average log fold changes across batch by cell type interactions
    sample.fit2 <- limma::contrasts.fit(sample.fit, contrasts=rowMeans(contrast.mat))
    sample.fit2 <- limma::treat(sample.fit2, lfc=1)
    sample.sum.res <- summary(limma::decideTests(sample.fit2))
    
    sample.sum.res
    
    sample.de.res <- limma::topTable(sample.fit2, n=Inf, sort="p", p=1.0)
    
    # set -inf to -308
    sample.de.res$adj.P.Val[sample.de.res$adj.P.Val == 0] <- .Machine$double.xmin
    
    sample.de.res$Sig <- 0
    sample.de.res$Sig[sample.de.res$adj.P.Val <= 0.05] <- 1
    sample.de.res$Sig <- as.factor(sample.de.res$Sig)
    
    sample.de.res$Diff <- 0
    sample.de.res$Diff[sample.de.res$logFC < 0 & sample.de.res$Sig == 1] <- -1
    sample.de.res$Diff[sample.de.res$logFC > 0 & sample.de.res$Sig == 1] <- 1
    
    sample.de.res$gene_id <- rownames(sample.de.res)
    
    de.list[[i]] <- list("GeneUp"=sample.de.res$gene_id[sample.de.res$Diff == 1],
                         "GeneDown"=sample.de.res$gene_id[sample.de.res$Diff == -1],
                         "TotalDE"=length(sample.de.res$gene_id[sample.de.res$Sig == 1]),
                         "Intersection"=length(intersect(sample.de.res$gene_id[sample.de.res$Sig == 1], de.symbols)))
    i <- i + 1
  }, warning = function(war){
    message(paste("WARNING: ", war))
  },
  error = function(err){
    message(paste("ERROR ", err))
    message("Trying a new sample")
  }, finally = {
  })
}

sample.de.res$gene_id <- rownames(sample.de.res)
sample.de.merge <- merge(sample.de.res,  gene_symbol,
                          by.x='gene_id', by.y='external_gene_name', all.x=TRUE)
sample.gene.top <- sample.de.merge$gene_id
sample.gene.top[sample.de.merge$logFC > -3 & sample.de.merge$logFC < 3] <- ""
sample.gene.top[is.na(sample.gene.top)] <- ""
sample.gene.top[sample.de.merge$Sig == 0] <- ""

ma.sample <- ggplot(sample.de.merge, aes(x=AveExpr, y=logFC, colour=as.factor(Sig))) +
  geom_point(size=1, alpha=0.8) + theme_mike() +
  scale_colour_Publication() +
  labs(x="Mean Cosine Normalized Expression",
       y=expression(paste("log"[2], " Fold Change")),
       title="MNN sub-sampled PP vs. Delta cells") +
  guides('colour'=guide_legend(title="Statistically significant")) +
  geom_text_repel(aes(label=sample.gene.top),
                  colour='black')

ma.sample

# what is the distribution of the number of DE genes and the intersection with the corrected DE genes?
# how many DEs
# check between each pair of iterations whether DE genes are in the same direction
concord.list <- list()
inter.vec <- c()
n.de <- c()
up.de <- c()
down.de <- c()

b <- 1
for(x in seq(1, n.iters-1)){
  comp.i <- de.list[[x]]
  comp.i$GeneUp
  inter.vec <- c(inter.vec, comp.i$Intersection)
  n.de <- c(n.de, comp.i$TotalDE)
  up.de <- c(up.de, length(comp.i$GeneUp))
  down.de <- c(down.de, length(comp.i$GeneDown))
  for(k in seq(1, n.iters-1)){
    comp.j <- de.list[[k]]
    if(x != k){
      concord.up <- length(intersect(comp.i$GeneUp, comp.j$GeneUp))/length(union(comp.i$GeneUp, comp.j$GeneUp))
      concord.down <- length(intersect(comp.i$GeneDown, comp.j$GeneDown))/length(union(comp.i$GeneDown, comp.j$GeneDown))
      discord.up <- length(intersect(comp.i$GeneUp, comp.j$GeneDown))/length(union(comp.i$GeneUp, comp.j$GeneDown))
      discord.down <- length(intersect(comp.i$GeneDown, comp.j$GeneUp))/length(union(comp.i$GeneDown, comp.j$GeneUp))
      concord.list[[b]] <- list("Up"=concord.up,
                               "Down"=concord.down,
                               "Discord.Up"=discord.up,
                               "Discord.Down"=discord.down)
      b <- b + 1
    }
    }
}

dir.cols <- c("#900967", "#f29e31")
names(dir.cols) <- c("Up", "Down")

de.df <- cbind.data.frame(inter.vec, n.de)
colnames(de.df) <- c("Intersection", "N.DE")

dir.df <- cbind.data.frame(up.de, down.de)
colnames(dir.df) <- c("N.Down", "N.Up")
dir.melt <- melt(dir.df)
dir.melt$variable <- gsub(dir.melt$variable, pattern="N.", replacement="")
colnames(dir.melt) <- c("Direction", "N")
dir.melt$Direction <- factor(dir.melt$Direction,
                             labels=c("Up", "Down"),
                             levels=c("Up", "Down"))

raw.v.correct <- ggplot(de.df, aes(x=Intersection, y=N.DE)) +
  geom_point() + theme_mike() +
  labs(x=expression("Corrected DE genes "*intersect("raw DE genes")),
       y=expression(paste("N raw DE genes"))) +
  theme(axis.title=element_text(size=12, face='bold'))

N.de.box <- ggplot(dir.melt, aes(x=Direction, y=N, fill=Direction)) +
  geom_boxplot() + theme_mike() +
  scale_fill_manual(values=dir.cols) +
  labs(x=expression("Differential Direction"), y=expression("N DE genes")) +
  guides(fill=FALSE)+
  theme(axis.title=element_text(face='bold', size=12))

concord.df <- do.call(rbind.data.frame, concord.list)
concord.melt <- melt(concord.df[, c("Up", "Down")])
colnames(concord.melt) <- c("Direction", "Concordant")

discord.melt <- melt(concord.df[, c("Discord.Up", "Discord.Down")])
colnames(discord.melt) <- c("Direction.D", "Discordant")
discord.melt$Direction.D <- gsub(discord.melt$Direction.D, pattern="Discord.", replacement="")

conc.df <- cbind.data.frame(concord.melt, discord.melt)
conc.df$Dir <- "Discordant"
conc.df$Dir[conc.df$Direction == "Up" & conc.df$Direction.D == "Up"] <- "Up"
conc.df$Dir[conc.df$Direction == "Down" & conc.df$Direction.D == "Down"] <- "Down"
conc.df$Dir <- factor(conc.df$Dir,
                      levels=c("Up", "Down"),
                      labels=c("Up", "Down"))

# there are no discordant direction DE genes
cond.box <- ggplot(conc.df, aes(x=Direction, y=Concordant, fill=Direction.D)) +
  geom_boxplot() + theme_mike() +
  scale_fill_manual(values=dir.cols) +
  labs(x=expression("Differential Direction"), y=expression("% overlap genes")) +
  guides(fill=FALSE) +
  theme(axis.title=element_text(face='bold', size=12))

grid.de <- plot_grid(raw.v.correct, N.de.box, cond.box, ncol=3,
          rel_widths=c(0.4, 0.3, 0.3))

ggsave(grid.de,
       filename="~/Dropbox/MNNfigures/DE_subsample.png",
       height=2.75, width=8.75, dpi=300)

