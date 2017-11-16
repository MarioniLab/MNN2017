## community detection on 78,000 cells, this is slow so will need to be executed on a cluster with ~50GB memory available

library(scran)
library(igraph)

pbmc <- read.table("Droplet/Results/PBMC_68k_corrected.tsv", sep="\t",
                   h=TRUE, stringsAsFactors=FALSE)
rownames(pnmc) <- pbmc$gene_id

tcell.hvg <- read.table("Droplet/Data/Tcell/Tcell_hvg.tsv",
                        h=T, sep="\t", stringsAsFactors=FALSE)

pbmc.hvg <- read.table("Droplet/Data/PBMC/PBMC_hvg.tsv",
                       h=T, sep="\t", stringsAsFactors=FALSE)

select.hvg <- intersect(tcell.hvg$HVG, pbmc.hvg$HVG)

# build an SNN graph
pbmc.snn <- buildSNNGraph(pbmc[select.hvg, 2:dim(pbmc)[2]], pc.approx=TRUE, d=30, k=10)
pbmc.comm <- cluster_walktrap(pbmc.snn, steps=5)

pbmc.cluster <- do.call(cbind.data.frame, list("Sample"=colnames(pbmc)[2:dim(pbmc)[2]],
                                               "Community"=as.character(membership(pbmc.comm))))
write.table(pbmc.cluster,
            file="Droplet/Results/PBMC_Corrected_communities.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)