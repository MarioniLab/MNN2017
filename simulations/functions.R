generateSamples <- function(ncells, ndims=2, ngenes=2000) 
# Generating simulated high-dimensional expression data containing
# cluster-based structure in a low-dimensional biological subspace.
{
    nbatches <- ncol(ncells) 
    nclusters <- nrow(ncells)

    # Cells are distributed according to a multivariate normal in a low-dimensional biological subspace. 
    # Each cell type has a different x/y center and a different SD.
    means <- matrix(rnorm(nclusters * ndims, sd=5), ncol=ndims)
    SDs <- matrix(rgamma(nclusters*ndims, 1, 1), ncol=ndims)

    all.dimensions <- all.clusters <- vector("list", nbatches)
    for (bdx in seq_len(nbatches)) {
        cur.ncells <- ncells[,bdx]
        cluster.id <- rep(seq_along(cur.ncells), cur.ncells)
        all.clusters[[bdx]] <- cluster.id

        total <- sum(cur.ncells)
        cur.dimensions <- matrix(0, total, ndims)
        for (idx in seq_len(ndims)) {
            cur.dimensions[,idx] <- rnorm(total, means[cluster.id,idx], sd=SDs[cluster.id,idx])
        }
        all.dimensions[[bdx]] <- cur.dimensions 
    }

    # Random projection to high-dimensional space.
    proj <- matrix(rnorm(ngenes*ndims), nrow=ngenes, ncol=ndims)
    for (bdx in seq_len(nbatches)) {
        mat <- tcrossprod(all.dimensions[[bdx]], proj)

        mat <- mat + rnorm(length(mat)) # Adding some random noise.
        mat <- t(mat)
        mat <- mat + rnorm(nrow(mat)) # Adding a random batch vector
        all.dimensions[[bdx]] <- mat 
    }

    # Add normally distributed noise.
    return(list(mat=all.dimensions, id=all.clusters))
}

library(scran)
library(sva)
library(limma)
library(Seurat)
runAllMethods <- function(...) 
# Running all batch correction methods.
{
	batches <- list(...)
    uncorrected <- do.call(cbind, batches)
	per.batch <- unlist(lapply(batches, FUN=ncol))
    batch.id <- rep(seq_along(batches), per.batch)

    Xmnn <- do.call(mnnCorrect2, c(batches, list(cos.norm.in=FALSE, approximate=TRUE)))
    Xlm <- removeBatchEffect(uncorrected, factor(batch.id))
    Xcom <- ComBat(uncorrected, factor(batch.id), mod=NULL, prior.plots = FALSE)

    mat <- list(uncorrected=t(uncorrected), MNN=Xmnn$corrected, limma=t(Xlm), ComBat=t(Xcom))

	if (length(batches)==2L) {
        colnames(batches[[1]]) <- paste0("Cell", seq_len(ncol(batches[[1]])), "-1")
        colnames(batches[[2]]) <- paste0("Cell", seq_len(ncol(batches[[2]])), "-2")
        rownames(batches[[1]]) <- rownames(batches[[2]]) <- paste0("Gene", seq_len(nrow(batches[[1]])))

        Se <- CreateSeuratObject(batches[[1]])
		Se2 <- CreateSeuratObject(batches[[2]])
	    Se@meta.data$group <- "group1"
	    Se2@meta.data$group <- "group2"

        Se <- ScaleData(Se)
        Se2 <- ScaleData(Se2)
        Y <- RunCCA(Se, Se2, genes.use=rownames(batches[[1]]), do.normalize=FALSE)
        suppressWarnings(Y <- AlignSubspace(Y, grouping.var="group", dims.align=1:20))
        mat$CCA <- Y@dr$cca.aligned@cell.embeddings
    } 

    return(list(mat=mat, batch=batch.id))
}

library(Rtsne)
plotResults <- function(mat, cluster.ids, batch.ids, 
    pch.choices=seq_len(max(batch.ids)),
    col.choices=rainbow(max(cluster.ids)),
    main="", ...) 
# Creating t-SNE plots of the (corrected) expression matrices.
{
    tout <- Rtsne(mat, ..., check_duplicates=FALSE)
    plot(tout$Y[,1], tout$Y[,2], xlab="t-SNE 1", ylab="t-SNE 2", 
        col=col.choices[cluster.ids],
        pch=pch.choices[batch.ids], main=main)
}

getVarExplained <- function(mat, cluster.ids, batch.ids) 
# Using the variance explained as a metric. Variance explained
# by cluster (while blocking on batch) should be high, variance 
# explained by batch (while blocking on cluster) should be low.
{
    v0 <- .GVE(mat, batch.ids)
    v1 <- .GVE(mat, cluster.ids)
    v2 <- .GVE(mat, paste0(cluster.ids, batch.ids))
    return(list(Cluster=1-v2/v0, Batch=1-v2/v1))
}

.GVE <- function(mat, block) {
    sum.var <- 0
    for (x in unique(block)) {
        chosen <- t(mat[x==block,])
        cur.var <- rowSums((chosen - rowMeans(chosen))^2)
        sum.var <- sum.var + cur.var 
    }
    return(sum(sum.var))
}
