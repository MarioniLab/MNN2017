generateSamples <- function(means, SDs, ncells, ngenes=2000) 
# Generating simulated high-dimensional expression data containing
# cluster-based structure in a low-dimensional biological subspace.
{
    ndims <- ncol(means)
    nbatches <- ncol(ncells) 

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
    proj <- matrix(rnorm(ngenes*total), nrow=ngenes, ncol=ndims)
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
    return(list(mat=list(uncorrected=t(uncorrected), MNN=Xmnn$corrected, limma=t(Xlm), ComBat=t(Xcom)), batch=batch.id))
}

library(Rtsne)
plotResults <- function(mat, cluster.ids, batch.ids, main="", ...) 
# Creating t-SNE plots of the (corrected) expression matrices.
{
    tout <- Rtsne(mat, ...)
    plot(tout$Y[,1], tout$Y[,2], xlab="t-SNE 1", ylab="t-SNE 2", 
        col=rainbow(length(unique(cluster.ids)))[cluster.ids],
        pch=batch.ids, main=main)
}

getVarExplained <- function(mat, cluster.ids, batch.ids) 
# Using the variance explained as a metric - 
{
    v0 <- .GVE(mat, batch.ids)
    v1 <- .GVE(mat, cluster.ids)
    v2 <- .GVE(mat, paste0(cluster.ids, batch.ids))
    return(list(Cluster=v2/v0, Batch=v2/v1))
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
