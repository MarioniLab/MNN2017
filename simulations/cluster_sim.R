# Loading in the functions.
source("functions.R")

########################################################################
# Our simulation involves three cell types/components.
# Cells are distributed according to a bivariate normal in a 2-D biological subspace. 
# Each cell type has a different x/y center and a different SD.

mus <- cbind(c(0,5,5),
             c(5,5,0))
sds <- cbind(c(1,0.1,1),
             c(1,0.1,0.5))

# Note that the different centers should not lie on the same `y=mx` line; this represents populations that differ only in library size. 
# Such differences should not be present in normalized data, and will be eliminated by the cosine normalization step.
# The centers above are chosen so as to guarantee good separation between the different components.

########################################################################
# Simulating an easy scenario where the composition is the same across batches.
    
prop <- c(300, 500, 200)
for (it in seq_len(10)) { 
    all.res <- generateSamples(mus, sds, cbind(prop, prop))
    cluster.id <- unlist(all.res$id)
    out <- runAllMethods(all.res$mat[[1]], all.res$mat[[2]])

    if (it==1L) {
        pdf("cluster_easy.pdf")
        for (type in names(out$mat)) {
            plotResults(out$mat[[type]], cluster.id, out$batch, main=type)
        }
        dev.off()
    }

    unc.stat <- getVarExplained(out$mat$uncorrected, cluster.id, out$batch)
    mnn.stat <- getVarExplained(out$mat$MNN, cluster.id, out$batch)
    lm.stat <- getVarExplained(out$mat$limma, cluster.id, out$batch)
    com.stat <- getVarExplained(out$mat$ComBat, cluster.id, out$batch)
    write.table(file="cluster_easy.tsv", 
        data.frame(Uncorrected=unc.stat, MNN=mnn.stat, limma=lm.stat, ComBat=com.stat),
        sep="\t", quote=FALSE, row.names=FALSE, append=(it>1L), col.names=(it==1L))
}

########################################################################
# Simulating a harder scenario where the composition changes across batches.
    
prop1 <- c(300, 500, 200)
prop2 <- c(100, 200, 800)
for (it in seq_len(10)) { 
    all.res <- generateSamples(mus, sds, cbind(prop1, prop2))
    cluster.id <- unlist(all.res$id)
    out <- runAllMethods(all.res$mat[[1]], all.res$mat[[2]])

    if (it==1L) {
        pdf("cluster_hard.pdf")
        for (type in names(out$mat)) {
            plotResults(out$mat[[type]], cluster.id, out$batch, main=type)
        }
        dev.off()
    }

    unc.stat <- getVarExplained(out$mat$uncorrected, cluster.id, out$batch)
    mnn.stat <- getVarExplained(out$mat$MNN, cluster.id, out$batch)
    lm.stat <- getVarExplained(out$mat$limma, cluster.id, out$batch)
    com.stat <- getVarExplained(out$mat$ComBat, cluster.id, out$batch)
    write.table(file="cluster_hard.tsv", 
        data.frame(Uncorrected=unc.stat, MNN=mnn.stat, limma=lm.stat, ComBat=com.stat),
        sep="\t", quote=FALSE, row.names=FALSE, append=(it>1L), col.names=(it==1L))
}

########################################################################
# Finally, simulating the hardest scenario where one population is just absent.
    
prop1 <- c(300, 500, 200)
prop2 <- c(800, 200, 0)
for (it in seq_len(10)) { 
    all.res <- generateSamples(mus, sds, cbind(prop1, prop2))
    cluster.id <- unlist(all.res$id)
    out <- runAllMethods(all.res$mat[[1]], all.res$mat[[2]])

    if (it==1L) {
        pdf("cluster_missing.pdf")
        for (type in names(out$mat)) {
            plotResults(out$mat[[type]], cluster.id, out$batch, main=type)
        }
        dev.off()
    }

    unc.stat <- getVarExplained(out$mat$uncorrected, cluster.id, out$batch)
    mnn.stat <- getVarExplained(out$mat$MNN, cluster.id, out$batch)
    lm.stat <- getVarExplained(out$mat$limma, cluster.id, out$batch)
    com.stat <- getVarExplained(out$mat$ComBat, cluster.id, out$batch)
    write.table(file="cluster_missing.tsv", 
        data.frame(Uncorrected=unc.stat, MNN=mnn.stat, limma=lm.stat, ComBat=com.stat),
        sep="\t", quote=FALSE, row.names=FALSE, append=(it>1L), col.names=(it==1L))
}

########################################################################
# End.

