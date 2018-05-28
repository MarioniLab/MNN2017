# This script quantifies the performance of the different batch correction methods
# on simulations involving clusters in low-dimensional space.

# Loading in the functions.
source("functions.R")

########################################################################
# Simulating an easy scenario where the composition is the same across batches.
    
prop <- c(300, 500, 200)
for (it in seq_len(10)) {
    set.seed(it*10000) # Avoid issues with seeds being reset within a function.
    all.res <- generateSamples(ncells=cbind(prop, prop))
    cluster.id <- unlist(all.res$id)
    out <- runAllMethods(all.res$mat[[1]], all.res$mat[[2]])

    if (it==1L) {
        pdf("cluster_easy.pdf")
        for (type in names(out$mat)) {
            plotResults(out$mat[[type]], cluster.id, out$batch, main=type, pch.choices=c(16, 2))
        }
        dev.off()
    }

    unc.stat <- getVarExplained(out$mat$uncorrected, cluster.id, out$batch)
    mnn.stat <- getVarExplained(out$mat$MNN, cluster.id, out$batch)
    lm.stat <- getVarExplained(out$mat$limma, cluster.id, out$batch)
    com.stat <- getVarExplained(out$mat$ComBat, cluster.id, out$batch)
    cca.stat <- getVarExplained(out$mat$CCA, cluster.id, out$batch)
    write.table(file="cluster_easy.tsv", 
        data.frame(Uncorrected=unc.stat, MNN=mnn.stat, limma=lm.stat, ComBat=com.stat, CCA=cca.stat),
        sep="\t", quote=FALSE, row.names=FALSE, append=(it>1L), col.names=(it==1L))
}

########################################################################
# Simulating a harder scenario where the composition changes across batches.
    
prop1 <- c(300, 500, 200)
prop2 <- c(100, 200, 800)
for (it in seq_len(10)) { 
    set.seed(it*10000) # Avoid issues with seeds being reset within a function.
    all.res <- generateSamples(ncells=cbind(prop1, prop2))
    cluster.id <- unlist(all.res$id)
    out <- runAllMethods(all.res$mat[[1]], all.res$mat[[2]])

    if (it==1L) {
        pdf("cluster_hard.pdf")
        for (type in names(out$mat)) {
            plotResults(out$mat[[type]], cluster.id, out$batch, main=type, pch.choices=c(16, 2))
        }
        dev.off()
    }

    unc.stat <- getVarExplained(out$mat$uncorrected, cluster.id, out$batch)
    mnn.stat <- getVarExplained(out$mat$MNN, cluster.id, out$batch)
    lm.stat <- getVarExplained(out$mat$limma, cluster.id, out$batch)
    com.stat <- getVarExplained(out$mat$ComBat, cluster.id, out$batch)
    cca.stat <- getVarExplained(out$mat$CCA, cluster.id, out$batch)
    write.table(file="cluster_hard.tsv", 
        data.frame(Uncorrected=unc.stat, MNN=mnn.stat, limma=lm.stat, ComBat=com.stat, CCA=cca.stat),
        sep="\t", quote=FALSE, row.names=FALSE, append=(it>1L), col.names=(it==1L))
}

########################################################################
# Finally, simulating the hardest scenario where one population is just absent.
    
prop1 <- c(300, 500, 200)
prop2 <- c(800, 200, 0)
for (it in seq_len(10)) { 
    set.seed(it*10000) # Avoid issues with seeds being reset within a function.
    all.res <- generateSamples(ncells=cbind(prop1, prop2))
    cluster.id <- unlist(all.res$id)
    out <- runAllMethods(all.res$mat[[1]], all.res$mat[[2]])

    if (it==1L) {
        pdf("cluster_missing.pdf")
        for (type in names(out$mat)) {
            plotResults(out$mat[[type]], cluster.id, out$batch, main=type, pch.choices=c(16, 2))
        }
        dev.off()
    }

    unc.stat <- getVarExplained(out$mat$uncorrected, cluster.id, out$batch)
    mnn.stat <- getVarExplained(out$mat$MNN, cluster.id, out$batch)
    lm.stat <- getVarExplained(out$mat$limma, cluster.id, out$batch)
    com.stat <- getVarExplained(out$mat$ComBat, cluster.id, out$batch)
    cca.stat <- getVarExplained(out$mat$CCA, cluster.id, out$batch)
    write.table(file="cluster_missing.tsv", 
        data.frame(Uncorrected=unc.stat, MNN=mnn.stat, limma=lm.stat, ComBat=com.stat, CCA=cca.stat),
        sep="\t", quote=FALSE, row.names=FALSE, append=(it>1L), col.names=(it==1L))
}

########################################################################
# End.

