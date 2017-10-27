# Here, we are trying to figure out where the remaining difference 
# between batches comes from, even after batch correction.

library(scran)
load("../logdataFandA_all.RData") 

# Keeping only MEPs for easier diagnostics. 
keepF <- colnames(logDataF3) %in% c("MEP")#, "GMP")
keepA <- colnames(logDataA3) %in% c("ERY")# "GMP")

logDataF3 <- logDataF3[,keepF]
logDataA3 <- logDataA3[,keepA]
first.batch <- rep(c(TRUE, FALSE), c(ncol(logDataF3), ncol(logDataA3)))

colF <- c(MEP="blue", CMP="orange", GMP="red")[colnames(logDataF3)]
colA <- c(ERY="lightblue", CMP="yellow", GMP="salmon")[colnames(logDataA3)]
colF[is.na(colF)] <- "grey50"
colA[is.na(colA)] <- "grey80"
allcolors <- c(colF, colA)

# Performing the correction and running PCA on the first 50 PCs.
mnn.out <- mnnCorrect(logDataF3, logDataA3, k=20, sigma=100, cos.norm=TRUE, svd.dim=NA)
X.mnn <- cbind(mnn.out$corrected[[1]], mnn.out$corrected[[2]])
t.mnn <- t(X.mnn)

pca.mnn <- prcomp(t.mnn, rank=50)
set.seed(100)
library(Rtsne)
tsne.mnn <- Rtsne(pca.mnn$x, perplexity=30, pca=FALSE)
Y <- tsne.mnn$Y
plot(Y[,1], Y[,2], cex=2, pch=21, col="black", bg=allcolors)

# Identifying divergent populations and testing for DE in the same batch.
out <- selectorPlot(Y[,1], Y[,2])

xout <- out
xout[[1]] <- xout[[1]] & first.batch # First cluster is the deviant MEP subpopulation in the first batch.
xout[[2]] <- !xout[[1]] & first.batch & allcolors=="blue" # Second cluster is the remaining MEPs in the first batch.
xout[[3]] <- !first.batch & allcolors=="lightblue" # Third cluster is the MEPs in the second batch.
to.test <- unlist(lapply(xout, which))

raw.all <- cbind(logDataF3, logDataA3)[,to.test]
cluster.id <- rep(seq_along(xout), unlist(lapply(xout, sum)))
res <- findMarkers(raw.all, cluster=cluster.id, pval.type="all")
head(res[[1]])

# Colouring by expression of certain genes.
library(RColorBrewer)
rbPal <- colorRampPalette(c('grey95', 'blue'))
curgene <- "ENSMUSG00000086503"
col <- rbPal(100)[cut(c(logDataF3[curgene,], logDataA3[curgene,]), breaks = 100)]
plot(Y[,1], Y[,2], cex=2, pch=21, col="black", bg=col)

# Checking the expression in the first batch.
tsne.F3 <- Rtsne(t(logDataF3), perplexity=30)
col <- rbPal(100)[cut(logDataF3[curgene,], breaks = 100)]
plot(tsne.F3$Y[,1], tsne.F3$Y[,2], pch=21, col="black", bg=col)

