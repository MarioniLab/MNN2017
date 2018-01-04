# This checks the use of a classifier on the simulated data,
# to distinguish between the two batches.

library(scran)
library(limma)
library(sva)
library(e1071)

fitSVM <- function(input, batch) {
    input <- prcomp(t(input), rank.=10)
    batch <- factor(batch)
    tuned <- tune(svm, train.x=input$x, train.y=batch, scale=FALSE,
        ranges = list(epsilon = seq(0,1,0.05), cost = 2^(-2:4)))
    return(tuned)
}

load("Sim.RData")
for (easy in c(FALSE, TRUE)) {
    if (easy) {
        B2 <- B2ii
        clust2 <- clust2ii
        prefix <- "easy_"
    } else {
        B2 <- B2i
        clust2 <- clust2i
        prefix <- ""
    }

    # Uncorrected.
    raw.all <- cbind(B1, B2)
    clust.cols <- c(clust1, clust2)
    batch.id <- rep(1:2, c(ncol(B1), ncol(B2)))
    raw.fitted <- fitSVM(raw.all, batch.id)   

    # MNN corrected (default parameters at time of testing).
    Xmnn <- mnnCorrect(B1, B2, k=20, sigma=1, cos.norm.in=FALSE, cos.norm.out=FALSE, var.adj=TRUE)
    corre <- cbind(Xmnn$corrected[[1]],Xmnn$corrected[[2]])
    mnn.fitted <- fitSVM(corre, batch.id)   

    # limma.
    Xlm <- removeBatchEffect(raw.all, factor(batch.id))
    lm.fitted <- fitSVM(Xlm, batch.id)

    # ComBat.
    cleandat <- ComBat(raw.all, factor(batch.id), mod=NULL, prior.plots = FALSE)
    com.fitted <- fitSVM(cleandat, batch.id)

    # Creating a plot of the classification errors.
    pdf(fig.path("figs", paste0(ifelse(easy, "easy_", ""), "classification.pdf")))
    barplot(c(Uncorrected=raw.fitted$best.performance,
              MNN=mnn.fitted$best.performance,
              limma=lm.fitted$best.performance,
              ComBat=com.fitted$best.performance),
            ylab="Classification error")
    dev.off()
}

