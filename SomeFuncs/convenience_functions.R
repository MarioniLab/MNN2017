# convenience functions for MNN2017 manuscript
################################################################################################
################################################################################################
## Define a couple of convenience functions for normalisation and highly variable gene detection

find_hvg <- function(dataframe, plot=FALSE, p.threshold=1e-2, return.ranks=FALSE){
  # define a set of highly variable gene for the GFP+ and GFP- separately
  require(MASS)
  require(limSolve)
  require(statmod)
  # assume gene names are in the rows
  # even if they aren't this will still get
  # the input row ordering
  gene.names <- rownames(dataframe)
  means <- rowMeans(dataframe, na.rm = T)
  vars <- apply(dataframe, 1, var, na.rm=T)
  cv2 <- vars/(means^2)
  
  minMeanForFit <- unname(quantile(means[which(cv2 > 0.2)], 0.8))
  
  # select genes with mean value greater than min value for fitting
  # remove values with 1/means == infinite
  recip.means <- 1/means
  recip.means[is.infinite(recip.means)] <- 0
  
  useForFit <- recip.means <= 1
  
  # fit with a gamma-distributed GLM
  fit <- glmgam.fit(cbind(a0 = 1, a1tilde=recip.means[!useForFit]), cv2[!useForFit])
  
  # calculate % variance explained by the model fit
  resid.var <- var(fitted.values(fit) - cv2[!useForFit])
  total.var <- var(cv2[!useForFit])
  
  # get fitted values and mean-dispersion dependence line
  a0 <- unname(fit$coefficients["a0"])
  a1 <- unname(fit$coefficients["a1tilde"])
  
  xg <- seq(0, max(means[means != Inf]), length.out=100000)
  vfit <- (a1/xg) + a0
  
  # add confidence intervals
  d.f <- ncol(dataframe) - 1
  
  # rank genes by the significance of their deviation from the fit
  # to call HVGs
  a.fit <- (a1/means) + a0
  varFitRatio <- vars/(a.fit * means^2)
  varOrder <- order(varFitRatio, decreasing=T)
  
  oed <- dataframe[varOrder, ]
  
  if(plot == TRUE){
    smoothScatter(means, cv2, xlab="Mean expression", ylab=expression("CV"^2))
    lines(xg, vfit, col="black", lwd=3 )
    lines(xg, vfit * qchisq(0.975, d.f)/d.f, lty=2, col="black")
    lines(xg, vfit * qchisq(0.025, d.f)/d.f,lty=2,col="black")
    # display the 100 most highly variable genes
    points(means[varOrder[1:100]], cv2[varOrder[1:100]], col='red')
  }
  
  pvals <- pchisq(varFitRatio * d.f, d.f, lower.tail = F)
  pvals[is.na(pvals)] <- 1.0
  adj.pvals <- p.adjust(pvals, method='fdr')
  HVG <- adj.pvals <= p.threshold
  
  if(return.ranks){
    # order p-values, then subset past a threshold
    rank.p <- adj.pvals[order(adj.pvals, decreasing=FALSE)]
    order.names <- gene.names[order(adj.pvals, decreasing=FALSE)]
    thr.p <- rank.p <= p.threshold
    HVG <- order.names[thr.p]
  }
  
  return(HVG)
}


size_factor_normalize <- function(dataframe, cell.sparse=0.95, gene.sparse=0.99,
                                  cluster.size=40){
  # take an input gene X cell (barcode) dataframe of
  # read counts/UMIs
  # output the cell-specific size factor normalized log2 expression values
  
  n.cells <- dim(dataframe)[2]
  n.genes <- dim(dataframe)[1]
  
  gene_sparsity <- (apply(dataframe == 0, MARGIN = 1, sum)/n.cells)
  keep_genes <- gene_sparsity < gene.sparse
  #dim(exprs(SIGAD8.10x.hg19)[keep_genes, ])
  data.nz <- dataframe[keep_genes, ]
  
  # remove a cell with very low counts, order of magnitude lower than all others
  cell_sparsity <- apply(data.nz == 0, MARGIN = 2, sum)/dim(data.nz)[1]
  keep_cells <- cell_sparsity < cell.sparse
  #dim(pdx.nz[, keep_cells])
  data.nz <- data.nz[, keep_cells]
  
  # let's try to use the deconvolution normalisation approach on these data
  # if its a dgTMatrix, need to convert to a dgCMatrix for beachmat
  if(class(data.nz)[1] == "dgTMatrix"){
    data.nz <- as(data.nz, "dgCMatrix")
  }
  
  sce <- SingleCellExperiment(list(counts=data.nz))
  
  clusters <- quickCluster(sce, min.size=cluster.size, 
                           max.size=3*cluster.size, method="igraph")
  max.size <- floor(cluster.size/2)
  
  # change the window size in 10% increments
  size.inc <- max.size * 0.1
  sce <- computeSumFactors(sce, 
                           sizes=c(max.size - (3 * size.inc),
                                   max.size - (2 * size.inc),
                                   max.size - size.inc,
                                   max.size),
                           positive=FALSE,
                           assay.type='counts', clusters=clusters)
  sce <- normalise(sce)
  
  sf.norm <- as.data.frame(as(exprs(sce), "matrix"))
  sf.norm$gene_id <- rownames(sf.norm)
  
  return(sf.norm)
}

tsne_wrapper <- function(dataframe, perplexity=30, is.dist=FALSE){
  require(Rtsne)
  set.seed(42)
  if(is.dist){
    data.tsne <- Rtsne(t(dataframe),
                       perplexity=perplexity, is_distance=is.dist, pca=FALSE)
    sample.names <- labels(dataframe)
    
  }
  else {
    dataframe <- dataframe[, !duplicated(t(dataframe))]
    data.tsne <- Rtsne(t(dataframe), perplexity=perplexity)
    sample.names <- colnames(dataframe)
  }
  data.map <- data.frame(data.tsne$Y)  
  colnames(data.map) <- c("Dim1", "Dim2")
  data.map$Sample <- sample.names
  
  return(data.map)
}
