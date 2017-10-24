# find highly variable genes in each data set separately a la Brennecke et al
# the assumption is that the normalizePancreas.R script has already been run so that the relevant files
# are in the appropriate directories.  If these files do not exist, then please run normalizePancreas.R
# first, or simply use source("Pancreas/normalizePancreas.R").  The objects from that script will NOT
# be retained in the environment if you use this approach, so files will be explicitly loaded.

# define a convenience function for HVG definition
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
  
  useForFit <- recip.means <= 0.1
  
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

##############
## GSE81076 ##
##############
gse81076.norm <- read.table("Pancreas/Data/GSE81076_SFnorm.tsv",
                            h=TRUE, sep="\t", stringsAsFactors=FALSE)

# remove gene id column for hvg discovery, last column
rownames(gse81076.norm) <- gse81076.norm$gene_id

# output is a boolean vector of the same order as the input matrix
gse81076.HVG <- find_hvg(dataframe=gse81076.norm[, 1:(dim(gse81076.norm)[2]-1)], 
                         p.threshold=0.05, plot=TRUE, return.ranks=FALSE)

# select the highly variable genes from the input dataframe column 'gene id'
gse81076.hvg_df <- cbind.data.frame(names(gse81076.HVG)[gse81076.HVG])
colnames(gse81076.hvg_df) <- "gene_id"

write.table(gse81076.hvg_df,
            file="Pancreas/Data/GSE81076-HVG.tsv",
            row.names=FALSE, quote=FALSE, sep="\t")

##############
## GSE85241 ##
##############
gse85241.norm <- read.table("Pancreas/Data/GSE85241_SFnorm.tsv",
                            h=TRUE, sep="\t", stringsAsFactors=FALSE)

# remove gene id column for hvg discovery, last column
rownames(gse85241.norm) <- gse85241.norm$gene_id


# output is a boolean vector of the same order as the input matrix
gse85241.HVG <- find_hvg(dataframe=gse85241.norm[, 1:(dim(gse85241.norm)[2]-1)], 
                         p.threshold=0.05, plot=FALSE, return.ranks=FALSE)

# select the highly variable genes from the input dataframe column 'gene id'
gse85241.hvg_df <- cbind.data.frame(names(gse85241.HVG)[gse85241.HVG])
colnames(gse85241.hvg_df) <- "gene_id"

write.table(gse85241.hvg_df,
            file="Pancreas/Data/GSE85241-HVG.tsv",
            row.names=FALSE, quote=FALSE, sep="\t")


##############
## GSE86473 ##
##############
gse86473.norm <- read.table("Pancreas/Data/GSE86473_SFnorm.tsv",
                            h=TRUE, sep="\t", stringsAsFactors=FALSE)

# remove gene id column for hvg discovery, last column
rownames(gse86473.norm) <- gse86473.norm$gene_id


# output is a boolean vector of the same order as the input matrix
gse86473.HVG <- find_hvg(dataframe=gse86473.norm[, 1:(dim(gse86473.norm)[2]-1)], 
                         p.threshold=0.05, plot=FALSE, return.ranks=FALSE)

# select the highly variable genes from the input dataframe column 'gene id'
gse86473.hvg_df <- cbind.data.frame(names(gse86473.HVG)[gse86473.HVG])
colnames(gse86473.hvg_df) <- "gene_id"

write.table(gse86473.hvg_df,
            file="Pancreas/Data/GSE86473-HVG.tsv",
            row.names=FALSE, quote=FALSE, sep="\t")

#################
## E-MTAB-5061 ##
#################
emtab5061.norm <- read.table("Pancreas/Data/E-MTAB-5061_SFnorm.tsv",
                            h=TRUE, sep="\t", stringsAsFactors=FALSE)

# remove gene id column for hvg discovery, last column
rownames(emtab5061.norm) <- emtab5061.norm$gene_id


# output is a boolean vector of the same order as the input matrix
emtab5061.HVG <- find_hvg(dataframe=emtab5061.norm[, 1:(dim(emtab5061.norm)[2]-1)], 
                         p.threshold=0.05, plot=FALSE, return.ranks=FALSE)

# select the highly variable genes from the input dataframe column 'gene id'
emtab5061.hvg_df <- cbind.data.frame(names(emtab5061.HVG)[emtab5061.HVG])
colnames(emtab5061.hvg_df) <- "gene_id"

write.table(emtab5061.hvg_df,
            file="Pancreas/Data/E-MTAB-5061-HVG.tsv",
            row.names=FALSE, quote=FALSE, sep="\t")

# clear the environment
rm(list=ls())
gc()


