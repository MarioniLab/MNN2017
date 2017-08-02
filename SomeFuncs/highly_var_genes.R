#counts: matrix with cells in columns and genes in rows (only endogenous genes, spike-ins to be removed)
#genes: vector with genes that are detected (have at least 1 count) in at least one cell



# result<-find.high.var.biol.genes(counts = counts, genes = genes, 
#                                  minMean=10,#only genes with normalized counts above this treshold will be used for fitting. This threshold should be such that it does not intersect the plateau at low expression levels (see plot in Brennecke's paper)
#                                  plot=TRUE,#option to show a figure of CV2 vs mean normalized counts
#                                  red.line = FALSE, #option to display a red vertical line in correspondence to the threshold chosen above
#                               p.adj.thr = 0.1 #threshold for significance (adjusted p-value, Benjamini&Hochberg correction)
#                               )
#The output "result" is a list with two elements:
#p.adj: adjusted p-value for each gene
#genes.high.var: genes with a p-value < p.adj.thr
##############################################

find.high.var.biol.genes<-function(counts, genes, red.line=TRUE, p.adj.thr,minMean=10, plot=FALSE){
  
  require(DESeq2)

  #compute size factors (with biological genes)
  #sf.genes<-estimateSizeFactorsForMatrix(raw.counts) #estimate size factor
  #data<-t( t(raw.counts) / sf.genes)
  #raw.counts<-raw.counts[genes,]
  data<-counts[genes,]
  
  #estimate the sample moments 
  means<-rowMeans(data)
  vars<-t(apply(data, 1, function(x) var((x))))
  
 
  names(vars)<-names(means)
  cv2<-vars/means^2
  names(cv2)<-names(vars)
  
  minMeanForFit <- minMean
  minMeanForFit
  
  
  if(plot==TRUE){
    # Next, we define a minimum mean value to exclude genes with low mean and hence high CV from fit as they would
    # otherwise skew it downwards.
    
    #plot
    par(cex=1.5)
    plot(means, cv2, 
         log='xy',  
         xlab="Average Normalized Read Counts", 
         ylab=expression(CV^2))
    abline(v=minMeanForFit, lwd=2, col="red", lty=2)
  }
  
  require(statmod)
  useForFit<-means >= minMeanForFit
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ), 
                     cv2[useForFit] ) 
  fit$coefficients
  
  
  xi <- mean( 1 / sf.genes)#sfHeLa )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"] - xi )
  c( a0, a1 )
  
  
  #find variable genes
  df <- ncol(data)-1
  
  chi2.values<-as.vector(df*(cv2/( (xi+a1)/means + a0 )))
  names(chi2.values)<-names(means)
  
  p.unadj<-pchisq(chi2.values, df, lower.tail=FALSE)
  
  p.adj <- p.adjust( p.unadj, "BH" )
  
  
  
  
  if(plot==TRUE){
    plot(means, cv2, 
         log='xy',  
         xlab=expression(bold("Average Normalized Read Counts")), 
         ylab=expression(bold(CV^2)),
         col = ifelse( p.adj < p.adj.thr, "red", "black" ))
    
    xg <- 10^seq( -2, 5, length.out=1000 )
    lines( xg, (xi+a1)/xg + a0, col="green", lwd=3 )
    
    if(red.line==TRUE){
      abline(v=minMeanForFit, lwd=2, col="red", lty=2)
    }
  }
  
  genes.high.var<-names(p.adj[p.adj<p.adj.thr])
  
  res<-list()
  
  res[[1]]<-p.adj
  res[[2]]<-genes.high.var
  
  names(res)<-c("p.adj", "genes.high.var")
  
  return(res)
  
  
}
