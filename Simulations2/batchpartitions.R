batchpartitions<-function(data,method,nclusters) { 
#batchpartitions<-function(data) { 
#this function partitions data to nclusters using kmeans on tsne   
#data G*n
# method  "DynamictTreeCut" or "kmeans  
# nclusters for kmeans on tsne, method=="DynamictTreeCut", does not need this input
  
  require(WGCNA)
  partitioned<-list()
  colidx<-list()
  count<-0
  
 
    all.dists <- (dist(t(data)))
  
    
    
    if (method=="kmeans") {
    require(Rtsne)
    set.seed(0)
    tsne<-Rtsne(all.dists, is_distance=TRUE)
    if(is.null(nclusters)) {
      nclusters <- 4 }
    ki<-kmeans(tsne$Y,nclusters)
    cluster<-ki$cluster
    }
    
    else if (method=="DynamictTreeCut") {
    set.seed(0)
    dendro<-hclust(all.dists)
    Ccut<-cutreeDynamic(dendro,cutHeight = NULL, minClusterSize = 100, method="tree")
    cluster<-Ccut+1
    nclusters<-max(cluster)
    }
    
    for (j in 1:nclusters){
      count<-count+1
     partitioned[[paste0("part",count)]]<-data[,cluster==j]
     colidx[[paste0("part",count)]]<-which(cluster==j)
    }
    
  #return(partitioned)
  return(colidx)
}