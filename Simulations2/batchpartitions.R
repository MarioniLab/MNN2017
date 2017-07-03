batchpartitions<-function(data,nclusters) { 
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
  
    require(Rtsne)
    set.seed(0)
    tsne<-Rtsne(all.dists, is_distance=TRUE)
    if(is.null(nclusters)) {
      nclusters <- 4 }
    ki<-kmeans(tsne$Y,nclusters)
    cluster<-ki$cluster

    for (j in 1:nclusters){
     count<-count+1
     partitioned[[paste0("part",count)]]<-data[,cluster==j]
     colidx[[paste0("part",count)]]<-which(cluster==j)
    }
    
  #return(partitioned)
  return(colidx)
}