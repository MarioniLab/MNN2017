matchgenenamestoENS<-function(genesnames,idknown=TRUE){
  library(biomaRt)
  mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org" )
  
  if (idknown==TRUE){
    out = getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), values = genesnames, mart = mart,filters = "ensembl_gene_id")#
  }
  else{
    out = getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), values = genesnames, mart = mart,filters = "mgi_symbol")  
  }
  
  return(out)
}