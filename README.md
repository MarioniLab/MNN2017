# MNN

This folder contains R codes used to generate figures 1-5 in the manuscript "Correcting batch effects in single-cell RNA sequencing data by matching mutual nearest neighbours", as well as 1-6 supplementary figures.

1. To generate the simulation figure in the main text, first run the source file `simulateBatches.R`, then run the source file `FourPlots_sim.R` in the Simulations folder.

2. To generate the simulation figure in the supplement (identical composition of cell types), first run the source file `easysimulateBatches.R`, then run the source file `FourPlots_easysim.R` in the Simulations folder  


3.To generate the haematopoietic data figures in the main text, first run the source file `preparedata.R`, then  run the source file `FourPlots_haem.R` in the Haematopoiesis folder.
4.Then, to generate the haematopoietic data figures in the supplement, run the source file `PCAs.R` in the Haematopoiesis folder.


5. First download gene expression data from the four public data sets (Gene X Cell matrices, meta data and highly variable gene lists) by running the bash script `DownloadData.sh`.  This will download a directory 'ftp2.cruk.cam.ac.uk', which contains two sub-directories; CELseq and Smartseq2.  These directories contain the data for the appropriate studies (denoted by GEO/ArrayExpress accession number detailed in the manuscript).
6. To generate any of the pancreas data sets results, you have to run the source file `preparedata.R` in the Pancreas folder first.
7. To correct the batch effects and generate t-SNE plots and the Silhouette boxplots for the pancreas data sets (Fig 4), run the source file `FourPlots_panc.R` in the Pancreas folder.
8. Then, to generate the pancreas PCA plots and the entropy of mixings boxplots in the supplement (Suppl. Fig.3), run the source file `entropyandPCAs.R` in the Pancrease folder.
9. To compare performance of MNN with locally variable batch effects versus a global batch effect settings (Suppl.Fig. 6), run the source file `local_global_batchvect.R` in the Pancrease folder.
10. Differential expression (Mike)


