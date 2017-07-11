#! /bin/bash

# this will be the download script to fetch the various normalised pancreas datasets from an FTP
# a directory ftp2.cruk.cam.ac.uk containing two sub-directories will be downloaded: CELseq and Smartseq2 which contain the following files:
# normalised log2 expression  - GSE81076-norm.tsv.gz, GSE85241-norm.tsv.gz, GSE86473-norm.tsv.gz, E-MTAB-5061-norm.tsv.gz
# highly variable genes - GSE81076-HVG.tsv, GSE85241-HVG.tsv, GSE86473-HVG.tsv, E-MTAB-5061-HVG.tsv
# meta data - GSE81076-marker_metadata.tsv, GSE85241-marker_metadata.tsv, GSE86473-marker_metadata.tsv, E-MTAB-5061-marker_metadata.tsv

wget -r --user jmlabftp --password HOBICAmeer6 ftp2.cruk.cam.ac.uk:/*