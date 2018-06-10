# Download the counts, metadata of Nestorowa et al. 2016

curl -L -o raw_data/GSE81682_HTSeq_counts.txt.gz https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81682&format=file&file=GSE81682%5FHTSeq%5Fcounts%2Etxt%2Egz
curl -L -o raw_data/metaF.txt http://blood.stemcells.cam.ac.uk/data/all_cell_types.txt

# Download the counts from Paul et al. 2015

curl -L -o raw_data/umitab_Amit.txt.gz https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE72857&format=file&file=GSE72857%5Fumitab%2Etxt%2Egz

# Download the HVG list.

curl -L -o raw_data/coordinates_gene_counts_flow_cytometry.txt.gz http://blood.stemcells.cam.ac.uk/data/coordinates_gene_counts_flow_cytometry.txt.gz
