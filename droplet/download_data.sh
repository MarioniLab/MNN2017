mkdir -p raw_data

# Downloading the 68K PBMC dataset.

curl -L -o raw_data/68k.tar.gz http://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz
cd raw_data
tar -xvf 68k.tar.gz
rm 68k.tar.gz
cd -
