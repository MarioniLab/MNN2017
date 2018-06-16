mkdir -p raw_data

# Downloading the 68K PBMC dataset.

curl -L -o raw_data/68k.tar.gz http://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz
cd raw_data
tar -xvf 68k.tar.gz
mv filtered_matrices_mex pbmc68k
rm 68k.tar.gz
cd -

# Downloading the 4K T-cell dataset.

curl -L -o raw_data/4k.tar.gz http://cf.10xgenomics.com/samples/cell-exp/2.1.0/t_4k/t_4k_filtered_gene_bc_matrices.tar.gz
cd raw_data
tar -xvf 4k.tar.gz
mv filtered_gene_bc_matrices t4k
rm 4k.tar.gz
cd -
