# this script is designed to run on linux and mac only, system calls to windows will fail(!)
# sourcing or executing this script in an open R session will generate the normalized gene expression
# matrices used for batch correction, along with the appropriate meta data for each study.

library(scran)
library(biomaRt)
library(limSolve)
library(scater)
library(SingleCellExperiment)

##############
## GSE81076 ##
##############
# download file from GEO
gse81076 <- 'GSE81076_D2_3_7_10_17.txt.gz'
if (!file.exists(gse81076)) { download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81076/suppl/GSE81076%5FD2%5F3%5F7%5F10%5F17%2Etxt%2Egz", 
                                            gse81076)}
gse81076.df <- read.table(gse81076, sep='\t', h=T, stringsAsFactors=F)
gse81076.ndim <- dim(gse81076.df)[2]

# construct the meta data from the cell names
donor.names <- unlist(regmatches(colnames(gse81076.df)[2:gse81076.ndim],
                                 gregexpr(pattern="D[0-9]{1,2}", 
                                          colnames(gse81076.df)[2:gse81076.ndim])))

plate.id <- unlist(lapply(strsplit(unlist(lapply(strsplit(head(colnames(gse81076.df)[2:gse81076.ndim]),
                                                          split="D[0-9]{1,2}", perl=T), 
                                                 FUN=function(x) paste0(x[2]))),
                                   fixed=T, split="_"),
                          FUN=function(c) paste0(c[1])))

protocol.id <- rep('CELseq', gse81076.ndim-1)
study.id <- rep('GSE81076', gse81076.ndim-1)
gse81076.meta <- data.frame(list('Donor' = donor.names,
                                 'Plate' = plate.id,
                                 'Protocol' = protocol.id,
                                 'Study' = study.id,
                                 'Sample' = colnames(gse81076.df)[2:gse81076.ndim]))
rownames(gse81076.meta) <- colnames(gse81076.df)[2:gse81076.ndim]
colnames(gse81076.df) <- gsub(colnames(gse81076.df), pattern='X', replacement='gene')

# remove superfluous suffixes from gene IDs
gse81076.df$gene <- gsub(gse81076.df$gene,
                         pattern="__chr[0-9]+", replacement="")
# remove the duplicated gene names
gse81076.df <- gse81076.df[!duplicated(gse81076.df$gene), ]

rownames(gse81076.df) <- gse81076.df$gene

# remove the gene ID column for downstream normalization
gse81076.df <- gse81076.df[, -1]

# remove cells and genes with all 0's
gene_sparsity <- (apply(gse81076.df == 0, MARGIN = 1, sum)/dim(gse81076.df)[2])
keep_genes <- gene_sparsity < 0.9
gse81076.nz <- gse81076.df[keep_genes, ]

cell_sparsity <- apply(gse81076.nz == 0, MARGIN = 2, sum)/dim(gse81076.nz)[1]
keep_cells <- cell_sparsity < 0.8
gse81076.nz <- gse81076.nz[, keep_cells]
gse81076.nz <- apply(gse81076.nz, 2, as.integer)

# use the spike in genes to estimate size factors for normalization
# all values show be non-negative integers
spikes <- grepl(rownames(gse81076.df[keep_genes, ]),
                pattern='ERCC')

sce <- SingleCellExperiment(list(counts = as.matrix(gse81076.nz)))
sce <- calculateQCMetrics(sce, feature_controls=list(Spikes=spikes))
isSpike(sce) <- spikes

clusters <- quickCluster(sce, get.spikes=TRUE, min.size=120)
sce <- computeSumFactors(sce, sizes=c(10, 20, 40, 60), positive=T,
                         assay.type='counts', clusters=clusters)
summary(sizeFactors(sce))
sce <- normalize(sce)

gse81076.norm <- data.frame(exprs(sce))

gse81076.norm$gene_id <- rownames(gse81076.df[keep_genes, ])
gse81076.norm$gene_id <- gsub(gse81076.norm$gene_id,
                              pattern="__chr[0-9X]+", replacement="")

write.table(gse81076.norm, sep='\t',
             file='Pancreas/Data/GSE81076_SFnorm.tsv',
             quote=FALSE, row.names=FALSE, col.names=TRUE)

write.table(gse81076.meta, sep="\t",
            file="Pancreas/Data/GSE81076_metadata.tsv",
            quote=FALSE, row.names=F, col.names=TRUE)


##############
## GSE85241 ##
##############
# clear environment and invoke garbage collector
rm(list=ls())
gc()

gse85241 <- 'GSE85241_cellsystems_dataset_4donors_updated.csv'
if (!file.exists(gse85241)) { download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85241/suppl/GSE85241%5Fcellsystems%5Fdataset%5F4donors%5Fupdated%2Ecsv%2Egz", 
                                            gse85241)}

gse85241.df <- read.table(gse85241, sep='\t', h=T, stringsAsFactors=F)
gse85241.df$gene_id <- rownames(gse85241.df)

# gene IDs are located in column X for these data
colnames(gse85241.df) <- gsub(colnames(gse85241.df), pattern="X",
                              replacement="gene_id")

donor.id <-unlist(lapply(strsplit(colnames(gse85241.df)[1:(dim(gse85241.df)[2]-1)],
                                  fixed=T, split="."), 
                         FUN=function(x) paste0(x[1])))

plate.id <- unlist(lapply(strsplit(colnames(gse85241.df)[1:(dim(gse85241.df)[2]-1)],
                                   fixed=T, split="."),
                          FUN=function(x) paste0(x[2])))

protocol.id <- rep('CELseq2', dim(gse85241.df)[2]-1)
study.id <- rep('GSE85241', dim(gse85241.df)[2]-1)

gse85241.meta <- data.frame(list('Donor' = donor.id,
                                 'Plate'= plate.id,
                                 'Protocol' = protocol.id,
                                 'Study' = study.id,
                                 'Sample' = colnames(gse85241.df)[1:(dim(gse85241.df)[2]-1)]))

rownames(gse85241.meta) <- colnames(gse85241.df)[1:(dim(gse85241.df)[2]-1)]

write.table(gse85241.meta,
            file="Pancreas/Data/GSE85241_metadata.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

# set gene IDs as rownames, remove gene ID column
# remove duplicated gene IDs

gse85241.df$gene_id <- gsub(gse85241.df$gene_id,
                            pattern="__chr[0-9X]+", replacement="")
gse85241.df <- gse85241.df[!duplicated(gse85241.df$gene_id), ]
rownames(gse85241.df) <- gse85241.df$gene_id
gse85241.df <- gse85241.df[, 1:(dim(gse85241.df)[2]-1)]

# remove cells and genes with all 0's
gene_sparsity <- (apply(gse85241.df == 0, MARGIN = 1, sum)/dim(gse85241.df)[2])
keep_genes <- gene_sparsity < 0.9
gse85241.nz <- gse85241.df[keep_genes, ]

cell_sparsity <- apply(gse85241.nz == 0, MARGIN = 2, sum)/dim(gse85241.nz)[1]
keep_cells <- cell_sparsity < 0.8
dim(gse85241.nz[, keep_cells])
gse85241.nz <- gse85241.nz[, keep_cells]
gse85241.nz <- apply(gse85241.nz, 2, as.integer)

spikes <- grepl(rownames(gse85241.df[keep_genes, ]),
                pattern='ERCC')
sce <- SingleCellExperiment(list(counts = as.matrix(gse85241.nz)))
sce <- calculateQCMetrics(sce, feature_controls=list(Spikes=spikes))
isSpike(sce) <- spikes

clusters <- quickCluster(sce, get.spikes=TRUE, min.size=120)
sce <- computeSumFactors(sce, sizes=c(10, 20, 40, 60), positive=T,
                         assay.type='counts', clusters=clusters)
summary(sizeFactors(sce))

sce <- normalize(sce)
gse85241.norm <- data.frame(exprs(sce))
gse85241.norm$gene_id <- rownames(gse85241.df[keep_genes, ])

write.table(gse85241.norm, sep='\t',
            file='Pancreas/Data/GSE85241_SFnorm.tsv',
            quote=F, row.names=F, col.names=T)

##############
## GSE86473 ##
##############
# clear environment and invoke garbage collector
rm(list=ls())
gc()

# the raw/processed counts table was not available for download from GEO for this data set
# these data were generated by mapping the original fastq's back to mm10, then using featureCounts to quantify
# against mm10 ensembl build 86
# data are contained in Pancreas/RawData for each cell type
alpha <- 'Pancreas/RawData/alpha-feature_counts.tsv.gz'
beta <- 'Pancreas/RawData/beta-feature_counts.tsv.gz'
delta <- 'Pancreas/RawData/delta-feature_counts.tsv.gz'
pp <- 'Pancreas/RawData/PP-feature_counts.tsv.gz'
# qc_out <- 'Pancreas/RawData/pancreas-smarseq2-qcout.tsv'

alpha.df <- read.table(alpha, sep='\t', h=T, stringsAsFactors=F)
beta.df <- read.table(beta, sep='\t', h=T, stringsAsFactors=F)
delta.df <- read.table(delta, sep='\t', h=T, stringsAsFactors=F)
pp.df <- read.table(pp, sep='\t', h=T, stringsAsFactors=F)

# qc.df <- read.table(qc_out, sep='\t', h=F, stringsAsFactors=F)
data.list <- list("alpha"=alpha.df, "beta"=beta.df, "delta"=delta.df,
                  "gamma"=pp.df)

gse86473.df <- Reduce(x=data.list,
                      f=function(x, y) merge(x, y, by='gene_id'))

# remove .dedup suffix
samp.names <- unlist(lapply(strsplit(colnames(gse86473.df), fixed=T,
                                     split="."), 
                            FUN=function(x) paste0(x[1])))
colnames(gse86473.df) <- tolower(samp.names)

gse86473.meta <- read.table('Pancreas/RawData/GSE86473_experimental_design.tsv',
                            sep='\t', h=T, stringsAsFactors=F)
gse86473.meta$Sample <- tolower(gse86473.meta$Sample)
gse86473.meta$Study <- "GSE86473"

# capitalize first letter of cell labels

gse86473.meta$CellType <- paste(toupper(substr(gse86473.meta$CellType, 1, 1)),
                                substr(gse86473.meta$CellType, 2, nchar(gse86473.meta$CellType)), sep="")

write.table(gse86473.meta,
            "Pancreas/Data/GSE86473_metadata.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

rownames(gse86473.df) <- gse86473.df$gene_id

# need to map ensemlb gene ids to hgnc symbols to match up with other Pancreas data sets
ensembl <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl', GRCh=37)
gene_symbol <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                     filters='ensembl_gene_id', mart=ensembl,
                     values=gse86473.df$gene_id)
dat <- merge(gse86473.df, gene_symbol,
             by.x='gene_id', by.y='ensembl_gene_id')
gene_symbols <- dat$external_gene_name
gse86473.df <- as.data.frame(append(dat[, -1],
                                    list(gene_id=gene_symbols),
                                    after=0))

# unset hgnc_symbol as a factor and remove duplicated IDs
gse86473.df$gene_id <- as.character(gse86473.df$gene_id)
gse86473.df <- gse86473.df[!duplicated(gse86473.df$gene_id), ]

# set the hgnc symbols as rownames and size factor normalize counts table
rownames(gse86473.df) <- gse86473.df$gene_id
gse86473.df <- gse86473.df[, 2:(dim(gse86473.df)[2]-1)]

# remove cells and genes with all 0's
gene_sparsity <- (apply(gse86473.df == 0, MARGIN = 1, sum)/dim(gse86473.df)[2])
keep_genes <- gene_sparsity < 0.9
dim(gse86473.df[keep_genes, ])
gse86473.nz <- gse86473.df[keep_genes, ]

cell_sparsity <- apply(gse86473.nz == 0, MARGIN = 2, sum)/dim(gse86473.nz)[1]
keep_cells <- cell_sparsity < 0.8
dim(gse86473.nz[, keep_cells])
gse86473.nz <- gse86473.nz[, keep_cells]
gse86473.nz <- apply(gse86473.nz, 2, as.integer)

sce <- SingleCellExperiment(list(counts = as.matrix(gse86473.nz)))
sce <- calculateQCMetrics(sce)
clusters <- quickCluster(sce, min.size=120)
sce <- computeSumFactors(sce, sizes=c(10, 20, 40, 60), positive=T,
                         assay.type='counts', clusters=clusters)
summary(sizeFactors((sce)))

sce <- normalize(sce)
gse86473.norm <- data.frame(exprs(sce))
gse86473.norm$gene_id <- rownames(gse86473.df[keep_genes, ])

write.table(gse86473.norm, sep='\t',
            file='Pancreas/Data/GSE86473_SFnorm.tsv',
            quote=F, row.names=F, col.names=T)

#################
## E-MTAB-5061 ##
#################
# clear environment and invoke garbage collector
rm(list=ls())
gc()

# the download file contains columns of RPKM & counts
# need to pull out just the integer gene count columns
emtab_combined <- "Pancreas/RawData/EMTAB5061_rpkm_counts.txt.zip"
if(!file.exists(emtab_combined)){ download.file("https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/files/E-MTAB-5061.processed.1.zip",
              emtab_combined)}

# need to unzip via a system call
unzip.file <- "Pancreas/RawData/pancreas_refseq_rpkms_counts_3514sc.txt"
if(!file.exists(unzip.file)) {system(paste0("unzip ", emtab_combined, " -d ", "Pancreas/RawData/"))}

emtab.df <- read.table(unzip.file,
                       h=FALSE, sep="\t", stringsAsFactors=F)

col.names <- unlist(read.table("Pancreas/RawData/pancreas_refseq_rpkms_counts_3514sc.txt",
                        h=FALSE, sep="\t", stringsAsFactors=F, comment.char="", nrows = 1))

# first 2 columns are gene symbol and NCBI ID
# because of the way that this table is constructed, there are two cells for each column, but they
# are not contiguous.  Therefore, the sample names need to be read in separately to the counts/rpkm table

emtab5061.df <- emtab.df[, c(1, 3517:dim(emtab.df)[2])]
colnames(emtab5061.df) <- gsub(col.names, pattern="#samples", replacement="gene_id")

# download sdrf file direct from arrayExpress
emtab.file <- "Pancreas/RawData/E-MTAB-5061.sdrf.txt"
if(!file.exists(emtab.file)) {download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/E-MTAB-5061.sdrf.txt",
              emtab.file)}

emtab.sdrf <- read.table("Pancreas/RawData/E-MTAB-5061.sdrf.txt",
                         h=TRUE, sep="\t", stringsAsFactors=FALSE)

# construct the appropriate meta data columns, i.e. donor, plate, protocol, study
emtab.meta <- emtab.sdrf[, c("Assay.Name", "Characteristics.cell.type.", "Characteristics.individual.")]
colnames(emtab.meta) <- c("Sample", "CellType", "Donor")
emtab.meta$Study <- "E-MTAB-5061"
emtab.meta$Protocol <- "SmartSeq2"

# remove the marked low quality cells
remove.cells <- unique(emtab.sdrf$Assay.Name[emtab.sdrf$Characteristics.single.cell.well.quality. == "low quality cell"])
emtab5061.df <- emtab5061.df[, !colnames(emtab5061.df) %in% remove.cells]
emtab.meta <- emtab.meta[!emtab.meta$Sample %in% remove.cells, ]
rownames(emtab.meta) <- emtab.meta$Sample

emtab.meta$CellType <- gsub(emtab.meta$CellType,
                            pattern=" cell", replacement="")

emtab.meta$CellType <- paste(toupper(substr(emtab.meta$CellType, 1, 1)),
                             substr(emtab.meta$CellType, 2, nchar(emtab.meta$CellType)), sep="")

write.table(emtab.meta,
            file="Pancreas/Data/E-MTAB-5061_metadata.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

# remove duplicated gene IDs and set to rownames
emtab5061.df <- emtab5061.df[!duplicated(emtab5061.df$gene_id), ]
rownames(emtab5061.df) <- emtab5061.df$gene_id

emtab5061.df <- emtab5061.df[, -1]

# remove cells and genes with all 0's
gene_sparsity <- (apply(emtab5061.df == 0, MARGIN = 1, sum)/dim(emtab5061.df)[2])
keep_genes <- gene_sparsity < 0.9
dim(emtab5061.df[keep_genes, ])
emtab5061.nz <- emtab5061.df[keep_genes, ]

cell_sparsity <- apply(emtab5061.nz == 0, MARGIN = 2, sum)/dim(emtab5061.nz)[1]
keep_cells <- cell_sparsity < 0.8
dim(emtab5061.nz[, keep_cells])
emtab5061.nz <- emtab5061.nz[, keep_cells]
emtab5061.nz <- apply(emtab5061.nz, 2, as.integer)

spikes <- grepl(x=rownames(emtab5061.df[keep_genes, ]), pattern="ERCC")
sce <- SingleCellExperiment(list(counts = as.matrix(emtab5061.nz)))
sce <- calculateQCMetrics(sce, feature_controls=list(Spikes=spikes))
isSpike(sce) <- spikes

clusters <- quickCluster(sce, min.size=120)
sce <- computeSumFactors(sce, sizes=c(10, 20, 40, 60), positive=T,
                         assay.type='counts', clusters=clusters)
summary(sizeFactors((sce)))
sce <- normalize(sce)
emtab.norm <- data.frame(exprs(sce))
emtab.norm$gene_id <- rownames(emtab5061.df[keep_genes, ])

write.table(emtab.norm,
            file="Pancreas/Data/E-MTAB-5061_SFnorm.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

# clear the R environment in case this script is directly sourced
rm(list=ls())
gc()

