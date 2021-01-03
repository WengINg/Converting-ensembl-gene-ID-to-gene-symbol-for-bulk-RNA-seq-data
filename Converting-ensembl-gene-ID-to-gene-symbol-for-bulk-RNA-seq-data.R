# Generating a ensembl gene ID to gene symbol and gene name table for human
library(biomaRt)
## Connects to the selected BioMart database and dataset
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")   ## use mmusculus_gene_ensembl for mouse
## Retrieves information from the BioMart database
ensg2symbol_human <- getBM(attributes = c('ensembl_transcript_id', 'ensembl_gene_id',
                                          'hgnc_symbol', 'wikigene_description'),mart = mart, useCache = F)  ## use mgi_symbol for mouse
## take a look at the first few rows of the table
head(ensg2symbol_human)
## save this table for future use
save(ensg2symbol_human, file = "ensg2symbol_human.RData")  

# Example
## Converting the ensembl gene ID rownames of featureCounts table from human RNA-seq to gene symbols
## load in the ensembl gene to gene symbol table generated above
load("~/ensg2symbol_human.RData")
## load in the data table
count_mat <- read.table('Your_featureCounts_results.txt',header=T,sep="\t", stringsAsFactors = F)  
## separate the gene information (first 6 columns) from the counts of genes in samples (the rest of the columns)
gene_info<-count_mat[,c(1,6)]
count_mat<-count_mat[,-c(1:6)]
## set ensembl gene id (Geneid) as rownames of count table
rownames(count_mat)<-gene_info$Geneid
## match the ensembl gene id (rownames of count_mat) to the ensembl gene id in ensg2symbol_human table
ind <- match(rownames(count_mat), ensg2symbol_human$ensembl_gene_id)         
## set the new rownames of count_mat using the matched gene symbol from the ensg2symbol_human table
rownames(count_mat)<-ensg2symbol_human$hgnc_symbol[ind]
