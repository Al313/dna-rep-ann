# BiocManager::install("biomaRt")

library(biomaRt)
library(stringr)
library(readxl)

`%notin%` <- Negate(`%in%`)


## This file was obtained from this article: https://doi.org/10.1038/s43018-020-0050-6
genes <- read_xlsx("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/43018_2020_50_MOESM3_ESM.xlsx", sheet = 7, skip = 1)



genes_of_interest <- unique(genes$Gene_somatic)
length(genes_of_interest)


# GRC37


#### Genes

grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
hg19 <- getBM(mart = grch37, attributes = c("ensembl_gene_id", "chromosome_name","start_position", "end_position", "hgnc_symbol", "hgnc_id"))

hg19_pruned <- hg19[nchar(hg19$chromosome_name) <=2,]
unique(hg19_pruned$chromosome_name)


hg19_interested <- hg19_pruned[hg19_pruned$hgnc_symbol %in% genes_of_interest,]
head(hg19_interested)

# 4 genes are not annotated
nrow(hg19_interested)
length(genes_of_interest)
genes_of_interest[genes_of_interest %notin% hg19_interested$hgnc_symbol]

# By pruning the extra chromosome atches we don't get any duplicates
hg19_interested[duplicated(hg19_interested$hgnc_symbol),]
hg19_interested[hg19_interested$hgnc_symbol == "MDC1",]


# Ordering column names
hg19_interested <- hg19_interested[,c(2,3,4,6,5,1)]
rownames(hg19_interested) <- 1:nrow(hg19_interested)

# toy_hg19_interested <- hg19_interested[hg19_interested$hgnc_symbol %in% c("MUTYH", "OGG1"),]


write.table(hg19_interested, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/genes.txt", sep = "\t", quote = F, row.names = F)




#### Exons



ensembl = useMart("ensembl", host = "grch37.ensembl.org")
ensembl = useDataset("hsapiens_gene_ensembl",mart = ensembl)
exon <- getBM(mart = ensembl, attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end","ensembl_exon_id", "ensembl_gene_id"))
exon_comb <- merge(exon, hg19_pruned[,c(1,3,4,5,6)], by = 'ensembl_gene_id', all.x = T)
head(exon)
head(hg19_pruned)
exon_interested <- exon_comb[exon_comb$hgnc_symbol %in% genes_of_interest,]


# Ordering column names
head(exon_interested)
exon_interested <- exon_interested[,c(2,3,4,5,8,1,6,7)]
rownames(exon_interested) <- 1:nrow(exon_interested)
colnames(exon_interested) <- c("chrom", "exon_start", "exon_end", "ensembl_exon_id", "hgnc_symbol", "ensembl_gene_id", "gene_start", "gene_end")

str(exon_interested)
# exon[exon$ensembl_gene_id == "ENSG00000114026",]
nrow(exon[exon$ensembl_gene_id == "ENSG00000114026",])


write.table(exon_interested, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/exons.txt", sep = "\t", quote = F, row.names = F)











# GRC38


mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
hsp <- getBM(mart = mart, attributes = c("ensembl_gene_id", "chromosome_name","start_position", "end_position", "hgnc_symbol", "hgnc_id"))
att <- listAttributes(mart,what = c("name","description","page"))
att$name[str_detect(att$name, pattern = "hgnc")]


mart2 = useMart('ensembl')
listDatasets(mart2)[str_detect(listDatasets(mart2)$dataset,pattern = "hsapiens"),]



ensembl38=useMart("ensembl")
ensembl38 = useDataset("hsapiens_gene_ensembl",mart=ensembl38)
filterlist <- list("3:9091628:10029903")
exon38 <- getBM(mart = ensembl38, filters = c("chromosomal_region"), values = filterlist, attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end","transcript_length","strand", "ensembl_gene_id", "ensembl_transcript_id","ensembl_exon_id"))
exon38[exon38$ensembl_gene_id == "ENSG00000114026",]
nrow(exon38[exon38$ensembl_gene_id == "ENSG00000114026",])

exon38.1 <- getBM(mart = ensembl38, filters = c("chromosomal_region"), values = filterlist, attributes = c("hgnc_symbol","ensembl_gene_id"))
results = merge(exon38,exon38.1,by='ensembl_gene_id',all.x=T)
?merge
