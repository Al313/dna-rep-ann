
devtools::load_all('/home/ali313/Documents/studies/genomic-data-analyst/projects/pcawg-hmf-comparison/clones/geneDriverAnnotator/')

# Loading mutSigExtractor package

library(geneDriverAnnotator)

if (dir.exists("/hpc/cuppen/")){
  
  base_dir <- "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/"
  devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mutSigExtractor/'))
  
} else {
  
  library(mutSigExtractor)
}

`%notin%` <- Negate(`%in%`)




if (dir.exists("/hpc/cuppen/")){
  wd <- "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/"
} else {
  wd <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/"
}


if (dir.exists("/hpc/cuppen/")){
  metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2-metadata_08072021.tsv"), sep = "\t", header = T, stringsAsFactors = F)
} else {
  metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2-metadata_08072021.tsv"), sep = "\t", header = T, stringsAsFactors = F)
}


metadata_included <- metadata[!(metadata$is_blacklisted),]


if (dir.exists("/hpc/cuppen/")){
  pcawg_meta <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/metadata-per-donor.Rds")
} else {
  pcawg_meta <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/metadata-per-donor.Rds")
}



if (dir.exists("/hpc/cuppen/")){
  hmf_meta <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj2/HMF/metadata.tsv", header = T,
                       sep = "\t", stringsAsFactors = F)
} else {
  hmf_meta <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj2/HMF/metadata.tsv", header = T,
                       sep = "\t", stringsAsFactors = F)
}


# =================================================================================================================================
# # Processing raw clinvar database downloaded from  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/ 
# # It was initially processed by remvoving the header section (lines starting with ##)
# 
# path_to_vcf <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/clinvar_20210710_processed.vcf.gz"
# 
# clinvar <- read.csv(file = path_to_vcf, sep = "\t", header = T, stringsAsFactors = F)
# 
# selelcted_info_fields <- getInfoValues(clinvar$INFO, keys = c("CLNSIG"))
# 
# 
# clinvar_processed <- cbind(clinvar[,1:5], selelcted_info_fields)
# 
# # reordering columns to be compatible with downstream functions according to  https://github.com/UMCUGenetics/geneDriverAnnotator/blob/master/R/preProcessSnvIndel.R
# clinvar_processed <- clinvar_processed[,c(1,2,4,5,6,3)]
# 
# order_chr <- c(1:22, "X", "Y")
# clinvar_processed <- clinvar_processed[order(factor(clinvar_processed$X.CHROM, levels = order_chr, ordered = T), clinvar_processed$POS),]
# 
# 
# 
# clinvar_processed$X.CHROM <- paste0("chr", clinvar_processed$X.CHROM)
# head(clinvar_processed)
# # colnames should not be included
# write.table(clinvar_processed, file = "/home/ali313/Desktop/clinvar_processed.bed",
#             sep = "\t", quote = F, row.names = F, col.names = F)



# =================================================================================================================================
# Annotating sample vcf files with clinvar

genes_of_interest <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/genes.txt.gz"
# exons_of_interest <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/exons.txt.gz"

genes_of_interest_df <- read.csv(file = genes_of_interest, sep = "\t", header = T, stringsAsFactors = F)


# sample_id <- "CPCT02010003T"
# 
# path <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/",
#                hmf_meta$setName[hmf_meta$sampleId == sample_id], "/purple/", sample_id, ".purple.somatic.vcf.gz")


if (dir.exists("/hpc/cuppen/")){
  manifest <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/manifest_HMF_PCAWG.gene_ann.txt.gz", sep = "\t",
                       header = T, stringsAsFactors = F)
} else {
  manifest <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/manifest_HMF_PCAWG.gene_ann.txt.gz", sep = "\t",
                       header = T, stringsAsFactors = F)
}

if (dir.exists("/hpc/cuppen/")){
  clinvar_db <- "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/clinvar_processed.bed.bgz"
} else {
  clinvar_db <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/clinvar_processed.bed.bgz"
}


if (dir.exists("/hpc/cuppen/")){
  hotspot_db <- "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/known_hotspots_hg19.txt.gz"
} else {
  hotspot_db <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/known_hotspots_hg19.txt.gz"
}


if (dir.exists("/hpc/cuppen/")){
  snpeff_path <- "/home/ali313/Documents/studies/genomic-data-analyst/projects/pcawg-hmf-comparison/clones/geneDriverAnnotator/dep/snpEff/snpEff.jar" # to be fixed
} else {
  snpeff_path <- "/home/ali313/Documents/studies/genomic-data-analyst/projects/pcawg-hmf-comparison/clones/geneDriverAnnotator/dep/snpEff/snpEff.jar"
}

if (dir.exists("/hpc/cuppen/")){
  snpsift_path <- "/home/ali313/Documents/studies/genomic-data-analyst/projects/pcawg-hmf-comparison/clones/geneDriverAnnotator/dep/snpEff/SnpSift.jar" # to be fixed
} else {
  snpsift_path <- "/home/ali313/Documents/studies/genomic-data-analyst/projects/pcawg-hmf-comparison/clones/geneDriverAnnotator/dep/snpEff/SnpSift.jar"
}



for (i in 1:1){
  if (dir.exists("/hpc/cuppen/")){
    
    # germVcf<- paste0(manifest[i, "dir"], manifest[i, "germ_vcf"])
    somVcf <- paste0(manifest[i, "dir"], manifest[i, "som_vcf"])
    
    sample_name <- manifest[i, "sample"]
    
    dir.create(file.path(paste0("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/diplotype-output/", sample_name, "/simple-mutation/som/")), recursive = T, showWarnings = FALSE)
    out_dir <- paste0("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/diplotype-output/", sample_name, "/simple-mutation/som/")
    
  } else {
    
    # germVcf<- paste0(local, manifest[i, "dir"], manifest[i, "germ_vcf"])
    somVcf <- paste0(local, manifest[i, "dir"], manifest[i, "som_vcf"])

    sample_name <- manifest[i, "sample"]
    
    dir.create(file.path(paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/diplotype-output/", sample_name, "/simple-mutation/som/")), recursive = T, showWarnings = FALSE)
    out_dir <- paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/diplotype-output/", sample_name, "/simple-mutation/som/")
  }
  
  
  filterVcf(vcf.file = somVcf, out.file = paste0(out_dir, "filteredVcf-genes.txt.gz"), bed.file = genes_of_interest, mode = "som")
  
  path_to_filtered_vcf_genes <- paste0(out_dir, "filteredVcf-genes.txt.gz")
            
                                       
   vcf <- variantsFromVcf(vcf.file = path_to_filtered_vcf_genes, 
                          merge.consecutive = T,
                          vcf.fields = c("CHROM", "POS", "REF", "ALT"), 
                          verbose = T)
   
   vcf$clinvar_sig <- geneDriverAnnotator::getClinSig(df = vcf, 
                                                      db.path = clinvar_db)
   
   vcf$is_hotspot_mut <- geneDriverAnnotator::detIsHotspotMut(df = vcf,
                                                              db.path = hotspot_db)
   
  
   geneDriverAnnotator::annotateVariantType(vcf.file = path_to_filtered_vcf_genes, out.file = paste0(out_dir, "snpeff.vcf.gz"), 
                                            snpeff.path = snpeff_path)
   
   
   
   geneDriverAnnotator::extractVcfFields(vcf.file = paste0(out_dir, "snpeff.vcf.gz"), out.file = paste0(out_dir, "snpeff-extracted.vcf.gz"),
                                         snpsift.path = snpsift_path)
   
  
   snpeff_info <- read.csv(file = paste0(out_dir, "snpeff-extracted.vcf.gz"), sep = "\t", header = T,
                           stringsAsFactors = F)
   
   snpeff_clinvar_info <- cbind(snpeff_info, vcf[,c("clinvar_sig", "is_hotspot_mut")])
   
   # snpeff_clinvar_info <- snpeff_clinvar_info[snpeff_clinvar_info$ensembl_gene_id %in% genes_of_interest_df$ensembl_gene_id,]
   
   mut_profile_som <- mkMutProfileSnvIndel(df.snv.indel = snpeff_clinvar_info)
   
   mut_profile_som$hgnc_symbol <- mut_profile_som$snpeff_gene
   mut_profile_som <- cbind(mut_profile_som[,c(1:8)], mut_profile_som[,16], mut_profile_som[,c(9:15)])
   colnames(mut_profile_som)[9] <- "hgnc_symbol"
   
   write.table(mut_profile_som, file = paste0(out_dir, "mut_profile_som_simp.txt"), sep = "\t", quote = F, row.names = F)
   
   
}









