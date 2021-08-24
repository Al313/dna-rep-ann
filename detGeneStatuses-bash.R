

args <- commandArgs(trailingOnly=T)


# library(seqminer)
library(devtools)
# library(geneDriverAnnotator)

# `%notin%` <- Negate(`%in%`)


if (dir.exists("/hpc/cuppen/")){
  devtools::load_all('/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/geneDriverAnnotator')
} else {
  devtools::load_all('/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/geneDriverAnnotator')
}


# Made a copy of it to
local <- "/home/ali313/Documents/studies/master/umc-project"


if (dir.exists("/hpc/cuppen/")){
  wd <- "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/"
} else {
  wd <- paste0(local,"/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/")
}

if (dir.exists("/hpc/cuppen/")){
  manifest <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/manifest_HMF_PCAWG.gene_ann.txt.gz", sep = "\t",
                       header = T, stringsAsFactors = F)
} else {
  manifest <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/manifest_HMF_PCAWG.gene_ann.txt.gz", sep = "\t",
                       header = T, stringsAsFactors = F)
}




# 
# # To check if all the necessary files exist
# 
# for (i in 1:nrow(manifest)) {
#   print(i)
#   somVcf <- paste0(manifest[i, "dir"], manifest[i, "som_vcf"])
#   cnv <- paste0(manifest[i, "dir"], manifest[i, "cnv_tsv"])
#   
#   sample_name <- manifest[i, "sample"]
#   
#   if (!file.exists(somVcf) | !file.exists(cnv)){
#     print(paste(sample_name, "is not processed"))
#   }
# }

# the clone on my pc is located at: /home/ali313/Documents/studies/genomic-data-analyst/projects/pcawg-hmf-comparison/clones/geneDriverAnnotator/dep/snpEff/snpEff.jar

if (dir.exists("/hpc/cuppen/")){
  snpeff_path <- "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/geneDriverAnnotator/dep/snpEff/snpEff.jar"
} else {
  snpeff_path <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/geneDriverAnnotator/dep/snpEff/snpEff.jar"
}

if (dir.exists("/hpc/cuppen/")){
  snpsift_path <- "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/geneDriverAnnotator/dep/snpEff/SnpSift.jar"
} else {
  snpsift_path <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/geneDriverAnnotator/dep/snpEff/SnpSift.jar"
}

#
# if (dir.exists("/hpc/cuppen/")){
#   clinvar_db <- "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/clinvar_processed.bed.bgz"
# } else {
#   clinvar_db <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/clinvar_processed.bed.bgz"
# }

if (dir.exists("/hpc/cuppen/")){
  JAVA <- "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/geneDriverAnnotator/dep/jre1.8.0_191/bin/java"
} else {
  JAVA <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/geneDriverAnnotator/dep/jre1.8.0_191/bin/java"
}




if (dir.exists("/hpc/cuppen/")){
  
  germVcf<- args[1]
  somVcf <- args[2]
  cnv <- args[3]
  sample_name <- args[4]
  
  out_dir <- paste0(wd,"diplotype-output/", sample_name)
  dir.create(file.path(out_dir), showWarnings = FALSE)
  
  
  genes_bed_file <- '/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/genes.txt.gz'
  exons_bed_file <- '/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/exons.txt.gz'
} else {
  
  germVcf<- paste0(local, manifest[i, "dir"], manifest[i, "germ_vcf"])
  somVcf <- paste0(local, manifest[i, "dir"], manifest[i, "som_vcf"])
  cnv <- paste0(local, manifest[i, "dir"], manifest[i, "cnv_tsv"])
  sample_name <- manifest[i, "sample"]
  
  out_dir <- paste0(wd,"diplotype-output/", sample_name)
  dir.create(file.path(out_dir), showWarnings = FALSE)
  
  
  genes_bed_file <- '/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/genes.txt.gz'
  exons_bed_file <- '/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/exons.txt.gz'
}



detGeneStatuses(
  out.dir = out_dir,
  input.file.paths = c(germ_vcf=NA, som_vcf=somVcf, cnv=cnv),
  sample.name = sample_name,
  java.path = JAVA,
  genes.bed.file = genes_bed_file,
  exons.bed.file = exons_bed_file,
  sel.cols.cnv=c(chrom='chromosome',start='start',end='end',total_cn='copyNumber',major_cn='majorAlleleCopyNumber',minor_cn='minorAlleleCopyNumber'),
  do.snpeff.ann = T,
  snpeff.path = snpeff_path,
  snpsift.path = snpsift_path
)









