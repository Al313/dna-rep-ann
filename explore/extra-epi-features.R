
### Purpose

# In this script I investigate the dynamics of mutational processes using mutational signature and mutational timing information.

### Loading required libraries

# library(ggplot2)
# library(ggrepel)
# library(ggpubr)
# library(magrittr)
# library(rcompanion)
# library(epitools)
library(stringr)
# library(tidyr)
# library(plyr)
# library(miscTools)
# library(effsize)
# library(gggibbous)




### Setting constants and reading in data


`%notin%` <- Negate(`%in%`)

local <- "/home/ali313/Documents/studies/master/umc-project"


if (dir.exists("/hpc/cuppen/")){
  wd <- "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/"
} else {
  wd <- paste0(local,"/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/")
}


if (dir.exists("/hpc/cuppen/")){
  hmf_meta <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/metadata_whitelisted.tsv", sep = "\t", header = T, stringsAsFactors = F)
} else {
  hmf_meta <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/metadata_whitelisted.tsv", sep = "\t", header = T, stringsAsFactors = F)
}


if (dir.exists("/hpc/cuppen/")){
  metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2-metadata_25082021.tsv"), sep = "\t", header = T, stringsAsFactors = F)
} else {
  metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2-metadata_25082021.tsv"), sep = "\t", header = T, stringsAsFactors = F)
}

metadata_included <- metadata[!(metadata$is_blacklisted),]

metadata_included$tmb <- rowSums(metadata_included[,22:23])


if (dir.exists("/hpc/cuppen/")){
  cosmic_sigs <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/cosmic-sigs.txt", sep = "\t", header = F)
} else {
  cosmic_sigs <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/cosmic-sigs.txt", sep = "\t", header = F)
}



# j <- 1
# sample_id <- metadata_included$sample_id[j]
# ann_vcf_test <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/updated/", sample_id, ".txt.gz"), header = T, stringsAsFactors = F, sep = "\t")
# 



#####################################################################3


# all.cols <- c("Active", "Active-2", "Heterochromatin", "Inactive", "Repressed", "TAD_bd")
# 
# num_row <- length(all.cols)
# 
# col.name <- cosmic_sigs$V1
# 
# mt <- matrix(nrow = num_row*nrow(metadata_included), ncol = (3 + length(col.name)))
# rownames(mt) <- as.character(1:(num_row*nrow(metadata_included)))
# colnames(mt) <- c("sample_id", "chromatin_tad", "info_exists", col.name)
# 
# not_include <- F
# first <- T
# # nrow(metadata_included)
# for (i in 1:nrow(metadata_included)){
#   print(i)
#   print("chromatin_tad")
#   sample_id <- metadata_included$sample_id[i]
#   # print(sample_id)
# 
# 
#   if (dir.exists("/hpc/cuppen/")){
#     ann_vcf <- tryCatch(read.csv(file = paste0("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/updated/", sample_id, ".txt.gz"), header = T, stringsAsFactors = F, sep = "\t"), error=function(e) NULL)
#   } else {
#     ann_vcf <- tryCatch(read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/updated/", sample_id, ".txt.gz"), header = T, stringsAsFactors = F, sep = "\t"), error=function(e) NULL)
#   }
# 
# 
#   if (!is.null(ann_vcf)){
#     exposure_threshold <- names(table(ann_vcf$mut_sign_updated)[100*table(ann_vcf$mut_sign_updated)/nrow(ann_vcf) >= 2.5]) ## Maybe when re-running next time set the threshold to 2 percent
#     if (length(exposure_threshold) >= 2){
# 
#       ann_vcf <- ann_vcf[ann_vcf$mut_sign_updated %in% exposure_threshold,]
# 
#       mm <- as.matrix(table(ann_vcf$mut_sign_updated, ann_vcf$chromatin_status_tads))
# 
#       mm <- mm[,!str_detect(colnames(mm), pattern = "\\/")]
# 
# 
#       m2 <- matrix(nrow = nrow(mm), ncol = length(colnames(mm)), dimnames = list(NULL, colnames(mm)))
#       m2[ , colnames(mm)] <- mm
#       rownames(m2) <- rownames(mm)
# 
#       mm <- m2
# 
#       cos_similarity_check <- as.numeric(sapply(str_split(rownames(mm), pattern = "\\."), "[[", 4)) >= 85
# 
# 
#       mm <- as.data.frame(mm)
#       mm <- mm[cos_similarity_check,]
# 
#       if (nrow(mm) >= 2){
# 
#         duplicated_check <- duplicated(sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2))
# 
# 
#         duplicates <- sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2)[duplicated_check]
# 
#         if (length(duplicates) > 0) {
#           for (j in 1:length(duplicates)){
#             duplicate_indeces <- which(sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2) == duplicates[j])
#             mm <- rbind(mm, colSums(mm[duplicate_indeces,]))
#             rownames(mm)[nrow(mm)] <- paste0(".", duplicates[j])
#             mm <- mm[-duplicate_indeces,]
#           }
#         }
# 
# 
#         rownames(mm) <- sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2)
# 
# 
#         mm <- t(mm)
# 
# 
#         if (!first){
#           i <- i + ((num_row-1)*(i-1))
#         }
# 
#         first <- F
#         mt[i:(i+num_row-1),"sample_id"] <- sample_id
#         mt[i:(i+num_row-1),"chromatin_tad"] <- rownames(mm)
#         mt[i:(i+num_row-1),"info_exists"] <- T
# 
#         mt[i:(i+num_row-1),colnames(mm)] <- mm
#       } else {not_include <- T}
#     } else {not_include <- T}
#   } else {not_include <- T}
# 
# 
#   if (not_include) {
#     if (!first){
#       i <- i + ((num_row-1)*(i-1))
#     }
#     first <- F
#     mt[i:(i+num_row-1),"sample_id"] <- sample_id
#     mt[i:(i+num_row-1),"chromatin_tad"] <- all.cols
#     mt[i:(i+num_row-1),"info_exists"] <- F
#   }
# 
#   not_include <- F
# }
# 
# mt <- as.data.frame(mt)
# mt$info_exists <- as.logical(mt$info_exists)
# mt[,4:ncol(mt)] <- sapply(mt[,4:ncol(mt)], as.numeric)
# 
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   write.table(mt, file = gzfile("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-chromatin-tad-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# } else {
#   write.table(mt, file = gzfile("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-chromatin-tad-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# }



if (dir.exists("/hpc/cuppen/")){
  ms_chromatin <- read.csv("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-chromatin-tad-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
} else {
  ms_chromatin <- read.csv("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-chromatin-tad-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
}


#####################################################################3

# 
# all.cols <- c("early", "late", "mid")
# 
# 
# num_row <- length(all.cols)
# 
# col.name <- cosmic_sigs$V1
# 
# mt <- matrix(nrow = num_row*nrow(metadata_included), ncol = (3 + length(col.name)))
# rownames(mt) <- as.character(1:(num_row*nrow(metadata_included)))
# colnames(mt) <- c("sample_id", "rep_timing", "info_exists", col.name)
# 
# not_include <- F
# first <- T
# # nrow(metadata_included)
# for (i in 7040:nrow(metadata_included)){
#   print(i)
#   print("rep_timing")
#   sample_id <- metadata_included$sample_id[i]
#   # print(sample_id)
#   
#   
#   if (dir.exists("/hpc/cuppen/")){
#     ann_vcf <- tryCatch(read.csv(file = paste0("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/updated/", sample_id, ".txt.gz"), header = T, stringsAsFactors = F, sep = "\t"), error=function(e) NULL)
#   } else {
#     ann_vcf <- tryCatch(read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/updated/", sample_id, ".txt.gz"), header = T, stringsAsFactors = F, sep = "\t"), error=function(e) NULL)
#   }
#   
#   
#   if (!is.null(ann_vcf)){
#     exposure_threshold <- names(table(ann_vcf$mut_sign_updated)[100*table(ann_vcf$mut_sign_updated)/nrow(ann_vcf) >= 2.5]) ## Maybe when re-running next time set the threshold to 2 percent
#     if (length(exposure_threshold) >= 2){
#       
#       ann_vcf <- ann_vcf[ann_vcf$mut_sign_updated %in% exposure_threshold,]
#       
#       mm <- as.matrix(table(ann_vcf$mut_sign_updated, ann_vcf$rep_timing))
#       
#       mm <- mm[,!str_detect(colnames(mm), pattern = "\\/")]
#       
#       
#       m2 <- matrix(nrow = nrow(mm), ncol = length(colnames(mm)), dimnames = list(NULL, colnames(mm)))
#       m2[ , colnames(mm)] <- mm
#       rownames(m2) <- rownames(mm)
#       
#       mm <- m2
#       
#       cos_similarity_check <- as.numeric(sapply(str_split(rownames(mm), pattern = "\\."), "[[", 4)) >= 85
#       
#       
#       mm <- as.data.frame(mm)
#       mm <- mm[cos_similarity_check,]
#       
#       if (nrow(mm) >= 2){
#         
#         duplicated_check <- duplicated(sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2))
#         
#         
#         duplicates <- sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2)[duplicated_check]
#         
#         if (length(duplicates) > 0) {
#           for (j in 1:length(duplicates)){
#             duplicate_indeces <- which(sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2) == duplicates[j])
#             mm <- rbind(mm, colSums(mm[duplicate_indeces,]))
#             rownames(mm)[nrow(mm)] <- paste0(".", duplicates[j])
#             mm <- mm[-duplicate_indeces,]
#           }
#         }
#         
#         
#         rownames(mm) <- sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2)
#         
#         
#         mm <- t(mm)
#         
#         
#         if (!first){
#           i <- i + ((num_row-1)*(i-1))
#         }
#         
#         first <- F
#         mt[i:(i+num_row-1),"sample_id"] <- sample_id
#         mt[i:(i+num_row-1),"rep_timing"] <- rownames(mm)
#         mt[i:(i+num_row-1),"info_exists"] <- T
#         
#         mt[i:(i+num_row-1),colnames(mm)] <- mm
#       } else {not_include <- T}
#     } else {not_include <- T}
#   } else {not_include <- T}
#   
#   
#   if (not_include) {
#     if (!first){
#       i <- i + ((num_row-1)*(i-1))
#     }
#     first <- F
#     mt[i:(i+num_row-1),"sample_id"] <- sample_id
#     mt[i:(i+num_row-1),"rep_timing"] <- all.cols
#     mt[i:(i+num_row-1),"info_exists"] <- F
#   }
#   
#   not_include <- F
# }
# 
# mt <- as.data.frame(mt)
# mt$info_exists <- as.logical(mt$info_exists)
# mt[,4:ncol(mt)] <- sapply(mt[,4:ncol(mt)], as.numeric)
# 
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   write.table(mt, file = gzfile("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-rep-timing-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# } else {
#   write.table(mt, file = gzfile("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-rep-timing-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# }
# 

if (dir.exists("/hpc/cuppen/")){
  ms_rep_timing <- read.csv("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-rep-timing-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
} else {
  ms_rep_timing <- read.csv("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-rep-timing-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
}


#####################################################################3



# all.cols <- c("-", "untranscribed", "transcribed")
# 
# num_row <- length(all.cols)
# 
# col.name <- cosmic_sigs$V1
# 
# mt <- matrix(nrow = num_row*nrow(metadata_included), ncol = (3 + length(col.name)))
# rownames(mt) <- as.character(1:(num_row*nrow(metadata_included)))
# colnames(mt) <- c("sample_id", "trp_str_ann", "info_exists", col.name)
# 
# not_include <- F
# first <- T
# # nrow(metadata_included)
# for (i in 1:nrow(metadata_included)){
#   print(i)
#   print("trp_str_ann")
#   sample_id <- metadata_included$sample_id[i]
#   # print(sample_id)
# 
# 
#   if (dir.exists("/hpc/cuppen/")){
#     ann_vcf <- tryCatch(read.csv(file = paste0("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/updated/", sample_id, ".txt.gz"), header = T, stringsAsFactors = F, sep = "\t"), error=function(e) NULL)
#   } else {
#     ann_vcf <- tryCatch(read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/updated/", sample_id, ".txt.gz"), header = T, stringsAsFactors = F, sep = "\t"), error=function(e) NULL)
#   }
# 
# 
#   if (!is.null(ann_vcf)){
#     exposure_threshold <- names(table(ann_vcf$mut_sign_updated)[100*table(ann_vcf$mut_sign_updated)/nrow(ann_vcf) >= 2.5]) ## Maybe when re-running next time set the threshold to 2 percent
#     if (length(exposure_threshold) >= 2){
# 
#       ann_vcf <- ann_vcf[ann_vcf$mut_sign_updated %in% exposure_threshold,]
# 
#       mm <- as.matrix(table(ann_vcf$mut_sign_updated, ann_vcf$trp_str_ann))
# 
#       mm <- mm[,!str_detect(colnames(mm), pattern = "\\/")]
# 
# 
#       m2 <- matrix(nrow = nrow(mm), ncol = length(colnames(mm)), dimnames = list(NULL, colnames(mm)))
#       m2[ , colnames(mm)] <- mm
#       rownames(m2) <- rownames(mm)
# 
#       mm <- m2
# 
#       cos_similarity_check <- as.numeric(sapply(str_split(rownames(mm), pattern = "\\."), "[[", 4)) >= 85
# 
# 
#       mm <- as.data.frame(mm)
#       mm <- mm[cos_similarity_check,]
# 
#       if (nrow(mm) >= 2){
# 
#         duplicated_check <- duplicated(sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2))
# 
# 
#         duplicates <- sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2)[duplicated_check]
# 
#         if (length(duplicates) > 0) {
#           for (j in 1:length(duplicates)){
#             duplicate_indeces <- which(sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2) == duplicates[j])
#             mm <- rbind(mm, colSums(mm[duplicate_indeces,]))
#             rownames(mm)[nrow(mm)] <- paste0(".", duplicates[j])
#             mm <- mm[-duplicate_indeces,]
#           }
#         }
# 
# 
#         rownames(mm) <- sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2)
# 
# 
#         mm <- t(mm)
# 
# 
#         if (!first){
#           i <- i + ((num_row-1)*(i-1))
#         }
# 
#         first <- F
#         mt[i:(i+num_row-1),"sample_id"] <- sample_id
#         mt[i:(i+num_row-1),"trp_str_ann"] <- rownames(mm)
#         mt[i:(i+num_row-1),"info_exists"] <- T
# 
#         mt[i:(i+num_row-1),colnames(mm)] <- mm
#       } else {not_include <- T}
#     } else {not_include <- T}
#   } else {not_include <- T}
# 
# 
#   if (not_include) {
#     if (!first){
#       i <- i + ((num_row-1)*(i-1))
#     }
#     first <- F
#     mt[i:(i+num_row-1),"sample_id"] <- sample_id
#     mt[i:(i+num_row-1),"trp_str_ann"] <- all.cols
#     mt[i:(i+num_row-1),"info_exists"] <- F
#   }
# 
#   not_include <- F
# }
# 
# mt <- as.data.frame(mt)
# mt$info_exists <- as.logical(mt$info_exists)
# mt[,4:ncol(mt)] <- sapply(mt[,4:ncol(mt)], as.numeric)
# 
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   write.table(mt, file = gzfile("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-trp-str-ann-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# } else {
#   write.table(mt, file = gzfile("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-trp-str-ann-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# }



if (dir.exists("/hpc/cuppen/")){
  ms_trp_strand <- read.csv("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-trp-str-ann-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
} else {
  ms_trp_strand <- read.csv("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-trp-str-ann-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
}



#####################################################################3






# 
# all.cols <- c("-", "left", "right")
# 
# num_row <- length(all.cols)
# 
# col.name <- cosmic_sigs$V1
# 
# mt <- matrix(nrow = num_row*nrow(metadata_included), ncol = (3 + length(col.name)))
# rownames(mt) <- as.character(1:(num_row*nrow(metadata_included)))
# colnames(mt) <- c("sample_id", "rep_str_ann", "info_exists", col.name)
# 
# not_include <- F
# first <- T
# # nrow(metadata_included)
# for (i in 1:nrow(metadata_included)){
#   print(i)
#   print("rep_str_ann")
#   sample_id <- metadata_included$sample_id[i]
#   # print(sample_id)
# 
# 
#   if (dir.exists("/hpc/cuppen/")){
#     ann_vcf <- tryCatch(read.csv(file = paste0("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/updated/", sample_id, ".txt.gz"), header = T, stringsAsFactors = F, sep = "\t"), error=function(e) NULL)
#   } else {
#     ann_vcf <- tryCatch(read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/updated/", sample_id, ".txt.gz"), header = T, stringsAsFactors = F, sep = "\t"), error=function(e) NULL)
#   }
# 
# 
#   if (!is.null(ann_vcf)){
#     exposure_threshold <- names(table(ann_vcf$mut_sign_updated)[100*table(ann_vcf$mut_sign_updated)/nrow(ann_vcf) >= 2.5]) ## Maybe when re-running next time set the threshold to 2 percent
#     if (length(exposure_threshold) >= 2){
# 
#       ann_vcf <- ann_vcf[ann_vcf$mut_sign_updated %in% exposure_threshold,]
# 
#       mm <- as.matrix(table(ann_vcf$mut_sign_updated, ann_vcf$rep_str_ann))
# 
# 
#       m2 <- matrix(nrow = nrow(mm), ncol = length(colnames(mm)), dimnames = list(NULL, colnames(mm)))
#       m2[ , colnames(mm)] <- mm
#       rownames(m2) <- rownames(mm)
# 
#       mm <- m2
# 
#       cos_similarity_check <- as.numeric(sapply(str_split(rownames(mm), pattern = "\\."), "[[", 4)) >= 85
# 
# 
#       mm <- as.data.frame(mm)
#       mm <- mm[cos_similarity_check,]
# 
#       if (nrow(mm) >= 2){
# 
#         duplicated_check <- duplicated(sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2))
# 
# 
#         duplicates <- sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2)[duplicated_check]
# 
#         if (length(duplicates) > 0) {
#           for (j in 1:length(duplicates)){
#             duplicate_indeces <- which(sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2) == duplicates[j])
#             mm <- rbind(mm, colSums(mm[duplicate_indeces,]))
#             rownames(mm)[nrow(mm)] <- paste0(".", duplicates[j])
#             mm <- mm[-duplicate_indeces,]
#           }
#         }
# 
# 
#         rownames(mm) <- sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2)
# 
# 
#         mm <- t(mm)
# 
# 
#         if (!first){
#           i <- i + ((num_row-1)*(i-1))
#         }
# 
#         first <- F
#         mt[i:(i+num_row-1),"sample_id"] <- sample_id
#         mt[i:(i+num_row-1),"rep_str_ann"] <- rownames(mm)
#         mt[i:(i+num_row-1),"info_exists"] <- T
# 
#         mt[i:(i+num_row-1),colnames(mm)] <- mm
#       } else {not_include <- T}
#     } else {not_include <- T}
#   } else {not_include <- T}
# 
# 
#   if (not_include) {
#     if (!first){
#       i <- i + ((num_row-1)*(i-1))
#     }
#     first <- F
#     mt[i:(i+num_row-1),"sample_id"] <- sample_id
#     mt[i:(i+num_row-1),"rep_str_ann"] <- all.cols
#     mt[i:(i+num_row-1),"info_exists"] <- F
#   }
# 
#   not_include <- F
# }
# 
# mt <- as.data.frame(mt)
# mt$info_exists <- as.logical(mt$info_exists)
# mt[,4:ncol(mt)] <- sapply(mt[,4:ncol(mt)], as.numeric)
# 
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   write.table(mt, file = gzfile("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-rep-str-ann-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# } else {
#   write.table(mt, file = gzfile("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-rep-str-ann-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# }


if (dir.exists("/hpc/cuppen/")){
  ms_rep_strand <- read.csv("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-rep-str-ann-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
} else {
  ms_rep_strand <- read.csv("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-rep-str-ann-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
}







#####################################################################3





# all.cols <- factor(c("coding", "fiveUTR", "intergenic", "intron" , "promoter", "threeUTR", "NA"), levels = c("coding", "fiveUTR", "intergenic", "intron", "promoter", "threeUTR", "NA"))
# 
# num_row <- length(all.cols)
# 
# col.name <- cosmic_sigs$V1
# 
# mt <- matrix(nrow = num_row*nrow(metadata_included), ncol = (3 + length(col.name)))
# rownames(mt) <- as.character(1:(num_row*nrow(metadata_included)))
# colnames(mt) <- c("sample_id", "genomic_func", "info_exists", col.name)
# 
# # i <- 1
# not_include <- F
# first <- T
# # nrow(metadata_included)
# for (i in 1:nrow(metadata_included)){
#   print(i)
#   print("genomic_func")
#   sample_id <- metadata_included$sample_id[i]
#   # print(sample_id)
# 
# 
#   if (dir.exists("/hpc/cuppen/")){
#     ann_vcf <- tryCatch(read.csv(file = paste0("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/updated/", sample_id, ".txt.gz"), header = T, stringsAsFactors = F, sep = "\t"), error=function(e) NULL)
#   } else {
#     ann_vcf <- tryCatch(read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/updated/", sample_id, ".txt.gz"), header = T, stringsAsFactors = F, sep = "\t"), error=function(e) NULL)
#   }
# 
# 
#   if (!is.null(ann_vcf)){
#     exposure_threshold <- names(table(ann_vcf$mut_sign_updated)[100*table(ann_vcf$mut_sign_updated)/nrow(ann_vcf) >= 2.5]) ## Maybe when re-running next time set the threshold to 2 percent
#     if (length(exposure_threshold) >= 2){
# 
#       ann_vcf <- ann_vcf[ann_vcf$mut_sign_updated %in% exposure_threshold,]
#       ann_vcf$genomic_func <- factor(ann_vcf$genomic_func, levels = levels(all.cols)[-length(all.cols)])
# 
#       mm <- as.matrix(table(ann_vcf$mut_sign_updated, ann_vcf$genomic_func, useNA= "always"))
#       mm <- as.matrix(mm[-nrow(mm),])
#       colnames(mm)[ncol(mm)] <- "NA"
# 
#       
#       
#       m2 <- matrix(nrow = nrow(mm), ncol = length(colnames(mm)), dimnames = list(NULL, colnames(mm)))
#       m2[ , colnames(mm)] <- mm
#       rownames(m2) <- rownames(mm)
# 
#       mm <- m2
# 
#       cos_similarity_check <- as.numeric(sapply(str_split(rownames(mm), pattern = "\\."), "[[", 4)) >= 85
# 
# 
#       mm <- as.data.frame(mm)
#       mm <- mm[cos_similarity_check,]
# 
#       if (nrow(mm) >= 2){
# 
#         duplicated_check <- duplicated(sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2))
# 
# 
#         duplicates <- sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2)[duplicated_check]
# 
#         if (length(duplicates) > 0) {
#           for (j in 1:length(duplicates)){
#             duplicate_indeces <- which(sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2) == duplicates[j])
#             mm <- rbind(mm, colSums(mm[duplicate_indeces,]))
#             rownames(mm)[nrow(mm)] <- paste0(".", duplicates[j])
#             mm <- mm[-duplicate_indeces,]
#           }
#         }
# 
# 
#         rownames(mm) <- sapply(str_split(rownames(mm), pattern = "\\."), "[[", 2)
# 
# 
#         mm <- t(mm)
# 
# 
#         if (!first){
#           i <- i + ((num_row-1)*(i-1))
#         }
# 
#         first <- F
#         mt[i:(i+num_row-1),"sample_id"] <- sample_id
#         mt[i:(i+num_row-1),"genomic_func"] <- rownames(mm)
#         mt[i:(i+num_row-1),"info_exists"] <- T
# 
#         mt[i:(i+num_row-1),colnames(mm)] <- mm
#       } else {not_include <- T}
#     } else {not_include <- T}
#   } else {not_include <- T}
# 
# 
#   if (not_include) {
#     if (!first){
#       i <- i + ((num_row-1)*(i-1))
#     }
#     first <- F
#     mt[i:(i+num_row-1),"sample_id"] <- sample_id
#     mt[i:(i+num_row-1),"genomic_func"] <- all.cols
#     mt[i:(i+num_row-1),"info_exists"] <- F
#   }
# 
#   not_include <- F
# }
# 
# mt <- as.data.frame(mt)
# mt$info_exists <- as.logical(mt$info_exists)
# mt[,4:ncol(mt)] <- sapply(mt[,4:ncol(mt)], as.numeric)
# 
# 
# ms_genomic_func$genomic_func[is.na(ms_genomic_func$genomic_func)] <- "NA"
# 
# if (dir.exists("/hpc/cuppen/")){
#   write.table(mt, file = gzfile("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-genomic-func-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# } else {
#   write.table(mt, file = gzfile("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-genomic-func-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# }



if (dir.exists("/hpc/cuppen/")){
  ms_genomic_func <- read.csv("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-genomic-func-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
} else {
  ms_genomic_func <- read.csv("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-genomic-func-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
}



#####################################################################
#####################################################################
#####################################################################







