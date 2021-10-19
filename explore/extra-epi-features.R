
### Purpose

# In this script I investigate the dynamics of mutational processes using mutational signature and mutational timing information.

### Loading required libraries

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(magrittr)
library(rcompanion)
library(epitools)
library(tidyr)
library(plyr)
library(miscTools)
library(effsize)
library(gggibbous)
library(stringr)




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

args <- commandArgs(trailingOnly=T)

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
#   write.table(mt, file = gzfile("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-chromatin-tad-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# } else {
#   write.table(mt, file = gzfile("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-chromatin-tad-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# }


# Preparing the data for visualization (fold change)

# if (args[1] == 0){
#   
#   if (dir.exists("/hpc/cuppen/")){
#     ms_chromatin <- read.csv("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-chromatin-tad-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   } else {
#     ms_chromatin <- read.csv("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-chromatin-tad-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   }
#   
#   
#   ms_chromatin_df <- data.frame(sample_id = character(7049*51), cosmic_sig = character(7049*51), fold_change_active_to_inactive = numeric(7049*51))
#   
#   # head(ms_chromatin_df[!is.na(ms_chromatin_df$fold_change_active_to_inactive),])
#   # ms_chromatin[ms_chromatin$sample_id == "CPCT02030502T",]
#   
#   
#   
#   # length(metadata_included$sample_id)
#   for (i in 1:length(metadata_included$sample_id)){
#     print("chromatin_tad")
#     print(i)
#     sample_id <- metadata_included$sample_id[i]
#     i <- i + (50*(i-1))
#     
#     for (j in 1:length(cosmic_sigs$V1)){
#       
#       ms_chromatin_df[i+j-1, "sample_id"] <- sample_id
#       ms_chromatin_df[i+j-1, "cosmic_sig"] <- cosmic_sigs$V1[j]
#       Active <-  ms_chromatin[ms_chromatin$sample_id == sample_id & ms_chromatin$chromatin_tad == "Active",cosmic_sigs$V1[j]]
#       Inactive <- ms_chromatin[ms_chromatin$sample_id == sample_id & ms_chromatin$chromatin_tad == "Inactive",cosmic_sigs$V1[j]]
#       
#       if (!is.na(Active) & !is.na(Inactive) & Active != 0 & Inactive != 0){
#         ms_chromatin_df[i+j-1, "fold_change_active_to_inactive"] <- Active/Inactive
#       } else {
#         ms_chromatin_df[i+j-1, "fold_change_active_to_inactive"] <- NA
#       }
#     }
#   }
#   
#   
#   if (dir.exists("/hpc/cuppen/")){
#     saveRDS(ms_chromatin_df, file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_chromatin_df.all.rds")
#   } else {
#     saveRDS(ms_chromatin_df, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_chromatin_df.all.rds"))
#   }
# 
# }



if (dir.exists("/hpc/cuppen/")){
  ms_chromatin_df <- readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_chromatin_df.all.rds")
} else {
  ms_chromatin_df <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_chromatin_df.all.rds"))
}




# # Preparing the data for visualization (for effect size)
# 
# if (args[1] == 0){
# 
#   if (dir.exists("/hpc/cuppen/")){
#     ms_chromatin <- read.csv("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-chromatin-tad-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   } else {
#     ms_chromatin <- read.csv("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-chromatin-tad-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   }
# 
# 
#   ms_chromatin_df2 <- data.frame(sample_id = character(7049*51), cosmic_sig = character(7049*51), active_count = numeric(7049*51), inactive_count = numeric(7049*51))
# 
#   # head(ms_chromatin_df2[!is.na(ms_chromatin_df2$fold_change_active_to_inactive),])
#   # ms_chromatin[ms_chromatin$sample_id == "CPCT02030502T",]
# 
# 
# 
#   # length(metadata_included$sample_id)
#   for (i in 1:length(metadata_included$sample_id)){
#     print("chromatin_tad")
#     print(i)
#     sample_id <- metadata_included$sample_id[i]
#     i <- i + (50*(i-1))
# 
#     for (j in 1:length(cosmic_sigs$V1)){
# 
#       ms_chromatin_df2[i+j-1, "sample_id"] <- sample_id
#       ms_chromatin_df2[i+j-1, "cosmic_sig"] <- cosmic_sigs$V1[j]
#       Active <-  ms_chromatin[ms_chromatin$sample_id == sample_id & ms_chromatin$chromatin_tad == "Active",cosmic_sigs$V1[j]]
#       Inactive <- ms_chromatin[ms_chromatin$sample_id == sample_id & ms_chromatin$chromatin_tad == "Inactive",cosmic_sigs$V1[j]]
# 
#       if (!is.na(Active) & !is.na(Inactive) & Active != 0 & Inactive != 0){
#         ms_chromatin_df2[i+j-1, c("active_count", "inactive_count")] <- c(Active, Inactive)
#       } else {
#         ms_chromatin_df2[i+j-1, c("active_count", "inactive_count")] <- NA
#       }
#     }
#   }
# 
# 
#   if (dir.exists("/hpc/cuppen/")){
#     saveRDS(ms_chromatin_df2, file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_chromatin_df2.all.rds")
#   } else {
#     saveRDS(ms_chromatin_df2, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_chromatin_df2.all.rds"))
#   }
# 
# }

if (dir.exists("/hpc/cuppen/")){
  ms_chromatin_df2 <- readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_chromatin_df2.all.rds")
} else {
  ms_chromatin_df2 <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_chromatin_df2.all.rds"))
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
#   write.table(mt, file = gzfile("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-rep-timing-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# } else {
#   write.table(mt, file = gzfile("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-rep-timing-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# }
# 



# Preparing the data for visualization (fold change)

# if (args[1] == 1){
#   
#   if (dir.exists("/hpc/cuppen/")){
#     ms_rep_timing <- read.csv("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-rep-timing-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   } else {
#     ms_rep_timing <- read.csv("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-rep-timing-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   }
#   
#   
#   
#   ms_rep_timing_df <- data.frame(sample_id = character(7049*51), cosmic_sig = character(7049*51), fold_change_late_to_early = numeric(7049*51))
#   
#   # head(ms_rep_timing_df[!is.na(ms_rep_timing_df$fold_change_late_to_early),])
#   # ms_rep_timing[ms_rep_timing$sample_id == "CPCT02030502T",]
#   
#   
#   
#   # length(metadata_included$sample_id)
#   for (i in 1:length(metadata_included$sample_id)){
#     print("rep_timing")
#     print(i)
#     sample_id <- metadata_included$sample_id[i]
#     i <- i + (50*(i-1))
#     
#     for (j in 1:length(cosmic_sigs$V1)){
#       
#       ms_rep_timing_df[i+j-1, "sample_id"] <- sample_id
#       ms_rep_timing_df[i+j-1, "cosmic_sig"] <- cosmic_sigs$V1[j]
#       late <-  ms_rep_timing[ms_rep_timing$sample_id == sample_id & ms_rep_timing$rep_timing == "late",cosmic_sigs$V1[j]]
#       early <- ms_rep_timing[ms_rep_timing$sample_id == sample_id & ms_rep_timing$rep_timing == "early",cosmic_sigs$V1[j]]
#       
#       if (!is.na(late) & !is.na(early) & late != 0 & early != 0){
#         ms_rep_timing_df[i+j-1, "fold_change_late_to_early"] <- late/early
#       } else {
#         ms_rep_timing_df[i+j-1, "fold_change_late_to_early"] <- NA
#       }
#     }
#   }
#   
#   
#   if (dir.exists("/hpc/cuppen/")){
#     saveRDS(ms_rep_timing_df, file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_rep_timing_df.all.rds")
#   } else {
#     saveRDS(ms_rep_timing_df, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_rep_timing_df.all.rds"))
#   }
# }


if (dir.exists("/hpc/cuppen/")){
  ms_rep_timing_df <- readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_rep_timing_df.all.rds")
} else {
  ms_rep_timing_df <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_rep_timing_df.all.rds"))
}



# # Preparing the data for visualization (for effect size)
# if (args[1] == 1){
# 
#   if (dir.exists("/hpc/cuppen/")){
#     ms_rep_timing <- read.csv("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-rep-timing-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   } else {
#     ms_rep_timing <- read.csv("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-rep-timing-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   }
# 
# 
# 
#   ms_rep_timing_df2 <- data.frame(sample_id = character(7049*51), cosmic_sig = character(7049*51), late_count = numeric(7049*51), early_count = numeric(7049*51))
# 
#   # head(ms_rep_timing_df2[!is.na(ms_rep_timing_df2$fold_change_late_to_early),])
#   # ms_rep_timing[ms_rep_timing$sample_id == "CPCT02030502T",]
# 
# 
# 
#   # length(metadata_included$sample_id)
#   for (i in 1:length(metadata_included$sample_id)){
#     print("rep_timing")
#     print(i)
#     sample_id <- metadata_included$sample_id[i]
#     i <- i + (50*(i-1))
# 
#     for (j in 1:length(cosmic_sigs$V1)){
# 
#       ms_rep_timing_df2[i+j-1, "sample_id"] <- sample_id
#       ms_rep_timing_df2[i+j-1, "cosmic_sig"] <- cosmic_sigs$V1[j]
#       late <-  ms_rep_timing[ms_rep_timing$sample_id == sample_id & ms_rep_timing$rep_timing == "late",cosmic_sigs$V1[j]]
#       early <- ms_rep_timing[ms_rep_timing$sample_id == sample_id & ms_rep_timing$rep_timing == "early",cosmic_sigs$V1[j]]
# 
#       if (!is.na(late) & !is.na(early) & late != 0 & early != 0){
#         ms_rep_timing_df2[i+j-1, c('late_count', "early_count")] <- c(late, early)
#       } else {
#         ms_rep_timing_df2[i+j-1, c('late_count', "early_count")] <- NA
#       }
#     }
#   }
# 
# 
#   if (dir.exists("/hpc/cuppen/")){
#     saveRDS(ms_rep_timing_df2, file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_rep_timing_df2.all.rds")
#   } else {
#     saveRDS(ms_rep_timing_df2, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_rep_timing_df2.all.rds"))
#   }
# }


if (dir.exists("/hpc/cuppen/")){
  ms_rep_timing_df2 <- readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_rep_timing_df2.all.rds")
} else {
  ms_rep_timing_df2 <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_rep_timing_df2.all.rds"))
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
#   write.table(mt, file = gzfile("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-trp-str-ann-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# } else {
#   write.table(mt, file = gzfile("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-trp-str-ann-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# }



# Preparing the data for visualization (fold change)

# if (args[1] == 2){
#   
#   if (dir.exists("/hpc/cuppen/")){
#     ms_trp_strand <- read.csv("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-trp-str-ann-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   } else {
#     ms_trp_strand <- read.csv("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-trp-str-ann-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   }
#   
#   
#   
#   
#   ms_trp_strand_df <- data.frame(sample_id = character(7049*51), cosmic_sig = character(7049*51), fold_change_trans_to_untrans = numeric(7049*51))
#   
#   # head(ms_trp_strand_df[!is.na(ms_trp_strand_df$fold_change_trans_to_untrans),])
#   # ms_trp_strand[ms_trp_strand$sample_id == "CPCT02030502T",]
#   
#   
#   
#   # length(metadata_included$sample_id)
#   for (i in 1:length(metadata_included$sample_id)){
#     print("trp_str_ann")
#     print(i)
#     sample_id <- metadata_included$sample_id[i]
#     i <- i + (50*(i-1))
#   
#     for (j in 1:length(cosmic_sigs$V1)){
#   
#       ms_trp_strand_df[i+j-1, "sample_id"] <- sample_id
#       ms_trp_strand_df[i+j-1, "cosmic_sig"] <- cosmic_sigs$V1[j]
#       transcribed <-  ms_trp_strand[ms_trp_strand$sample_id == sample_id & ms_trp_strand$trp_str_ann == "transcribed",cosmic_sigs$V1[j]]
#       untranscribed <- ms_trp_strand[ms_trp_strand$sample_id == sample_id & ms_trp_strand$trp_str_ann == "untranscribed",cosmic_sigs$V1[j]]
#   
#       if (!is.na(transcribed) & !is.na(untranscribed) & transcribed != 0 & untranscribed != 0){
#         ms_trp_strand_df[i+j-1, "fold_change_trans_to_untrans"] <- transcribed/untranscribed
#       } else {
#         ms_trp_strand_df[i+j-1, "fold_change_trans_to_untrans"] <- NA
#       }
#     }
#   }
#   
#   if (dir.exists("/hpc/cuppen/")){
#     saveRDS(ms_trp_strand_df, file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_trp_strand_df.all.rds")
#   } else {
#     saveRDS(ms_trp_strand_df, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_trp_strand_df.all.rds"))
#   }
#   
# }


if (dir.exists("/hpc/cuppen/")){
  ms_trp_strand_df <- readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_trp_strand_df.all.rds")
} else {
  ms_trp_strand_df <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_trp_strand_df.all.rds"))
}





# Preparing the data for visualization (for effect size)

# if (args[1] == 2){
# 
#   if (dir.exists("/hpc/cuppen/")){
#     ms_trp_strand <- read.csv("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-trp-str-ann-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   } else {
#     ms_trp_strand <- read.csv("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-trp-str-ann-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   }
# 
# 
# 
# 
#   ms_trp_strand_df2 <- data.frame(sample_id = character(7049*51), cosmic_sig = character(7049*51), trans_count = numeric(7049*51), untrans_count = numeric(7049*51))
# 
#   # head(ms_trp_strand_df2[!is.na(ms_trp_strand_df2$fold_change_trans_to_untrans),])
#   # ms_trp_strand[ms_trp_strand$sample_id == "CPCT02030502T",]
# 
# 
# 
#   # length(metadata_included$sample_id)
#   for (i in 1:length(metadata_included$sample_id)){
#     print("trp_str_ann")
#     print(i)
#     sample_id <- metadata_included$sample_id[i]
#     i <- i + (50*(i-1))
# 
#     for (j in 1:length(cosmic_sigs$V1)){
# 
#       ms_trp_strand_df2[i+j-1, "sample_id"] <- sample_id
#       ms_trp_strand_df2[i+j-1, "cosmic_sig"] <- cosmic_sigs$V1[j]
#       transcribed <-  ms_trp_strand[ms_trp_strand$sample_id == sample_id & ms_trp_strand$trp_str_ann == "transcribed",cosmic_sigs$V1[j]]
#       untranscribed <- ms_trp_strand[ms_trp_strand$sample_id == sample_id & ms_trp_strand$trp_str_ann == "untranscribed",cosmic_sigs$V1[j]]
# 
#       if (!is.na(transcribed) & !is.na(untranscribed) & transcribed != 0 & untranscribed != 0){
#         ms_trp_strand_df2[i+j-1, c("trans_count", "untrans_count")] <- c(transcribed, untranscribed)
#       } else {
#         ms_trp_strand_df2[i+j-1, c("trans_count", "untrans_count")] <- NA
#       }
#     }
#   }
# 
#   if (dir.exists("/hpc/cuppen/")){
#     saveRDS(ms_trp_strand_df2, file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_trp_strand_df2.all.rds")
#   } else {
#     saveRDS(ms_trp_strand_df2, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_trp_strand_df2.all.rds"))
#   }
# 
# }


if (dir.exists("/hpc/cuppen/")){
  ms_trp_strand_df2 <- readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_trp_strand_df2.all.rds")
} else {
  ms_trp_strand_df2 <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_trp_strand_df2.all.rds"))
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
#   write.table(mt, file = gzfile("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-rep-str-ann-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# } else {
#   write.table(mt, file = gzfile("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-rep-str-ann-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# }


# Preparing the data for visualization (fold change)

# if (args[1] == 3){
#   
#   if (dir.exists("/hpc/cuppen/")){
#     ms_rep_strand <- read.csv("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-rep-str-ann-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   } else {
#     ms_rep_strand <- read.csv("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-rep-str-ann-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   }
#   
#   
#   
#   ms_rep_strand_df <- data.frame(sample_id = character(7049*51), cosmic_sig = character(7049*51), fold_change_right_to_left = numeric(7049*51))
#   
#   
#   # head(ms_rep_strand_df[!is.na(ms_rep_strand_df$fold_change_right_to_left),])
#   # ms_rep_strand[ms_rep_strand$sample_id == "CPCT02030502T",]
#   
#   
#   # length(metadata_included$sample_id)
#   for (i in 1:length(metadata_included$sample_id)){
#     print("rep_str_ann")
#     print(i)
#     sample_id <- metadata_included$sample_id[i]
#     i <- i + (50*(i-1))
#   
#     for (j in 1:length(cosmic_sigs$V1)){
#   
#       ms_rep_strand_df[i+j-1, "sample_id"] <- sample_id
#       ms_rep_strand_df[i+j-1, "cosmic_sig"] <- cosmic_sigs$V1[j]
#       right <-  ms_rep_strand[ms_rep_strand$sample_id == sample_id & ms_rep_strand$rep_str_ann == "right",cosmic_sigs$V1[j]]
#       left <- ms_rep_strand[ms_rep_strand$sample_id == sample_id & ms_rep_strand$rep_str_ann == "left",cosmic_sigs$V1[j]]
#   
#       if (!is.na(right) & !is.na(left) & right != 0 & left != 0){
#         ms_rep_strand_df[i+j-1, "fold_change_right_to_left"] <- right/left
#       } else {
#         ms_rep_strand_df[i+j-1, "fold_change_right_to_left"] <- NA
#       }
#     }
#   }
#   
#   
#   if (dir.exists("/hpc/cuppen/")){
#     saveRDS(ms_rep_strand_df, file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_rep_strand_df.all.rds")
#   } else {
#     saveRDS(ms_rep_strand_df, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_rep_strand_df.all.rds"))
#   }
# }


if (dir.exists("/hpc/cuppen/")){
  ms_rep_strand_df <- readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_rep_strand_df.all.rds")
} else {
  ms_rep_strand_df <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_rep_strand_df.all.rds"))
}





# Preparing the data for visualization (for effet size)

# if (args[1] == 3){
# 
#   if (dir.exists("/hpc/cuppen/")){
#     ms_rep_strand <- read.csv("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-rep-str-ann-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   } else {
#     ms_rep_strand <- read.csv("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-rep-str-ann-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   }
# 
# 
# 
#   ms_rep_strand_df2 <- data.frame(sample_id = character(7049*51), cosmic_sig = character(7049*51), right_replicating_count = numeric(7049*51), left_replicating_count = numeric(7049*51))
# 
# 
#   # head(ms_rep_strand_df2[!is.na(ms_rep_strand_df2$fold_change_right_to_left),])
#   # ms_rep_strand[ms_rep_strand$sample_id == "CPCT02030502T",]
# 
# 
#   # length(metadata_included$sample_id)
#   for (i in 1:length(metadata_included$sample_id)){
#     print("rep_str_ann")
#     print(i)
#     sample_id <- metadata_included$sample_id[i]
#     i <- i + (50*(i-1))
# 
#     for (j in 1:length(cosmic_sigs$V1)){
# 
#       ms_rep_strand_df2[i+j-1, "sample_id"] <- sample_id
#       ms_rep_strand_df2[i+j-1, "cosmic_sig"] <- cosmic_sigs$V1[j]
#       right <-  ms_rep_strand[ms_rep_strand$sample_id == sample_id & ms_rep_strand$rep_str_ann == "right",cosmic_sigs$V1[j]]
#       left <- ms_rep_strand[ms_rep_strand$sample_id == sample_id & ms_rep_strand$rep_str_ann == "left",cosmic_sigs$V1[j]]
# 
#       if (!is.na(right) & !is.na(left) & right != 0 & left != 0){
#         ms_rep_strand_df2[i+j-1, c("right_replicating_count", "left_replicating_count")] <- c(right, left)
#       } else {
#         ms_rep_strand_df2[i+j-1, c("right_replicating_count", "left_replicating_count")] <- NA
#       }
#     }
#   }
# 
# 
#   if (dir.exists("/hpc/cuppen/")){
#     saveRDS(ms_rep_strand_df2, file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_rep_strand_df2.all.rds")
#   } else {
#     saveRDS(ms_rep_strand_df2, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_rep_strand_df2.all.rds"))
#   }
# }


if (dir.exists("/hpc/cuppen/")){
  ms_rep_strand_df2 <- readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_rep_strand_df2.all.rds")
} else {
  ms_rep_strand_df2 <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_rep_strand_df2.all.rds"))
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
# 
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
#     mt[i:(i+num_row-1),"genomic_func"] <- as.character(levels(all.cols))
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
# mt[is.na(mt$genomic_func),"genomic_func"] <- "NA"
# 
# if (dir.exists("/hpc/cuppen/")){
#   write.table(mt, file = gzfile("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-genomic-func-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# } else {
#   write.table(mt, file = gzfile("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-genomic-func-2.5-percent-threshold.txt.gz"), sep = "\t", quote = F, row.names = F)
# }
# 
# 
# # Preparing the data for visualization (fold change)
# ms_genomic_func_df[ms_genomic_func_df$sample_id != "",]
# ms_genomic_func[ms_genomic_func$sample_id == "CPCT02230109T",]
# m <- 1
# i <- 16
# if (args[1] == 4){
# 
#   if (dir.exists("/hpc/cuppen/")){
#     ms_genomic_func <- read.csv("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-genomic-func-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   } else {
#     ms_genomic_func <- read.csv("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-genomic-func-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   }
#   ms_genomic_func[is.na(ms_genomic_func$genomic_func),"genomic_func"] <- "NA"
# 
# 
#   ms_genomic_func_df <- data.frame(sample_id = character(7049*51), cosmic_sig = character(7049*51), fold_change_coding_to_intron = numeric(7049*51), fold_change_coding_to_intergenic = numeric(7049*51))
# 
# 
#   # head(ms_genomic_func_df[!is.na(ms_genomic_func_df$fold_change_coding_to_intron),])
#   # ms_genomic_func[ms_genomic_func$sample_id == "CPCT02030502T",]
# 
# 
#   # length(metadata_included$sample_id)
#   for (i in 1:length(metadata_included$sample_id)){
#     # print("genomic_func")
#     print(i)
#     sample_id <- metadata_included$sample_id[i]
#     i <- i + (50*(i-1))
# 
#     for (j in 1:length(cosmic_sigs$V1)){
#       ms_genomic_func_df[i+j-1, "sample_id"] <- sample_id
#       ms_genomic_func_df[i+j-1, "cosmic_sig"] <- cosmic_sigs$V1[j]
#       coding <-  ms_genomic_func[ms_genomic_func$sample_id == sample_id & ms_genomic_func$genomic_func == "coding",cosmic_sigs$V1[j]]
#       intron <-  ms_genomic_func[ms_genomic_func$sample_id == sample_id & ms_genomic_func$genomic_func == "intron",cosmic_sigs$V1[j]]
#       intergenic <- ms_genomic_func[ms_genomic_func$sample_id == sample_id & ms_genomic_func$genomic_func == "intergenic",cosmic_sigs$V1[j]]
# 
#       if (length(coding) !=0 & length(intron) !=0){
#         if(!is.na(coding) & !is.na(intron) & coding != 0 & intron != 0){
#           ms_genomic_func_df[i+j-1, "fold_change_coding_to_intron"] <- coding/intron
#         } else {
#           ms_genomic_func_df[i+j-1, "fold_change_coding_to_intron"] <- NA
#         }
#       } else {
#         ms_genomic_func_df[i+j-1, "fold_change_coding_to_intron"] <- NA
#         }
# 
# 
#       if (length(coding) !=0 & length(intergenic) !=0){
#         if(!is.na(coding) & !is.na(intergenic) & coding != 0 & intergenic != 0) {
#           m <- m + 1
#           print(paste0("m is ",m))
#           ms_genomic_func_df[i+j-1, "fold_change_coding_to_intergenic"] <- coding/intergenic
#         } else {
#           ms_genomic_func_df[i+j-1, "fold_change_coding_to_intergenic"] <- NA
#         }
#       } else {
#         ms_genomic_func_df[i+j-1, "fold_change_coding_to_intergenic"] <- NA
#         }
# 
# 
#     }
#   }
# 
#   if (dir.exists("/hpc/cuppen/")){
#     saveRDS(ms_genomic_func_df, file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_genomic_func_df.all.rds")
#   } else {
#     saveRDS(ms_genomic_func_df, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_genomic_func_df.all.rds"))
#   }
# 
# }

if (dir.exists("/hpc/cuppen/")){
  ms_genomic_func_df <- readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_genomic_func_df.all.rds")
} else {
  ms_genomic_func_df <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_genomic_func_df.all.rds"))
}





# Preparing the data for visualization (fold change)

# if (args[1] == 4){
# 
#   if (dir.exists("/hpc/cuppen/")){
#     ms_genomic_func <- read.csv("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-genomic-func-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   } else {
#     ms_genomic_func <- read.csv("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms-genomic-func-2.5-percent-threshold.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
#   }
#   ms_genomic_func[is.na(ms_genomic_func$genomic_func),"genomic_func"] <- "NA"
# 
# 
#   ms_genomic_func_df2 <- data.frame(sample_id = character(7049*51), cosmic_sig = character(7049*51), coding_count = numeric(7049*51), intron_count = numeric(7049*51), intergenic_count = numeric(7049*51))
# 
# 
#   # head(ms_genomic_func_df2[!is.na(ms_genomic_func_df2$fold_change_coding_to_intron),])
#   # ms_genomic_func[ms_genomic_func$sample_id == "CPCT02030502T",]
# 
# 
#   # length(metadata_included$sample_id)
#   for (i in 1:length(metadata_included$sample_id)){
#     print("genomic_func")
#     print(i)
#     sample_id <- metadata_included$sample_id[i]
#     i <- i + (50*(i-1))
# 
#     for (j in 1:length(cosmic_sigs$V1)){
#       ms_genomic_func_df2[i+j-1, "sample_id"] <- sample_id
#       ms_genomic_func_df2[i+j-1, "cosmic_sig"] <- cosmic_sigs$V1[j]
#       coding <-  ms_genomic_func[ms_genomic_func$sample_id == sample_id & ms_genomic_func$genomic_func == "coding",cosmic_sigs$V1[j]]
#       intron <-  ms_genomic_func[ms_genomic_func$sample_id == sample_id & ms_genomic_func$genomic_func == "intron",cosmic_sigs$V1[j]]
#       intergenic <- ms_genomic_func[ms_genomic_func$sample_id == sample_id & ms_genomic_func$genomic_func == "intergenic",cosmic_sigs$V1[j]]
# 
# 
#       if (length(coding) !=0 & length(intergenic) !=0){
#         if(!is.na(coding) & !is.na(intergenic) & coding != 0 & intergenic != 0) {
#           ms_genomic_func_df2[i+j-1, c("coding_count", "intron_count", "intergenic_count")] <- c(coding, intron, intergenic)
#         } else {
#           ms_genomic_func_df2[i+j-1, c("coding_count", "intron_count", "intergenic_count")] <- NA
#         }
#       } else {
#         ms_genomic_func_df2[i+j-1, c("coding_count", "intron_count", "intergenic_count")] <- NA }
# 
# 
#     }
#   }
# 
#   if (dir.exists("/hpc/cuppen/")){
#     saveRDS(ms_genomic_func_df2, file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_genomic_func_df2.all.rds")
#   } else {
#     saveRDS(ms_genomic_func_df2, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_genomic_func_df2.all.rds"))
#   }
# 
# }

if (dir.exists("/hpc/cuppen/")){
  ms_genomic_func_df2 <- readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_genomic_func_df2.all.rds")
} else {
  ms_genomic_func_df2 <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/epi-genomic/ms_genomic_func_df2.all.rds"))
}





#####################################################################
# ms_chromatin      Fig19
#####################################################################
head(ms_chromatin_df)

df2 <- ms_chromatin_df

df2 <- merge(df2, metadata_included[,c("sample_id", "cancer_type", "cancer_type_code" ,"is_metastatic", "cohort")], by = "sample_id")

sigs <- unique(df2$cosmic_sig)
cancerTypes <- unique(df2$cancer_type)
cancerTypesCodes <- unique(df2$cancer_type_code)
cohort <- unique(df2$cohort)
nrow <- length(sigs)*length(cancerTypes)*length(cohort)



dff <- data.frame(cosmic_sig = character(nrow), cancer_types = character(nrow), cancer_types_code = character(nrow), cohort = character(nrow), median_foldChange = numeric(nrow), mean_foldChange = numeric(nrow), number_of_samples = integer(nrow))



for (i in 1:length(sigs)){
  for (j in 1:length(cancerTypes)){
    jj <- j + ((length(cohort)-1)*(j-1))
    for (k in 1:length(cohort)){
      tmp_df <- df2[df2$cosmic_sig == as.character(sigs[i]) & df2$cancer_type == cancerTypes[j] & df2$cohort == cohort[k] & !is.na(df2$fold_change_active_to_inactive),]
      dff[length(cancerTypes)*length(cohort)*(i-1)  + jj+k-1,1:4] <- c(as.character(sigs[i]), cancerTypes[j], cancerTypesCodes[j], cohort[k])
      if (nrow(tmp_df) >= 10){
        dff[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:7] <- c(median(tmp_df$fold_change_active_to_inactive, na.rm = T), mean(tmp_df$fold_change_active_to_inactive, na.rm = T), nrow(tmp_df))
      } else {
        dff[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:7] <- c(NA, NA, nrow(tmp_df))
      }
    }
  }
}

dff <- dff[!is.na(dff$median_foldChange),]



dff$cohort <- factor(dff$cohort)


# saveRDS(dff, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/dff.fig19.rds"))


dff <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/dff.fig19.rds"))

dff2 <- dff

summary(dff2$median_foldChange)
dff2$median_foldChange <- dff2$median_foldChange*(1/0.25679)

# dff2$is_metastatic <- as.logical(dff2$is_metastatic)

dff2_prim <- dff2[dff2$cohort == "PCAWG",]
dff2_metas <- dff2[dff2$cohort == "HMF",]




str(dff2_prim)
dot_plot_prim <- dff2_prim %>%
  ggplot(aes(x=cosmic_sig, y = cancer_types_code, fill=log2(median_foldChange)), size = 5) +
  theme_bw() +
  gggibbous::geom_moon(ratio=1, stroke = 0.3) +
  scale_fill_distiller(name = "log2(median of active/inactive regions)", palette='Spectral') +
  # scale_size_continuous(name="Effect size (Cohen's d)", range=c(2,12)) +
  # scale_linetype((name='Significance (< 0.05)')) +
  theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, size = 10),
        axis.text.y=element_text(vjust=0.5, size = 10),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 30, hjust = 0.5),
        plot.caption = element_text(size = 10)) +
  labs(caption = "10 sample + 2.5% threshold", title = "Mutational Signature Activity in Chromatin context in PCAWG Tumors\n") +
  xlab("\nCosmic signature") +
  ylab ("Cancer types")



for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig19-ms-chromatin-median-fraction-per-cancer-hmf-dotplot.png", height = 960, width = 1440)
    print(dot_plot_metas)
    dev.off()

    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig19-ms-chromatin-median-fraction-per-cancer-pcawg-dotplot.png", height = 960, width = 1440)
    print(dot_plot_prim)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig19-ms-chromatin-median-fraction-per-cancer-hmf-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_metas)
    dev.off()

    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig19-ms-chromatin-median-fraction-per-cancer-pcawg-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_prim)
    dev.off()
  }
}


### dotplot with comparison

df3 <- ms_chromatin_df2

df3 <- merge(df3, metadata_included[,c("sample_id", "cancer_type", "cancer_type_code" ,"is_metastatic", "cohort")], by = "sample_id")


df3 <- df3[!is.na(df3$active_count) | !is.na(df3$inactive_count),]
summary(df3$active_count)
summary(df3$inactive_count)
df3$active_count <- df3$active_count*(608/147)

sigs <- unique(df3$cosmic_sig)
cancerTypes <- unique(df3$cancer_type)
cancerTypesCodes <- unique(df3$cancer_type_code)
cohort <- unique(df3$cohort)
nrow <- length(sigs)*length(cancerTypes)*length(cohort)


dfff3 <- data.frame(cosmic_sig = character(nrow), cancer_types = character(nrow), cancer_types_code = character(nrow), cohort = character(nrow), median_active =  numeric(nrow), median_inactive = numeric(nrow), p_val = numeric(nrow), eff_size = numeric(nrow), number_of_samples = integer(nrow))


for (i in 1:length(sigs)){
  for (j in 1:length(cancerTypes)){
    jj <- j + ((length(cohort)-1)*(j-1))
    for (k in 1:length(cohort)){
      tmp_df <- df3[df3$cosmic_sig == as.character(sigs[i]) & df3$cancer_type == cancerTypes[j] & df3$cohort == cohort[k] & !is.na(df3$active_count),]
      active <- tmp_df$active_count
      inactive <- tmp_df$inactive_count
      dfff3[length(cancerTypes)*length(cohort)*(i-1)  + jj+k-1,1:4] <- c(as.character(sigs[i]), cancerTypes[j], cancerTypesCodes[j], cohort[k])
      if (nrow(tmp_df) >= 10){
        res <- wilcox.test(active, inactive)
        effsize_res <- cohen.d(active, inactive, hedges.correction=T)
        dfff3[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:8] <- c(median(active, na.rm = T), median(inactive, na.rm = T), res$p.value, abs(effsize_res$estimate))
        dfff3[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,9] <- c(nrow(tmp_df))
      } else {
        dfff3[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:8] <- rep(NA, times = 4)
        dfff3[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,9] <- c(nrow(tmp_df))
      }
    }
  }
}

dfff3 <- dfff3[!is.na(dfff3$p_val),]
dfff3$ratio_of_medians_active_to_inactive <- dfff3$median_active/dfff3$median_inactive
summary(dfff3$ratio_of_medians_active_to_inactive)

# dfff3$ratio_of_medians_active_to_inactive <- dfff3$ratio_of_medians_active_to_inactive*(1/1.0435)

dfff3 <- tidyr::gather(dfff3, key="active_inactive", value="median", c(median_active, median_inactive))
dfff3$active_inactive[dfff3$active_inactive == "median_active"] <- T
dfff3$active_inactive[dfff3$active_inactive == "median_inactive"] <- F
dfff3$active_inactive <- as.logical(dfff3$active_inactive)
nrow(dfff3[dfff3$active_inactive,])

str(dfff3)

dfff3$p_val <- dfff3$p_val < 0.05
dfff3$p_val <- factor(dfff3$p_val, levels = c("FALSE", "TRUE"))






# borderline color

dfff3$compare <- NA

sigs <- unique(dfff3$cosmic_sig)



for (i in 1:length(sigs)){
  tmp_df <- dfff3[dfff3$cosmic_sig == sigs[i],]
  cancerTypes <- unique(tmp_df$cancer_types)
  for (j in 1:length(cancerTypes)){
    tmp_df <- dfff3[dfff3$cosmic_sig == sigs[i] & dfff3$cancer_types == cancerTypes[j],]
    cohort <- unique(tmp_df$cohort)
    for (k in 1:length(cohort)){
      # print(sigs[i])
      # print(cancerTypes[j])
      # print(cohort[k])
      tmp_df <- dfff3[dfff3$cosmic_sig == sigs[i] & dfff3$cancer_types == cancerTypes[j] & dfff3$cohort == cohort[k],]
      if (tmp_df$median[tmp_df$active_inactive] >= tmp_df$median[!(tmp_df$active_inactive)]){
        dfff3[dfff3$cosmic_sig == sigs[i] & dfff3$cancer_types == cancerTypes[j] & dfff3$cohort == cohort[k],"compare"] <- "active"
      } else {
        dfff3[dfff3$cosmic_sig == sigs[i] & dfff3$cancer_types == cancerTypes[j] & dfff3$cohort == cohort[k],"compare"] <- "inactive"
      }
    }
  }
}


dfff3$compare <- factor(dfff3$compare)






dfff3_prim <- dfff3[dfff3$cohort == "PCAWG",]
dfff3_metas <- dfff3[dfff3$cohort == "HMF",]


# abs(log10(p_val))
dot_plot_prim <- dfff3_prim %>%
  ggplot(aes(x=cosmic_sig, y = cancer_types_code, size = eff_size)) +
  theme_bw() +
  gggibbous::geom_moon(aes(ratio=0.5, right=active_inactive, fill=log2(median), linetype = p_val, color = compare), stroke = 0.3) +
  scale_fill_distiller(name = "log2(median)", palette='Spectral') +
  scale_size_continuous(name="Effect size (Cohen's d)", range=c(2,12)) +
  scale_linetype((name='Significance (< 0.05)')) +
  theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, size = 10),
        axis.text.y=element_text(vjust=0.5, size = 10),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 30, hjust = 0.5),
        plot.caption = element_text(size = 10)) +
  labs(caption = "10 sample + 2.5% threshold \n left: inactive - right: active", title = "Mutational Signature Activity in Chromatin context in PCAWG Tumors\n") +
  xlab("\nCosmic signature") +
  ylab ("Cancer types")





for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig19-1-ms-chromatin-median-fraction-per-cancer-hmf-dotplot.png", height = 960, width = 1440)
    print(dot_plot_metas)
    dev.off()

    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig19-1-ms-chromatin-median-fraction-per-cancer-pcawg-dotplot.png", height = 960, width = 1440)
    print(dot_plot_prim)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig19-1-ms-chromatin-median-fraction-per-cancer-hmf-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_metas)
    dev.off()

    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig19-1-ms-chromatin-median-fraction-per-cancer-pcawg-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_prim)
    dev.off()
  }
}




#####################################################################
# ms_rep_timing    Fig20
#####################################################################
head(ms_rep_timing)


df2 <- ms_rep_timing_df

df2 <- merge(df2, metadata_included[,c("sample_id", "cancer_type", "cancer_type_code" ,"is_metastatic", "cohort")], by = "sample_id")

sigs <- unique(df2$cosmic_sig)
cancerTypes <- unique(df2$cancer_type)
cancerTypesCodes <- unique(df2$cancer_type_code)
cohort <- unique(df2$cohort)
nrow <- length(sigs)*length(cancerTypes)*length(cohort)



dff <- data.frame(cosmic_sig = character(nrow), cancer_types = character(nrow), cancer_types_code = character(nrow), cohort = character(nrow), median_foldChange = numeric(nrow), mean_foldChange = numeric(nrow), number_of_samples = integer(nrow))

# df2 <- df2[!is.na(df2$fold_change_late_to_early),]

for (i in 1:length(sigs)){
  for (j in 1:length(cancerTypes)){
    jj <- j + ((length(cohort)-1)*(j-1))
    for (k in 1:length(cohort)){
      tmp_df <- df2[df2$cosmic_sig == as.character(sigs[i]) & df2$cancer_type == cancerTypes[j] & df2$cohort == cohort[k] & !is.na(df2$fold_change_late_to_early),]
      dff[length(cancerTypes)*length(cohort)*(i-1)  + jj+k-1,1:4] <- c(as.character(sigs[i]), cancerTypes[j], cancerTypesCodes[j], cohort[k])
      if (nrow(tmp_df) >= 10){
        dff[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:7] <- c(median(tmp_df$fold_change_late_to_early, na.rm = T), mean(tmp_df$fold_change_late_to_early, na.rm = T), nrow(tmp_df))
      } else {
        dff[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:7] <- c(NA, NA, nrow(tmp_df))
      }
    }
  }
}

dff <- dff[!is.na(dff$median_foldChange),]



dff$cohort <- factor(dff$cohort)


# saveRDS(dff, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/dff.fig20.rds"))


dff <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/dff.fig20.rds"))

dff2 <- dff

summary(dff2$median_foldChange)
dff2$median_foldChange <- dff2$median_foldChange*(1/3.454)

# dff2$is_metastatic <- as.logical(dff2$is_metastatic)

dff2_prim <- dff2[dff2$cohort == "PCAWG",]
dff2_metas <- dff2[dff2$cohort == "HMF",]




str(dff2_prim)
dot_plot_metas <- dff2_metas %>%
  ggplot(aes(x=cosmic_sig, y = cancer_types_code, fill=log2(median_foldChange)), size = 5) +
  theme_bw() +
  gggibbous::geom_moon(ratio=1, stroke = 0.3) +
  scale_fill_distiller(name = "log2(median of late/early replicating regions)", palette='Spectral') +
  # scale_size_continuous(name="Effect size (Cohen's d)", range=c(2,12)) +
  # scale_linetype((name='Significance (< 0.05)')) +
  theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, size = 10),
        axis.text.y=element_text(vjust=0.5, size = 10),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 30, hjust = 0.5),
        plot.caption = element_text(size = 10)) +
  labs(caption = "10 sample + 2.5% threshold", title = "Mutational Signature Activity in Replication Timing context in HMF Tumors\n") +
  xlab("\nCosmic signature") +
  ylab ("Cancer types")



for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig20-ms-replication-timing-median-fraction-per-cancer-hmf-dotplot.png", height = 960, width = 1440)
    print(dot_plot_metas)
    dev.off()

    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig20-ms-replication-timing-median-fraction-per-cancer-pcawg-dotplot.png", height = 960, width = 1440)
    print(dot_plot_prim)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig20-ms-replication-timing-median-fraction-per-cancer-hmf-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_metas)
    dev.off()

    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig20-ms-replication-timing-median-fraction-per-cancer-pcawg-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_prim)
    dev.off()
  }
}






### dotplot with comparison

df3 <- ms_rep_timing_df2

df3 <- merge(df3, metadata_included[,c("sample_id", "cancer_type", "cancer_type_code" ,"is_metastatic", "cohort")], by = "sample_id")


df3 <- df3[!is.na(df3$late_count) | !is.na(df3$early_count),]
summary(df3$late_count)
summary(df3$early_count)
df3$early_count <- df3$early_count*(457/123)

sigs <- unique(df3$cosmic_sig)
cancerTypes <- unique(df3$cancer_type)
cancerTypesCodes <- unique(df3$cancer_type_code)
cohort <- unique(df3$cohort)
nrow <- length(sigs)*length(cancerTypes)*length(cohort)


dfff3 <- data.frame(cosmic_sig = character(nrow), cancer_types = character(nrow), cancer_types_code = character(nrow), cohort = character(nrow), median_late =  numeric(nrow), median_early = numeric(nrow), p_val = numeric(nrow), eff_size = numeric(nrow), number_of_samples = integer(nrow))


for (i in 1:length(sigs)){
  for (j in 1:length(cancerTypes)){
    jj <- j + ((length(cohort)-1)*(j-1))
    for (k in 1:length(cohort)){
      tmp_df <- df3[df3$cosmic_sig == as.character(sigs[i]) & df3$cancer_type == cancerTypes[j] & df3$cohort == cohort[k] & !is.na(df3$late_count),]
      late <- tmp_df$late_count
      early <- tmp_df$early_count
      dfff3[length(cancerTypes)*length(cohort)*(i-1)  + jj+k-1,1:4] <- c(as.character(sigs[i]), cancerTypes[j], cancerTypesCodes[j], cohort[k])
      if (nrow(tmp_df) >= 10){
        res <- wilcox.test(late, early)
        effsize_res <- cohen.d(late, early, hedges.correction=T)
        dfff3[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:8] <- c(median(late, na.rm = T), median(early, na.rm = T), res$p.value, abs(effsize_res$estimate))
        dfff3[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,9] <- c(nrow(tmp_df))
      } else {
        dfff3[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:8] <- rep(NA, times = 4)
        dfff3[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,9] <- c(nrow(tmp_df))
      }
    }
  }
}

dfff3 <- dfff3[!is.na(dfff3$p_val),]
dfff3$ratio_of_medians_late_to_early <- dfff3$median_late/dfff3$median_early
summary(dfff3$ratio_of_medians_late_to_early)

# dfff3$ratio_of_medians_late_to_early <- dfff3$ratio_of_medians_late_to_early*(1/0.25229)

dfff3 <- tidyr::gather(dfff3, key="late_early", value="median", c(median_late, median_early))
dfff3$late_early[dfff3$late_early == "median_late"] <- T
dfff3$late_early[dfff3$late_early == "median_early"] <- F
dfff3$late_early <- as.logical(dfff3$late_early)
nrow(dfff3[dfff3$late_early,])

str(dfff3)

dfff3$p_val <- dfff3$p_val < 0.05
dfff3$p_val <- factor(dfff3$p_val, levels = c("FALSE", "TRUE"))





# borderline color

dfff3$compare <- NA

sigs <- unique(dfff3$cosmic_sig)



for (i in 1:length(sigs)){
  tmp_df <- dfff3[dfff3$cosmic_sig == sigs[i],]
  cancerTypes <- unique(tmp_df$cancer_types)
  for (j in 1:length(cancerTypes)){
    tmp_df <- dfff3[dfff3$cosmic_sig == sigs[i] & dfff3$cancer_types == cancerTypes[j],]
    cohort <- unique(tmp_df$cohort)
    for (k in 1:length(cohort)){
      # print(sigs[i])
      # print(cancerTypes[j])
      # print(cohort[k])
      tmp_df <- dfff3[dfff3$cosmic_sig == sigs[i] & dfff3$cancer_types == cancerTypes[j] & dfff3$cohort == cohort[k],]
      if (tmp_df$median[tmp_df$late_early] >= tmp_df$median[!(tmp_df$late_early)]){
        dfff3[dfff3$cosmic_sig == sigs[i] & dfff3$cancer_types == cancerTypes[j] & dfff3$cohort == cohort[k],"compare"] <- "late"
      } else {
        dfff3[dfff3$cosmic_sig == sigs[i] & dfff3$cancer_types == cancerTypes[j] & dfff3$cohort == cohort[k],"compare"] <- "early"
      }
    }
  }
}


dfff3$compare <- factor(dfff3$compare)








dfff3_prim <- dfff3[dfff3$cohort == "PCAWG",]
dfff3_metas <- dfff3[dfff3$cohort == "HMF",]


# abs(log10(p_val))
dot_plot_prim <- dfff3_prim %>%
  ggplot(aes(x=cosmic_sig, y = cancer_types_code, size = eff_size)) +
  theme_bw() +
  gggibbous::geom_moon(aes(ratio=0.5, right=late_early, fill=log2(median), linetype = p_val, color = compare), stroke = 0.3) +
  scale_fill_distiller(name = "log2(median)", palette='Spectral') +
  scale_size_continuous(name="Effect size (Cohen's d)", range=c(2,12)) +
  scale_linetype((name='Significance (< 0.05)')) +
  theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, size = 10),
        axis.text.y=element_text(vjust=0.5, size = 10),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 30, hjust = 0.5),
        plot.caption = element_text(size = 10)) +
  labs(caption = "10 sample + 2.5% threshold \n left: early - right: late", title = "Mutational Signature Activity in Replication Timing context in PCAWG Tumors\n") +
  xlab("\nCosmic signature") +
  ylab ("Cancer types")





for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig20-1-ms-replication-timing-median-fraction-per-cancer-hmf-dotplot.png", height = 960, width = 1440)
    print(dot_plot_metas)
    dev.off()

    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig20-1-ms-replication-timing-median-fraction-per-cancer-pcawg-dotplot.png", height = 960, width = 1440)
    print(dot_plot_prim)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig20-1-ms-replication-timing-median-fraction-per-cancer-hmf-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_metas)
    dev.off()

    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig20-1-ms-replication-timing-median-fraction-per-cancer-pcawg-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_prim)
    dev.off()
  }
}





#####################################################################
# ms_trp_strand     Fig21
#####################################################################
head(ms_trp_strand_df)



df2 <- ms_trp_strand_df

df2 <- merge(df2, metadata_included[,c("sample_id", "cancer_type", "cancer_type_code" ,"is_metastatic", "cohort")], by = "sample_id")

sigs <- unique(df2$cosmic_sig)
cancerTypes <- unique(df2$cancer_type)
cancerTypesCodes <- unique(df2$cancer_type_code)
cohort <- unique(df2$cohort)
nrow <- length(sigs)*length(cancerTypes)*length(cohort)



dff <- data.frame(cosmic_sig = character(nrow), cancer_types = character(nrow), cancer_types_code = character(nrow), cohort = character(nrow), median_foldChange = numeric(nrow), mean_foldChange = numeric(nrow), number_of_samples = integer(nrow))

# df2 <- df2[!is.na(df2$fold_change_trans_to_untrans),]

for (i in 1:length(sigs)){
  for (j in 1:length(cancerTypes)){
    jj <- j + ((length(cohort)-1)*(j-1))
    for (k in 1:length(cohort)){
      tmp_df <- df2[df2$cosmic_sig == as.character(sigs[i]) & df2$cancer_type == cancerTypes[j] & df2$cohort == cohort[k] & !is.na(df2$fold_change_trans_to_untrans),]
      dff[length(cancerTypes)*length(cohort)*(i-1)  + jj+k-1,1:4] <- c(as.character(sigs[i]), cancerTypes[j], cancerTypesCodes[j], cohort[k])
      if (nrow(tmp_df) >= 10){
        dff[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:7] <- c(median(tmp_df$fold_change_trans_to_untrans, na.rm = T), mean(tmp_df$fold_change_trans_to_untrans, na.rm = T), nrow(tmp_df))
      } else {
        dff[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:7] <- c(NA, NA, nrow(tmp_df))
      }
    }
  }
}

dff <- dff[!is.na(dff$median_foldChange),]



dff$cohort <- factor(dff$cohort)


# saveRDS(dff, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/dff.fig21.rds"))


dff <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/dff.fig21.rds"))

dff2 <- dff

summary(dff2$median_foldChange)
dff2$median_foldChange <- dff2$median_foldChange*(1/0.9673)

# dff2$is_metastatic <- as.logical(dff2$is_metastatic)

dff2_prim <- dff2[dff2$cohort == "PCAWG",]
dff2_metas <- dff2[dff2$cohort == "HMF",]




str(dff2_prim)
dot_plot_prim <- dff2_prim %>%
  ggplot(aes(x=cosmic_sig, y = cancer_types_code, fill=log2(median_foldChange)), size = 5) +
  theme_bw() +
  gggibbous::geom_moon(ratio=1, stroke = 0.3) +
  scale_fill_distiller(name = "log2(median of transcribed/untranscribed strand regions)", palette='Spectral') +
  # scale_size_continuous(name="Effect size (Cohen's d)", range=c(2,12)) +
  # scale_linetype((name='Significance (< 0.05)')) +
  theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, size = 10),
        axis.text.y=element_text(vjust=0.5, size = 10),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 30, hjust = 0.5),
        plot.caption = element_text(size = 10)) +
  labs(caption = "10 sample + 2.5% threshold", title = "Mutational Signature Activity in Transcription Strand context in PCAWG Tumors\n") +
  xlab("\nCosmic signature") +
  ylab ("Cancer types")



for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig21-ms-transcription-strand-median-fraction-per-cancer-hmf-dotplot.png", height = 960, width = 1440)
    print(dot_plot_metas)
    dev.off()

    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig21-ms-transcription-strand-median-fraction-per-cancer-pcawg-dotplot.png", height = 960, width = 1440)
    print(dot_plot_prim)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig21-ms-transcription-strand-median-fraction-per-cancer-hmf-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_metas)
    dev.off()

    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig21-ms-transcription-strand-median-fraction-per-cancer-pcawg-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_prim)
    dev.off()
  }
}









### dotplot with comparison

df3 <- ms_trp_strand_df2

df3 <- merge(df3, metadata_included[,c("sample_id", "cancer_type", "cancer_type_code" ,"is_metastatic", "cohort")], by = "sample_id")


df3 <- df3[!is.na(df3$trans_count) | !is.na(df3$untrans_count),]
summary(df3$trans_count)
summary(df3$untrans_count)
df3$untrans_count <- df3$untrans_count*(196/202)

sigs <- unique(df3$cosmic_sig)
cancerTypes <- unique(df3$cancer_type)
cancerTypesCodes <- unique(df3$cancer_type_code)
cohort <- unique(df3$cohort)
nrow <- length(sigs)*length(cancerTypes)*length(cohort)


dfff3 <- data.frame(cosmic_sig = character(nrow), cancer_types = character(nrow), cancer_types_code = character(nrow), cohort = character(nrow), median_trans =  numeric(nrow), median_untrans = numeric(nrow), p_val = numeric(nrow), eff_size = numeric(nrow), number_of_samples = integer(nrow))


for (i in 1:length(sigs)){
  for (j in 1:length(cancerTypes)){
    jj <- j + ((length(cohort)-1)*(j-1))
    for (k in 1:length(cohort)){
      tmp_df <- df3[df3$cosmic_sig == as.character(sigs[i]) & df3$cancer_type == cancerTypes[j] & df3$cohort == cohort[k] & !is.na(df3$trans_count),]
      trans <- tmp_df$trans_count
      untrans <- tmp_df$untrans_count
      dfff3[length(cancerTypes)*length(cohort)*(i-1)  + jj+k-1,1:4] <- c(as.character(sigs[i]), cancerTypes[j], cancerTypesCodes[j], cohort[k])
      if (nrow(tmp_df) >= 10){
        res <- wilcox.test(trans, untrans)
        effsize_res <- cohen.d(trans, untrans, hedges.correction=T)
        dfff3[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:8] <- c(median(trans, na.rm = T), median(untrans, na.rm = T), res$p.value, abs(effsize_res$estimate))
        dfff3[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,9] <- c(nrow(tmp_df))
      } else {
        dfff3[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:8] <- rep(NA, times = 4)
        dfff3[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,9] <- c(nrow(tmp_df))
      }
    }
  }
}

dfff3 <- dfff3[!is.na(dfff3$p_val),]
dfff3$ratio_of_medians_trans_to_untrans <- dfff3$median_trans/dfff3$median_untrans
summary(dfff3$ratio_of_medians_trans_to_untrans)

dfff3$ratio_of_medians_trans_to_untrans <- dfff3$ratio_of_medians_trans_to_untrans*(1/0.25229)

dfff3 <- tidyr::gather(dfff3, key="trans_untrans", value="median", c(median_trans, median_untrans))
dfff3$trans_untrans[dfff3$trans_untrans == "median_trans"] <- T
dfff3$trans_untrans[dfff3$trans_untrans == "median_untrans"] <- F
dfff3$trans_untrans <- as.logical(dfff3$trans_untrans)
nrow(dfff3[dfff3$trans_untrans,])

str(dfff3)

dfff3$p_val <- dfff3$p_val < 0.05
dfff3$p_val <- factor(dfff3$p_val, levels = c("FALSE", "TRUE"))




dfff3_prim <- dfff3[dfff3$cohort == "PCAWG",]
dfff3_metas <- dfff3[dfff3$cohort == "HMF",]


# abs(log10(p_val))
dot_plot_prim <- dfff3_prim %>%
  ggplot(aes(x=cosmic_sig, y = cancer_types_code, size = eff_size)) +
  theme_bw() +
  gggibbous::geom_moon(aes(ratio=0.5, right=trans_untrans, fill=log2(median), linetype = p_val), stroke = 0.3) +
  scale_fill_distiller(name = "log2(median)", palette='Spectral') +
  scale_size_continuous(name="Effect size (Cohen's d)", range=c(2,12)) +
  scale_linetype((name='Significance (< 0.05)')) +
  theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, size = 10),
        axis.text.y=element_text(vjust=0.5, size = 10),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 30, hjust = 0.5),
        plot.caption = element_text(size = 10)) +
  labs(caption = "10 sample + 2.5% threshold \n left: untrans - right: trans", title = "Mutational Signature Activity in Transcription Strand context in PCAWG Tumors\n") +
  xlab("\nCosmic signature") +
  ylab ("Cancer types")





for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig21-1-ms-transcription-strand-median-fraction-per-cancer-hmf-dotplot.png", height = 960, width = 1440)
    print(dot_plot_metas)
    dev.off()

    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig21-1-ms-transcription-strand-median-fraction-per-cancer-pcawg-dotplot.png", height = 960, width = 1440)
    print(dot_plot_prim)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig21-1-ms-transcription-strand-median-fraction-per-cancer-hmf-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_metas)
    dev.off()

    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig21-1-ms-transcription-strand-median-fraction-per-cancer-pcawg-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_prim)
    dev.off()
  }
}








#####################################################################
# ms_rep_strand       Fig22
#####################################################################
head(ms_rep_strand_df)



df2 <- ms_rep_strand_df

df2 <- merge(df2, metadata_included[,c("sample_id", "cancer_type", "cancer_type_code" ,"is_metastatic", "cohort")], by = "sample_id")

sigs <- unique(df2$cosmic_sig)
cancerTypes <- unique(df2$cancer_type)
cancerTypesCodes <- unique(df2$cancer_type_code)
cohort <- unique(df2$cohort)
nrow <- length(sigs)*length(cancerTypes)*length(cohort)



dff <- data.frame(cosmic_sig = character(nrow), cancer_types = character(nrow), cancer_types_code = character(nrow), cohort = character(nrow), median_foldChange = numeric(nrow), mean_foldChange = numeric(nrow), number_of_samples = integer(nrow))

# df2 <- df2[!is.na(df2$fold_change_right_to_left),]

for (i in 1:length(sigs)){
  for (j in 1:length(cancerTypes)){
    jj <- j + ((length(cohort)-1)*(j-1))
    for (k in 1:length(cohort)){
      tmp_df <- df2[df2$cosmic_sig == as.character(sigs[i]) & df2$cancer_type == cancerTypes[j] & df2$cohort == cohort[k] & !is.na(df2$fold_change_right_to_left),]
      dff[length(cancerTypes)*length(cohort)*(i-1)  + jj+k-1,1:4] <- c(as.character(sigs[i]), cancerTypes[j], cancerTypesCodes[j], cohort[k])
      if (nrow(tmp_df) >= 10){
        dff[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:7] <- c(median(tmp_df$fold_change_right_to_left, na.rm = T), mean(tmp_df$fold_change_right_to_left, na.rm = T), nrow(tmp_df))
      } else {
        dff[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:7] <- c(NA, NA, nrow(tmp_df))
      }
    }
  }
}

dff <- dff[!is.na(dff$median_foldChange),]



dff$cohort <- factor(dff$cohort)


# saveRDS(dff, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/dff.fig22.rds"))


dff <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/dff.fig22.rds"))

dff2 <- dff

summary(dff2$median_foldChange)
dff2$median_foldChange <- dff2$median_foldChange*(1/1.0040)

# dff2$is_metastatic <- as.logical(dff2$is_metastatic)

dff2_prim <- dff2[dff2$cohort == "PCAWG",]
dff2_metas <- dff2[dff2$cohort == "HMF",]




str(dff2_prim)
dot_plot_prim <- dff2_prim %>%
  ggplot(aes(x=cosmic_sig, y = cancer_types_code, fill=log2(median_foldChange)), size = 5) +
  theme_bw() +
  gggibbous::geom_moon(ratio=1, stroke = 0.3) +
  scale_fill_distiller(name = "log2(median of right/left replicating regions)", palette='Spectral') +
  # scale_size_continuous(name="Effect size (Cohen's d)", range=c(2,12)) +
  # scale_linetype((name='Significance (< 0.05)')) +
  theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, size = 10),
        axis.text.y=element_text(vjust=0.5, size = 10),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 30, hjust = 0.5),
        plot.caption = element_text(size = 10)) +
  labs(caption = "10 sample + 2.5% threshold", title = "Mutational Signature Activity in Replication Orientation context in PCAWG Tumors\n") +
  xlab("\nCosmic signature") +
  ylab ("Cancer types")



for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig22-ms-replication-strand-median-fraction-per-cancer-hmf-dotplot.png", height = 960, width = 1440)
    print(dot_plot_metas)
    dev.off()

    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig22-ms-replication-strand-median-fraction-per-cancer-pcawg-dotplot.png", height = 960, width = 1440)
    print(dot_plot_prim)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig22-ms-replication-strand-median-fraction-per-cancer-hmf-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_metas)
    dev.off()

    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig22-ms-replication-strand-median-fraction-per-cancer-pcawg-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_prim)
    dev.off()
  }
}









### dotplot with comparison

df3 <- ms_rep_strand_df2

df3 <- merge(df3, metadata_included[,c("sample_id", "cancer_type", "cancer_type_code" ,"is_metastatic", "cohort")], by = "sample_id")


df3 <- df3[!is.na(df3$right_replicating_count) | !is.na(df3$left_replicating_count),]
summary(df3$right_replicating_count)
summary(df3$left_replicating_count)
df3$left_replicating_count <- df3$left_replicating_count*(208/206)

sigs <- unique(df3$cosmic_sig)
cancerTypes <- unique(df3$cancer_type)
cancerTypesCodes <- unique(df3$cancer_type_code)
cohort <- unique(df3$cohort)
nrow <- length(sigs)*length(cancerTypes)*length(cohort)


dfff3 <- data.frame(cosmic_sig = character(nrow), cancer_types = character(nrow), cancer_types_code = character(nrow), cohort = character(nrow), median_right_replicating =  numeric(nrow), median_left_replicating = numeric(nrow), p_val = numeric(nrow), eff_size = numeric(nrow), number_of_samples = integer(nrow))


for (i in 1:length(sigs)){
  for (j in 1:length(cancerTypes)){
    jj <- j + ((length(cohort)-1)*(j-1))
    for (k in 1:length(cohort)){
      tmp_df <- df3[df3$cosmic_sig == as.character(sigs[i]) & df3$cancer_type == cancerTypes[j] & df3$cohort == cohort[k] & !is.na(df3$right_replicating_count),]
      right_replicating <- tmp_df$right_replicating_count
      left_replicating <- tmp_df$left_replicating_count
      dfff3[length(cancerTypes)*length(cohort)*(i-1)  + jj+k-1,1:4] <- c(as.character(sigs[i]), cancerTypes[j], cancerTypesCodes[j], cohort[k])
      if (nrow(tmp_df) >= 10){
        res <- wilcox.test(right_replicating, left_replicating)
        effsize_res <- cohen.d(right_replicating, left_replicating, hedges.correction=T)
        dfff3[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:8] <- c(median(right_replicating, na.rm = T), median(left_replicating, na.rm = T), res$p.value, abs(effsize_res$estimate))
        dfff3[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,9] <- c(nrow(tmp_df))
      } else {
        dfff3[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:8] <- rep(NA, times = 4)
        dfff3[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,9] <- c(nrow(tmp_df))
      }
    }
  }
}

dfff3 <- dfff3[!is.na(dfff3$p_val),]
dfff3$ratio_of_medians_right_replicating_to_left_replicating <- dfff3$median_right_replicating/dfff3$median_left_replicating
summary(dfff3$ratio_of_medians_right_replicating_to_left_replicating)

dfff3$ratio_of_medians_right_replicating_to_left_replicating <- dfff3$ratio_of_medians_right_replicating_to_left_replicating*(1/0.25229)

dfff3 <- tidyr::gather(dfff3, key="right_replicating_left_replicating", value="median", c(median_right_replicating, median_left_replicating))
dfff3$right_replicating_left_replicating[dfff3$right_replicating_left_replicating == "median_right_replicating"] <- T
dfff3$right_replicating_left_replicating[dfff3$right_replicating_left_replicating == "median_left_replicating"] <- F
dfff3$right_replicating_left_replicating <- as.logical(dfff3$right_replicating_left_replicating)
nrow(dfff3[dfff3$right_replicating_left_replicating,])

str(dfff3)

dfff3$p_val <- dfff3$p_val < 0.05
dfff3$p_val <- factor(dfff3$p_val, levels = c("FALSE", "TRUE"))




dfff3_prim <- dfff3[dfff3$cohort == "PCAWG",]
dfff3_metas <- dfff3[dfff3$cohort == "HMF",]


# abs(log10(p_val))
dot_plot_metas <- dfff3_metas %>%
  ggplot(aes(x=cosmic_sig, y = cancer_types_code, size = eff_size)) +
  theme_bw() +
  gggibbous::geom_moon(aes(ratio=0.5, right=right_replicating_left_replicating, fill=log2(median), linetype = p_val), stroke = 0.3) +
  scale_fill_distiller(name = "log2(median)", palette='Spectral') +
  scale_size_continuous(name="Effect size (Cohen's d)", range=c(2,12)) +
  scale_linetype((name='Significance (< 0.05)')) +
  theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, size = 10),
        axis.text.y=element_text(vjust=0.5, size = 10),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 30, hjust = 0.5),
        plot.caption = element_text(size = 10)) +
  labs(caption = "10 sample + 2.5% threshold \n left: left_replicating - right: right_replicating", title = "Mutational Signature Activity in Replication Strand context in HMF Tumors\n") +
  xlab("\nCosmic signature") +
  ylab ("Cancer types")





for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig22-1-ms-replication-strand-median-fraction-per-cancer-hmf-dotplot.png", height = 960, width = 1440)
    print(dot_plot_metas)
    dev.off()

    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig22-1-ms-replication-strand-median-fraction-per-cancer-pcawg-dotplot.png", height = 960, width = 1440)
    print(dot_plot_prim)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig22-1-ms-replication-strand-median-fraction-per-cancer-hmf-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_metas)
    dev.off()

    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig22-1-ms-replication-strand-median-fraction-per-cancer-pcawg-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_prim)
    dev.off()
  }
}







#####################################################################
# ms_genomic_func        Fig23
#####################################################################
head(ms_genomic_func_df)

nrow(ms_genomic_func_df[!is.na(ms_genomic_func_df$fold_change_coding_to_intergenic) & ms_genomic_func_df$fold_change_coding_to_intergenic != 0,])
nrow(df2[!is.na(df2$fold_change_coding_to_intergenic),])


df2 <- ms_genomic_func_df

df2 <- merge(df2, metadata_included[,c("sample_id", "cancer_type", "cancer_type_code" ,"is_metastatic", "cohort")], by = "sample_id")


sigs <- unique(df2$cosmic_sig)
cancerTypes <- unique(df2$cancer_type)
cancerTypesCodes <- unique(df2$cancer_type_code)
cohort <- unique(df2$cohort)
nrow <- length(sigs)*length(cancerTypes)*length(cohort)



dff <- data.frame(cosmic_sig = character(nrow), cancer_types = character(nrow), cancer_types_code = character(nrow), cohort = character(nrow), median_foldChange = numeric(nrow), mean_foldChange = numeric(nrow), number_of_samples = integer(nrow))

df2 <- df2[!is.na(df2$fold_change_coding_to_intergenic),]
i <- 1
j <- 2
k <- 1

for (i in 1:length(sigs)){
  for (j in 1:length(cancerTypes)){
    jj <- j + ((length(cohort)-1)*(j-1))
    for (k in 1:length(cohort)){
      tmp_df <- df2[df2$cosmic_sig == as.character(sigs[i]) & df2$cancer_type == cancerTypes[j] & df2$cohort == cohort[k] & !is.na(df2$fold_change_coding_to_intergenic),]
      dff[length(cancerTypes)*length(cohort)*(i-1)  + jj+k-1,1:4] <- c(as.character(sigs[i]), cancerTypes[j], cancerTypesCodes[j], cohort[k])
      if (nrow(tmp_df) >= 10){
        dff[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:7] <- c(median(tmp_df$fold_change_coding_to_intergenic, na.rm = T), mean(tmp_df$fold_change_coding_to_intergenic, na.rm = T), nrow(tmp_df))
      } else {
        dff[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:7] <- c(NA, NA, nrow(tmp_df))
      }
    }
  }
}

dff <- dff[!is.na(dff$median_foldChange),]



dff$cohort <- factor(dff$cohort)


# saveRDS(dff, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/dff.fig23.rds"))


dff <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/dff.fig23.rds"))

dff2 <- dff

summary(dff2$median_foldChange)
dff2$median_foldChange <- dff2$median_foldChange*(1/1.0040)

# dff2$is_metastatic <- as.logical(dff2$is_metastatic)

dff2_prim <- dff2[dff2$cohort == "PCAWG",]
dff2_metas <- dff2[dff2$cohort == "HMF",]




str(dff2_prim)
dot_plot_prim <- dff2_prim %>%
  ggplot(aes(x=cosmic_sig, y = cancer_types_code, fill=log2(median_foldChange)), size = 5) +
  theme_bw() +
  gggibbous::geom_moon(ratio=1, stroke = 0.3) +
  scale_fill_distiller(name = "log2(median of coding/intergenic regions)", palette='Spectral') +
  # scale_size_continuous(name="Effect size (Cohen's d)", range=c(2,12)) +
  # scale_linetype((name='Significance (< 0.05)')) +
  theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, size = 10),
        axis.text.y=element_text(vjust=0.5, size = 10),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 30, hjust = 0.5),
        plot.caption = element_text(size = 10)) +
  labs(caption = "10 sample + 2.5% threshold", title = "Mutational Signature Activity in Genomic Functional Units context in PCAWG Tumors\n") +
  xlab("\nCosmic signature") +
  ylab ("Cancer types")



for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig23-ms-genomic-function-coding-to-intergenic-median-fraction-per-cancer-hmf-dotplot.png", height = 960, width = 1440)
    print(dot_plot_metas)
    dev.off()

    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig23-ms-genomic-function-coding-to-intergenic-median-fraction-per-cancer-pcawg-dotplot.png", height = 960, width = 1440)
    print(dot_plot_prim)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig23-ms-genomic-function-coding-to-intergenic-median-fraction-per-cancer-hmf-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_metas)
    dev.off()

    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig23-ms-genomic-function-coding-to-intergenic-median-fraction-per-cancer-pcawg-dotplot.pdf", height = 14, width = 21)
    print(dot_plot_prim)
    dev.off()
  }
}



