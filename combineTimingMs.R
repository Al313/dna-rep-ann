

# Loading required libraries

library(epitools)
library(rcompanion)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(gridExtra)
library(gghighlight)
library(stringr)
library(dplyr)


# args <- commandArgs(trailingOnly=T)


# Assigning path directories

local <- "/home/ali313/Documents/studies/master/umc-project"





if (dir.exists("/hpc/cuppen/")){
  wd <- "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/"
} else {
  wd <- paste0(local,"/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/")
}




# Reading in the data

if (dir.exists("/hpc/cuppen/")){
  hmf_meta <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/metadata_whitelisted.tsv", sep = "\t", header = T, stringsAsFactors = F)
} else {
  hmf_meta <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/metadata_whitelisted.tsv", sep = "\t", header = T, stringsAsFactors = F)
}




if (dir.exists("/hpc/cuppen/")){
  metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2-metadata_23072021.tsv"), sep = "\t", header = T, stringsAsFactors = F)
} else {
  metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2-metadata_23072021.tsv"), sep = "\t", header = T, stringsAsFactors = F)
}

metadata_included <- metadata[!(metadata$is_blacklisted),]

metadata_included$tmb <- rowSums(metadata_included[,22:23])



### VCF annotated with MS and timing information:



# mut_tim <- read.csv(file = paste0(local, "/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/timing/DO36159/DO36159.mutationaltiming.tsv.gz"), header = T, sep = "\t", stringsAsFactors = F)
# mut_tim2 <- read.csv(file = paste0(local, "/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/timing/151111_HMFreg0017_HMF0206_HMF0207_CPCT02010240/CPCT02010240T.mutationaltiming.tsv.gz"), header = T, sep = "\t", stringsAsFactors = F)
# mut_tim$CHROM <- as.character(mut_tim$CHROM)
# # 
# # 
# mut_sig <- read.csv(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/muts_assigned/DO36159.txt.gz"), header = T, sep = "\t", stringsAsFactors = F)
# mut_sig2 <- read.csv(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/muts_assigned/CPCT02010240T.txt.gz"), header = T, sep = "\t", stringsAsFactors = F)
# # 
# colnames(mut_sig) <- toupper(colnames(mut_sig))
# # 
# # 
# str(mut_tim)
# str(mut_sig)
# mut_sig$CHROM <- sapply(str_split(mut_sig$CHROM, pattern = "r"), "[[", 2)
# # 
# shared_df <- intersect(mut_tim[,1:4], mut_sig[,1:4])
# 
# colnames(mut_tim)
# combined_df <- merge(merge(shared_df, mut_sig, by = c("CHROM", "POS", "REF", "ALT")), mut_tim[,c(1:4,17:18,24:26)], by = c("CHROM", "POS", "REF", "ALT"))
# 
# 
# 
# 
# write.table(combined_df, file = gzfile(paste0(wd, "timing-ms-combined/", sample_id, ".txt.gz")), quote = F, row.names = F, sep = "\t")
# 
# 
# metadata_included[metadata_included$cohort == "PCAWG",]


# pcawg : c(4786, 4787)
# hmf : 1:nrow(metadata_included[1:2,])
# CPCT02330041T

# timing_absence_counter <- 0
# signature_absence_counter <- 0
# timing_absent <- character()
# signature_absent <- character()
# 
# for (i in args[1]:args[2]){
#   sample_id <- metadata_included$sample_id[i]
#   # print(i)
#   # print(sample_id)
#   
#   if (metadata_included$cohort[i] == "HMF"){
#     
#     
#     set_name <- hmf_meta$setName[hmf_meta$sampleId == sample_id]
#     
#     if (dir.exists("/hpc/cuppen/")){
#       mut_timing_path <- paste0("/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/timing/", set_name, "/", sample_id, ".mutationaltiming.tsv.gz")
#     } else {
#       mut_timing_path <- paste0(local, "/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/timing/", set_name, "/", sample_id, ".mutationaltiming.tsv.gz")
#     }
#     
#     if (file.exists(mut_timing_path)){
#       mut_timing_file <- T
#       mut_tim <- tryCatch(read.csv(file = mut_timing_path, header = T, sep = "\t", stringsAsFactors = F), error=function(e) NULL)
#       if (is.null(mut_tim)) {
#         mut_timing_file <- F
#         timing_absence_counter <- timing_absence_counter + 1
#         timing_absent <- append(timing_absent, sample_id)
#       }
#     } else {
#       mut_timing_file <- F
#       timing_absence_counter <- timing_absence_counter + 1
#       timing_absent <- append(timing_absent, sample_id)
#       # print(paste0("Mutational timing file for ", sample_id, " does not exist!"))
#     }
#     
#     
#     if (dir.exists("/hpc/cuppen/")){
#       mut_signature_path <- paste0("/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/muts_assigned/", sample_id, ".txt.gz")
#     } else {
#       mut_signature_path <- paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/muts_assigned/", sample_id, ".txt.gz")
#     }
#     
#     if (file.exists(mut_signature_path)){
#       mut_signature_file <- T
#       # mut_sig <- read.csv(file = mut_signature_path, header = T, sep = "\t", stringsAsFactors = F)
#       mut_sig <- tryCatch(read.csv(file = mut_signature_path, header = T, sep = "\t", stringsAsFactors = F), error=function(e) NULL)
#       if (is.null(mut_sig)) {
#         mut_signature_file <- F
#         signature_absence_counter <- signature_absence_counter + 1
#         signature_absent <- append(signature_absent, sample_id)
#       }
#     } else {
#       mut_signature_file <- F
#       signature_absence_counter <- signature_absence_counter + 1
#       signature_absent <- append(signature_absent, sample_id)
#       # print(paste0("Mutational signature file for ", sample_id, " does not exist!"))
#     }
#     
#     
#   } else if (metadata_included$cohort[i] == "PCAWG"){
#     
#     
#     if (dir.exists("/hpc/cuppen/")){
#       mut_timing_path <- paste0("/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/timing/", sample_id, "/", sample_id, ".mutationaltiming.tsv.gz")
#     } else {
#       mut_timing_path <- paste0(local, "/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/timing/", sample_id, "/", sample_id, ".mutationaltiming.tsv.gz")
#     }
#     
#     if (file.exists(mut_timing_path)){
#       mut_timing_file <- T
#       mut_tim <- tryCatch(read.csv(file = mut_timing_path, header = T, sep = "\t", stringsAsFactors = F), error=function(e) NULL)
#       if (is.null(mut_tim)) {
#         mut_timing_file <- F
#         timing_absence_counter <- timing_absence_counter + 1
#         timing_absent <- append(timing_absent, sample_id)
#       }
#     } else {
#       mut_timing_file <- F
#       timing_absence_counter <- timing_absence_counter + 1
#       timing_absent <- append(timing_absent, sample_id)
#       # print(paste0("Mutational timing file for ", sample_id, " does not exist!"))
#     }
#     
#     
#     if (dir.exists("/hpc/cuppen/")){
#       mut_signature_path <- paste0("/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/muts_assigned/", sample_id, ".txt.gz")
#     } else {
#       mut_signature_path <- paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/muts_assigned/", sample_id, ".txt.gz")
#     }
#     
#     if (file.exists(mut_signature_path)){
#       mut_signature_file <- T
#       mut_sig <- tryCatch(read.csv(file = mut_signature_path, header = T, sep = "\t", stringsAsFactors = F), error=function(e) NULL)
#       if (is.null(mut_sig)) {
#         mut_signature_file <- F
#         signature_absence_counter <- signature_absence_counter + 1
#         signature_absent <- append(signature_absent, sample_id)
#       }
#     } else {
#       mut_signature_file <- F
#       signature_absence_counter <- signature_absence_counter + 1
#       signature_absent <- append(signature_absent, sample_id)
#       # print(paste0("Mutational signature file for ", sample_id, " does not exist!"))
#     }
#     
#     
#   }
#   
#   if (mut_timing_file & mut_signature_file){
#     
#     mut_tim$CHROM <- as.character(mut_tim$CHROM)
#     
#     colnames(mut_sig) <- toupper(colnames(mut_sig))
#     
#     mut_sig$CHROM <- sapply(str_split(mut_sig$CHROM, pattern = "r"), "[[", 2)
#     
#     shared_df <- intersect(mut_tim[,1:4], mut_sig[,1:4])
#     
#     
#     combined_df <- merge(merge(shared_df, mut_sig, by = c("CHROM", "POS", "REF", "ALT")), mut_tim[,c(1:4,17:18,24:26)], by = c("CHROM", "POS", "REF", "ALT"))
#     
#     
#     
#     write.table(combined_df, file = gzfile(paste0(wd, "timing-ms-combined/", sample_id, ".txt.gz")), quote = F, row.names = F, sep = "\t")
#   }
# }
# 
# 
# cat('timing_absence_counter equals to \n', timing_absence_counter)
# cat("\n timing absent file list: \n", paste(timing_absent, collapse = " "))
# cat("\n signature_absence_counter equals to \n", signature_absence_counter)
# cat("\n signature absent file list: \n", paste(signature_absent, collapse = " "))
# 
# 




### WGD molecular timing information

# wgd_timing_df <- data.frame(sample_id = character(), isWGD = logical(), molecular_timing = numeric(), cohort = character(), is_metastatic = logical())
# 
# 
# for (i in 1:nrow(metadata_included)){
#   
#   sample_id <- metadata_included$sample_id[i]
#   print(i)
#   print(sample_id)
#   
#   if (metadata_included$cohort[i] == "HMF"){
#     
#     cohort = "HMF"
#     is_metastatic = metadata_included$is_metastatic[i]
#     
#     set_name <- hmf_meta$setName[hmf_meta$sampleId == sample_id]
#     # print(set_name)
#     
#     if (dir.exists("/hpc/cuppen/")){
#       wgd_timing_path <- paste0("/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/timing/", set_name, "/", sample_id, ".WGD_time.tsv")
#     } else {
#       wgd_timing_path <- paste0(local, "/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/timing/", set_name, "/", sample_id, ".WGD_time.tsv")
#     }
#     
#     
#     if (file.exists(wgd_timing_path)){
#       wgd_timing_file <- T
#     } else {
#       wgd_timing_file <- F
#     }
#     
#   }
#     
#     
#   if (metadata_included$cohort[i] == "PCAWG"){
#     
#     cohort = "PCAWG"
#     is_metastatic = metadata_included$is_metastatic[i]
#     
#     if (dir.exists("/hpc/cuppen/")){
#       wgd_timing_path <- paste0("/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/timing/", sample_id, "/", sample_id, ".WGD_time.tsv")
#     } else {
#       wgd_timing_path <- paste0(local, "/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/timing/", sample_id, "/", sample_id, ".WGD_time.tsv")
#     }
#     
#     
#     if (file.exists(wgd_timing_path)){
#       wgd_timing_file <- T
#     } else {
#       wgd_timing_file <- F
#     }
#     
#   }
#     
#     
#     
#     
#     
#     
#     
#     
#     if (wgd_timing_file){
#       wgd_timing <- read.csv(file = wgd_timing_path, header = T, sep = "\t", stringsAsFactors = F)
#       if (nrow(wgd_timing) > 1){
#         print(paste0("Number of rows for ", sample_id, " is: ", nrow(wgd_timing)))
#       }
#       
#       wgd_timing_df <- rbind(wgd_timing_df, cbind(wgd_timing[1,c(2,3,1)], cohort, is_metastatic))
#       
#     } else {
#       wgd_timing_df <- rbind(wgd_timing_df, c(sample_id, NA, NA, cohort, is_metastatic))
#     }
#     
#   
#   
# }
# 
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   write.table(wgd_timing_df, file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/wgd-timing.txt", quote = F, row.names = F, sep = "\t")
# } else {
#   write.table(wgd_timing_df, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/wgd-timing.txt", quote = F, row.names = F, sep = "\t")
# }


## Those samples that their molecular timing is not determined (is NA in the file) have more than one row. Specifically they have 100 rows which are duplicated.





if (dir.exists("/hpc/cuppen/")){
  wgd_timing_df <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/wgd-timing.txt", header = T, stringsAsFactors = F, sep = "\t")
} else {
  wgd_timing_df <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/wgd-timing.txt", header = T, stringsAsFactors = F, sep = "\t")
}

wgd_timing_df <- merge(wgd_timing_df, metadata_included[,c("sample_id", "cancer_type")], by.x = "sample", by.y = "sample_id")
# str(wgd_timing_df)
# sum(as.numeric(table(wgd_timing_df$cohort, useNA = "always")))
# sum(as.numeric(table(wgd_timing_df$is_metastatic, useNA = "always")))
# table(wgd_timing_df$isWGD, useNA = "always")                  # for 897 samples the wgd timing file is missing
# table(wgd_timing_df$molecular_timing, useNA = "always")       # for 1410 samples molecular timing information is missing





#### Summarize mutational timing


# for (i in 1:nrow(metadata_included)){
#   print(i)
#   sample_id <- metadata_included$sample_id[i]
#   file_path <- paste0(wd, "timing-ms-combined/", sample_id, ".txt.gz")
#   if (file.exists(file_path)){
#     combined_df <- read.csv(file = file_path, header = T, sep = "\t", stringsAsFactors = F)
#     uu <- as.data.frame(table(combined_df$timing_class), stringsAsFactors = F)
# 
#     for (category in c("clonal [early]", "clonal [late]", "clonal [NA]", "subclonal")){
#       if (category %notin% uu$Var1){
#         uu <- rbind(uu, c(category, 0))
#       }
#     }
#     uu$Freq <- as.integer(uu$Freq)
#     global_timing_info[i,1:2] <- c(sample_id, metadata_included$cohort[i])
#     global_timing_info[i,3:7] <- c(uu$Freq[uu$Var1 == "clonal [early]"], uu$Freq[uu$Var1 == "clonal [late]"], uu$Freq[uu$Var1 == "clonal [NA]"], sum(uu$Freq[uu$Var1 != "subclonal"]), uu$Freq[uu$Var1 == "subclonal"])
#   } else {
#     global_timing_info[i,1:2] <- c(sample_id, metadata_included$cohort[i])
#     global_timing_info[i,3:7] <- rep(NA, times = 5)
#   }
# }
# 
# 
# write.table(global_timing_info, file = gzfile(paste0(wd, "r-objects/global-timing-df.txt.gz")), quote = F, row.names = F, sep = "\t")
# 
# 
# global_timing_info <- read.csv(file = paste0(wd, "r-objects/global-timing-df.txt.gz"), header = T, sep = "\t", stringsAsFactors = F)
# 
# 
# global_timing_info <- merge(global_timing_info, metadata_included[,c(3,6,8,9,14,18,19,20,21,22,23,24)], by = "sample_id")
# 
# write.table(global_timing_info, file = gzfile(paste0(wd, "r-objects/global-timing-with-metadata.txt.gz")), quote = F, row.names = F, sep = "\t")
