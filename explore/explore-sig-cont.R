

## Loading packages


if (dir.exists("/hpc/cuppen/")){
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
} else {
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(magrittr)
}

## Defining functions



`%notin%` <- Negate(`%in%`)


## Defining constants

local <- "/home/ali313/Documents/studies/master/umc-project"


if (dir.exists("/hpc/cuppen/")){
  wd <- "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/"
} else {
  wd <- paste0(local,"/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/")
}


## Loading in data

if (dir.exists("/hpc/cuppen/")){
  metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2-metadata_23072021.tsv"), sep = "\t", header = T, stringsAsFactors = F)
} else {
  metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2-metadata_23072021.tsv"), sep = "\t", header = T, stringsAsFactors = F)
}

metadata_included <- metadata[!(metadata$is_blacklisted),]

metadata_included$tmb <- rowSums(metadata_included[,22:23])




if (dir.exists("/hpc/cuppen/")){
  genes <- read.csv("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/43018_2020_50_MOESM3_ESM.txt", sep = "\t", header = T, stringsAsFactors = F)
} else {
  genes <- read.csv("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/43018_2020_50_MOESM3_ESM.txt", sep = "\t", header = T, stringsAsFactors = F)
}

# genes <- read_xlsx("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/43018_2020_50_MOESM3_ESM.xlsx", sheet = 7, skip = 1)
genes <- genes[,-5]


# Remove genes that are not annotated
pruned_genes <- genes[genes$Gene_somatic %notin% c("BABAM2", "FAAP100", "FAAP20", "FAAP24"),]

pruned_genes <- pruned_genes %>% group_by(Type, Gene_somatic) %>% summarise(Type = Type, Pathway = paste(Pathway, collapse = "|"), Pathway_abb = paste(Pathway_abb, collapse = "|")) %>% distinct() %>% relocate(Gene_somatic, .before = Type)


pruned_genes[str_detect(pruned_genes$Gene_somatic, pattern = "POL"),"Pathway_abb"] <- paste0(pull(pruned_genes[str_detect(pruned_genes$Gene_somatic, pattern = "POL"),"Pathway_abb"]), "|POL")
pruned_genes[str_detect(pruned_genes$Gene_somatic, pattern = "POL"),"Pathway"] <- paste0(pull(pruned_genes[str_detect(pruned_genes$Gene_somatic, pattern = "POL"),"Pathway"]), "|DNA Polymerases")




if (dir.exists("/hpc/cuppen/")){
  sig_cont <- readRDS("/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/contribs.fit_lsq.raw_merged.rds")
} else {
  sig_cont <- readRDS("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/contribs.fit_lsq.raw_merged.rds")
}


if (dir.exists("/hpc/cuppen/")){
  sig_cont_by_tissue <- readRDS("/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/fit_lsq.by_tissue_type/denovo_contribs.lsq.by_tissue_type.rds")
} else {
  sig_cont_by_tissue <- readRDS("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/fit_lsq.by_tissue_type/denovo_contribs.lsq.by_tissue_type.rds")
}




ddr_sig_df <- data.frame(pathway_abb = c(unique(genes$Pathway_abb), "POL"), associated_sbs_sig = c("SBS18;SBS30;SBS36",
                   NA,
                   NA,
                   "SBS3",
                   NA,
                   NA,
                   "SBS6;SBS14;SBS15;SBS20;SBS21;SBS26;SBS44",
                   NA,
                   "SBS8;SBS16",
                   NA,
                   NA,
                   NA,
                   "SBS9;SBS10a;SBS10b;SBS10c;SBS10d;SBS14;SBS20"), associated_dbs_sig = c(NA,
                                               NA,
                                               NA,
                                               NA,
                                               NA,
                                               NA,
                                               "DBS6;DBS10",
                                               NA,
                                               NA,
                                               NA,
                                               NA,
                                               NA,
                                               "DBS3"), associated_id_sig = c(NA,
                                                                           NA,
                                                                           NA,
                                                                           "ID6",
                                                                           NA,
                                                                           NA,
                                                                           "ID7",
                                                                           "ID8",
                                                                           NA,
                                                                           NA,
                                                                           NA,
                                                                           NA,
                                                                           "ID1;ID2"))



# bi-allelic annotation                        
# metadata_included_ann_bi <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/metadata-included-ddrd-annotated-biallelic.txt",
#                                      sep = "\t", header = T, stringsAsFactors = F)


if (dir.exists("/hpc/cuppen/")){
  metadata_included_ann_bi <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/metadata-included-ddrd-annotated-biallelic.txt",
                                       sep = "\t", header = T, stringsAsFactors = F)
} else {
  metadata_included_ann_bi <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/metadata-included-ddrd-annotated-biallelic.txt",
                                       sep = "\t", header = T, stringsAsFactors = F)
}

pathways <- c("No deficiency", unique(genes$Pathway_abb), "POL")

for (pathway in pathways) {
  metadata_included_ann_bi[,pathway] <- NA
  metadata_included_ann_bi[,pathway] <- str_detect(metadata_included_ann_bi$dna_repair_deficiency, pattern = pathway)
}



metadata_included_ann_bi_tibb <- tidyr::gather(metadata_included_ann_bi, colnames(metadata_included_ann_bi[,27:40]), key = "pathway_abb", value = "is_dysfunctional", )
metadata_included_ann_bi_tibb <- metadata_included_ann_bi_tibb[metadata_included_ann_bi_tibb$is_dysfunctional,]


metadata_included_ann_bi_tibb$pathway_abb <- factor(metadata_included_ann_bi_tibb$pathway_abb, levels = pathways)
metadata_included_ann_bi_tibb$cohort <- factor(metadata_included_ann_bi_tibb$cohort, levels = c("PCAWG", "HMF"))
metadata_included_ann_bi_tibb$whole_genome_duplication <- factor(metadata_included_ann_bi_tibb$whole_genome_duplication)
metadata_included_ann_bi_tibb$stage[metadata_included_ann_bi_tibb$is_metastatic] <- "Metastatic"
metadata_included_ann_bi_tibb$stage[!metadata_included_ann_bi_tibb$is_metastatic] <- "Primary"
metadata_included_ann_bi_tibb$stage <- factor(metadata_included_ann_bi_tibb$stage, levels = c("Primary", "Metastatic"))

# Remove hypermutated samples
metadata_included_ann_bi_tibb2 <- metadata_included_ann_bi_tibb[!(metadata_included_ann_bi_tibb$is_hypermutated),]




# mono-allelic annotation

# metadata_included_ann_mono <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/metadata-included-ddrd-annotated-monoallelic.txt",
#                                        sep = "\t", header = T, stringsAsFactors = F)



if (dir.exists("/hpc/cuppen/")){
  metadata_included_ann_mono <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/metadata-included-ddrd-annotated-monoallelic.txt",
                                         sep = "\t", header = T, stringsAsFactors = F)
} else {
  metadata_included_ann_mono <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/metadata-included-ddrd-annotated-monoallelic.txt",
                                         sep = "\t", header = T, stringsAsFactors = F)
}

pathways <- c("No deficiency", unique(genes$Pathway_abb), "POL")

for (pathway in pathways) {
  metadata_included_ann_mono[,pathway] <- NA
  metadata_included_ann_mono[,pathway] <- str_detect(metadata_included_ann_mono$dna_repair_deficiency, pattern = pathway)
}

metadata_included_ann_mono_tibb <- tidyr::gather(metadata_included_ann_mono, colnames(metadata_included_ann_mono[,27:40]), key = "pathway_abb", value = "is_dysfunctional", )
metadata_included_ann_mono_tibb <- metadata_included_ann_mono_tibb[metadata_included_ann_mono_tibb$is_dysfunctional,]


metadata_included_ann_mono_tibb$pathway_abb <- factor(metadata_included_ann_mono_tibb$pathway_abb, levels = pathways)
metadata_included_ann_mono_tibb$cohort <- factor(metadata_included_ann_mono_tibb$cohort, levels = c("PCAWG", "HMF"))
metadata_included_ann_mono_tibb$whole_genome_duplication <- factor(metadata_included_ann_mono_tibb$whole_genome_duplication)
metadata_included_ann_mono_tibb$stage[metadata_included_ann_mono_tibb$is_metastatic] <- "Metastatic"
metadata_included_ann_mono_tibb$stage[!metadata_included_ann_mono_tibb$is_metastatic] <- "Primary"
metadata_included_ann_mono_tibb$stage <- factor(metadata_included_ann_mono_tibb$stage, levels = c("Primary", "Metastatic"))


# Remove hypermutated samples
metadata_included_ann_mono_tibb2 <- metadata_included_ann_mono_tibb[!(metadata_included_ann_mono_tibb$is_hypermutated),]




# ========================================================================================================================================
# 
# pathway_sigs_bi <- list()
# 
# for (pathway in ddr_sig_df$pathway_abb){
#   sigs <- unlist(str_split(ddr_sig_df[ddr_sig_df$pathway_abb == pathway, 2:4], pattern = ";"))
#   sigs <- sigs[!is.na(sigs)]
#   c <- sig_cont[metadata_included_ann_bi_tibb$sample_id[metadata_included_ann_bi_tibb$pathway_abb == pathway]]
#   if (length(sigs) > 0){
#     pathway_sigs_bi[[pathway]] <- list()
#     for (i in 1:length(c)){
#       u <- F
#       for (sig in sigs){
#         u <- u | lapply(lapply(c, FUN = colnames), FUN = str_detect, pattern = paste0(sig,"\\."))[[i]]
#       }
#       pathway_sigs_bi[[pathway]][[names(c[i])]] <- c[[i]][u][as.logical(c[[i]][u] != 0)]
#       if (ncol(pathway_sigs_bi[[pathway]][[names(c[i])]]) == 0) {
#         pathway_sigs_bi[[pathway]][[names(c[i])]] <- NA
#       }
#     }
#   }
# }
# 
# 
# 
# pathway_sigs_mono <- list()
# pathway <- "POL"
# 
# for (pathway in ddr_sig_df$pathway_abb){
#   sigs <- unlist(str_split(ddr_sig_df[ddr_sig_df$pathway_abb == pathway, 2:4], pattern = ";"))
#   sigs <- sigs[!is.na(sigs)]
#   c <- sig_cont[metadata_included_ann_mono_tibb$sample_id[metadata_included_ann_mono_tibb$pathway_abb == pathway]]
#   if (length(sigs) > 0){
#     pathway_sigs_mono[[pathway]] <- list()
#     for (i in 1:length(c)){
#       u <- F
#       for (sig in sigs){
#         u <- u | lapply(lapply(c, FUN = colnames), FUN = str_detect, pattern = paste0(sig,"\\."))[[i]]
#       }
#       pathway_sigs_mono[[pathway]][[names(c[i])]] <- c[[i]][u][as.logical(c[[i]][u] != 0)]
#       if (ncol(pathway_sigs_mono[[pathway]][[names(c[i])]]) == 0) {
#         pathway_sigs_mono[[pathway]][[names(c[i])]] <- NA
#       }
#     }
#   }
# }
# 
# 
# 
# names(pathway_sigs_mono)
# 
# 
# ## BER
# sigs <- unlist(str_split(ddr_sig_df[ddr_sig_df$pathway_abb == "BER", 2:4], pattern = ";"))
# sigs <- sigs[!is.na(sigs)]
# 
# sig <- "SBS18"
# pathway_sigs_mono[["BER"]][["CPCT02290056T"]][lapply(lapply(pathway_sigs_mono[["BER"]], FUN = colnames), FUN = str_detect, pattern = paste0(sig,"\\."))[["CPCT02290056T"]]]
# 
# sbs18 <- NA
# mean(sbs18, na.rm = T)
# sum(sbs18, na.rm = T)/length(pathway_sigs_mono[["BER"]])
# summary(sbs18)
# for (sig in sigs[1]){
#   for (i in 1:length(pathway_sigs_mono[["BER"]])){
#     sbs18 <- append(sbs18, unlist(pathway_sigs_mono[["BER"]][[names(pathway_sigs_mono[["BER"]])[i]]][lapply(lapply(pathway_sigs_mono[["BER"]], FUN = colnames), FUN = str_detect, pattern = paste0(sig,"\\."))[[i]]]))
#   }
# }
# 
# head(sig_cont)
# sbs18_tot <- NA
# 
# for (sig in sigs[1]){
#   for (i in 1:length(sig_cont)){
#     sbs18_tot <- append(sbs18_tot, unlist(sig_cont[[names(sig_cont)[i]]][lapply(lapply(sig_cont, FUN = colnames), FUN = str_detect, pattern = paste0(sig,"\\."))[[i]]]))
#   }
# }
# 
# 
# 
# 
# pathway_sigs_mono[["BER"]]
# 
# 
# length(names(pathway_sigs_bi[[1]][as.logical(!is.na(pathway_sigs_bi[[1]]))])) / nrow(metadata_included_ann_bi_tibb[metadata_included_ann_bi_tibb$pathway_abb == "BER",])
# 
# length(names(pathway_sigs_mono[[1]][as.logical(!is.na(pathway_sigs_mono[[1]]))])) / nrow(metadata_included_ann_mono_tibb[metadata_included_ann_mono_tibb$pathway_abb == "BER",])
# 
# 
# 
# 
# 
# length(names(pathway_sigs_bi[[2]][as.logical(!is.na(pathway_sigs_bi[[2]]))])) / nrow(metadata_included_ann_bi_tibb[metadata_included_ann_bi_tibb$pathway_abb == "HDR",])
# 
# 
# length(names(pathway_sigs_mono[[2]][as.logical(!is.na(pathway_sigs_mono[[2]]))])) / nrow(metadata_included_ann_mono_tibb[metadata_included_ann_mono_tibb$pathway_abb == "HDR",])
# 
# 
# 
# 
# length(names(pathway_sigs_bi[[3]][as.logical(!is.na(pathway_sigs_bi[[3]]))])) / nrow(metadata_included_ann_bi_tibb[metadata_included_ann_bi_tibb$pathway_abb == "MMR",])
# 
# 
# length(names(pathway_sigs_mono[[3]][as.logical(!is.na(pathway_sigs_mono[[3]]))])) / nrow(metadata_included_ann_mono_tibb[metadata_included_ann_mono_tibb$pathway_abb == "MMR",])
# 
# 
# 
# 
# 
# length(names(pathway_sigs_bi[[4]][as.logical(!is.na(pathway_sigs_bi[[4]]))])) / nrow(metadata_included_ann_bi_tibb[metadata_included_ann_bi_tibb$pathway_abb == "NHEJ",])
# 
# 
# length(names(pathway_sigs_mono[[4]][as.logical(!is.na(pathway_sigs_mono[[4]]))])) / nrow(metadata_included_ann_mono_tibb[metadata_included_ann_mono_tibb$pathway_abb == "NHEJ",])
# 
# 
# 
# 
# 
# length(names(pathway_sigs_bi[[5]][as.logical(!is.na(pathway_sigs_bi[[5]]))])) / nrow(metadata_included_ann_bi_tibb[metadata_included_ann_bi_tibb$pathway_abb == "NER",])
# 
# 
# length(names(pathway_sigs_mono[[5]][as.logical(!is.na(pathway_sigs_mono[[5]]))])) / nrow(metadata_included_ann_mono_tibb[metadata_included_ann_mono_tibb$pathway_abb == "NER",])
# 
# 
# 
# 
# 
# length(names(pathway_sigs_bi[[6]][as.logical(!is.na(pathway_sigs_bi[[6]]))])) / nrow(metadata_included_ann_bi_tibb[metadata_included_ann_bi_tibb$pathway_abb == "POL",])
# 
# 
# length(names(pathway_sigs_mono[[6]][as.logical(!is.na(pathway_sigs_mono[[6]]))])) / nrow(metadata_included_ann_mono_tibb[metadata_included_ann_mono_tibb$pathway_abb == "POL",])
# 



# ========================================================================================================================================


# sig_cont_by_tissue <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/fit_lsq.by_tissue_type/denovo_contribs.lsq.by_tissue_type.rds")

if (dir.exists("/hpc/cuppen/")){
  sig_cont_by_tissue <- readRDS("/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/fit_lsq.by_tissue_type/denovo_contribs.lsq.by_tissue_type.rds")
} else {
  sig_cont_by_tissue <- readRDS("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/fit_lsq.by_tissue_type/denovo_contribs.lsq.by_tissue_type.rds")
}

# colnames(sig_cont_by_tissue[["Lymphoid"]])
# sig_cont_by_tissue[["Lymphoid"]] <- sig_cont_by_tissue[["Lymphoid"]][,-12]
colnames(sig_cont_by_tissue[["Lymphoid"]])[12] <- paste0(colnames(sig_cont_by_tissue[["Lymphoid"]])[12], ".NA.NA.NA")

# colnames(sig_cont_by_tissue[["Pancreas"]])
colnames(sig_cont_by_tissue[["Pancreas"]])[10] <- paste0(colnames(sig_cont_by_tissue[["Pancreas"]])[10], ".NA.NA.NA")
# colnames(sig_cont_by_tissue[["Pancreas"]]["CPCT02210054T",])

# colnames(sig_cont_by_tissue[["Pancreas_NET"]])
colnames(sig_cont_by_tissue[["Pancreas_NET"]])[8] <- paste0(colnames(sig_cont_by_tissue[["Pancreas_NET"]])[8], ".NA.NA.NA")
# sig_cont_by_tissue[["Pancreas_NET"]] <- sig_cont_by_tissue[["Pancreas_NET"]][,-8]
# colnames(sig_cont_by_tissue[["Pancreas_NET"]]["CPCT02010688T",])


# rownames(sig_cont_by_tissue[["Pancreas"]]) %in% rownames(sig_cont_by_tissue[["Pancreas_NET"]])




merged_cont <- data.frame(Sample_id = NA, Tissue = NA, Sum = NA, Signature = NA, Contribution = NA)
merged_cont <- merged_cont[-1,]

for (name in names(sig_cont_by_tissue)){
  # tissue_abb <- substr(name, start = 1, stop = 3)
  colnames(sig_cont_by_tissue[[name]]) <- paste0(name, "_", colnames(sig_cont_by_tissue[[name]]))
  
  mmm <- data.frame(sig_cont_by_tissue[[name]])
  
  mmm <- mmm %>% mutate(Sum = rowSums(mmm, na.rm = T))
  
  mmm <- mmm %>% mutate(Sample_id = rownames(mmm), Tissue = name) %>% relocate(where(is.character), .before = where(is.numeric)) %>% gather(colnames(mmm)[1:(ncol(mmm)-1)], key = "Signature", value = "Contribution")
  
  merged_cont <- rbind(merged_cont, mmm)
}


merged_cont <- merged_cont %>% mutate(Cosmic_signature = unlist(lapply(strsplit(merged_cont$Signature, "\\."), FUN = "[[", 2)), Cosine_similarity = unlist(lapply(strsplit(merged_cont$Signature, "\\."), FUN = "[[", 4)))
merged_cont$Cosine_similarity <- as.numeric(merged_cont$Cosine_similarity)

merged_cont <- merged_cont[merged_cont$Sample_id %in% metadata_included$sample_id,]

# This is to sum de-novo signatures that are attributed to the same cosmic signature within the same cohort

merged_cont3 <- aggregate(merged_cont['Contribution'], by=merged_cont[c('Sample_id', "Sum", "Cosmic_signature")], sum)
merged_cont4 <- aggregate(merged_cont['Cosine_similarity'], by=merged_cont[c('Sample_id', "Sum", "Cosmic_signature")], FUN = mean)
# all(merged_cont4$Sample_id == merged_cont3$Sample_id)
merged_cont3$Cosine_similarity <- merged_cont4$Cosine_similarity


# merged_cont3[merged_cont3$Sample_id == "WIDE01010227T",]
# merged_cont4[merged_cont4$Sample_id == "WIDE01010227T",]
# merged_cont3[merged_cont3$Sample_id == "WIDE01010227T",]


merged_cont <- merged_cont3




# ========================================================================================================================================
# ========================================================================================================================================
# merged_cont
# head(merged_cont)
# str(merged_cont)
# head(merged_cont2)
merged_cont2 <- NA
merged_cont2 <- merged_cont[merged_cont$Cosine_similarity > 90,]
# nrow(merged_cont2[!is.na(merged_cont2$Sample_id),])
merged_cont2 <- merged_cont2[!is.na(merged_cont2$Sample_id),]
# nrow(merged_cont2)
# ddr_sig_df
# 
# mean(merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_bi$sample_id[metadata_included_ann_bi$BER] & merged_cont$Cosmic_signature =="SBS36"])
# length(merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_bi$sample_id[metadata_included_ann_bi$BER] & merged_cont$Cosmic_signature =="SBS36"])
# mean(merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_bi$sample_id[!metadata_included_ann_bi$BER] & merged_cont$Cosmic_signature =="SBS36"])
# length(merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_bi$sample_id[!metadata_included_ann_bi$BER] & merged_cont$Cosmic_signature =="SBS36"])
# boxplot(merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_mono$sample_id[metadata_included_ann_mono$BER] & merged_cont$Cosmic_signature =="SBS36"], merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_bi$sample_id[metadata_included_ann_bi$BER] & merged_cont$Cosmic_signature =="SBS36"], merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_bi$sample_id[!metadata_included_ann_bi$BER] & merged_cont$Cosmic_signature =="SBS36"])
# 
# 
# mean(merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_mono$sample_id[metadata_included_ann_mono$BER] & merged_cont$Cosmic_signature =="SBS36"])
# mean(merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_mono$sample_id[!metadata_included_ann_mono$BER] & merged_cont$Cosmic_signature =="SBS36"])
# 
# 
# 
# mean(merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_bi$sample_id[metadata_included_ann_bi$HDR] & merged_cont$Cosmic_signature =="ID6"])
# mean(merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_bi$sample_id[!metadata_included_ann_bi$HDR] & merged_cont$Cosmic_signature =="ID6"])
# 
# mean(merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_mono$sample_id[metadata_included_ann_mono$HDR] & merged_cont$Cosmic_signature =="ID6"])
# mean(merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_mono$sample_id[!metadata_included_ann_mono$HDR] & merged_cont$Cosmic_signature =="ID6"])
# 
# 
# mean(merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_bi$sample_id[metadata_included_ann_bi$MMR] & merged_cont$Cosmic_signature =="SBS44"])
# mean(merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_bi$sample_id[!metadata_included_ann_bi$MMR] & merged_cont$Cosmic_signature =="SBS44"])
# 
# mean(merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_mono$sample_id[metadata_included_ann_mono$MMR] & merged_cont$Cosmic_signature =="SBS44"])
# mean(merged_cont$Contribution[merged_cont$Sample_id %in% metadata_included_ann_mono$sample_id[!metadata_included_ann_mono$MMR] & merged_cont$Cosmic_signature =="SBS44"])
# 



# plot_df2 <- data.frame(Sample_id = character(60000), Pathway = character(60000), Cosmic_signature = character(60000), Annotation = character(60000), 
#                       Absolute_contribution = numeric(60000), Relative_contribution = numeric(60000), 
#                       Number_of_samples_with_the_signature = numeric(60000), Number_of_samples_deficient = numeric(60000))
# 
# i <- 0
# for (pathway in ddr_sig_df$pathway_abb[!is.na(ddr_sig_df$associated_sbs_sig) | !is.na(ddr_sig_df$associated_sbs_sig) | !is.na(ddr_sig_df$associated_sbs_sig)]) {
#   sigs <- unlist(str_split(ddr_sig_df[ddr_sig_df$pathway_abb == pathway, 2:4], pattern = ";"))
#   sigs <- sigs[!is.na(sigs)]
#   for (sig in sigs) {
#     for (ann in c("Bi", "Mono", "None")){
#       if (ann == "Bi") {
#         for (sample in metadata_included_ann_bi$sample_id[metadata_included_ann_bi[,pathway]]){
#           if (length(merged_cont2$Contribution[merged_cont2$Sample_id == sample & merged_cont2$Cosmic_signature == sig]) > 0){
#             i <- i + 1
#             print(i)
#             plot_df2$Sample_id[i] <- sample
#             plot_df2$Pathway[i] <- pathway
#             plot_df2$Cosmic_signature[i] <- sig
#             plot_df2$Annotation[i] <- ann
#             plot_df2$Absolute_contribution[i] <- merged_cont2$Contribution[merged_cont2$Sample_id == sample & merged_cont2$Cosmic_signature == sig]
#             plot_df2$Relative_contribution[i] <- merged_cont2$Contribution[merged_cont2$Sample_id == sample & merged_cont2$Cosmic_signature == sig]/merged_cont2$Sum[merged_cont2$Sample_id == sample & merged_cont2$Cosmic_signature == sig]
#             plot_df2$Number_of_samples_with_the_signature[i] <- length(merged_cont2$Contribution[merged_cont2$Sample_id %in% metadata_included_ann_bi$sample_id[metadata_included_ann_bi[,pathway]] & merged_cont2$Cosmic_signature == sig])
#             plot_df2$Number_of_samples_deficient[i] <- length(metadata_included_ann_bi$sample_id[metadata_included_ann_bi[,pathway]])
#           }
#         }
#       } else if (ann == "Mono"){
#         for (sample in metadata_included_ann_mono$sample_id[metadata_included_ann_mono[,pathway]]){
#           if (sample %in% merged_cont2$Sample_id & length(merged_cont2$Contribution[merged_cont2$Sample_id == sample & merged_cont2$Cosmic_signature == sig]) > 0){
#             i <- i + 1
#             print(i)
#             plot_df2$Sample_id[i] <- sample
#             plot_df2$Pathway[i] <- pathway
#             plot_df2$Cosmic_signature[i] <- sig
#             plot_df2$Annotation[i] <- ann
#             plot_df2$Absolute_contribution[i] <- merged_cont2$Contribution[merged_cont2$Sample_id == sample & merged_cont2$Cosmic_signature == sig]
#             plot_df2$Relative_contribution[i] <- merged_cont2$Contribution[merged_cont2$Sample_id == sample & merged_cont2$Cosmic_signature == sig]/merged_cont2$Sum[merged_cont2$Sample_id == sample & merged_cont2$Cosmic_signature == sig]
#             plot_df2$Number_of_samples_with_the_signature[i] <- length(merged_cont2$Contribution[merged_cont2$Sample_id %in% metadata_included_ann_mono$sample_id[metadata_included_ann_mono[,pathway]] & merged_cont2$Cosmic_signature == sig])
#             plot_df2$Number_of_samples_deficient[i] <- length(metadata_included_ann_mono$sample_id[metadata_included_ann_mono[,pathway]])
#           }
#         }
#       } else if (ann == "None") {
#         for (sample in metadata_included_ann_mono$sample_id[!metadata_included_ann_mono[,pathway]]){
#           if (sample %in% merged_cont2$Sample_id & length(merged_cont2$Contribution[merged_cont2$Sample_id == sample & merged_cont2$Cosmic_signature == sig]) > 0){
#             i <- i + 1
#             print(i)
#             plot_df2$Sample_id[i] <- sample
#             plot_df2$Pathway[i] <- pathway
#             plot_df2$Cosmic_signature[i] <- sig
#             plot_df2$Annotation[i] <- ann
#             plot_df2$Absolute_contribution[i] <- merged_cont2$Contribution[merged_cont2$Sample_id == sample & merged_cont2$Cosmic_signature == sig]
#             plot_df2$Relative_contribution[i] <- merged_cont2$Contribution[merged_cont2$Sample_id == sample & merged_cont2$Cosmic_signature == sig]/merged_cont2$Sum[merged_cont2$Sample_id == sample & merged_cont2$Cosmic_signature == sig]
#             plot_df2$Number_of_samples_with_the_signature[i] <- length(merged_cont2$Contribution[merged_cont2$Sample_id %in% metadata_included_ann_mono$sample_id[!metadata_included_ann_mono[,pathway]] & merged_cont2$Cosmic_signature == sig])
#             plot_df2$Number_of_samples_deficient[i] <- length(metadata_included_ann_mono$sample_id[!metadata_included_ann_mono[,pathway]])
#           }
#         }
#       }
#     }
#   }
# }
# 
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(plot_df, file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/sig-all-cosine.rds")
# } else {
#   saveRDS(plot_df, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/sig-all-cosine.rds")
# }


# Plot_df is when we consider all the signatures regardless of their cosine similarity

if (dir.exists("/hpc/cuppen/")){
  plot_df <- readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/sig-all-cosine.rds")
} else {
  plot_df <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/sig-all-cosine.rds")
}
plot_df <- plot_df[-(47707:nrow(plot_df)),]

plot_df$Pathway <- factor(plot_df$Pathway)
plot_df$Cosmic_signature <- factor(plot_df$Cosmic_signature)
plot_df$Annotation <- factor(plot_df$Annotation)


# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(plot_df2, file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/sig-above-90-cosine.rds")
# } else {
#   saveRDS(plot_df2, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/sig-above-90-cosine.rds")
# }


# Plot_df2 is when we consider only the signatures  with cosine similarity of 90 percent or higher

if (dir.exists("/hpc/cuppen/")){
  plot_df2 <- readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/sig-above-90-cosine.rds")
} else {
  plot_df2 <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/sig-above-90-cosine.rds")
}

# nrow(plot_df2[plot_df2$Sample_id != "",])
plot_df2 <- plot_df2[-(26834:nrow(plot_df2)),]

plot_df2$Pathway <- factor(plot_df2$Pathway)
plot_df2$Cosmic_signature <- factor(plot_df2$Cosmic_signature)
plot_df2$Annotation <- factor(plot_df2$Annotation)
# 
# 
# 
# 
# 
# pathway <- "POL"
# sig <- "ID1"
# nrow(merged_cont2[merged_cont2$Sample_id %in% metadata_included_ann_mono$sample_id[!metadata_included_ann_mono[,pathway]] & merged_cont2$Cosmic_signature == sig,])
# 
# nrow(merged_cont2)
# nrow(merged_cont2[merged_cont2$Sample_id %in% metadata_included_ann_mono$sample_id[!metadata_included_ann_mono[,pathway]]& merged_cont2$Cosmic_signature == sig,])
# merged_cont2[merged_cont2$Sample_id %in% metadata_included_ann_mono$sample_id[!metadata_included_ann_mono[,pathway]]& merged_cont2$Cosmic_signature == sig,"Sample_id"][duplicated(merged_cont2[merged_cont2$Sample_id %in% metadata_included_ann_mono$sample_id[!metadata_included_ann_mono[,pathway]]& merged_cont2$Cosmic_signature == sig,"Sample_id"])]
# nrow(metadata_included_ann_mono[!(metadata_included_ann_mono[,pathway]),])
# 
# 
# merged_cont2[merged_cont2$Sample_id == "WIDE01010227T",]
# 
# x <- merged_cont[,c(4,6)]
# 
# unique(merged_cont[merged_cont$Cosmic_signature == "ID1","Signature"])


# # reating geom_text data
number_labels_plot_df <- unique(plot_df[,c(2:4,7:8)])
number_labels_plot_df2 <- unique(plot_df2[,c(2:4,7:8)])



ddr_sig_plot_abs <- plot_df %>% ggplot(aes(x = Cosmic_signature, y = log2(Absolute_contribution), color = Annotation)) +
  facet_wrap(~ Pathway, scales = "free") +
  geom_boxplot(outlier.shape = NA) +
  geom_text(data = number_labels_plot_df, aes(x = Cosmic_signature, y = -5, label = paste0(Number_of_samples_with_the_signature, "/",Number_of_samples_deficient), color = Annotation), size = 3, position = position_dodge(width = 0.9), angle = 90)



for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig15-ddrd-sig-rel-contribution-boxplot-with-hypermutated-abs.png", width = 960, height = 960)
    print(ddr_sig_plot_abs)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig15-ddrd-sig-rel-contribution-boxplot-with-hypermutated-abs.pdf", width = 14, height = 14)
    print(ddr_sig_plot_abs)
    dev.off()
  }
}



ddr_sig_plot_abs_90 <- plot_df2 %>% ggplot(aes(x = Cosmic_signature, y = log2(Absolute_contribution), color = Annotation)) + facet_wrap(~ Pathway, scales = "free") +
  geom_boxplot(outlier.shape = NA) +
  geom_text(data = number_labels_plot_df2, aes(x = Cosmic_signature, y = -5, label = paste0(Number_of_samples_with_the_signature, "/",Number_of_samples_deficient), color = Annotation), size = 3, position = position_dodge(width = 0.9), angle = 90)



for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig15-ddrd-sig-rel-contribution-boxplot-with-hypermutated-abs-90.png", width = 960, height = 960)
    print(ddr_sig_plot_abs_90)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig15-ddrd-sig-rel-contribution-boxplot-with-hypermutated-abs-90.pdf", width = 14, height = 14)
    print(ddr_sig_plot_abs_90)
    dev.off()
  }
}



ddr_sig_plot_rel <- plot_df %>% ggplot(aes(x = Cosmic_signature, y = Relative_contribution*100, color = Annotation)) + facet_wrap(~ Pathway, scales = "free") +
  geom_boxplot(outlier.shape = NA) +
  geom_text(data = number_labels_plot_df, aes(x = Cosmic_signature, y = -5, label = paste0(Number_of_samples_with_the_signature, "/",Number_of_samples_deficient), color = Annotation), size = 3, position = position_dodge(width = 0.9), angle = 90)



for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig15-ddrd-sig-rel-contribution-boxplot-with-hypermutated-rel.png", width = 960, height = 960)
    print(ddr_sig_plot_rel)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig15-ddrd-sig-rel-contribution-boxplot-with-hypermutated-rel.pdf", width = 14, height = 14)
    print(ddr_sig_plot_rel)
    dev.off()
  }
}



ddr_sig_plot_rel_90 <- plot_df2 %>% ggplot(aes(x = Cosmic_signature, y = Relative_contribution*100, color = Annotation)) + facet_wrap(~ Pathway, scales = "free") +
  geom_boxplot(outlier.shape = NA) +
  geom_text(data = number_labels_plot_df2, aes(x = Cosmic_signature, y = -5, label = paste0(Number_of_samples_with_the_signature, "/",Number_of_samples_deficient), color = Annotation), size = 3, position = position_dodge(width = 0.9), angle = 90)



for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig15-ddrd-sig-rel-contribution-boxplot-with-hypermutated-rel-90.png", width = 960, height = 960)
    print(ddr_sig_plot_rel_90)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig15-ddrd-sig-rel-contribution-boxplot-with-hypermutated-rel-90.pdf", width = 14, height = 14)
    print(ddr_sig_plot_rel_90)
    dev.off()
  }
}
# 
# 
# refit_sig_cont <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/annotation-beds/refit-signatures-contribution.txt.gz", sep = "\t", header = T, stringsAsFactors = F)
# nrow(refit_sig_cont)
# nrow(refit_sig_cont)
# refit_sig_cont <- refit_sig_cont(row)





# ========================================================================================================================================
# Here I make the same plots but this time primary and metastatic are plotted separately


head(plot_df2)
plot_df2$cohort <- NA

for (i in 1:nrow(plot_df2)) {
  plot_df2[i,"cohort"] <- metadata_included$is_metastatic[metadata_included$sample_id == plot_df2$Sample_id[i]]
}


plot_df2[which(plot_df2$cohort),"cohort"] <- "Metastatic"
plot_df2[plot_df2$cohort == "FALSE","cohort"] <- "Primary"




ddr_sig_plot_rel_90_prim <- plot_df2[plot_df2$cohort == "Primary",] %>% ggplot(aes(x = Cosmic_signature, y = Relative_contribution*100, color = Annotation)) + facet_wrap(~ Pathway, scales = "free") +
  geom_boxplot(outlier.shape = NA) 
  # geom_text(data = number_labels_plot_df2_prim, aes(x = Cosmic_signature, y = -5, label = paste0(Number_of_samples_with_the_signature, "/",Number_of_samples_deficient), color = Annotation), size = 3, position = position_dodge(width = 0.9), angle = 90)


ddr_sig_plot_rel_90_metas <- plot_df2[plot_df2$cohort == "Metastatic",] %>% ggplot(aes(x = Cosmic_signature, y = Relative_contribution*100, color = Annotation)) + facet_wrap(~ Pathway, scales = "free") +
  geom_boxplot(outlier.shape = NA) 
  # geom_text(data = number_labels_plot_df2_metas, aes(x = Cosmic_signature, y = -5, label = paste0(Number_of_samples_with_the_signature, "/",Number_of_samples_deficient), color = Annotation), size = 3, position = position_dodge(width = 0.9), angle = 90)



for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig15-ddrd-sig-rel-contribution-boxplot-with-hypermutated-rel-90-metas.png", width = 960, height = 960)
    print(ddr_sig_plot_rel_90_metas)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig15-ddrd-sig-rel-contribution-boxplot-with-hypermutated-rel-90-metas.pdf", width = 14, height = 14)
    print(ddr_sig_plot_rel_90_metas)
    dev.off()
  }
}







# ========================================================================================================================================
# ========================================================================================================================================
# Validation


diplotype_final <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/all-diplotypes-metadata-merged.txt.gz",
                            sep = "\t", header = T, stringsAsFactors = F)


# Internship R-objects
processed_refit_sig_cont <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-contribution-refit/processed-refit-sig-contribution.rds")

diplotypes <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/annotation-beds/diplotypes_germPonFilt.txt.gz",
                       sep = "\t", header = T, stringsAsFactors = F)




nrow(metadata_included[metadata_included$cohort == "HMF" & metadata_included$sample_id %in% metadata_included_ann_bi[metadata_included_ann_bi$BER,"sample_id"],])
metadata_included_ann_bi[metadata_included_ann_bi$BER,"sample_id"]


head(metadata)
nrow(metadata)
prev_extra <- metadata$sampleId[metadata$pathway == "BER"][metadata$sampleId[metadata$pathway == "BER"] %in% metadata_included[metadata_included$cohort == "HMF" & metadata_included$sample_id %notin% metadata_included_ann_bi[metadata_included_ann_bi$BER,"sample_id"],"sample_id"]]
now_extra <- metadata_included[metadata_included$cohort == "HMF" & metadata_included$sample_id %in% metadata_included_ann_bi[metadata_included_ann_bi$BER,"sample_id"],"sample_id"][metadata_included[metadata_included$cohort == "HMF" & metadata_included$sample_id %in% metadata_included_ann_bi[metadata_included_ann_bi$BER,"sample_id"],"sample_id"] %notin% metadata$sampleId[metadata$pathway == "BER"]]

metadata_included[metadata_included$cohort == "HMF" & metadata_included$sample_id %in% metadata_included_ann_bi[metadata_included_ann_bi$BER,"sample_id"],"sample_id"]



head(diplotype_final)
diplotype_final[diplotype_final$sample == "CPCT02010469T" & diplotype_final$hgnc_symbol == "MUTYH",]
diplotypes[diplotypes$sample == "CPCT02010469T",]


diplotype_final[diplotype_final$sample == "CPCT02020664T" & diplotype_final$hgnc_symbol == "MBD4",]
diplotypes[diplotypes$sample == "CPCT02020664T",]


head(processed_refit_sig_cont)
mean(processed_refit_sig_cont$SBS36[processed_refit_sig_cont$sampleId %in% prev_extra])
mean(processed_refit_sig_cont$SBS36[processed_refit_sig_cont$sampleId %in% now_extra])


mean(processed_refit_sig_cont$SBS36[processed_refit_sig_cont$sampleId %in% metadata$sampleId[metadata$pathway == "BER"]])
mean(processed_refit_sig_cont$SBS18[processed_refit_sig_cont$sampleId %in% metadata_included[metadata_included$cohort == "HMF" & metadata_included$sample_id %in% metadata_included_ann_bi[!metadata_included_ann_bi$BER,"sample_id"],"sample_id"]])

mean(processed_refit_sig_cont$SBS44[processed_refit_sig_cont$sampleId %in% metadata$sampleId[metadata$pathway == "MMR"]])
mean(processed_refit_sig_cont$SBS44[processed_refit_sig_cont$sampleId %in% metadata_included[metadata_included$cohort == "HMF" & metadata_included$sample_id %in% metadata_included_ann_bi[!(metadata_included_ann_bi$MMR),"sample_id"],"sample_id"]])


head(merged_cont)
mean(merged_cont[merged_cont$Sample_id %in% metadata$sampleId[metadata$pathway == "MMR"] & merged_cont$Cosmic_signature == "SBS44","Contribution"]/merged_cont[merged_cont$Sample_id %in% metadata$sampleId[metadata$pathway == "MMR"] & merged_cont$Cosmic_signature == "SBS44","Sum"])
mean(merged_cont[merged_cont$Sample_id %in% metadata_included[metadata_included$cohort == "HMF" & metadata_included$sample_id %in% metadata_included_ann_bi[(metadata_included_ann_bi$MMR),"sample_id"],"sample_id"] & merged_cont$Cosmic_signature == "SBS44","Contribution"]/merged_cont[merged_cont$Sample_id %in% metadata_included[metadata_included$cohort == "HMF" & metadata_included$sample_id %in% metadata_included_ann_bi[(metadata_included_ann_bi$MMR),"sample_id"],"sample_id"] & merged_cont$Cosmic_signature == "SBS44","Sum"])




mean(merged_cont[merged_cont$Sample_id %in% now_extra & merged_cont$Cosmic_signature == "SBS18","Contribution"])
