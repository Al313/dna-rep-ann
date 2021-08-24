

## Loading packages
# library(rcompanion)
# library(readxl)
library(dplyr)
library(tibble)
library(stringr)
# library(magrittr)
# library(ggplot2)
# library(gridExtra)
# library(ggpubr)
# library(grid)
library(tibble)
library(tidyr)
# library(epitools)


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

pruned_genes %<>% group_by(Type, Gene_somatic) %>% summarise(Type = Type, Pathway = paste(Pathway, collapse = "|"), Pathway_abb = paste(Pathway_abb, collapse = "|")) %>% distinct() %>% relocate(Gene_somatic, .before = Type)


pruned_genes[str_detect(pruned_genes$Gene_somatic, pattern = "POL"),"Pathway_abb"] <- paste0(pull(pruned_genes[str_detect(pruned_genes$Gene_somatic, pattern = "POL"),"Pathway_abb"]), "|POL")
pruned_genes[str_detect(pruned_genes$Gene_somatic, pattern = "POL"),"Pathway"] <- paste0(pull(pruned_genes[str_detect(pruned_genes$Gene_somatic, pattern = "POL"),"Pathway"]), "|DNA Polymerases")

# tibble::view(pruned_genes[str_detect(pruned_genes$Pathway_abb, pattern =  "HDR"),])
# pruned_genes$Gene_somatic[str_detect(pruned_genes$Pathway_abb, pattern =  "HDR") & pruned_genes$Type == "core"]

# =====================================================================================================================================


### Merging diplotype info, pathway info, and clinical metadata: 
# !!!(Updated by cancer_types_HMF_PCAWG_2-metadata_23072021.tsv)

# ## This file was obtained from this article: https://doi.org/10.1038/s43018-020-0050-6
# 
# 
# 
# pruned_core_genes <- pruned_genes[pruned_genes$Type == "core",]
# pruned_accessory_genes <- pruned_genes[pruned_genes$Type == "accessory",]
# pruned_accessory_negreg_genes <- pruned_genes[pruned_genes$Type == "accessory (negative regulator)",]
# 
# 
# 
# 
# diplotypes <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/all_diplotypes.txt.gz",
#                        sep = "\t", header = T, stringsAsFactors = F)
# 
# 
# diplotypes2 <- merge(diplotypes, pruned_core_genes[,c(1,4)], by.x = "hgnc_symbol", by.y = "Gene_somatic", sort = F, all.x = T)
# colnames(diplotypes2)[ncol(diplotypes2)] <- "core_pathway"
# 
# diplotypes3 <- merge(diplotypes2, pruned_accessory_genes[,c(1,4)], by.x = "hgnc_symbol", by.y = "Gene_somatic", sort = F, all.x = T)
# colnames(diplotypes3)[ncol(diplotypes3)] <- "accessory_pathway"
# 
# diplotypes4 <- merge(diplotypes3, pruned_accessory_negreg_genes[,c(1,4)], by.x = "hgnc_symbol", by.y = "Gene_somatic", sort = F, all.x = T)
# colnames(diplotypes4)[ncol(diplotypes4)] <- "accessory_neg_reg_pathway"
# 
# 
# 
# diplotype_final <- merge(diplotypes4, metadata_included[,c(1,3,6:10,13:24)], by.x = "sample", by.y = "sample_id")
# 
# write.table(diplotype_final, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/all-diplotypes-metadata-merged.txt",
#             sep = "\t", quote = F, row.names = F)





if (dir.exists("/hpc/cuppen/")){
  diplotype_final <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/all-diplotypes-metadata-merged.txt.gz",
                              sep = "\t", header = T, stringsAsFactors = F)
  
} else {
  diplotype_final <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/all-diplotypes-metadata-merged.txt.gz",
                              sep = "\t", header = T, stringsAsFactors = F)
  
}

# =====================================================================================================================================

### Merging sample level ddrd metadata (used by cancer_types_HMF_PCAWG_2-metadata_23072021.tsv) _ biallelic hit:


metadata_included$dna_repair_deficiency <- NA
metadata_included$dna_repair_genes <- NA
i <- 0


# !!! A: In this for loop I used my own sampling annotation criteria which would find the pathway with absolute max count of mutated core genes.

for (sample_id in unique(diplotype_final$sample)) {
  i <- i + 1
  print(i)
  bam <- diplotype_final$core_pathway[diplotype_final$sample == sample_id & (diplotype_final$a1.max_score >= 4 | diplotype_final$a2.max_score >= 4) & diplotype_final$biall_status != "loh_only" & !is.na(diplotype_final$core_pathway)]
  # bam[!is.na(bam)]
  decoy <- paste(bam[!is.na(bam)], collapse = "_")


  if (nchar(decoy) < 1) {
    nr_max = 0
  } else {

    decoy_table <- table(unlist(str_split(decoy, pattern = "_")[[1]] %>% str_split(pattern = "\\|")))

    if (length(decoy_table) > 1 & "CCR" %in% names(decoy_table)){
      decoy_table["CCR"] <- 0
    }


    max_repeat <- max(decoy_table)
    nr_max <- sum(as.vector(decoy_table) == max_repeat)
  }

  sam <- diplotype_final$hgnc_symbol[diplotype_final$sample == sample_id & diplotype_final$a1.max_score >= 3 & diplotype_final$a2.max_score >= 3 & !is.na(diplotype_final$core_pathway)]

  if (nr_max == 1) {
    consensus <- names(which.max(decoy_table))
    # print(consensus)

    repair_genes <- paste(sam[str_detect(consensus, pruned_genes$Pathway_abb[pruned_genes$Gene_somatic %in% sam & pruned_genes$Type == "core"])], collapse = "_")
    # print(repair_genes)
  } else if (nr_max > 1) {
    # print("Too complicated to be considered!")
    consensus <- "No consensus"
    repair_genes <- NA
  } else if (nr_max == 0) {
    # print("No DNA repair deficiency")
    consensus <- "No deficiency"
    repair_genes <- NA
  }

  metadata_included$dna_repair_deficiency[metadata_included$sample_id == sample_id] <- consensus
  metadata_included$dna_repair_genes[metadata_included$sample_id == sample_id] <- repair_genes
}

write.table(metadata_included, file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/metadata-included-ddrd-annotated-monoallelic-parsimony.txt",
            sep = "\t", quote = F, row.names = F)

metadata_included_ann_mono_parsimony <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/metadata-included-ddrd-annotated-monoallelic-parsimony.txt",
                                     sep = "\t", header = T, stringsAsFactors = F)


pathways <- c("No deficiency", unique(genes$Pathway_abb), "POL")

for (pathway in pathways) {
  metadata_included_ann_mono_parsimony[,pathway] <- NA
  metadata_included_ann_mono_parsimony[,pathway] <- str_detect(metadata_included_ann_mono_parsimony$dna_repair_deficiency, pattern = pathway)
}


metadata_included_ann_mono_parsimony_tibb <- tidyr::gather(metadata_included_ann_mono_parsimony, colnames(metadata_included_ann_mono_parsimony[,27:40]), key = "pathway_abb", value = "is_dysfunctional", )
metadata_included_ann_mono_parsimony_tibb <- metadata_included_ann_mono_parsimony_tibb[metadata_included_ann_mono_parsimony_tibb$is_dysfunctional,]


metadata_included_ann_mono_parsimony_tibb$pathway_abb <- factor(metadata_included_ann_mono_parsimony_tibb$pathway_abb, levels = pathways)
metadata_included_ann_mono_parsimony_tibb$cohort <- factor(metadata_included_ann_mono_parsimony_tibb$cohort, levels = c("PCAWG", "HMF"))
metadata_included_ann_mono_parsimony_tibb$whole_genome_duplication <- factor(metadata_included_ann_mono_parsimony_tibb$whole_genome_duplication)
metadata_included_ann_mono_parsimony_tibb$stage[metadata_included_ann_mono_parsimony_tibb$is_metastatic] <- "Metastatic"
metadata_included_ann_mono_parsimony_tibb$stage[!metadata_included_ann_mono_parsimony_tibb$is_metastatic] <- "Primary"
metadata_included_ann_mono_parsimony_tibb$stage <- factor(metadata_included_ann_mono_parsimony_tibb$stage, levels = c("Primary", "Metastatic"))


# Arne's questions
table(metadata_included_ann_mono_parsimony_tibb$pathway_abb[metadata_included_ann_mono_parsimony_tibb$is_metastatic & !(metadata_included_ann_mono_parsimony_tibb$is_hypermutated)])
nrow(metadata_included_ann_mono_parsimony_tibb[metadata_included_ann_mono_parsimony_tibb$is_metastatic & !(metadata_included_ann_mono_parsimony_tibb$is_hypermutated),])
table(metadata_included_ann_mono_parsimony$hr_status[metadata_included_ann_mono_parsimony$HDR  & !(metadata_included_ann_mono_parsimony$is_hypermutated)])
table(metadata_included_ann_mono_parsimony$hr_status[!(metadata_included_ann_mono_parsimony$HDR) & !(metadata_included_ann_mono_parsimony$is_hypermutated)])
nrow(metadata_included_ann_mono_parsimony[metadata_included_ann_mono_parsimony$is_metastatic & !(metadata_included_ann_mono_parsimony$is_hypermutated),])



# Remove hypermutated samples
metadata_included_ann_mono_parsimony_tibb2 <- metadata_included_ann_mono_parsimony_tibb[!(metadata_included_ann_mono_parsimony_tibb$is_hypermutated),]




# =====================================================================================================================================

# !!! A: In this for loop I did not through away any information even for samples in which two or more pathways have equal number of mutated core genes.


# for (sample_id in unique(diplotype_final$sample)) {
#   i <- i + 1
#   print(i)
#   bam <- diplotype_final$core_pathway[diplotype_final$sample == sample_id & diplotype_final$a1.max_score >= 3 & diplotype_final$a2.max_score >= 3 & !is.na(diplotype_final$core_pathway)]
#   # bam[!is.na(bam)]
#   decoy <- paste(bam, collapse = "_")
# 
# 
#   if (nchar(decoy) < 1) {
#     nr_max = 0
# 
#     } else {
#     nr_max = 1
#     sam <- diplotype_final$hgnc_symbol[diplotype_final$sample == sample_id & diplotype_final$a1.max_score >= 3 & diplotype_final$a2.max_score >= 3 & !is.na(diplotype_final$core_pathway)]
# 
#   }
# 
# 
#   if (nr_max >= 1) {
# 
#     repair_genes <- repair_genes <- paste(sam, collapse = "_")
#     # print(repair_genes)
#   } else if (nr_max == 0) {
#     # print("No DNA repair deficiency")
#     decoy <- "No deficiency"
#     repair_genes <- NA
#   }
# 
#   metadata_included$dna_repair_deficiency[metadata_included$sample_id == sample_id] <- decoy
#   metadata_included$dna_repair_genes[metadata_included$sample_id == sample_id] <- repair_genes
# }
# 
# 
# 
# write.table(metadata_included, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/metadata-included-ddrd-annotated-biallelic.txt",
#             sep = "\t", quote = F, row.names = F)
# 





# metadata_included_ann_bi <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/metadata-included-ddrd-annotated-biallelic.txt",
#             sep = "\t", header = T, stringsAsFactors = F)
# 
# pathways <- c("No deficiency", unique(genes$Pathway_abb), "POL")
# 
# for (pathway in pathways) {
#   metadata_included_ann_bi[,pathway] <- NA
#   metadata_included_ann_bi[,pathway] <- str_detect(metadata_included_ann_bi$dna_repair_deficiency, pattern = pathway)
# }
# 
# 
# 
# metadata_included_ann_bi_tibb <- tidyr::gather(metadata_included_ann_bi, colnames(metadata_included_ann_bi[,27:40]), key = "pathway_abb", value = "is_dysfunctional", )
# metadata_included_ann_bi_tibb <- metadata_included_ann_bi_tibb[metadata_included_ann_bi_tibb$is_dysfunctional,]
# 
# 
# metadata_included_ann_bi_tibb$pathway_abb <- factor(metadata_included_ann_bi_tibb$pathway_abb, levels = pathways)
# metadata_included_ann_bi_tibb$cohort <- factor(metadata_included_ann_bi_tibb$cohort, levels = c("PCAWG", "HMF"))
# metadata_included_ann_bi_tibb$whole_genome_duplication <- factor(metadata_included_ann_bi_tibb$whole_genome_duplication)
# metadata_included_ann_bi_tibb$stage[metadata_included_ann_bi_tibb$is_metastatic] <- "Metastatic"
# metadata_included_ann_bi_tibb$stage[!metadata_included_ann_bi_tibb$is_metastatic] <- "Primary"
# metadata_included_ann_bi_tibb$stage <- factor(metadata_included_ann_bi_tibb$stage, levels = c("Primary", "Metastatic"))
# 
# # Remove hypermutated samples
# metadata_included_ann_bi_tibb2 <- metadata_included_ann_bi_tibb[!(metadata_included_ann_bi_tibb$is_hypermutated),]
# 

# =====================================================================================================================================

### Merging sample level ddrd metadata (used by cancer_types_HMF_PCAWG_2-metadata_23072021.tsv) _ mono-allelic hit:



# metadata_included$dna_repair_deficiency <- NA
# metadata_included$dna_repair_genes <- NA
# i <- 0



# for (sample_id in unique(diplotype_final$sample)) {
#   i <- i + 1
#   print(i)
#   bam <- diplotype_final$core_pathway[diplotype_final$sample == sample_id & (diplotype_final$a1.max_score >= 4 | diplotype_final$a2.max_score >= 4) & diplotype_final$biall_status != "loh_only" & !is.na(diplotype_final$core_pathway)]
#   # bam[!is.na(bam)]
#   decoy <- paste(bam, collapse = "_")
# 
# 
#   if (nchar(decoy) < 1) {
#     nr_max = 0
# 
#   } else {
#     nr_max = 1
#     sam <- diplotype_final$hgnc_symbol[diplotype_final$sample == sample_id & (diplotype_final$a1.max_score >= 4 | diplotype_final$a2.max_score >= 4) & diplotype_final$biall_status != "loh_only" & !is.na(diplotype_final$core_pathway)]
# 
#   }
# 
# 
#   if (nr_max >= 1) {
# 
#     repair_genes <- repair_genes <- paste(sam, collapse = "_")
#     # print(repair_genes)
#   } else if (nr_max == 0) {
#     # print("No DNA repair deficiency")
#     decoy <- "No deficiency"
#     repair_genes <- NA
#   }
# 
#   metadata_included$dna_repair_deficiency[metadata_included$sample_id == sample_id] <- decoy
#   metadata_included$dna_repair_genes[metadata_included$sample_id == sample_id] <- repair_genes
# }
# 
# 
# write.table(metadata_included, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/metadata-included-ddrd-annotated-monoallelic.txt",
#             sep = "\t", quote = F, row.names = F)





# metadata_included_ann_mono <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/metadata-included-ddrd-annotated-monoallelic.txt",
#                                      sep = "\t", header = T, stringsAsFactors = F)
# 
# pathways <- c("No deficiency", unique(genes$Pathway_abb), "POL")
# 
# for (pathway in pathways) {
#   metadata_included_ann_mono[,pathway] <- NA
#   metadata_included_ann_mono[,pathway] <- str_detect(metadata_included_ann_mono$dna_repair_deficiency, pattern = pathway)
# }
# 
# 
# metadata_included_ann_mono_tibb <- tidyr::gather(metadata_included_ann_mono, colnames(metadata_included_ann_mono[,27:40]), key = "pathway_abb", value = "is_dysfunctional", )
# metadata_included_ann_mono_tibb <- metadata_included_ann_mono_tibb[metadata_included_ann_mono_tibb$is_dysfunctional,]
# 
# 
# metadata_included_ann_mono_tibb$pathway_abb <- factor(metadata_included_ann_mono_tibb$pathway_abb, levels = pathways)
# metadata_included_ann_mono_tibb$cohort <- factor(metadata_included_ann_mono_tibb$cohort, levels = c("PCAWG", "HMF"))
# metadata_included_ann_mono_tibb$whole_genome_duplication <- factor(metadata_included_ann_mono_tibb$whole_genome_duplication)
# metadata_included_ann_mono_tibb$stage[metadata_included_ann_mono_tibb$is_metastatic] <- "Metastatic"
# metadata_included_ann_mono_tibb$stage[!metadata_included_ann_mono_tibb$is_metastatic] <- "Primary"
# metadata_included_ann_mono_tibb$stage <- factor(metadata_included_ann_mono_tibb$stage, levels = c("Primary", "Metastatic"))
# 
# 
# # Arne's questions
# table(metadata_included_ann_mono_tibb$pathway_abb[metadata_included_ann_mono_tibb$is_metastatic & !(metadata_included_ann_mono_tibb$is_hypermutated)])
# nrow(metadata_included_ann_mono_tibb[metadata_included_ann_mono_tibb$is_metastatic & !(metadata_included_ann_mono_tibb$is_hypermutated),])
# table(metadata_included_ann_mono$hr_status[metadata_included_ann_mono$HDR  & !(metadata_included_ann_mono$is_hypermutated)])
# table(metadata_included_ann_mono$hr_status[!(metadata_included_ann_mono$HDR) & !(metadata_included_ann_mono$is_hypermutated)])
# nrow(metadata_included_ann_mono[metadata_included_ann_mono$is_metastatic & !(metadata_included_ann_mono$is_hypermutated),])
# 
# 
# 
# # Remove hypermutated samples
# metadata_included_ann_mono_tibb2 <- metadata_included_ann_mono_tibb[!(metadata_included_ann_mono_tibb$is_hypermutated),]

# =====================================================================================================================================




# # Fig 1: using biallelic definition
# 
# d <- paste0(names(table(metadata_included_ann_bi_tibb2$pathway_abb)), " (n=", as.vector(table(metadata_included_ann_bi_tibb2$pathway_abb)), ")")
# 
# 
# tmb_plot_bi <- metadata_included_ann_bi_tibb %>% ggplot(aes(x = pathway_abb, y = log2(tmb), color = stage)) +
#   geom_boxplot(outlier.shape = NA) +
#   # stat_compare_means(method = "t.test") +
#   # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins,
#   #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
#   # scale_color_manual(values = c("blue", "green")) +
#   # geom_jitter(color="black", size=0.4, alpha=0.5) +
#   geom_point(position=position_jitterdodge(), size = 0.4, aes(color = stage), alpha = 0.5) +
#   scale_color_manual(values = c("#ff0101", "#010dff")) +
#   # stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = F) +
#   theme_bw() +
#   ggtitle("DDR deficiency Effects on TMB (bi-allelic criteria) \n") +
#   ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "wilcox.test", label = "p.signif") +
#   # ylim(0, 12) +
#   xlab("DNA repair pathway mutated \n") +
#   scale_x_discrete(labels = d) +
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# 
# 
# for (i in 1:2){
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig1-ddrd-tmb-boxplot-bi-no-hypermutated.png", width = 960, height = 960)
#     print(tmb_plot_bi)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig1-ddrd-tmb-boxplot-bi-no-hypermutated.pdf", width = 14, height = 14)
#     print(tmb_plot_bi)
#     dev.off()
#   }
# }
# 
# 
# 
# 
# # =====================================================================================================================================
# 
# 
# # Fig 2: Using monoallelic definition (primary and metastatic together)
# 
# d <- paste0(names(table(metadata_included_ann_mono_tibb2$pathway_abb)), " (n=", as.vector(table(metadata_included_ann_mono_tibb2$pathway_abb)), ")")
# 
# # normalizing_ratio <- mean(metadata_included_ann_mono_tibb2[metadata_included_ann_mono_tibb2$cohort == "HMF" & metadata_included_ann_mono_tibb2$pathway_abb == "no_deficiency", "tmb"], na.rm = T)/
# #   mean(metadata_included_ann_mono_tibb2[metadata_included_ann_mono_tibb2$cohort == "PCAWG" & metadata_included_ann_mono_tibb2$pathway_abb == "no_deficiency", "tmb"], na.rm = T)
# # 
# # metadata_included_ann_mono_tibb2[metadata_included_ann_mono_tibb2$cohort == "PCAWG", "tmb"] <- normalizing_ratio*metadata_included_ann_mono_tibb2[metadata_included_ann_mono_tibb2$cohort == "PCAWG", "tmb"]
# 
# 
# 
# tmb_plot_mono <- metadata_included_ann_mono_tibb2 %>% ggplot(aes(x = pathway_abb, y = log2(tmb), color = stage)) +
#   geom_boxplot(outlier.shape = NA) +
#   # stat_compare_means(method = "t.test") +
#   # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins,
#   #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
#   # scale_color_manual(values = c("blue", "green")) +
#   # geom_jitter(color="black", size=0.4, alpha=0.5) +
#   geom_point(position=position_jitterdodge(), size = 0.4, aes(color = stage), alpha = 0.5) +
#   scale_color_manual(values = c("#ff0101", "#010dff")) +
#   # stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = F) +
#   theme_bw() +
#   ggtitle("DDR deficiency Effects on TMB  (mono-allelic criteria)\n") +
#   ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   # ylim(0, 12) +
#   xlab("DNA repair pathway mutated \n") +
#   scale_x_discrete(labels = d) +
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# 
# 
# for (i in 1:2){
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig2-ddrd-tmb-boxplot-mono-no-hypermutated.png", width = 960, height = 960)
#     print(tmb_plot_mono)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig2-ddrd-tmb-boxplot-mono-no-hypermutated.pdf", width = 14, height = 14)
#     print(tmb_plot_mono)
#     dev.off()
#   }
# }
# 
# 
# # =====================================================================================================================================
# 
# # Fig 3: mono-allelic primary only
# 
# metadata_included_ann_mono_tibb2_primary <- metadata_included_ann_mono_tibb2[metadata_included_ann_mono_tibb2$stage == "Primary",]
# 
# d <- paste0(names(table(metadata_included_ann_mono_tibb2_primary$pathway_abb)), " (n=", as.vector(table(metadata_included_ann_mono_tibb2_primary$pathway_abb)), ")")
# 
# 
# tmb_plot_mono_primary <- metadata_included_ann_mono_tibb2_primary %>% ggplot(aes(x = pathway_abb, y = log2(tmb))) +
#   geom_boxplot(outlier.shape = NA) +
#   # stat_compare_means(method = "t.test") +
#   # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins,
#   #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
#   # scale_color_manual(values = c("blue", "green")) +
#   # geom_jitter(color="black", size=0.4, alpha=0.5) +
#   # geom_point(position=position_jitterdodge(), size = 0.4, alpha = 0.5) +
#   # scale_color_manual(values = c("#ff0101", "#010dff")) +
#   stat_compare_means(method = "anova") +
#   theme_bw() +
#   ggtitle("DDR deficiency Effects on TMB  in Primary Cancers (mono-allelic criteria)\n") +
#   ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   # ylim(0, 12) +
#   xlab("DNA repair pathway mutated \n") +
#   scale_x_discrete(labels = d) +
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# for (i in 1:2){
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig3-ddrd-tmb-boxplot-mono-no-hypermutated-primary.png", width = 960, height = 960)
#     print(tmb_plot_mono_primary)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig3-ddrd-tmb-boxplot-mono-no-hypermutated-primary.pdf", width = 14, height = 14)
#     print(tmb_plot_mono_primary)
#     dev.off()
#   }
# }
# 
# 
# 
# 
# 
# # =====================================================================================================================================
# 
# # Fig 4: mono-allelic metastatic only
# 
# metadata_included_ann_mono_tibb2_metastatic <- metadata_included_ann_mono_tibb2[metadata_included_ann_mono_tibb2$stage == "Metastatic",]
# 
# d <- paste0(names(table(metadata_included_ann_mono_tibb2_metastatic$pathway_abb)), " (n=", as.vector(table(metadata_included_ann_mono_tibb2_metastatic$pathway_abb)), ")")
# 
# 
# tmb_plot_mono_metastatic <- metadata_included_ann_mono_tibb2_metastatic %>% ggplot(aes(x = pathway_abb, y = log2(tmb))) +
#   geom_boxplot(outlier.shape = NA) +
#   # stat_compare_means(method = "t.test") +
#   # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins,
#   #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
#   # scale_color_manual(values = c("blue", "green")) +
#   # geom_jitter(color="black", size=0.4, alpha=0.5) +
#   # geom_point(position=position_jitterdodge(), size = 0.4, alpha = 0.5) +
#   # scale_color_manual(values = c("#ff0101", "#010dff")) +
#   stat_compare_means(method = "anova") +
#   theme_bw() +
#   ggtitle("DDR deficiency Effects on TMB  in Metastatic Cancers (mono-allelic criteria)\n") +
#   ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   # ylim(0, 12) +
#   xlab("DNA repair pathway mutated \n") +
#   scale_x_discrete(labels = d) +
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# 
# 
# for (i in 1:2){
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig4-ddrd-tmb-boxplot-mono-no-hypermutated-metastatic.png", width = 960, height = 960)
#     print(tmb_plot_mono_metastatic)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig4-ddrd-tmb-boxplot-mono-no-hypermutated-metastatic.pdf", width = 14, height = 14)
#     print(tmb_plot_mono_metastatic)
#     dev.off()
#   }
# }
# 
# 
# 
# # =====================================================================================================================================
# 
# 
# 
# 
# # Fig 5: mono-allelic primary only divided by wgd status
# 
# d <- paste0(names(table(metadata_included_ann_mono_tibb2_primary$pathway_abb)), " (n=", as.vector(table(metadata_included_ann_mono_tibb2_primary$pathway_abb)), ")")
# 
# 
# tmb_plot_mono_primary_wgd <- metadata_included_ann_mono_tibb2_primary %>% ggplot(aes(x = pathway_abb, y = log2(tmb), color = whole_genome_duplication)) +
#   geom_boxplot(outlier.shape = NA) +
#   # stat_compare_means(method = "t.test") +
#   # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins,
#   #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
#   # scale_color_manual(values = c("blue", "green")) +
#   # geom_jitter(color="black", size=0.4, alpha=0.5) +
#   # geom_point(position=position_jitterdodge(), size = 0.4, alpha = 0.5) +
#   # scale_color_manual(values = c("#ff0101", "#010dff")) +
#   theme_bw() +
#   stat_compare_means(method = "wilcox.test", label = "p.signif") +
#   ggtitle("DDR deficiency Effects on TMB in Primary Cancers based on WGD status (mono-allelic criteria)\n") +
#   ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   # ylim(0, 12) +
#   xlab("DNA repair pathway mutated \n") +
#   scale_x_discrete(labels = d) +
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# for (i in 1:2){
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig5-ddrd-tmb-boxplot-mono-no-hypermutated-primary-wgd.png", width = 960, height = 960)
#     print(tmb_plot_mono_primary_wgd)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig5-ddrd-tmb-boxplot-mono-no-hypermutated-primary-wgd.pdf", width = 14, height = 14)
#     print(tmb_plot_mono_primary_wgd)
#     dev.off()
#   }
# }
# 
# 
# # =====================================================================================================================================
# 
# 
# 
# 
# # Fig 6: mono-allelic metastatic only divided by wgd status
# 
# d <- paste0(names(table(metadata_included_ann_mono_tibb2_metastatic$pathway_abb)), " (n=", as.vector(table(metadata_included_ann_mono_tibb2_metastatic$pathway_abb)), ")")
# 
# 
# tmb_plot_mono_metastatic_wgd <- metadata_included_ann_mono_tibb2_metastatic %>% ggplot(aes(x = pathway_abb, y = log2(tmb), color = whole_genome_duplication)) +
#   geom_boxplot(outlier.shape = NA) +
#   # stat_compare_means(method = "t.test") +
#   # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins,
#   #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
#   # scale_color_manual(values = c("blue", "green")) +
#   # geom_jitter(color="black", size=0.4, alpha=0.5) +
#   # geom_point(position=position_jitterdodge(), size = 0.4, alpha = 0.5) +
#   # scale_color_manual(values = c("#ff0101", "#010dff")) +
#   # stat_compare_means(method = "anova") +
#   theme_bw() +
#   stat_compare_means(method = "wilcox.test", label = "p.signif") +
#   ggtitle("DDR deficiency Effects on TMB in Metastatic Cancers based on WGD status (mono-allelic criteria)\n") +
#   ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", aes(label = "p.signif")) +
#   # ylim(0, 12) +
#   xlab("DNA repair pathway mutated \n") +
#   scale_x_discrete(labels = d) +
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# 
# for (i in 1:2){
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig6-ddrd-tmb-boxplot-mono-no-hypermutated-metastatic-wgd.png", width = 960, height = 960)
#     print(tmb_plot_mono_metastatic_wgd)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig6-ddrd-tmb-boxplot-mono-no-hypermutated-metastatic-wgd.pdf", width = 14, height = 14)
#     print(tmb_plot_mono_metastatic_wgd)
#     dev.off()
#   }
# }
# 
# 
# # =====================================================================================================================================
# 
# # Fig 7: bi-allelic primary only
# 
# 
# metadata_included_ann_bi_tibb2_primary <- metadata_included_ann_bi_tibb2[metadata_included_ann_bi_tibb2$stage == "Primary",]
# 
# d <- paste0(names(table(metadata_included_ann_bi_tibb2_primary$pathway_abb)), " (n=", as.vector(table(metadata_included_ann_bi_tibb2_primary$pathway_abb)), ")")
# 
# 
# tmb_plot_bi_primary <- metadata_included_ann_bi_tibb2_primary %>% ggplot(aes(x = pathway_abb, y = log2(tmb))) +
#   geom_boxplot(outlier.shape = NA) +
#   # stat_compare_means(method = "t.test") +
#   # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins,
#   #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
#   # scale_color_manual(values = c("blue", "green")) +
#   # geom_jitter(color="black", size=0.4, alpha=0.5) +
#   # geom_point(position=position_jitterdodge(), size = 0.4, alpha = 0.5) +
#   # scale_color_manual(values = c("#ff0101", "#010dff")) +
#   stat_compare_means(method = "anova") +
#   theme_bw() +
#   ggtitle("DDR deficiency Effects on TMB  in Primary Cancers (bi-allelic criteria)\n") +
#   ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   # ylim(0, 12) +
#   xlab("DNA repair pathway mutated \n") +
#   scale_x_discrete(labels = d) +
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# for (i in 1:2){
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig7-ddrd-tmb-boxplot-bi-no-hypermutated-primary.png", width = 960, height = 960)
#     print(tmb_plot_bi_primary)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig7-ddrd-tmb-boxplot-bi-no-hypermutated-primary.pdf", width = 14, height = 14)
#     print(tmb_plot_bi_primary)
#     dev.off()
#   }
# }
# 
# 
# # =====================================================================================================================================
# 
# # Fig 8: bi-allelic metastatic only
# 
# 
# metadata_included_ann_bi_tibb2_metastatic <- metadata_included_ann_bi_tibb2[metadata_included_ann_bi_tibb2$stage == "Metastatic",]
# 
# d <- paste0(names(table(metadata_included_ann_bi_tibb2_metastatic$pathway_abb)), " (n=", as.vector(table(metadata_included_ann_bi_tibb2_metastatic$pathway_abb)), ")")
# 
# 
# tmb_plot_bi_metastatic <- metadata_included_ann_bi_tibb2_metastatic %>% ggplot(aes(x = pathway_abb, y = log2(tmb))) +
#   geom_boxplot(outlier.shape = NA) +
#   # stat_compare_means(method = "t.test") +
#   # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins,
#   #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
#   # scale_color_manual(values = c("blue", "green")) +
#   # geom_jitter(color="black", size=0.4, alpha=0.5) +
#   # geom_point(position=position_jitterdodge(), size = 0.4, alpha = 0.5) +
#   # scale_color_manual(values = c("#ff0101", "#010dff")) +
#   stat_compare_means(method = "anova") +
#   theme_bw() +
#   ggtitle("DDR deficiency Effects on TMB  in Metastatic Cancers (bi-allelic criteria)\n") +
#   ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   # ylim(0, 12) +
#   xlab("DNA repair pathway mutated \n") +
#   scale_x_discrete(labels = d) +
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# 
# 
# for (i in 1:2){
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig8-ddrd-tmb-boxplot-bi-no-hypermutated-metastatic.png", width = 960, height = 960)
#     print(tmb_plot_bi_metastatic)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig8-ddrd-tmb-boxplot-bi-no-hypermutated-metastatic.pdf", width = 14, height = 14)
#     print(tmb_plot_bi_metastatic)
#     dev.off()
#   }
# }
# 
# # =====================================================================================================================================
# 
# # Fig 9: bi-allelic primary only cancer type specific
# 
# metadata_included_ann_bi_tibb_primary <- metadata_included_ann_bi_tibb[metadata_included_ann_bi_tibb$stage == "Primary",]
# 
# 
# 
# 
# m <- as.data.frame(table(metadata_included_ann_bi_tibb_primary$cancer_type, metadata_included_ann_bi_tibb_primary$pathway_abb))
# y <- m[m$Freq >= 5,]
# 
# 
# metadata_included_ann_bi_tibb_primary$tti <- FALSE
# for (i in 1:nrow(y)) {
#   metadata_included_ann_bi_tibb_primary$tti[metadata_included_ann_bi_tibb_primary$cancer_type == y[i,"Var1"] & metadata_included_ann_bi_tibb_primary$pathway_abb == y[i,"Var2"]] <- TRUE
# }
# 
# 
# tmb_plot_bi_primary <- metadata_included_ann_bi_tibb_primary[metadata_included_ann_bi_tibb_primary$tti,] %>% ggplot(aes(x = pathway_abb, y = log2(tmb))) + facet_wrap(~cancer_type) +
#   geom_boxplot(outlier.shape = NA) +
#   # stat_compare_means(method = "t.test") +
#   # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins,
#   #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
#   # scale_color_manual(values = c("blue", "green")) +
#   # geom_jitter(color="black", size=0.4, alpha=0.5) +
#   # geom_point(position=position_jitterdodge(), size = 0.4, alpha = 0.5) +
#   # scale_color_manual(values = c("#ff0101", "#010dff")) +
#   stat_compare_means(method = "anova") +
#   theme_bw() +
#   ggtitle("DDR deficiency Effects on TMB  in Primary Cancers (bi-allelic criteria) _ facetted by cancer type\n") +
#   ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   # ylim(0, 12) +
#   xlab("DNA repair pathway mutated \n") +
#   # scale_x_discrete(labels = d) +
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# 
# 
# 
# for (i in 1:2){
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig9-ddrd-tmb-boxplot-bi-with-hypermutated-primary-cancer-type.png", width = 960, height = 960)
#     print(tmb_plot_bi_primary)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig9-ddrd-tmb-boxplot-bi-with-hypermutated-primary-cancer-type.pdf", width = 14, height = 14)
#     print(tmb_plot_bi_primary)
#     dev.off()
#   }
# }
# 
# 
# # =====================================================================================================================================
# 
# # Fig 10: bi-allelic metastatic only cancer type specific
# 
# 
# metadata_included_ann_bi_tibb_metastatic <- metadata_included_ann_bi_tibb[metadata_included_ann_bi_tibb$stage == "Metastatic",]
# 
# 
# m <- as.data.frame(table(metadata_included_ann_bi_tibb_metastatic$cancer_type, metadata_included_ann_bi_tibb_metastatic$pathway_abb))
# y <- m[m$Freq >= 5,]
# 
# 
# metadata_included_ann_bi_tibb_metastatic$tti <- FALSE
# for (i in 1:nrow(y)) {
#   metadata_included_ann_bi_tibb_metastatic$tti[metadata_included_ann_bi_tibb_metastatic$cancer_type == y[i,"Var1"] & metadata_included_ann_bi_tibb_metastatic$pathway_abb == y[i,"Var2"]] <- TRUE
# }
# 
# 
# tmb_plot_bi_metastatic <- metadata_included_ann_bi_tibb_metastatic[metadata_included_ann_bi_tibb_metastatic$tti,] %>% ggplot(aes(x = pathway_abb, y = log2(tmb))) + facet_wrap(~cancer_type) +
#   geom_boxplot(outlier.shape = NA) +
#   # stat_compare_means(method = "t.test") +
#   # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins,
#   #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
#   # scale_color_manual(values = c("blue", "green")) +
#   # geom_jitter(color="black", size=0.4, alpha=0.5) +
#   # geom_point(position=position_jitterdodge(), size = 0.4, alpha = 0.5) +
#   # scale_color_manual(values = c("#ff0101", "#010dff")) +
#   stat_compare_means(method = "anova") +
#   theme_bw() +
#   ggtitle("DDR deficiency Effects on TMB  in Metastatic Cancers (bi-allelic criteria) _ facetted by cancer type\n") +
#   ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   # ylim(0, 12) +
#   xlab("DNA repair pathway mutated \n") +
#   scale_x_discrete(labels = d) +
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# 
# 
# for (i in 1:2){
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig10-ddrd-tmb-boxplot-bi-with-hypermutated-metastatic-cancer-type.png", width = 960, height = 960)
#     print(tmb_plot_bi_metastatic)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig10-ddrd-tmb-boxplot-bi-with-hypermutated-metastatic-cancer-type.pdf", width = 14, height = 14)
#     print(tmb_plot_bi_metastatic)
#     dev.off()
#   }
# }
# 
# 
# 
# # =====================================================================================================================================
# 
# # Fig 11: mono-allelic primary only cancer type specific
# 
# metadata_included_ann_mono_tibb_primary <- metadata_included_ann_mono_tibb[metadata_included_ann_mono_tibb$stage == "Primary",]
# 
# 
# m <- as.data.frame(table(metadata_included_ann_mono_tibb_primary$cancer_type, metadata_included_ann_mono_tibb_primary$pathway_abb))
# y <- m[m$Freq >= 5,]
# 
# 
# metadata_included_ann_mono_tibb_primary$tti <- FALSE
# for (i in 1:nrow(y)) {
#   metadata_included_ann_mono_tibb_primary$tti[metadata_included_ann_mono_tibb_primary$cancer_type == y[i,"Var1"] & metadata_included_ann_mono_tibb_primary$pathway_abb == y[i,"Var2"]] <- TRUE
# }
# 
# 
# tmb_plot_mono_primary <- metadata_included_ann_mono_tibb_primary[metadata_included_ann_mono_tibb_primary$tti,] %>% ggplot(aes(x = pathway_abb, y = log2(tmb))) + facet_wrap(~cancer_type) +
#   geom_boxplot(outlier.shape = NA) +
#   # stat_compare_means(method = "t.test") +
#   # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins,
#   #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
#   # scale_color_manual(values = c("blue", "green")) +
#   # geom_jitter(color="black", size=0.4, alpha=0.5) +
#   # geom_point(position=position_jitterdodge(), size = 0.4, alpha = 0.5) +
#   # scale_color_manual(values = c("#ff0101", "#010dff")) +
#   stat_compare_means(method = "anova") +
#   theme_bw() +
#   ggtitle("DDR deficiency Effects on TMB  in Primary Cancers (mono-allelic criteria) _ facetted by cancer type\n") +
#   ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   # ylim(0, 12) +
#   xlab("DNA repair pathway mutated \n") +
#   # scale_x_discrete(labels = d) +
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# 
# 
# 
# for (i in 1:2){
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig11-ddrd-tmb-boxplot-mono-with-hypermutated-primary-cancer-type.png", width = 960, height = 960)
#     print(tmb_plot_mono_primary)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig11-ddrd-tmb-boxplot-mono-with-hypermutated-primary-cancer-type.pdf", width = 14, height = 14)
#     print(tmb_plot_mono_primary)
#     dev.off()
#   }
# }
# 
# 
# # =====================================================================================================================================
# 
# # Fig 12: mono-allelic metastatic only cancer type specific
# 
# 
# 
# 
# metadata_included_ann_mono_tibb_metastatic <- metadata_included_ann_mono_tibb[metadata_included_ann_mono_tibb$stage == "Metastatic",]
# 
# 
# 
# 
# m <- as.data.frame(table(metadata_included_ann_mono_tibb_metastatic$cancer_type, metadata_included_ann_mono_tibb_metastatic$pathway_abb))
# y <- m[m$Freq >= 5,]
# 
# 
# metadata_included_ann_mono_tibb_metastatic$tti <- FALSE
# for (i in 1:nrow(y)) {
#   metadata_included_ann_mono_tibb_metastatic$tti[metadata_included_ann_mono_tibb_metastatic$cancer_type == y[i,"Var1"] & metadata_included_ann_mono_tibb_metastatic$pathway_abb == y[i,"Var2"]] <- TRUE
# }
# 
# 
# tmb_plot_mono_metastatic <- metadata_included_ann_mono_tibb_metastatic[metadata_included_ann_mono_tibb_metastatic$tti,] %>% ggplot(aes(x = pathway_abb, y = log2(tmb))) + facet_wrap(~cancer_type) +
#   geom_boxplot(outlier.shape = NA) +
#   # stat_compare_means(method = "t.test") +
#   # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins,
#   #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
#   # scale_color_manual(values = c("blue", "green")) +
#   # geom_jitter(color="black", size=0.4, alpha=0.5) +
#   # geom_point(position=position_jitterdodge(), size = 0.4, alpha = 0.5) +
#   # scale_color_manual(values = c("#ff0101", "#010dff")) +
#   stat_compare_means(method = "anova") +
#   theme_bw() +
#   ggtitle("DDR deficiency Effects on TMB  in Metastatic Cancers (mono-allelic criteria) _ facetted by cancer type \n") +
#   ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   # ylim(0, 12) +
#   xlab("DNA repair pathway mutated \n") +
#   scale_x_discrete(labels = d) +
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# 
# 
# for (i in 1:2){
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig12-ddrd-tmb-boxplot-mono-with-hypermutated-metastatic-cancer-type.png", width = 960, height = 960)
#     print(tmb_plot_mono_metastatic)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig12-ddrd-tmb-boxplot-mono-with-hypermutated-metastatic-cancer-type.pdf", width = 14, height = 14)
#     print(tmb_plot_mono_metastatic)
#     dev.off()
#   }
# }
# 
# 
# # =====================================================================================================================================
# 
# # Fig 14: wgd frequencies
# 
# 
# wgd_true_prim <- table(metadata_included_ann_mono_tibb2_primary[metadata_included_ann_mono_tibb2_primary$whole_genome_duplication == "TRUE","pathway_abb"])
# wgd_false_prim <- table(metadata_included_ann_mono_tibb2_primary[metadata_included_ann_mono_tibb2_primary$whole_genome_duplication == "FALSE","pathway_abb"])
# 
# wgd_info_prim <- data.frame(pathway = names(wgd_true_prim), with_WGD = as.vector(wgd_true_prim), without_WGD = as.vector(wgd_false_prim))
# wgd_info_prim[,"WGD(+)/WGD(-)"] <- wgd_info_prim$with_WGD/wgd_info_prim$without_WGD
# wgd_info_prim[,"Cohort"] <- "Primary"
# rownames(wgd_info_prim) <- wgd_info_prim$pathway
# wgd_info_prim <- wgd_info_prim[,-1]
# 
# 
# 
# j <- 1
# odds_df_prim <- data.frame(pathway = NA,
#                            cohort = NA,
#                            odds_ratio = NA,
#                             lower_confidence_level = NA,
#                             upper_confidence_level = NA,
#                             effect_size = NA,
#                             p_val = NA)
# i <- 2
# for (i in 2:nrow(wgd_info_prim)) {
#   ff <- oddsratio.wald(as.matrix(wgd_info_prim[c(j,i),1:2])[c(2,1),])
#   odds_df_prim[i-1,] <- c(rownames(wgd_info_prim)[i], "Primary", ff$measure[2,1], ff$measure[2,2], ff$measure[2,3], as.vector(cramerV(as.matrix(wgd_info_prim[c(j,i),1:2]))), ff$p.value[2,3])
# }
# 
# 
# 
# wgd_true_metas <- table(metadata_included_ann_mono_tibb2_metastatic[metadata_included_ann_mono_tibb2_metastatic$whole_genome_duplication == "TRUE","pathway_abb"])
# wgd_false_metas <- table(metadata_included_ann_mono_tibb2_metastatic[metadata_included_ann_mono_tibb2_metastatic$whole_genome_duplication == "FALSE","pathway_abb"])
# 
# wgd_info_metas <- data.frame(pathway = names(wgd_true_metas), with_WGD = as.vector(wgd_true_metas), without_WGD = as.vector(wgd_false_metas))
# wgd_info_metas[,"WGD(+)/WGD(-)"] <- wgd_info_metas$with_WGD/wgd_info_metas$without_WGD
# wgd_info_metas[,"Cohort"] <- "Metastatic"
# rownames(wgd_info_metas) <- wgd_info_metas$pathway
# wgd_info_metas <- wgd_info_metas[,-1]
# 
# 
# j <- 1
# odds_df_metas <- data.frame(pathway = NA,
#                             cohort = NA,
#                             odds_ratio = NA,
#                             lower_confidence_level = NA,
#                             upper_confidence_level = NA,
#                             effect_size = NA,
#                             p_val = NA)
# 
# for (i in 2:nrow(wgd_info_metas)) {
#   ff <- oddsratio.wald(as.matrix(wgd_info_metas[c(j,i),1:2])[c(2,1),])
#   odds_df_metas[i-1,] <- c(rownames(wgd_info_metas)[i], "Metastatic", ff$measure[2,1], ff$measure[2,2], ff$measure[2,3], as.vector(cramerV(as.matrix(wgd_info_metas[c(j,i),1:2]))), ff$p.value[2,3])
# }
# 
# 
# odds_wgd_merged <- rbind(odds_df_prim, odds_df_metas)
# 
# odd_numbers <- seq(1,14,2)
# 
# str(odds_wgd_merged)
# 
# odds_wgd_merged[,c(3,4,5,6,7)] <- sapply(odds_wgd_merged[,c(3,4,5,6,7)], as.numeric)
# 
# 
# odds_wgd_plot <- odds_wgd_merged %>% ggplot(aes(x = pathway, y = odds_ratio, group = cohort)) +
#   geom_rect(data = odds_wgd_merged[odd_numbers, ], xmin = odd_numbers - 0.5, xmax = odd_numbers + 
#               0.5, ymin = -Inf, ymax = Inf, fill = 'grey', alpha = 0.5) +
#   geom_point(position = position_dodge(width = 1), aes(color = cohort, shape = cohort, size = effect_size)) +
#   geom_errorbar(position = position_dodge(width = 1), aes(ymin=lower_confidence_level, ymax=upper_confidence_level, color = cohort), lwd = 0.2) +
#   scale_color_manual(values = c("blue", "red")) +
#   scale_shape_manual(values=c(16, 18)) + 
#   # geom_line(aes(color = cohort)) +
#   theme_classic() +
#   ggtitle("Whole Genome Duplication Prevalence\n based on DDR Status \n") +
#   # ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   # ylim(0, 12) +
#   coord_cartesian(ylim = c(0,10))+
#   xlab("DNA repair pathway mutated \n") +
#   ylab("Odds Ratio of WGD(+)/WGD(-)") +
#   # scale_x_discrete(labels = d) +
#   scale_size_continuous(name = "Effect Size")+
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   # labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# 
# 
# for (i in 1:2){
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig14-ddrd-odds-wgd-mono-no-hypermutated.png", width = 960, height = 480)
#     print(odds_wgd_plot)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig14-ddrd-odds-wgd-mono-no-hypermutated.pdf", width = 14, height = 7)
#     print(odds_wgd_plot)
#     dev.off()
#   }
# }
# 
# 
# wgd_info_merged <- rbind(wgd_info_prim[,c(1,4,5)], wgd_info_metas[,c(1,4,5)])
# 
# 
# wgd_info_merged <- wgd_info_merged[order(match(wgd_info_merged$pathway, pathways)),]
# 
# wgd_info_merged$pathway <- factor(wgd_info_merged$pathway, levels = pathways)
# wgd_info_merged$Cohort <- factor(wgd_info_merged$Cohort, levels = c("Primary", "Metastatic"))
# 
# 
# my_plot <- wgd_info_merged %>% ggplot(aes(x = pathway, y = `WGD(+)/WGD(-)`, color = Cohort, shape = Cohort)) +
#   geom_dotplot(binaxis = "y", stackdir = "center") +
#   scale_color_manual(values = c("blue", "green"))
#   
# 
# # geom_point
# wgd_plot <- wgd_info_merged %>% ggplot(aes(x = pathway, y = `WGD(+)/WGD(-)`, group= Cohort)) +
#   geom_point(aes(color = Cohort, shape=Cohort), size = 5) +
#   scale_color_manual(values = c("blue", "red")) +
#   scale_shape_manual(values=c(16, 18)) + 
#   geom_line(aes(color = Cohort)) +
#   theme_classic() +
#   ggtitle("Whole Genome Duplication Prevalence based on DDR Status (mono-allelic criteria) \n") +
#   # ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   # ylim(0, 12) +
#   xlab("DNA repair pathway mutated \n") +
#   # scale_x_discrete(labels = d) +
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   # labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# # geom_bar
# wgd_plot <- wgd_info_merged %>% ggplot(aes(x = pathway, y = `WGD(+)/WGD(-)`)) +
#   geom_bar(position=position_dodge(0.5), stat = "identity",width=0.5,aes(fill = Cohort)) +
#   scale_fill_manual(values = alpha(c("blue", "red"), 0.8)) +
#   theme_classic() +
#   ggtitle("Whole Genome Duplication Prevalence based on DDR Status (mono-allelic criteria) \n") +
#   # ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   # ylim(0, 12) +
#   xlab("DNA repair pathway mutated \n") +
#   # scale_x_discrete(labels = d) +
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   # labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# for (i in 1:2){
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig13-ddrd-wgd-mono-no-hypermutated.png", width = 960, height = 960)
#     print(wgd_plot)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig13-ddrd-wgd-mono-no-hypermutated.pdf", width = 14, height = 14)
#     print(wgd_plot)
#     dev.off()
#   }
# }
# 
# 
# 
# 
# 
# 
# 
# 
# # =====================================================================================================================================
# # Fig 16: wgd frequencies per cancer type
# 
# colnames(table(metadata_included_ann_mono_tibb2_primary$whole_genome_duplication[metadata_included_ann_mono_tibb2_primary$pathway_abb == "CCR"], metadata_included_ann_mono_tibb2_primary$cancer_type[metadata_included_ann_mono_tibb2_primary$pathway_abb == "CCR"])) %in%
# colnames(table(metadata_included_ann_mono_tibb2_primary$whole_genome_duplication[metadata_included_ann_mono_tibb2_primary$pathway_abb == "No deficiency"], metadata_included_ann_mono_tibb2_primary$cancer_type[metadata_included_ann_mono_tibb2_primary$pathway_abb == "No deficiency"]))
# 
# metadata_included_ann_mono_tibb2_primary$whole_genome_duplication[metadata_included_ann_mono_tibb2_primary$pathway_abb == "CCR" & !(str_detect(metadata_included_ann_mono_tibb2_primary$dna_repair_genes, pattern = "TP53"))]
# 
# ccr <- table(metadata_included_ann_mono_tibb2_primary$whole_genome_duplication[metadata_included_ann_mono_tibb2_primary$pathway_abb == "CCR"], metadata_included_ann_mono_tibb2_primary$cancer_type[metadata_included_ann_mono_tibb2_primary$pathway_abb == "CCR"])
# no_def <- table(metadata_included_ann_mono_tibb2_primary$whole_genome_duplication[metadata_included_ann_mono_tibb2_primary$pathway_abb == "No deficiency"], metadata_included_ann_mono_tibb2_primary$cancer_type[metadata_included_ann_mono_tibb2_primary$pathway_abb == "No deficiency"])[,colnames(table(metadata_included_ann_mono_tibb2_primary$whole_genome_duplication[metadata_included_ann_mono_tibb2_primary$pathway_abb == "CCR"], metadata_included_ann_mono_tibb2_primary$cancer_type[metadata_included_ann_mono_tibb2_primary$pathway_abb == "CCR"]))]
# 
# 
# ccr <- ccr +1
# no_def <- no_def + 1
# 
# options(scipen=999)
# 
# per_cancer_type_wgd_ccr_odds_prim <- data.frame(tissue = character(26), odds_ratio = numeric(26), lower_ci_level = numeric(26), upper_ci_level = numeric(26), effect_size = numeric(26), p_value = numeric(26))
# 
# for (i in 1:length(colnames(ccr))){
#   
#   working_mat <- t(as.matrix(data.frame(ccr = ccr[,i], no_def = no_def[,i])))[c(2,1),]
#   gg <- oddsratio.wald(working_mat)
#   per_cancer_type_wgd_ccr_odds_prim[i,] <- c(colnames(ccr)[i], gg$measure[2,1], gg$measure[2,2], gg$measure[2,3], as.vector(cramerV(working_mat)), gg$p.value[2,3])
# }
# 
# 
# per_cancer_type_wgd_ccr_odds_prim <- per_cancer_type_wgd_ccr_odds_prim[per_cancer_type_wgd_ccr_odds_prim$p_value < 0.01,]
# 
# per_cancer_type_wgd_ccr_odds_prim$tissue <- factor(per_cancer_type_wgd_ccr_odds_prim$tissue)
# per_cancer_type_wgd_ccr_odds_prim[,2] <- as.numeric(per_cancer_type_wgd_ccr_odds_prim[,2])
# per_cancer_type_wgd_ccr_odds_prim[,3] <- as.numeric(per_cancer_type_wgd_ccr_odds_prim[,3])
# per_cancer_type_wgd_ccr_odds_prim[,4] <- as.numeric(per_cancer_type_wgd_ccr_odds_prim[,4])
# per_cancer_type_wgd_ccr_odds_prim[,5] <- as.numeric(per_cancer_type_wgd_ccr_odds_prim[,5])
# per_cancer_type_wgd_ccr_odds_prim[,6] <- as.numeric(per_cancer_type_wgd_ccr_odds_prim[,6])
# 
# 
# 
# str(per_cancer_type_wgd_ccr_odds_prim)
# 
# 
# odds_wgd_ccr_plot_primary <- per_cancer_type_wgd_ccr_odds_prim %>% ggplot(aes(x = tissue, y = odds_ratio)) +
#   geom_point(aes(size = effect_size)) +
#   # geom_errorbar( aes(ymin=lower_ci_level, ymax=upper_ci_level), lwd = 0.2) +
#   theme_classic() +
#   ggtitle("Whole Genome Duplication Prevalence\n based on DDR Status \n") +
#   # ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   # ylim(0, 12) +
#   # coord_cartesian(ylim = c(0,10))+
#   xlab("DNA repair pathway mutated \n") +
#   ylab("Odds Ratio of WGD(+)/WGD(-)") +
#   # scale_x_discrete(labels = d) +
#   scale_size_continuous(name = "Effect Size")+
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   # labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# 
# 
# 
# 
# metadata_included_ann_mono_tibb_primary <- metadata_included_ann_mono_tibb[metadata_included_ann_mono_tibb$stage == "Primary",]
# 
# 
# metadata_included_ann_mono_tibb_metastatic <- metadata_included_ann_mono_tibb[metadata_included_ann_mono_tibb$stage == "Metastatic",]
# 
# metadata_included_ann_mono_tibb2_primary <- metadata_included_ann_mono_tibb2[metadata_included_ann_mono_tibb2$stage == "Primary",]
# 
# metadata_included_ann_mono_tibb2_metastatic <- metadata_included_ann_mono_tibb2[metadata_included_ann_mono_tibb2$stage == "Metastatic",]
# 
# 
# 
# ccr <- table(metadata_included_ann_mono_tibb2_metastatic$whole_genome_duplication[metadata_included_ann_mono_tibb2_metastatic$pathway_abb == "CCR"], metadata_included_ann_mono_tibb2_metastatic$cancer_type[metadata_included_ann_mono_tibb2_metastatic$pathway_abb == "CCR"])
# ccr[,"Gastrointestinal tract unknown"]
# 
# # colnames(table(metadata_included_ann_mono_tibb2_metastatic$whole_genome_duplication[metadata_included_ann_mono_tibb2_metastatic$pathway_abb == "CCR"], metadata_included_ann_mono_tibb2_metastatic$cancer_type[metadata_included_ann_mono_tibb2_metastatic$pathway_abb == "CCR"]))[!colnames(table(metadata_included_ann_mono_tibb2_metastatic$whole_genome_duplication[metadata_included_ann_mono_tibb2_metastatic$pathway_abb == "CCR"], metadata_included_ann_mono_tibb2_metastatic$cancer_type[metadata_included_ann_mono_tibb2_metastatic$pathway_abb == "CCR"])) %in%
# # colnames(table(metadata_included_ann_mono_tibb2_metastatic$whole_genome_duplication[metadata_included_ann_mono_tibb2_metastatic$pathway_abb == "No deficiency"], metadata_included_ann_mono_tibb2_metastatic$cancer_type[metadata_included_ann_mono_tibb2_metastatic$pathway_abb == "No deficiency"]))]
# 
# drops <- c("Gastrointestinal tract unknown", "Lung cancer mixed")
# ccr <- table(metadata_included_ann_mono_tibb2_metastatic$whole_genome_duplication[metadata_included_ann_mono_tibb2_metastatic$pathway_abb == "CCR"], metadata_included_ann_mono_tibb2_metastatic$cancer_type[metadata_included_ann_mono_tibb2_metastatic$pathway_abb == "CCR"])[,!(colnames(ccr) %in% drops)]
# ccr <- ccr +1
# 
# no_def <- table(metadata_included_ann_mono_tibb2_metastatic$whole_genome_duplication[metadata_included_ann_mono_tibb2_metastatic$pathway_abb == "No deficiency"], metadata_included_ann_mono_tibb2_metastatic$cancer_type[metadata_included_ann_mono_tibb2_metastatic$pathway_abb == "No deficiency"])[,colnames(ccr)]
# no_def <- no_def +1
# 
# # colnames(ccr) == colnames(no_def)
# 
# per_cancer_type_wgd_ccr_odds_metas <- data.frame(tissue = character(43), odds_ratio = numeric(43), lower_ci_level = numeric(43), upper_ci_level = numeric(43), effect_size = numeric(43), p_value = numeric(43))
# 
# 
# for (i in 1:length(colnames(ccr))){
#   
#   working_mat <- t(as.matrix(data.frame(ccr = ccr[,i], no_def = no_def[,i])))[c(2,1),]
#   gg <- oddsratio.wald(working_mat)
#   per_cancer_type_wgd_ccr_odds_metas[i,] <- c(colnames(ccr)[i], gg$measure[2,1], gg$measure[2,2], gg$measure[2,3], as.vector(cramerV(working_mat)), gg$p.value[2,3])
# }
# 
# 
# per_cancer_type_wgd_ccr_odds_metas <- per_cancer_type_wgd_ccr_odds_metas[per_cancer_type_wgd_ccr_odds_metas$p_value < 0.01,]
# 
# 
# per_cancer_type_wgd_ccr_odds_metas$tissue <- factor(per_cancer_type_wgd_ccr_odds_metas$tissue)
# per_cancer_type_wgd_ccr_odds_metas[,2] <- as.numeric(per_cancer_type_wgd_ccr_odds_metas[,2])
# per_cancer_type_wgd_ccr_odds_metas[,3] <- as.numeric(per_cancer_type_wgd_ccr_odds_metas[,3])
# per_cancer_type_wgd_ccr_odds_metas[,4] <- as.numeric(per_cancer_type_wgd_ccr_odds_metas[,4])
# per_cancer_type_wgd_ccr_odds_metas[,5] <- as.numeric(per_cancer_type_wgd_ccr_odds_metas[,5])
# per_cancer_type_wgd_ccr_odds_metas[,6] <- as.numeric(per_cancer_type_wgd_ccr_odds_metas[,6])
# 
# 
# 
# 
# odds_wgd_ccr_plot_metastatic2
# 
# 
# odds_wgd_ccr_plot_metastatic <- per_cancer_type_wgd_ccr_odds_metas %>% ggplot(aes(x = tissue, y = odds_ratio)) +
#   geom_point(aes(size = effect_size)) +
#   # geom_errorbar( aes(ymin=lower_ci_level, ymax=upper_ci_level), lwd = 0.2) +
#   theme_classic() +
#   ggtitle("Whole Genome Duplication Prevalence\n based on DDR Status \n") +
#   # ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   # ylim(0, 12) +
#   # coord_cartesian(ylim = c(0,10))+
#   xlab("DNA repair pathway mutated \n") +
#   ylab("Odds Ratio of WGD(+)/WGD(-)") +
#   # scale_x_discrete(labels = d) +
#   scale_size_continuous(name = "Effect Size")+
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   # labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# 
# per_cancer_type_wgd_ccr_odds_prim$cohort <- "Primary"
# 
# per_cancer_type_wgd_ccr_odds_metas$cohort <- "Metastatic"
# 
# 
# per_cancer_type_wgd_ccr_odds_combined <- rbind(per_cancer_type_wgd_ccr_odds_prim, per_cancer_type_wgd_ccr_odds_metas)
# per_cancer_type_wgd_ccr_odds_combined$cohort <- factor(per_cancer_type_wgd_ccr_odds_combined$cohort)
# 
# 
# odds_wgd_ccr_plot_combined <- per_cancer_type_wgd_ccr_odds_combined %>% ggplot(aes(x = tissue, y = odds_ratio)) +
#   geom_point(aes(size = effect_size, color = cohort)) +
#   # geom_errorbar( aes(ymin=lower_ci_level, ymax=upper_ci_level), lwd = 0.2) +
#   theme_classic() +
#   geom_text(x = 11, y = 50, label = "P.val < 0.01") +
#   ggtitle("CCR-deficient samples show propensity \n to have WGD across cancer types \n") +
#   # ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   # ylim(0, 12) +
#   # coord_cartesian(ylim = c(0,10))+
#   xlab("Cancer Types \n") +
#   geom_hline(yintercept = 1, linetype = "dashed") +
#   ylab("Odds Ratio of WGD(+)/WGD(-) _  \nCCRd compared to DDR-proficient") +
#   # scale_x_discrete(labels = d) +
#   scale_size_continuous(name = "Effect Size")+
#   scale_color_discrete(name = "Cohort")+
#   theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 13)) +
#   theme(axis.title.y = element_text(face = "italic", size = 13)) +
#   # labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# 
# 
# for (i in 1:2){
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig16-odds-wgd-ccr-plot-combined.png")
#     print(odds_wgd_ccr_plot_combined)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig16-odds-wgd-ccr-plot-combined.pdf")
#     print(odds_wgd_ccr_plot_combined)
#     dev.off()
#   }
# }
# 
# 
# # wgd frequencies ccr per cancer type (no tp53)
# #
# 
# ccr <- table(metadata_included_ann_mono_tibb2_primary$whole_genome_duplication[metadata_included_ann_mono_tibb2_primary$pathway_abb == "CCR" & !(str_detect(metadata_included_ann_mono_tibb2_primary$dna_repair_genes, pattern = "TP53"))],
#              metadata_included_ann_mono_tibb2_primary$cancer_type[metadata_included_ann_mono_tibb2_primary$pathway_abb == "CCR" & !(str_detect(metadata_included_ann_mono_tibb2_primary$dna_repair_genes, pattern = "TP53"))])
# 
# no_def <- table(metadata_included_ann_mono_tibb2_primary$whole_genome_duplication[metadata_included_ann_mono_tibb2_primary$pathway_abb == "No deficiency"], metadata_included_ann_mono_tibb2_primary$cancer_type[metadata_included_ann_mono_tibb2_primary$pathway_abb == "No deficiency"])[,colnames(table(metadata_included_ann_mono_tibb2_primary$whole_genome_duplication[metadata_included_ann_mono_tibb2_primary$pathway_abb == "CCR"], metadata_included_ann_mono_tibb2_primary$cancer_type[metadata_included_ann_mono_tibb2_primary$pathway_abb == "CCR"]))]
# 
# # colnames(ccr) %in% colnames(no_def)
# ccr <- ccr +1
# no_def <- no_def + 1
# 
# 
# 
# 
# options(scipen=999)
# 
# per_cancer_type_wgd_ccr_no_tp53_odds_prim <- data.frame(tissue = character(19), odds_ratio = numeric(19), lower_ci_level = numeric(19), upper_ci_level = numeric(19), effect_size = numeric(19), p_value = numeric(19))
# 
# for (i in 1:length(colnames(ccr))){
#   
#   working_mat <- t(as.matrix(data.frame(ccr = ccr[,i], no_def = no_def[,i])))[c(2,1),]
#   gg <- oddsratio.wald(working_mat)
#   per_cancer_type_wgd_ccr_no_tp53_odds_prim[i,] <- c(colnames(ccr)[i], gg$measure[2,1], gg$measure[2,2], gg$measure[2,3], as.vector(cramerV(working_mat)), gg$p.value[2,3])
# }
# 
# 
# 
# per_cancer_type_wgd_ccr_no_tp53_odds_prim <- per_cancer_type_wgd_ccr_no_tp53_odds_prim[per_cancer_type_wgd_ccr_no_tp53_odds_prim$p_value < 0.01,]
# 
# per_cancer_type_wgd_ccr_no_tp53_odds_prim$tissue <- factor(per_cancer_type_wgd_ccr_no_tp53_odds_prim$tissue)
# per_cancer_type_wgd_ccr_no_tp53_odds_prim[,2] <- as.numeric(per_cancer_type_wgd_ccr_no_tp53_odds_prim[,2])
# per_cancer_type_wgd_ccr_no_tp53_odds_prim[,3] <- as.numeric(per_cancer_type_wgd_ccr_no_tp53_odds_prim[,3])
# per_cancer_type_wgd_ccr_no_tp53_odds_prim[,4] <- as.numeric(per_cancer_type_wgd_ccr_no_tp53_odds_prim[,4])
# per_cancer_type_wgd_ccr_no_tp53_odds_prim[,5] <- as.numeric(per_cancer_type_wgd_ccr_no_tp53_odds_prim[,5])
# per_cancer_type_wgd_ccr_no_tp53_odds_prim[,6] <- as.numeric(per_cancer_type_wgd_ccr_no_tp53_odds_prim[,6])
# 
# 
# 
# 
# 
# 
# 
# 
# ccr <- table(metadata_included_ann_mono_tibb2_metastatic$whole_genome_duplication[metadata_included_ann_mono_tibb2_metastatic$pathway_abb == "CCR" & !(str_detect(metadata_included_ann_mono_tibb2_metastatic$dna_repair_genes, pattern = "TP53"))], metadata_included_ann_mono_tibb2_metastatic$cancer_type[metadata_included_ann_mono_tibb2_metastatic$pathway_abb == "CCR" & !(str_detect(metadata_included_ann_mono_tibb2_metastatic$dna_repair_genes, pattern = "TP53"))])
# 
# no_def <- table(metadata_included_ann_mono_tibb2_metastatic$whole_genome_duplication[metadata_included_ann_mono_tibb2_metastatic$pathway_abb == "No deficiency"], metadata_included_ann_mono_tibb2_metastatic$cancer_type[metadata_included_ann_mono_tibb2_metastatic$pathway_abb == "No deficiency"])[,colnames(ccr)]
# 
# colnames(ccr) %in% colnames(no_def)
# 
# 
# per_cancer_type_wgd_ccr__no_tp53_odds_metas <- data.frame(tissue = character(24), odds_ratio = numeric(24), lower_ci_level = numeric(24), upper_ci_level = numeric(24), effect_size = numeric(24), p_value = numeric(24))
# 
# 
# for (i in 1:length(colnames(ccr))){
#   
#   working_mat <- t(as.matrix(data.frame(ccr = ccr[,i], no_def = no_def[,i])))[c(2,1),]
#   gg <- oddsratio.wald(working_mat)
#   per_cancer_type_wgd_ccr__no_tp53_odds_metas[i,] <- c(colnames(ccr)[i], gg$measure[2,1], gg$measure[2,2], gg$measure[2,3], as.vector(cramerV(working_mat)), gg$p.value[2,3])
# }
# 
# 
# per_cancer_type_wgd_ccr__no_tp53_odds_metas <- per_cancer_type_wgd_ccr__no_tp53_odds_metas[per_cancer_type_wgd_ccr__no_tp53_odds_metas$p_value < 0.01,]
# 
# 
# per_cancer_type_wgd_ccr__no_tp53_odds_metas$tissue <- factor(per_cancer_type_wgd_ccr__no_tp53_odds_metas$tissue)
# per_cancer_type_wgd_ccr__no_tp53_odds_metas[,2] <- as.numeric(per_cancer_type_wgd_ccr__no_tp53_odds_metas[,2])
# per_cancer_type_wgd_ccr__no_tp53_odds_metas[,3] <- as.numeric(per_cancer_type_wgd_ccr__no_tp53_odds_metas[,3])
# per_cancer_type_wgd_ccr__no_tp53_odds_metas[,4] <- as.numeric(per_cancer_type_wgd_ccr__no_tp53_odds_metas[,4])
# per_cancer_type_wgd_ccr__no_tp53_odds_metas[,5] <- as.numeric(per_cancer_type_wgd_ccr__no_tp53_odds_metas[,5])
# per_cancer_type_wgd_ccr__no_tp53_odds_metas[,6] <- as.numeric(per_cancer_type_wgd_ccr__no_tp53_odds_metas[,6])
# 
# 
# 
# 
# 
# 
# # =====================================================================================================================================
# # Fig 17: ddr deficiency frequencies per cancer type
# 
# 
# sum(metadata_included_ann_mono$BER)
# 
# length(unique(metadata_included_ann_mono_tibb$pathway_abb))
# 
# 
# ddr_prevalence_met_vs_pri_odds <- data.frame(pathway = character(13), odds_ratio = numeric(13), lower_ci_level = numeric(13), upper_ci_level = numeric(13), effect_size = numeric(13), p_value = numeric(13))
# 
# 
# 
# for (i in 1:(length(unique(metadata_included_ann_mono_tibb$pathway_abb)))-1) {
#   working_mat <- as.matrix(table(metadata_included_ann_mono_tibb$cohort, metadata_included_ann_mono_tibb$pathway_abb)[,c(1,i+1)])
#   hh <- oddsratio.wald(working_mat)
#   ddr_prevalence_met_vs_pri_odds[i,] <- c(as.vector(unique(metadata_included_ann_mono_tibb$pathway_abb)[i+1]), hh$measure[2,1], hh$measure[2,2], hh$measure[2,3], as.vector(cramerV(working_mat)), hh$p.value[2,3])
# }
# 
# ddr_prevalence_met_vs_pri_odds
# 
# 
# 
# 
# 
# ddr_prevalence_met_vs_pri_odds$pathway <- factor(ddr_prevalence_met_vs_pri_odds$pathway)
# ddr_prevalence_met_vs_pri_odds[,2] <- as.numeric(ddr_prevalence_met_vs_pri_odds[,2])
# ddr_prevalence_met_vs_pri_odds[,3] <- as.numeric(ddr_prevalence_met_vs_pri_odds[,3])
# ddr_prevalence_met_vs_pri_odds[,4] <- as.numeric(ddr_prevalence_met_vs_pri_odds[,4])
# ddr_prevalence_met_vs_pri_odds[,5] <- as.numeric(ddr_prevalence_met_vs_pri_odds[,5])
# ddr_prevalence_met_vs_pri_odds[,6] <- as.numeric(ddr_prevalence_met_vs_pri_odds[,6])
# 
# str(ddr_prevalence_met_vs_pri_odds)
# 
# 
# 
# 
# ddr_prevalence_met_vs_pri_odds_plot <- ddr_prevalence_met_vs_pri_odds %>% ggplot(aes(x = log2(odds_ratio), y = -log10(p_value))) +
#   geom_point(aes( size = effect_size)) +
#   # geom_text(x = 2, y = 9, label = "P.val < 0.01") +
#   # geom_errorbar( aes(ymin=lower_ci_level, ymax=upper_ci_level), lwd = 0.2) +
#   theme_classic() +
#   ggrepel::geom_text_repel(aes(label = pathway)) +
#   ggtitle("DDR Deficiency Prevalence \n") +
#   # ylab("-log(P.value) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   xlim(c(-3.5, 3.5)) +
#   scale_y_continuous(name = "-log10(P.value)", breaks = c(0,5,10,20,30,40,50), labels = c("0","5","10","20","30","40","50")) +
#   # coord_cartesian(ylim = c(0,10))+
#   xlab("log2(Odds ratio (Metastatic vs. Primary))\n") +
#   # ylab("Odds Ratio (Metastatic vs. Primary)") +
#   # scale_x_discrete(labels = d) +
#   geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
#   # guides(colour = guide_colorbar(reverse=T)) +
#   scale_size_continuous(name = c("Effect Size", "Odds ratio"))+
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   # labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# 
# 
# 
# 
# for (i in 1:2) {
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig17-odds-ratio-metastatic-vs-primary.png")
#     print(ddr_prevalence_met_vs_pri_odds_plot)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig17-odds-ratio-metastatic-vs-primary.pdf")
#     print(ddr_prevalence_met_vs_pri_odds_plot)
#     dev.off()
#   }
# }
# 
# 
# # Fig 18: DDR DEFICIENCY prevalence per cancer type
# 
# cancer_types <- unique(metadata_included_ann_mono_tibb$cancer_type)
# cancer_type <- "Urothelial cancer"
# 
# tmp_df <- metadata_included_ann_mono_tibb[metadata_included_ann_mono_tibb$cancer_type == cancer_type,]
# 
# 
# tmp_tabl <- as.data.frame.matrix(table(tmp_df$cohort, tmp_df$pathway_abb))
# 
# tmp_tabl[tmp_tabl < 4] <- NA 
# 
# tmp_tabl <- tmp_tabl[,colSums(is.na(tmp_tabl)) == 0]
# 
# class(tmp_tabl)
# is.data.frame(tmp_tabl)
# 
# oddsratio.wald(as.matrix(tmp_tabl[,c(1,3)]))
# 
# ddr_prevalence_met_vs_pri_odds_per_cancer <- data.frame(cancer_type = character(62), pathway = character(62), odds_ratio = numeric(62), lower_ci_level = numeric(62), upper_ci_level = numeric(62), effect_size = numeric(62), p_value = numeric(62))
# 
# 
# i <- 0
# for (cancer_type in unique(metadata_included_ann_mono_tibb$cancer_type)){
#   tmp_df <- metadata_included_ann_mono_tibb[metadata_included_ann_mono_tibb$cancer_type == cancer_type,]
#   tmp_tabl <- as.data.frame.matrix(table(tmp_df$cohort, tmp_df$pathway_abb))
#   tmp_tabl[tmp_tabl < 4] <- NA 
#   tmp_tabl <- tmp_tabl[,colSums(is.na(tmp_tabl)) == 0]
#   if (is.data.frame(tmp_tabl)){
#     if (ncol(tmp_tabl) >= 2 & colnames(tmp_tabl)[1] == "No deficiency"){
#       print(cancer_type)
#       print(tmp_tabl)
#       for (pathway in colnames(tmp_tabl)[-1]) {
#         i <- i + 1
#         working_mat <- as.matrix(tmp_tabl[,c("No deficiency",pathway)])
#         jj <- oddsratio.wald(working_mat)
#         ddr_prevalence_met_vs_pri_odds_per_cancer[i,] <- c(cancer_type, pathway, jj$measure[2,1], jj$measure[2,2], jj$measure[2,3], as.vector(cramerV(working_mat)), jj$p.value[2,3])
#       }
#     }
#   }
# }
# 
# 
# str(ddr_prevalence_met_vs_pri_odds_per_cancer)
# 
# ddr_prevalence_met_vs_pri_odds_per_cancer$cancer_type <- factor(ddr_prevalence_met_vs_pri_odds_per_cancer$cancer_type)
# ddr_prevalence_met_vs_pri_odds_per_cancer$pathway <- factor(ddr_prevalence_met_vs_pri_odds_per_cancer$pathway)
# ddr_prevalence_met_vs_pri_odds_per_cancer[,3] <- as.numeric(ddr_prevalence_met_vs_pri_odds_per_cancer[,3])
# ddr_prevalence_met_vs_pri_odds_per_cancer[,4] <- as.numeric(ddr_prevalence_met_vs_pri_odds_per_cancer[,4])
# ddr_prevalence_met_vs_pri_odds_per_cancer[,5] <- as.numeric(ddr_prevalence_met_vs_pri_odds_per_cancer[,5])
# ddr_prevalence_met_vs_pri_odds_per_cancer[,6] <- as.numeric(ddr_prevalence_met_vs_pri_odds_per_cancer[,6])
# ddr_prevalence_met_vs_pri_odds_per_cancer[,7] <- as.numeric(ddr_prevalence_met_vs_pri_odds_per_cancer[,7])
# 
# 
# 
# 
# 
# 
# ddr_prevalence_met_vs_pri_odds_plot_per_cancer <- ddr_prevalence_met_vs_pri_odds_per_cancer %>% ggplot(aes(x = log2(odds_ratio), y = -log10(p_value))) + facet_wrap(~cancer_type) +
#   geom_point(aes( size = effect_size, color = pathway)) +
#   ggrepel::geom_text_repel(aes(label = pathway), size = 2, fontface = "bold") +
#   # geom_text(x = 2, y = 9, label = "P.val < 0.01") +
#   # geom_errorbar( aes(ymin=lower_ci_level, ymax=upper_ci_level), lwd = 0.2) +
#   theme_classic() +
#   ggtitle("DDR Deficiency Prevalence \n") +
#   # ylab("-log(P.value) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   xlim(c(-4, 4)) +
#   scale_y_continuous(name = "-log10(P.value)", breaks = c(0,1,2,3,4,5,6), labels = c("0","1","2","3","4","5","6"), limits = c(0,10)) +
#   # coord_cartesian(ylim = c(0,10))+
#   xlab("log2(Odds ratio (Metastatic vs. Primary))\n") +
#   # ylab("Odds Ratio (Metastatic vs. Primary)") +
#   # scale_x_discrete(labels = d) +
#   geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
#   # guides(colour = guide_colorbar(reverse=T)) +
#   scale_size_continuous(name = c("Effect Size", "Odds ratio"))+
#   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 15)) +
#   theme(axis.title.y = element_text(face = "italic", size = 15)) +
#   # labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# 
# 
# 
# for (i in 1:2) {
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig18-odds-ratio-metastatic-vs-primary-per-cancer.png", height = 960, width = 960)
#     print(ddr_prevalence_met_vs_pri_odds_plot_per_cancer)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig18-odds-ratio-metastatic-vs-primary-per-cancer.pdf", height = 14, width = 14)
#     print(ddr_prevalence_met_vs_pri_odds_plot_per_cancer)
#     dev.off()
#   }
# }
# 
# 
# 
# # =====================================================================================================================================
# 
# # producing ddr-wgd plot per cancer type for other pathways other than ccr
# 
# 
# ddr_wgd_odds_per_cancer <- data.frame(cancer_type = character(75), pathway = character(75), odds_ratio = numeric(75), lower_ci_level = numeric(75), upper_ci_level = numeric(75), effect_size = numeric(75), p_value = numeric(75))
# 
# 
# i <- 0
# metadata_included_ann_mono_tibb_primary
# for (cancer_type in unique(metadata_included_ann_mono_tibb_primary$cancer_type)){
#   tmp_df <- metadata_included_ann_mono_tibb_primary[metadata_included_ann_mono_tibb_primary$cancer_type == cancer_type,]
#   tmp_tabl <- as.data.frame.matrix(table(tmp_df$whole_genome_duplication, tmp_df$pathway_abb))
#   tmp_tabl[tmp_tabl < 4] <- NA 
#   tmp_tabl <- tmp_tabl[,colSums(is.na(tmp_tabl)) == 0]
#   if (is.data.frame(tmp_tabl)){
#     if (ncol(tmp_tabl) >= 2 & colnames(tmp_tabl)[1] == "No deficiency"){
#       print(cancer_type)
#       print(tmp_tabl)
#       for (pathway in colnames(tmp_tabl)[-1]) {
#         i <- i + 1
#         working_mat <- as.matrix(tmp_tabl[,c("No deficiency",pathway)])
#         jj <- oddsratio.wald(working_mat)
#         ddr_wgd_odds_per_cancer[i,] <- c(cancer_type, pathway, jj$measure[2,1], jj$measure[2,2], jj$measure[2,3], as.vector(cramerV(working_mat)), jj$p.value[2,3])
#       }
#     }
#   }
# }
# 
# 
# 
# ddr_wgd_odds_per_cancer <- ddr_wgd_odds_per_cancer[ddr_wgd_odds_per_cancer$pathway != "CCR" & ddr_wgd_odds_per_cancer$p_value < 0.01,]
# 
# 
# ddr_wgd_odds_per_cancer$cancer_type <- factor(ddr_wgd_odds_per_cancer$cancer_type)
# ddr_wgd_odds_per_cancer$pathway <- factor(ddr_wgd_odds_per_cancer$pathway)
# ddr_wgd_odds_per_cancer[,3] <- as.numeric(ddr_wgd_odds_per_cancer[,3])
# ddr_wgd_odds_per_cancer[,4] <- as.numeric(ddr_wgd_odds_per_cancer[,4])
# ddr_wgd_odds_per_cancer[,5] <- as.numeric(ddr_wgd_odds_per_cancer[,5])
# ddr_wgd_odds_per_cancer[,6] <- as.numeric(ddr_wgd_odds_per_cancer[,6])
# ddr_wgd_odds_per_cancer[,7] <- as.numeric(ddr_wgd_odds_per_cancer[,7])
# 
# 
# 
# 
# text_dat <- data.frame(x_pos = c(1.2,5.2,1.2,1.2), y_pos = c(6,6,6,6), lab = c("P.val < 0.01", "P.val < 0.01", "P.val < 0.01", "P.val < 0.01"), cancer_type = c("Breast cancer", "Non small cell lung cancer", "Ovarian cancer", "Skin melanoma"))
# 
# 
# odds_wgd_all_plot_combined <- ddr_wgd_odds_per_cancer %>% ggplot(aes(x = pathway, y = odds_ratio)) + facet_wrap(~ cancer_type, scales = "free") +
#   geom_point(aes(size = effect_size)) +
#   # geom_errorbar( aes(ymin=lower_ci_level, ymax=upper_ci_level), lwd = 0.2) +
#   theme_classic() +
#   geom_text(data = text_dat, aes(x = x_pos, y = y_pos, label = lab)) +
#   ggtitle("DDR-deficiency and WGD across cancer types (Metastatic) \n") +
#   # ylab("Genomic Mutations/TMB (log2) \n") +
#   # stat_compare_means(method = "t.test", label = "p.signif") +
#   ylim(0, 7) +
#   # coord_cartesian(ylim = c(0,10))+
#   xlab("DDR Pathway Mutated \n") +
#   geom_hline(yintercept = 1, linetype = "dashed") +
#   ylab("Odds Ratio of WGD(+)/WGD(-) _  \n (compared to DDR-proficient)") +
#   # scale_x_discrete(labels = d) +
#   scale_size_continuous(name = "Effect Size")+
#   scale_color_discrete(name = "Cohort")+
#   theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) +
#   theme(axis.title.x = element_text(face = "italic", size = 13)) +
#   theme(axis.title.y = element_text(face = "italic", size = 13)) +
#   # labs(fill = "Cancer Type") +
#   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# 
# 
# 
# 
# 
# 
# for (i in 1:2) {
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig16-odds-wgd-ddr-plot-metastatic.png")
#     print(odds_wgd_all_plot_combined)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig16-odds-wgd-ddr-plot-metastatic.pdf")
#     print(odds_wgd_all_plot_combined)
#     dev.off()
#   }
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##### Produce overview plots/figures
# 
# # metadata_included_ann$dna_repair_deficiency <- factor(metadata_included_ann$dna_repair_deficiency, levels = names(table(metadata_included$dna_repair_deficiency))[c(11,10,1:9,12,13)])
# # metadata_included_ann$cohort <- factor(metadata_included_ann$cohort, levels = c("PCAWG", "HMF"))
# # metadata_included_ann$stage[metadata_included_ann$is_metastatic] <- "Metastatic"
# # metadata_included_ann$stage[!metadata_included_ann$is_metastatic] <- "Primary"
# # metadata_included_ann$stage <- factor(metadata_included_ann$stage, levels = c("Primary", "Metastatic"))
# # 
# # 
# # d <- paste0(names(table(metadata_included_ann$dna_repair_deficiency)), " (n=", as.vector(table(metadata_included_ann$dna_repair_deficiency)), ")")
# # 
# # my_comparisons <- list(c("Primary", "Metastatic"))
# # 
# # tmb_plot <- metadata_included_ann %>% ggplot(aes(x = dna_repair_deficiency, y = log2(tmb), color = stage)) +
# #   geom_boxplot(outlier.shape = NA) +
# #   # stat_compare_means(method = "t.test") +
# #   # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins, 
# #   #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
# #   # scale_color_manual(values = c("blue", "green")) +
# #   # geom_jitter(color="black", size=0.4, alpha=0.5) +
# #   geom_point(position=position_jitterdodge(), size = 0.4, aes(color = stage), alpha = 0.5) +
# #   scale_color_manual(values = c("#ff0101", "#010dff")) +
# #   # stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = F) +
# #   theme_bw() +
# #   ggtitle("DDR deficiency Effects on TMB  \n") +
# #   ylab("Genomic Mutations/TMB (log2) \n") +
# #   stat_compare_means(method = "wilcox.test", label = "p.signif") +
# #   # ylim(0, 12) +
# #   xlab("DNA repair pathway mutated \n") +
# #   scale_x_discrete(labels = d) +
# #   theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
# #   theme(axis.title.x = element_text(face = "italic", size = 15)) +
# #   theme(axis.title.y = element_text(face = "italic", size = 15)) +
# #   labs(fill = "Cancer Type") +
# #   theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
# #   theme(plot.margin = unit(c(1,1,1,1), "cm")) 
# # 
# # 
# # wilcox.test(metadata_included_ann[metadata_included_ann$stage == "Primary" & metadata_included_ann$dna_repair_deficiency == "NHEJ","tmb"],
# #             metadata_included_ann[metadata_included_ann$stage == "Metastatic" & metadata_included_ann$dna_repair_deficiency == "NHEJ","tmb"])
# # 
# # 
# # 
# # for (i in 1:2){
# #   if (i == 1) {
# #     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/ddrd-tmb-boxplot.png", width = 960, height = 960)
# #     print(tmb_plot)
# #     dev.off()
# #   }
# #   if (i == 2) {
# #     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/ddrd-tmb-boxplot.pdf", width = 14, height = 14)
# #     print(tmb_plot)
# #     dev.off()
# #   }
# # }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
