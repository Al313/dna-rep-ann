# In this script I'll get the SV counts for the samples



# library(vcfR)
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(gridExtra)



if (dir.exists("/hpc/cuppen/")){
  hmf_meta <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/metadata_whitelisted.tsv", sep = "\t", header = T, stringsAsFactors = F)
} else {
  hmf_meta <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/metadata_whitelisted.tsv", sep = "\t", header = T, stringsAsFactors = F)
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



(dir.exists("/hpc/cuppen/")){
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





metadata_included_ann_mono <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/metadata-included-ddrd-annotated-monoallelic.txt",
                                       sep = "\t", header = T, stringsAsFactors = F)

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




# 
# 
# simple_sv_events <- c("Deletion (> 100kb)", "Deletion (< 100kb)", "Duplication (> 100kb)", "Duplication (< 100kb)", "Inversion")
# sv_count_matrix <- matrix(nrow = nrow(metadata_included), ncol = length(simple_sv_events))
# colnames(sv_count_matrix) <- simple_sv_events
# rownames(sv_count_matrix) <- metadata_included$sample_id
# 
# 
# 
# for (i in 1:nrow(metadata_included)) {
#   print(i)
#   file_there <- T
#   if (metadata_included$cohort[i] == "HMF"){
#     
#     
#     sample_id <- metadata_included$sample_id[i]
#     print(sample_id)
#     set_name <- hmf_meta$setName[hmf_meta$sampleId == sample_id]
# 
#     if (dir.exists("/hpc/cuppen/")){
#       sv_count_file <- paste0("/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/",
#                               set_name, "/linx14/", sample_id, ".linx.vis_sv_data.tsv")
#       if (file.exists(sv_count_file)){
#         vcf <- read.csv(file = sv_count_file, header = T, sep = "\t", stringsAsFactors = F)
#       } else {
#         print(paste0(sv_count_file, " does not exist!!"))
#         file_there <- F
#         }
#     } else {
#       sv_count_file <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/",
#                               set_name, "/linx14/", sample_id, ".linx.vis_sv_data.tsv")
#       
#       if (file.exists(sv_count_file)){
#         vcf <- read.csv(file = sv_count_file, header = T, sep = "\t", stringsAsFactors = F)
#       } else {
#         print(paste0(sv_count_file, " does not exist!!"))
#         file_there <- F
#       }
#     }
# 
#     if (file_there){
#       
#       vcf <- mutate(vcf, length = abs(PosEnd - PosStart))
#       
#       
#       
#       vcf1 <- vcf[vcf$ResolvedType == "DEL" & vcf$length > 100000,]
#       sv_count_matrix[i,1] <- nrow(vcf1)
#       
#       vcf2 <- vcf[vcf$ResolvedType == "DEL" & vcf$length < 100000,]
#       sv_count_matrix[i,2] <- nrow(vcf2)
#       
#       vcf3 <- vcf[vcf$ResolvedType == "DUP" & vcf$length > 100000,]
#       sv_count_matrix[i,3] <- nrow(vcf3)
#       
#       vcf4 <- vcf[vcf$ResolvedType == "DUP" & vcf$length < 100000,]
#       sv_count_matrix[i,4] <- nrow(vcf4)
#       
#       vcf5 <- vcf[vcf$ResolvedType == "INV",]
#       sv_count_matrix[i,5] <- nrow(vcf5)
#       
#     } else {
#       sv_count_matrix[i,] <- NA
#     }
#     
# 
#   } else if (metadata_included$cohort[i] == "PCAWG") {
#     
#     print(metadata_included$sample_id[i])
#     if (dir.exists("/hpc/cuppen/")){
#       sv_count_file <- paste0("/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
#                               metadata_included$patient_id[i], "-from-jar/linx14/", metadata_included$patient_id[i], "T.linx.vis_sv_data.tsv")
#       
#       if (file.exists(sv_count_file)){
#         vcf <- read.csv(file = sv_count_file, header = T, sep = "\t", stringsAsFactors = F)
#       } else {
#         print(paste0(sv_count_file, " does not exist!!"))
#         file_there <- F
#       }
#       
#       } else {
#         sv_count_file <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
#                                 metadata_included$patient_id[i], "-from-jar/linx14/", metadata_included$patient_id[i], "T.linx.vis_sv_data.tsv")
#         if (file.exists(sv_count_file)){
#           vcf <- read.csv(file = sv_count_file, header = T, sep = "\t", stringsAsFactors = F)
#         } else {
#           print(paste0(sv_count_file, " does not exist!!"))
#           file_there <- F
#         }
#       }
# 
#         if (file_there){
#           vcf <- mutate(vcf, length = abs(PosEnd - PosStart))
#   
#   
#   
#           vcf1 <- vcf[vcf$ResolvedType == "DEL" & vcf$length > 100000,]
#           sv_count_matrix[i,1] <- nrow(vcf1)
#   
#           vcf2 <- vcf[vcf$ResolvedType == "DEL" & vcf$length < 100000,]
#           sv_count_matrix[i,2] <- nrow(vcf2)
#   
#           vcf3 <- vcf[vcf$ResolvedType == "DUP" & vcf$length > 100000,]
#           sv_count_matrix[i,3] <- nrow(vcf3)
#   
#           vcf4 <- vcf[vcf$ResolvedType == "DUP" & vcf$length < 100000,]
#           sv_count_matrix[i,4] <- nrow(vcf4)
#   
#           vcf5 <- vcf[vcf$ResolvedType == "INV",]
#           sv_count_matrix[i,5] <- nrow(vcf5)
#         
#         } else {
#           sv_count_matrix[i,] <- NA
#         }
# 
# 
#   }
# }
# 
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(sv_count_matrix, file = paste0(wd,"r-objects/simple_sv_count.rds"))
# } else {
#   saveRDS(sv_count_matrix, file = paste0(wd,"r-objects/simple_sv_count.rds"))
# }
# 
# 



if (dir.exists("/hpc/cuppen/")){
  sv_count_matrix <- readRDS(file = paste0(wd,"r-objects/simple_sv_count.rds"))
} else {
  sv_count_matrix <- readRDS(file = paste0(wd,"r-objects/simple_sv_count.rds"))
}

sv_count_df <- as.data.frame(sv_count_matrix)

sv_count_df %<>% mutate(sample_id = rownames(sv_count_df)) %>% relocate(where(is.character), .before = where(is.numeric))

rownames(sv_count_df) <- 1:nrow(sv_count_df)

metadata_included_ann_mono_tibb <- merge(metadata_included_ann_mono_tibb, sv_count_df, by = "sample_id")

sv_event_level <- colnames(metadata_included_ann_mono_tibb)[30:33][c(2,1,4,3)]
metadata_included_ann_mono_tibb <- metadata_included_ann_mono_tibb[,-34]

metadata_included_ann_mono_tibb_tibb <- tidyr::gather(metadata_included_ann_mono_tibb, colnames(metadata_included_ann_mono_tibb[,30:33]), key = "sv_event", value = "sv_counts")
metadata_included_ann_mono_tibb_tibb$sv_event <- factor(metadata_included_ann_mono_tibb_tibb$sv_event, level = sv_event_level)


metadata_included_ann_mono_tibb_tibb_prim <- metadata_included_ann_mono_tibb_tibb[!(metadata_included_ann_mono_tibb_tibb$is_metastatic),]
metadata_included_ann_mono_tibb_tibb_metas <- metadata_included_ann_mono_tibb_tibb[metadata_included_ann_mono_tibb_tibb$is_metastatic,]



## 
metadata_included_ann_mono_tibb2 <- merge(metadata_included_ann_mono_tibb2, sv_count_df, by = "sample_id")

sv_event_level <- colnames(metadata_included_ann_mono_tibb2)[30:33][c(2,1,4,3)]
metadata_included_ann_mono_tibb2 <- metadata_included_ann_mono_tibb2[,-34]

metadata_included_ann_mono_tibb_tibb2 <- tidyr::gather(metadata_included_ann_mono_tibb2, colnames(metadata_included_ann_mono_tibb2[,30:33]), key = "sv_event", value = "sv_counts")
metadata_included_ann_mono_tibb_tibb2$sv_event <- factor(metadata_included_ann_mono_tibb_tibb2$sv_event, level = sv_event_level)


metadata_included_ann_mono_tibb_tibb2_prim <- metadata_included_ann_mono_tibb_tibb2[!(metadata_included_ann_mono_tibb_tibb2$is_metastatic),]
metadata_included_ann_mono_tibb_tibb2_metas <- metadata_included_ann_mono_tibb_tibb2[metadata_included_ann_mono_tibb_tibb2$is_metastatic,]

##
sv_plot_mono_primary <- metadata_included_ann_mono_tibb_tibb2_prim %>% ggplot(aes(x = sv_event, y = log2(sv_counts), fill = pathway_abb)) +
  geom_boxplot(outlier.shape = NA) +
  # stat_compare_means(method = "t.test") +
  # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins,
  #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
  # scale_color_manual(values = c("blue", "green")) +
  # geom_jitter(color="black", size=0.4, alpha=0.5) +
  # geom_point(position=position_jitterdodge(), size = 0.4, alpha = 0.5) +
  # scale_color_manual(values = c("#ff0101", "#010dff")) +
  stat_compare_means(method = "anova") +
  theme_bw() +
  ggtitle("DDR deficiency Effects on SV Event Frequency in Primary Cancers \n") +
  ylab("SV Count (log2) \n") +
  # stat_compare_means(method = "t.test", label = "p.signif") +
  ylim(0, 12) +
  xlab("SV Event\n") +
  # scale_x_discrete(labels = d) +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.title.x = element_text(face = "italic", size = 15)) +
  theme(axis.title.y = element_text(face = "italic", size = 15)) +
  labs(fill = "DDR Mutated Pathway") +
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))





for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig19-ddrd-sv_plot_mono_primary.png", width = 960, height = 960)
    print(sv_plot_mono_primary)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig19-ddrd-sv_plot_mono_primary.pdf", width = 14, height = 14)
    print(sv_plot_mono_primary)
    dev.off()
  }
}




m <- as.data.frame(table(metadata_included_ann_mono_tibb_tibb2_prim$cancer_type, metadata_included_ann_mono_tibb_tibb2_prim$pathway_abb))
y <- m[m$Freq >= 5,]



metadata_included_ann_mono_tibb_tibb2_prim$tti <- FALSE
for (i in 1:nrow(y)) {
  metadata_included_ann_mono_tibb_tibb2_prim$tti[metadata_included_ann_mono_tibb_tibb2_prim$cancer_type == y[i,"Var1"] & metadata_included_ann_mono_tibb_tibb2_prim$pathway_abb == y[i,"Var2"]] <- TRUE
}




sv_plot_mono_primary_per_cancer_type <- metadata_included_ann_mono_tibb_tibb2_prim[metadata_included_ann_mono_tibb_tibb2_prim$tti,] %>% ggplot(aes(x = sv_event, y = log2(sv_counts), fill = pathway_abb)) + facet_wrap(~cancer_type) +
  geom_boxplot(outlier.shape = NA) +
  # stat_compare_means(method = "t.test") +
  # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins,
  #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
  # scale_color_manual(values = c("blue", "green")) +
  # geom_jitter(color="black", size=0.4, alpha=0.5) +
  # geom_point(position=position_jitterdodge(), size = 0.4, alpha = 0.5) +
  # scale_color_manual(values = c("#ff0101", "#010dff")) +
  stat_compare_means(method = "anova") +
  theme_bw() +
  ggtitle("DDR deficiency Effects on SV Event Frequency in Primary Cancers Per Cancer Type \n") +
  ylab("SV Count (log2) \n") +
  # stat_compare_means(method = "t.test", label = "p.signif") +
  ylim(0, 12) +
  xlab("SV Event \n") +
  # scale_x_discrete(labels = d) +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.title.x = element_text(face = "italic", size = 15)) +
  theme(axis.title.y = element_text(face = "italic", size = 15)) +
  labs(fill = "DDR Mutated Pathway") +
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))



for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig19-ddrd-sv_plot_mono_primary_per_cancer_type.png", width = 960, height = 960)
    print(sv_plot_mono_primary_per_cancer_type)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig19-ddrd-sv_plot_mono_primary_per_cancer_type.pdf", width = 21, height = 21)
    print(sv_plot_mono_primary_per_cancer_type)
    dev.off()
  }
}









sv_plot_mono_metastatic <- metadata_included_ann_mono_tibb_tibb2_metas %>% ggplot(aes(x = sv_event, y = log2(sv_counts), fill = pathway_abb)) +
  geom_boxplot(outlier.shape = NA) +
  # stat_compare_means(method = "t.test") +
  # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins,
  #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
  # scale_color_manual(values = c("blue", "green")) +
  # geom_jitter(color="black", size=0.4, alpha=0.5) +
  # geom_point(position=position_jitterdodge(), size = 0.4, alpha = 0.5) +
  # scale_color_manual(values = c("#ff0101", "#010dff")) +
  stat_compare_means(method = "anova") +
  theme_bw() +
  ggtitle("DDR deficiency Effects on SV Event Frequency in Metastatic Cancers\n") +
  ylab("SV Count (log2) \n") +
  # stat_compare_means(method = "t.test", label = "p.signif") +
  ylim(0, 12) +
  xlab("SV Event \n") +
  # scale_x_discrete(labels = d) +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.title.x = element_text(face = "italic", size = 15)) +
  theme(axis.title.y = element_text(face = "italic", size = 15)) +
  labs(fill = "DDR Mutated Pathway") +
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))





for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig20-ddrd-sv_plot_mono_metastatic.png", width = 960, height = 960)
    print(sv_plot_mono_metastatic)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig20-ddrd-sv_plot_mono_metastatic.pdf", width = 14, height = 14)
    print(sv_plot_mono_metastatic)
    dev.off()
  }
}





m <- as.data.frame(table(metadata_included_ann_mono_tibb_tibb2_metas$cancer_type, metadata_included_ann_mono_tibb_tibb2_metas$pathway_abb))
y <- m[m$Freq >= 5,]



metadata_included_ann_mono_tibb_tibb2_metas$tti <- FALSE
for (i in 1:nrow(y)) {
  metadata_included_ann_mono_tibb_tibb2_metas$tti[metadata_included_ann_mono_tibb_tibb2_metas$cancer_type == y[i,"Var1"] & metadata_included_ann_mono_tibb_tibb2_metas$pathway_abb == y[i,"Var2"]] <- TRUE
}




sv_plot_mono_metastatic_per_cancer_type <- metadata_included_ann_mono_tibb_tibb2_metas[metadata_included_ann_mono_tibb_tibb2_metas$tti,] %>% ggplot(aes(x = sv_event, y = log2(sv_counts), fill = pathway_abb)) + facet_wrap(~cancer_type) +
  geom_boxplot(outlier.shape = NA) +
  # stat_compare_means(method = "t.test") +
  # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins,
  #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
  # scale_color_manual(values = c("blue", "green")) +
  # geom_jitter(color="black", size=0.4, alpha=0.5) +
  # geom_point(position=position_jitterdodge(), size = 0.4, alpha = 0.5) +
  # scale_color_manual(values = c("#ff0101", "#010dff")) +
  # stat_compare_means(method = "anova") +
  theme_bw() +
  ggtitle("DDR deficiency Effects on SV Event Frequency in Metastatic Cancers Per Cancer Type\n") +
  ylab("SV Count (log2) \n") +
  # stat_compare_means(method = "t.test", label = "p.signif") +
  ylim(0, 12) +
  xlab("SV Event \n") +
  # scale_x_discrete(labels = d) +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.title.x = element_text(face = "italic", size = 15)) +
  theme(axis.title.y = element_text(face = "italic", size = 15)) +
  labs(fill = "DDR Mutated Pathway") +
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))



for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/png/fig20-ddrd-sv_plot_mono_metastatic_per_cancer_type.png", width = 960, height = 960)
    print(sv_plot_mono_metastatic_per_cancer_type)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs/pdf/fig20-ddrd-sv_plot_mono_metastatic_per_cancer_type.pdf", width = 21, height = 21)
    print(sv_plot_mono_metastatic_per_cancer_type)
    dev.off()
  }
}


