
### Purpose

# In this script I investigate the dynamics of mutational processes using mutational signature and mutational timing information.

### Loading required libraries

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(magrittr)
library(rcompanion)
library(epitools)
library(stringr)
library(tidyr)
library(plyr)
library(miscTools)
library(effsize)
library(gggibbous)




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


# if (dir.exists("/hpc/cuppen/")){
#   sig_cont <- readRDS("/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/contribs.fit_lsq.raw_merged.rds")
# } else {
#   sig_cont <- readRDS("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/contribs.fit_lsq.raw_merged.rds")
# }
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   sig_cont_by_tissue <- readRDS("/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/fit_lsq.by_tissue_type/denovo_contribs.lsq.by_tissue_type.rds")
# } else {
#   sig_cont_by_tissue <- readRDS("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/fit_lsq.by_tissue_type/denovo_contribs.lsq.by_tissue_type.rds")
# }


global_timing_info_com <- read.csv(file = paste0(wd, "r-objects/global-timing-with-metadata.txt.gz"), header = T, sep = "\t", stringsAsFactors = F)


global_timing_info_com <- global_timing_info_com[global_timing_info_com$sample_id != "DO217817",]


if (dir.exists("/hpc/cuppen/")){
  wgd_timing_df <- read.csv(file = paste0(wd, "r-objects/wgd-timing.txt"), sep = "\t", header = T, stringsAsFactors = F)
} else {
  wgd_timing_df <- read.csv(file = paste0(wd, "r-objects/wgd-timing.txt"), sep = "\t", header = T, stringsAsFactors = F)
}

wgd_timing_df <- merge(wgd_timing_df, metadata_included[,c("sample_id", "cancer_type", "whole_genome_duplication", "tmb", "is_hypermutated")], by.x = "sample", by.y = "sample_id")
# ================================================================================================================================================
# ================================================================================================================================================
# WGD

# ************************************************************************************************************************************************
# Fig1: WGD frequency between metastatic and primary per cancer type

cancer_types <- unique(wgd_timing_df$cancer_type)
length(cancer_types)

per_cancer_type_wgd_metas_vs_prim_timing_annot <- data.frame(cancaer_type = character(59), odds_ratio = numeric(59), lower_ci_level = numeric(59), upper_ci_level = numeric(59), effect_size = numeric(59), p_value = numeric(59), annotation = "MutationTimingR")


for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  working_df <- wgd_timing_df[wgd_timing_df$cancer_type == cancer_type,]
  working_mat <- as.matrix(table(working_df$is_metastatic, working_df$isWGD))
  
  if (nrow(working_mat) > 1 & ncol(working_mat) > 1){
    if (cancer_type != "Thyroid cancer"){
      ff <- oddsratio.wald(working_mat)
      per_cancer_type_wgd_metas_vs_prim_timing_annot[i,1:6] <- c(cancer_type, ff$measure[2,1], ff$measure[2,2], ff$measure[2,3], as.vector(cramerV(working_mat)), ff$p.value[2,3])
    } else {
      working_mat <- working_mat + 1
      ff <- oddsratio.wald(working_mat)
      per_cancer_type_wgd_metas_vs_prim_timing_annot[i,1:6] <- c(cancer_type, ff$measure[2,1], ff$measure[2,2], ff$measure[2,3], as.vector(cramerV(working_mat)), ff$p.value[2,3])
    }
  } else {
    per_cancer_type_wgd_metas_vs_prim_timing_annot[i,1:6] <- c(cancer_type, NA, NA, NA, NA, NA)
  }
}



per_cancer_type_wgd_metas_vs_prim_metadata_annot <- data.frame(cancaer_type = character(59), odds_ratio = numeric(59), lower_ci_level = numeric(59), upper_ci_level = numeric(59), effect_size = numeric(59), p_value = numeric(59), annotation = "Metadata")


for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  working_df <- metadata_included[metadata_included$cancer_type == cancer_type,]
  working_mat <- as.matrix(table(working_df$is_metastatic, working_df$whole_genome_duplication))
  print(cancer_type)
  print(working_mat)
  
  if (nrow(working_mat) > 1 & ncol(working_mat) > 1){
    ff <- oddsratio.wald(working_mat)
    per_cancer_type_wgd_metas_vs_prim_metadata_annot[i,1:6] <- c(cancer_type, ff$measure[2,1], ff$measure[2,2], ff$measure[2,3], as.vector(cramerV(working_mat)), ff$p.value[2,3])
  } else {
    per_cancer_type_wgd_metas_vs_prim_metadata_annot[i,1:6] <- c(cancer_type, NA, NA, NA, NA, NA)
  }
}


per_cancer_type_wgd_metas_vs_prim_combined <- rbind(per_cancer_type_wgd_metas_vs_prim_timing_annot, per_cancer_type_wgd_metas_vs_prim_metadata_annot)
str(per_cancer_type_wgd_metas_vs_prim_combined)

per_cancer_type_wgd_metas_vs_prim_combined$cancaer_type <- factor(per_cancer_type_wgd_metas_vs_prim_combined$cancaer_type, levels = cancer_types)
per_cancer_type_wgd_metas_vs_prim_combined$odds_ratio <- as.numeric(per_cancer_type_wgd_metas_vs_prim_combined$odds_ratio)
per_cancer_type_wgd_metas_vs_prim_combined$lower_ci_level <- as.numeric(per_cancer_type_wgd_metas_vs_prim_combined$lower_ci_level)
per_cancer_type_wgd_metas_vs_prim_combined$upper_ci_level <- as.numeric(per_cancer_type_wgd_metas_vs_prim_combined$upper_ci_level)
per_cancer_type_wgd_metas_vs_prim_combined$effect_size <- as.numeric(per_cancer_type_wgd_metas_vs_prim_combined$effect_size)
per_cancer_type_wgd_metas_vs_prim_combined$p_value <- as.numeric(per_cancer_type_wgd_metas_vs_prim_combined$p_value)
per_cancer_type_wgd_metas_vs_prim_combined$annotation <- factor(per_cancer_type_wgd_metas_vs_prim_combined$annotation)

per_cancer_type_wgd_metas_vs_prim_combined <- per_cancer_type_wgd_metas_vs_prim_combined[!is.na(per_cancer_type_wgd_metas_vs_prim_combined$odds_ratio),]


s




wgd_prevalence_met_vs_pri_odds_plot <- per_cancer_type_wgd_metas_vs_prim_combined %>% ggplot(aes(x = log2(odds_ratio), y = -log10(p_value))) +
  geom_point(aes( size = effect_size, color = annotation, shape = annotation)) +
  # geom_text(x = 2, y = 9, label = "P.val < 0.01") +
  # geom_errorbar( aes(ymin=lower_ci_level, ymax=upper_ci_level), lwd = 0.2) +
  theme_classic() +
  gghighlight() +
  facet_wrap(~ annotation) +
  ggrepel::geom_text_repel(data = per_cancer_type_wgd_metas_vs_prim_combined[-log10(per_cancer_type_wgd_metas_vs_prim_combined$p_value) > 2,], aes(label = cancaer_type), size = 3, max.overlaps = 20) +
  ggtitle("WGD Frequency \n") +
  # ylab("-log(P.value) \n") +
  # stat_compare_means(method = "t.test", label = "p.signif") +
  xlim(c(-3.5, 8)) +
  scale_y_continuous(name = "-log10(P.value)", breaks = c(0,2,10,20,30,40), labels = c("0","2","10","20","30","40")) +
  # coord_cartesian(ylim = c(0,10))+
  xlab("log2(Odds ratio (Metastatic vs. Primary))\n") +
  # ylab("Odds Ratio (Metastatic vs. Primary)") +
  # scale_x_discrete(labels = d) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  # guides(colour = guide_colorbar(reverse=T)) +
  scale_size_continuous(name = c("Effect Size", "Odds ratio"))+
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.title.x = element_text(face = "italic", size = 15)) +
  theme(axis.title.y = element_text(face = "italic", size = 15)) +
  # labs(fill = "Cancer Type") +
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.6, face = "bold")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))


for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig1-odds-wgd-metas-vs-prim.png", height = 920, width = 1160)
    print(wgd_prevalence_met_vs_pri_odds_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig1-odds-wgd-metas-vs-prim.pdf", height = 14, width = 17.5)
    print(wgd_prevalence_met_vs_pri_odds_plot)
    dev.off()
  }
}






shared <- intersect(per_cancer_type_wgd_metas_vs_prim_timing_annot[!(is.na(per_cancer_type_wgd_metas_vs_prim_timing_annot$effect_size)),1], per_cancer_type_wgd_metas_vs_prim_metadata_annot[!(is.na(per_cancer_type_wgd_metas_vs_prim_metadata_annot$effect_size)),1])



for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig1-odds-wgd-metas-vs-prim-effectsize.png")
    plot(per_cancer_type_wgd_metas_vs_prim_timing_annot$effect_size[per_cancer_type_wgd_metas_vs_prim_timing_annot$cancaer_type %in% shared], per_cancer_type_wgd_metas_vs_prim_metadata_annot$effect_size[per_cancer_type_wgd_metas_vs_prim_metadata_annot$cancaer_type %in% shared], xlab = "MutationTimingR annotation (effect size)", ylab = "Metadata annotation (effect size)")
    lines(x = c(0,0.5), y = c(0,0.5), color = "red")
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig1-odds-wgd-metas-vs-prim-effectsize.pdf")
    plot(per_cancer_type_wgd_metas_vs_prim_timing_annot$effect_size[per_cancer_type_wgd_metas_vs_prim_timing_annot$cancaer_type %in% shared], per_cancer_type_wgd_metas_vs_prim_metadata_annot$effect_size[per_cancer_type_wgd_metas_vs_prim_metadata_annot$cancaer_type %in% shared], xlab = "MutationTimingR annotation (effect size)", ylab = "Metadata annotation (effect size)")
    lines(x = c(0,0.5), y = c(0,0.5), color = "red")
    dev.off()
  }
}


# ************************************************************************************************************************************************
# Fig2: density plot of wgd timing



for (cancer_type in cancer_types) {
  # print(cancer_type)
  if (nrow(wgd_timing_df[wgd_timing_df$cancer_type == cancer_type & wgd_timing_df$isWGD & !is.na(wgd_timing_df$molecular_timing),]) < 5){
    wgd_timing_df <- wgd_timing_df[wgd_timing_df$cancer_type != cancer_type,]
  }
}




wgd_timing_df$cancer_type <- factor(wgd_timing_df$cancer_type, levels = unique(wgd_timing_df$cancer_type))

density_timing_plot <- wgd_timing_df[!(is.na(wgd_timing_df$molecular_timing)) & wgd_timing_df$isWGD,] %>% ggplot(aes(x = molecular_timing, fill = is_metastatic, ..density..)) + #facet_wrap(~ cancer_type, scale = "free_y") +
  geom_density(alpha = 0.5, adjust = 0.5) +
  ggtitle("density")

# count_timing_per_plot <- wgd_timing_df[!(is.na(wgd_timing_df$molecular_timing)) & wgd_timing_df$isWGD,] %>% ggplot(aes(x = molecular_timing, fill = is_metastatic, ..count..)) + #facet_wrap(~ cancer_type, scale = "free_y") +
#   geom_density(alpha = 0.5, adjust = 0.5, stat = "density", position = "identity") +
#   ggtitle("count")


count_timing_plot <- wgd_timing_df[!(is.na(wgd_timing_df$molecular_timing)) & wgd_timing_df$isWGD,]  %>% ggplot(aes(x = molecular_timing, fill = is_metastatic)) + #facet_wrap(~ cancer_type, scale = "free_y") +
  geom_histogram(alpha = 0.5, bins = 20, position="identity")+
  ggtitle("count")




scaled_timing_per_plot <- wgd_timing_df[!(is.na(wgd_timing_df$molecular_timing)) & wgd_timing_df$isWGD,] %>% ggplot(aes(x = molecular_timing, fill = is_metastatic, ..scaled..)) + #facet_wrap(~ cancer_type, scale = "free_y") +
  geom_density(alpha = 0.5, adjust = 0.5) +
  ggtitle("scaled")





for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig2-wgd-timing-plot-combined.png", height = 920, width = 460)
    print(grid.arrange(count_timing_plot, density_timing_plot, scaled_timing_per_plot))
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig2-wgd-timing-plot-combined.pdf", height = 14, width = 7)
    print(grid.arrange(count_timing_plot, density_timing_plot, scaled_timing_per_plot))
    dev.off()
  }
}












density_timing_per_cancer_type_plot <- wgd_timing_df[!(is.na(wgd_timing_df$molecular_timing)) & wgd_timing_df$isWGD,] %>% ggplot(aes(x = molecular_timing, fill = is_metastatic, ..density..)) + facet_wrap(~ cancer_type, scale = "free_y") +
  geom_density(alpha = 0.5, adjust = 0.2) +
  ggtitle("density")

# count_timing_per_cancer_type_plot <- wgd_timing_df[!(is.na(wgd_timing_df$molecular_timing)) & wgd_timing_df$isWGD,] %>% ggplot(aes(x = molecular_timing, fill = is_metastatic, ..count..)) + facet_wrap(~ cancer_type, scale = "free_y") +
#   geom_density(alpha = 0.5, adjust = 0.2) +
#   ggtitle("count")


count_timing_per_cancer_type_plot <- wgd_timing_df[!(is.na(wgd_timing_df$molecular_timing)) & wgd_timing_df$isWGD,]  %>% ggplot(aes(x = molecular_timing, fill = is_metastatic)) + facet_wrap(~ cancer_type, scale = "free_y") +
  geom_histogram(alpha = 0.5, bins = 20, position="identity")+
  ggtitle("count")






options(scipen=999) 
wgd_timing_df$Mann_whittney_p_value <- NA

for (cancer_type in unique(wgd_timing_df$cancer_type)) {
  
  print(cancer_type)
  dd <- tryCatch(wilcox.test(wgd_timing_df$molecular_timing[!(is.na(wgd_timing_df$molecular_timing)) & wgd_timing_df$isWGD & !(is.na(wgd_timing_df$molecular_timing)) & wgd_timing_df$cancer_type == cancer_type & wgd_timing_df$is_metastatic],  wgd_timing_df$molecular_timing[!(is.na(wgd_timing_df$molecular_timing)) & wgd_timing_df$isWGD & !(is.na(wgd_timing_df$molecular_timing)) & wgd_timing_df$cancer_type == cancer_type & !(wgd_timing_df$is_metastatic)]), error=function(e) NULL)
  if (!is.null(dd)){
    print(dd$p.value)
    wgd_timing_df$Mann_whittney_p_value[wgd_timing_df$cancer_type == cancer_type] <- round(dd$p.value, digits = 5)
  } else {
    wgd_timing_df$Mann_whittney_p_value[wgd_timing_df$cancer_type == cancer_type] <- "Not applicable"
  }
  
}


wgd_timing_df$Mann_whittney_p_value <- as.character(wgd_timing_df$Mann_whittney_p_value)

text_df <- distinct(wgd_timing_df[!(is.na(wgd_timing_df$molecular_timing)) & wgd_timing_df$isWGD,c("cancer_type", "Mann_whittney_p_value")])
text_df <- text_df[order(match(text_df$cancer_type, unique(wgd_timing_df$cancer_type))),]

# text_df_n <- text_df[text_df$Mann_whittney_p_value != "Not applicable",]
# text_df_n$Mann_whittney_p_value <- as.numeric(text_df_n$Mann_whittney_p_value)
# text_df_n[text_df_n$Mann_whittney_p_value < 0.05,]

scaled_timing_per_cancer_type_plot <- wgd_timing_df[!(is.na(wgd_timing_df$molecular_timing)) & wgd_timing_df$isWGD,] %>% ggplot(aes(x = molecular_timing, fill = is_metastatic, ..scaled..)) + facet_wrap(~ cancer_type) +
  geom_density(alpha = 0.5, adjust = 0.2) +
  ylim(c(0,1.5)) +
  geom_text(data = text_df, mapping = aes(label = Mann_whittney_p_value, x = 0.2, y = 1.15), inherit.aes = F) +
  ggtitle("scaled")


table(wgd_timing_df$molecular_timing[wgd_timing_df$cancer_type == "Osteosarcoma" & !(is.na(wgd_timing_df$molecular_timing)) &  wgd_timing_df$isWGD], wgd_timing_df$is_metastatic[wgd_timing_df$cancer_type == "Osteosarcoma" & !(is.na(wgd_timing_df$molecular_timing)) &  wgd_timing_df$isWGD])



for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig2-wgd-timing-per-cancer-type-plot-combined.png", height = 1380, width = 1160)
    print(grid.arrange(count_timing_per_cancer_type_plot, density_timing_per_cancer_type_plot, scaled_timing_per_cancer_type_plot))
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig2-wgd-timing-per-cancer-type-plot-combined.pdf", height = 21, width = 17.5)
    print(grid.arrange(count_timing_per_cancer_type_plot, density_timing_per_cancer_type_plot, scaled_timing_per_cancer_type_plot))
    dev.off()
  }
}


# ************************************************************************************************************************************************
# Fig3: violin plot of wgd timing + mann-whittney test

violin_timing_per_cancer_type_plot <- wgd_timing_df[!(is.na(wgd_timing_df$molecular_timing)) & wgd_timing_df$isWGD,] %>% ggplot(aes(x = is_metastatic, y = molecular_timing, fill = is_metastatic)) + facet_wrap(~ cancer_type, scale = "free_y") +
  geom_violin(alpha = 0.5) +
  ylim(c(0,1.5)) +
  geom_text(data = text_df, mapping = aes(label = Mann_whittney_p_value, x = 1.5, y = 1.25), inherit.aes = F) +
  ggtitle("WGD Timing \n") +
  ylab("Molecular Timing \n") +
  xlab("Cohort \n") +
  guides(fill=guide_legend(title="Cohort")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) +
  theme(axis.title.x = element_text(face = "italic", size = 12)) +
  theme(axis.title.y = element_text(face = "italic", size = 12)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))


for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig3-wgd-timing-per-cancer-type-violin-plot.png", height = 920, width = 920)
    print(violin_timing_per_cancer_type_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig3-wgd-timing-per-cancer-type-violin-plot.pdf", height = 14, width = 14)
    print(violin_timing_per_cancer_type_plot)
    dev.off()
  }
}







# ================================================================================================================================================
# ================================================================================================================================================
# Clonality

# !!! This is to investigate the discrepancy between clonality data from hmf pipeline and that of muttiontimerR package


sample_id <- "CPCT02010003T"
set_name <- hmf_meta$setName[hmf_meta$sampleId == sample_id]

if (dir.exists("/hpc/cuppen/")){
  mut_timing_path <- paste0("/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/timing/", set_name, "/", sample_id, ".mutationaltiming.tsv.gz")
} else {
  mut_timing_path <- paste0(local, "/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/timing/", set_name, "/", sample_id, ".mutationaltiming.tsv.gz")
}

mut_tim <- tryCatch(read.csv(file = mut_timing_path, header = T, sep = "\t", stringsAsFactors = F), error=function(e) NULL)
table(mut_tim$timing_class)


metadata_included[metadata_included$sample_id == "CPCT02010003T",]





for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/clonality-discrepancy-plot.png")
    plot(log2(global_timing_info_com$subclonal + 1), log2(global_timing_info_com$subclonal_tmb + 1))
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/clonality-discrepancy-plot.pdf")
    plot(log2(global_timing_info_com$subclonal + 1), log2(global_timing_info_com$subclonal_tmb + 1))
    dev.off()
  }
}




## 
global_timing_info_com$clonal_to_subclonal_metadata <- log2(global_timing_info_com$clonal_tmb/global_timing_info_com$subclonal_tmb)
global_timing_info_com[,6] <- global_timing_info_com[,6] - global_timing_info_com[,4]
global_timing_info_com[,7] <- global_timing_info_com[,7] + global_timing_info_com[,4]



cancer_types <- unique(global_timing_info_com$cancer_type)


# ************************************************************************************************************************************************
# fig 4 clonality per cancer type between met and prim

ff <- data.frame(cancer_type = character(59), cancer_type_abb = character(59), median_clona_met = numeric(59), median_clona_prim = numeric(59))

for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included$cancer_type_code[metadata_included$cancer_type == cancer_type])
  tmp_df <- global_timing_info_com[global_timing_info_com$cancer_type == cancer_type,]
  
  ff[i,1:2] <- c(cancer_type, cancer_type_abb)
  if (nrow(tmp_df[tmp_df$is_metastatic & tmp_df$subclonal_tmb != 0,]) >= 5 & nrow(tmp_df[!(tmp_df$is_metastatic) & tmp_df$subclonal_tmb != 0,]) >= 5){
    res <- wilcox.test(tmp_df$clonal_tmb[tmp_df$subclonal_tmb != 0 & tmp_df$is_metastatic]/tmp_df$tmb[tmp_df$subclonal_tmb != 0 & tmp_df$is_metastatic], tmp_df$clonal_tmb[tmp_df$subclonal_tmb != 0 & !(tmp_df$is_metastatic)]/tmp_df$tmb[tmp_df$subclonal_tmb != 0 & !(tmp_df$is_metastatic)])
    if (res$p.value > 0.05){
      ff[i,2] <- NA
    }
    ff[i,3:4] <- c(median(tmp_df$clonal_tmb[tmp_df$is_metastatic & tmp_df$subclonal_tmb != 0]/tmp_df$tmb[tmp_df$is_metastatic & tmp_df$subclonal_tmb != 0]), median(tmp_df$clonal_tmb[!(tmp_df$is_metastatic) & tmp_df$subclonal_tmb != 0]/tmp_df$tmb[!(tmp_df$is_metastatic) & tmp_df$subclonal_tmb != 0]))
  } else {
    ff[i,3:4] <- rep(NA, times = 2)
  }
}




ff$cancer_type <- factor(ff$cancer_type)
ff$cancer_type_abb <- factor(ff$cancer_type_abb)



plloott_median <- ff[!(is.na(ff$median_clona_met)) & !(is.na(ff$median_clona_prim)),] %>% ggplot(aes(x = 100*median_clona_prim, y = 100*median_clona_met, color = cancer_type)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb), max.overlaps = 20) +
  xlim(c(60,100)) +
  xlab("primary median %") +
  ylim(c(60,100)) +
  ylab("Metastatsis median %") +
  ggtitle("Percentage of Clonal mutations (HMF Pipeline)") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw()


global_timing_info_com_cp <- global_timing_info_com


for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  tmp_df <- global_timing_info_com_cp[global_timing_info_com_cp$cancer_type == cancer_type,]
  if (nrow(tmp_df[tmp_df$is_metastatic & tmp_df$subclonal_tmb != 0,]) >= 5 & nrow(tmp_df[!(tmp_df$is_metastatic) & tmp_df$subclonal_tmb != 0,]) >= 5){
    
  } else {
    global_timing_info_com_cp <- global_timing_info_com_cp[global_timing_info_com_cp$cancer_type != cancer_type,]
  }
}


global_timing_info_com_cp$cancer_type <- factor(global_timing_info_com_cp$cancer_type)




plloott_box <- global_timing_info_com_cp[global_timing_info_com_cp$subclonal_tmb != 0,] %>% ggplot(aes(x = is_metastatic, y = clonal_tmb/tmb, color = cancer_type)) + facet_wrap(~ cancer_type) +
  geom_boxplot() +
  # ggrepel::geom_text_repel(aes(label = cancer_type_abb)) +
  # xlim(c(50,100)) +
  xlab("Metastatic") +
  stat_compare_means(comparisons = list(c('FALSE', 'TRUE')), method = "wilcox.test", label = "p.signif") +
  ylab("Percentage of Clonal mutations (HMF Pipeline)") +
  ylim(c(0,1.5)) +
  ggtitle("Distribution of clonality percentage per cancer type (HMF Pipeline)") 



com_plot <- ggarrange(plloott_mean, plloott_median, plloott_box,  ncol=3, common.legend = TRUE, legend="bottom")

for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig4-clonality-prim-vs-metas-per-cancer.png", height = 920, width = 1380)
    print(com_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig4-clonality-prim-vs-metas-per-cancer.pdf", height = 14, width = 21)
    print(com_plot)
    dev.off()
  }
}









gg <- data.frame(cancer_type = character(59), cancer_type_abb = character(59), mean_clona_met = numeric(59), mean_clona_prim = numeric(59))

for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included$cancer_type_code[metadata_included$cancer_type == cancer_type])
  tmp_df <- global_timing_info_com[global_timing_info_com$cancer_type == cancer_type,]
  gg[i,1:2] <- c(cancer_type, cancer_type_abb)
  if (nrow(tmp_df[tmp_df$is_metastatic & !is.na(tmp_df$total_clonal) & tmp_df$subclonal != 0,]) >= 5 & nrow(tmp_df[!(tmp_df$is_metastatic) & !is.na(tmp_df$total_clonal) & tmp_df$subclonal != 0,]) >= 5){
    gg[i,3:4] <- c(mean(tmp_df$total_clonal[tmp_df$is_metastatic & tmp_df$subclonal != 0]/tmp_df$tmb[tmp_df$is_metastatic & tmp_df$subclonal != 0], na.rm = T), mean(tmp_df$total_clonal[!(tmp_df$is_metastatic) & tmp_df$subclonal != 0]/tmp_df$tmb[!(tmp_df$is_metastatic) & tmp_df$subclonal != 0], na.rm = T))
  } else {
    gg[i,3:4] <- rep(NA, times = 2)
  }
}



gg$cancer_type <- factor(gg$cancer_type)
gg$cancer_type_abb <- factor(gg$cancer_type_abb)

plloott_mean <- gg[!(is.na(gg$mean_clona_met)) & !(is.na(gg$mean_clona_prim)),] %>% ggplot(aes(x = 100*mean_clona_prim, y = 100*mean_clona_met, color = cancer_type)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb), max.overlaps = 20) +
  xlim(c(50,100)) +
  xlab("primary mean %") +
  ylim(c(50,100)) +
  ylab("Metastatsis mean %") +
  ggtitle("Percentage of Clonal mutations (mutationTimerR)") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw()



global_timing_info_com_cp <- global_timing_info_com


for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  tmp_df <- global_timing_info_com_cp[global_timing_info_com_cp$cancer_type == cancer_type,]
  if (nrow(tmp_df[tmp_df$is_metastatic & !is.na(tmp_df$total_clonal) & tmp_df$subclonal != 0,]) >= 5 & nrow(tmp_df[!(tmp_df$is_metastatic) & !is.na(tmp_df$total_clonal) & tmp_df$subclonal != 0,]) >= 5){
    
  } else {
    global_timing_info_com_cp <- global_timing_info_com_cp[global_timing_info_com_cp$cancer_type != cancer_type,]
  }
}


global_timing_info_com_cp$cancer_type <- factor(global_timing_info_com_cp$cancer_type)


global_timing_info_com_cp <- global_timing_info_com_cp[!is.na(global_timing_info_com_cp$subclonal),]


plloott_box <- global_timing_info_com_cp[global_timing_info_com_cp$subclonal != 0,] %>% ggplot(aes(x = is_metastatic, y = total_clonal/tmb, color = cancer_type)) + facet_wrap(~ cancer_type) +
  geom_boxplot() +
  # ggrepel::geom_text_repel(aes(label = cancer_type_abb)) +
  # xlim(c(50,100)) +
  xlab("Metastatic") +
  stat_compare_means(comparisons = list(c('FALSE', 'TRUE')), method = "wilcox.test", label = "p.signif") +
  ylab("Percentage of Clonal mutations (HMF Pipeline)") +
  ylim(c(0,1.5)) +
  ggtitle("Distribution of clonality percentage per cancer type (mutationTimerR)") 






com_plot <- ggarrange(plloott_mean, plloott_median, plloott_box, ncol=3, common.legend = TRUE, legend="bottom")

for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig4-clonality-prim-vs-metas-per-cancer-mutationTimerR-2.png", height = 920, width = 1380)
    print(com_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig4-clonality-prim-vs-metas-per-cancer-mutationTimerR-2.pdf", height = 14, width = 21)
    print(com_plot)
    dev.off()
  }
}



#### new timing with more lenient settings for pcawg samples



global_timing_info_com <- read.csv(file = paste0(wd, "r-objects/global-timing-with-metadata.txt.gz"), header = T, sep = "\t", stringsAsFactors = F)


global_timing_info_com <- global_timing_info_com[global_timing_info_com$sample_id != "DO217817",]




global_new_timing_info_com <- read.csv(file = paste0(wd, "r-objects/global-new-timing-with-metadata.txt.gz"), header = T, sep = "\t", stringsAsFactors = F)

global_new_timing_info_com <- global_new_timing_info_com[global_new_timing_info_com$sample_id != "DO217817",]

global_new_timing_info_com$tmb_based_on_mutationtimer <- as.integer(rowSums(global_new_timing_info_com[,c("total_clonal", "subclonal")], na.rm = T))
global_new_timing_info_com <- merge(global_new_timing_info_com, metadata_included[,c("sample_id", "cancer_type_code")], by = "sample_id")
## 
# global_new_timing_info_com$clonal_to_subclonal_metadata <- log2(global_new_timing_info_com$clonal_tmb/global_new_timing_info_com$subclonal_tmb)
# global_new_timing_info_com[,6] <- global_new_timing_info_com[,6] - global_new_timing_info_com[,4]
# global_new_timing_info_com[,7] <- global_new_timing_info_com[,7] + global_new_timing_info_com[,4]



cancer_types <- unique(global_new_timing_info_com$cancer_type)





gg <- data.frame(cancer_type = character(59), cancer_type_abb = character(59), mean_clona_met = numeric(59), mean_clona_prim = numeric(59))

for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included$cancer_type_code[metadata_included$cancer_type == cancer_type])
  tmp_df <- global_new_timing_info_com[global_new_timing_info_com$cancer_type == cancer_type,]
  gg[i,1:2] <- c(cancer_type, cancer_type_abb)
  if (nrow(tmp_df[tmp_df$is_metastatic & !is.na(tmp_df$total_clonal) & tmp_df$subclonal != 0,]) >= 5 & nrow(tmp_df[!(tmp_df$is_metastatic) & !is.na(tmp_df$total_clonal) & tmp_df$subclonal != 0,]) >= 5){
    gg[i,3:4] <- c(mean(tmp_df$total_clonal[tmp_df$is_metastatic & tmp_df$subclonal != 0]/tmp_df$tmb_based_on_mutationtimer[tmp_df$is_metastatic & tmp_df$subclonal != 0], na.rm = T), mean(tmp_df$total_clonal[!(tmp_df$is_metastatic) & tmp_df$subclonal != 0]/tmp_df$tmb_based_on_mutationtimer[!(tmp_df$is_metastatic) & tmp_df$subclonal != 0], na.rm = T))
  } else {
    gg[i,3:4] <- rep(NA, times = 2)
  }
}


gg$cancer_type <- paste0(gg$cancer_type, " (", gg$cancer_type_abb , ")")

gg$cancer_type <- factor(gg$cancer_type)
gg$cancer_type_abb <- factor(gg$cancer_type_abb)

plloott_mean <- gg[!(is.na(gg$mean_clona_met)) & !(is.na(gg$mean_clona_prim)),] %>% ggplot(aes(x = 100*mean_clona_prim, y = 100*mean_clona_met, color = cancer_type)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb), max.overlaps = 20) +
  xlim(c(50,100)) +
  xlab("primary mean %") +
  ylim(c(50,100)) +
  ylab("Metastatsis mean %") +
  ggtitle("Percentage of Clonal mutations (new mutationTimerR)") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw()



global_new_timing_info_com_cp <- global_new_timing_info_com


for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  tmp_df <- global_new_timing_info_com_cp[global_new_timing_info_com_cp$cancer_type == cancer_type,]
  if (nrow(tmp_df[tmp_df$is_metastatic & !is.na(tmp_df$total_clonal) & tmp_df$subclonal != 0,]) >= 5 & nrow(tmp_df[!(tmp_df$is_metastatic) & !is.na(tmp_df$total_clonal) & tmp_df$subclonal != 0,]) >= 5){
    
  } else {
    global_new_timing_info_com_cp <- global_new_timing_info_com_cp[global_new_timing_info_com_cp$cancer_type != cancer_type,]
  }
}


global_new_timing_info_com_cp$cancer_type <- paste0(global_new_timing_info_com_cp$cancer_type, " (", global_new_timing_info_com_cp$cancer_type_code , ")")

global_new_timing_info_com_cp$cancer_type <- factor(global_new_timing_info_com_cp$cancer_type)


global_new_timing_info_com_cp <- global_new_timing_info_com_cp[!is.na(global_new_timing_info_com_cp$subclonal),]


plloott_box <- global_new_timing_info_com_cp[global_new_timing_info_com_cp$subclonal != 0,] %>% ggplot(aes(x = is_metastatic, y = total_clonal/tmb_based_on_mutationtimer, color = cancer_type)) + facet_wrap(~ cancer_type) +
  geom_boxplot() +
  # ggrepel::geom_text_repel(aes(label = cancer_type_abb)) +
  # xlim(c(50,100)) +
  xlab("Metastatic") +
  stat_compare_means(comparisons = list(c('FALSE', 'TRUE')), method = "wilcox.test", label = "p.signif") +
  ylab("Percentage of Clonal mutations (HMF Pipeline)") +
  ylim(c(0,1.5)) +
  ggtitle("Distribution of clonality percentage per cancer type (new mutationTimerR)") 






com_plot <- ggarrange(plloott_mean, plloott_median, plloott_box, ncol=3, common.legend = TRUE, legend="bottom")

for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig4-new-timing-clonality-prim-vs-metas-per-cancer-mutationTimerR.png", height = 920, width = 1380)
    print(com_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig4-new-timing-clonality-prim-vs-metas-per-cancer-mutationTimerR.pdf", height = 14, width = 21)
    print(com_plot)
    dev.off()
  }
}



## comparing

nrow(global_timing_info_com[global_timing_info_com$cohort == "PCAWG" & !is.na(global_timing_info_com$total_clonal) & global_timing_info_com$subclonal != 0,])
nrow(global_new_timing_info_com[global_new_timing_info_com$cohort == "PCAWG" & !is.na(global_new_timing_info_com$total_clonal) & global_new_timing_info_com$subclonal != 0,])

int_sam <- intersect(global_timing_info_com[global_timing_info_com$cohort == "PCAWG" & !is.na(global_timing_info_com$total_clonal) & global_timing_info_com$subclonal != 0,"sample_id"], 
          global_new_timing_info_com[global_new_timing_info_com$cohort == "PCAWG" & !is.na(global_new_timing_info_com$total_clonal) & global_new_timing_info_com$subclonal != 0,"sample_id"])

oo <- global_timing_info_com[global_timing_info_com$sample_id %in% int_sam, c("sample_id", "total_clonal", "tmb")]
colnames(oo)[2:3] <- c("total_clonal_old", "tmb_old")

nn <- global_new_timing_info_com[global_new_timing_info_com$sample_id %in% int_sam, c("sample_id", "total_clonal", "tmb")]
colnames(nn)[2:3] <- c("total_clonal_new", "tmb_new")


mm <- merge(oo, nn, by = "sample_id")

wd
png(filename = paste0(wd,"explore/figs-wgd/png/old-new-clonal.png"))
plot(mm$total_clonal_old[mm$sample_id %in% int_sam]/mm$tmb_old[mm$sample_id %in% int_sam], mm$total_clonal_new[mm$sample_id %in% int_sam]/mm$tmb_new[mm$sample_id %in% int_sam], ylab = "old_timing clonal percentage", xlab = "new_timing clonal percentage", main = "PCAWG Clonality")

abline(coef = c(0,1), col = "red")

dev.off()

# ************************************************************************************************************************************************
# fig 5 clonality per cancer type between wgd and non-wgd


clonality_wgd <- global_timing_info_com[!is.infinite(global_timing_info_com$clonal_to_subclonal_metadata),] %>% ggplot(aes(x = whole_genome_duplication, y = clonal_to_subclonal_metadata)) +
  geom_boxplot()


hh <- data.frame(cancer_type = character(59), cancer_type_abb = character(59), mean_clona_wgd = numeric(59), mean_clona_non_wgd = numeric(59))

for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included$cancer_type_code[metadata_included$cancer_type == cancer_type])
  tmp_df <- global_timing_info_com[global_timing_info_com$cancer_type == cancer_type,]
  hh[i,1:2] <- c(cancer_type, cancer_type_abb)
  if (nrow(tmp_df[tmp_df$whole_genome_duplication & tmp_df$subclonal_tmb != 0,]) >= 5 & nrow(tmp_df[!(tmp_df$whole_genome_duplication) & tmp_df$subclonal_tmb != 0,]) >= 5){
    res <- wilcox.test(tmp_df$clonal_tmb[tmp_df$subclonal_tmb != 0 & tmp_df$whole_genome_duplication]/tmp_df$tmb[tmp_df$subclonal_tmb != 0 & tmp_df$whole_genome_duplication], tmp_df$clonal_tmb[tmp_df$subclonal_tmb != 0 & !(tmp_df$whole_genome_duplication)]/tmp_df$tmb[tmp_df$subclonal_tmb != 0 & !(tmp_df$whole_genome_duplication)])
    if (res$p.value > 0.05){
      hh[i,2] <- NA
    }
    hh[i,3:4] <- c(mean(tmp_df$clonal_tmb[tmp_df$whole_genome_duplication & tmp_df$subclonal_tmb != 0]/tmp_df$tmb[tmp_df$whole_genome_duplication & tmp_df$subclonal_tmb != 0], na.rm = T), mean(tmp_df$clonal_tmb[!(tmp_df$whole_genome_duplication) & tmp_df$subclonal_tmb != 0]/tmp_df$tmb[!(tmp_df$whole_genome_duplication) & tmp_df$subclonal_tmb != 0], na.rm = T))
  } else {
    hh[i,3:4] <- rep(NA, times = 2)
  }
}


hh$cancer_type <- factor(hh$cancer_type)
hh$cancer_type_abb <- factor(hh$cancer_type_abb)

plot_mean <- hh[!(is.na(hh$mean_clona_non_wgd)) & !(is.na(hh$mean_clona_wgd)),] %>% ggplot(aes(x = 100*mean_clona_wgd, y = 100*mean_clona_non_wgd, color = cancer_type)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb), max.overlaps = 20) +
  xlim(c(60,100)) +
  xlab("wgd mean %") +
  ylim(c(60,100)) +
  ylab("non-wgd mean %") +
  ggtitle("Percentage of Clonal mutations (HMF Pipeline)") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw()




global_timing_info_com_cp <- global_timing_info_com


for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  tmp_df <- global_timing_info_com_cp[global_timing_info_com_cp$cancer_type == cancer_type,]
  if (nrow(tmp_df[tmp_df$whole_genome_duplication & tmp_df$subclonal_tmb != 0,]) >= 5 & nrow(tmp_df[!(tmp_df$whole_genome_duplication) & tmp_df$subclonal_tmb != 0,]) >= 5){
    
  } else {
    global_timing_info_com_cp <- global_timing_info_com_cp[global_timing_info_com_cp$cancer_type != cancer_type,]
  }
}

global_timing_info_com_cp$cancer_type <- factor(global_timing_info_com_cp$cancer_type)




plloott_box <- global_timing_info_com_cp[global_timing_info_com_cp$subclonal_tmb != 0,] %>% ggplot(aes(x = whole_genome_duplication, y = clonal_tmb/tmb, color = cancer_type)) + facet_wrap(~ cancer_type) +
  geom_boxplot() +
  # ggrepel::geom_text_repel(aes(label = cancer_type_abb)) +
  # xlim(c(50,100)) +
  xlab("WGD") +
  stat_compare_means(comparisons = list(c('FALSE', 'TRUE')), method = "wilcox.test", label = "p.signif") +
  ylab("Percentage of Clonal mutations (HMF Pipeline)") +
  ylim(c(0,1.5)) +
  ggtitle("Distribution of clonality percentage per cancer type (HMF Pipeline)") 





com_plot <- ggarrange(plot_mean, plot_median, plloott_box, ncol=3, common.legend = TRUE, legend="bottom")


for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig5-clonality-wgd-vs-nonwgd-per-cancer.png", height = 920, width = 1380)
    print(com_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig5-clonality-wgd-vs-nonwgd-per-cancer.pdf", height = 14, width = 21)
    print(com_plot)
    dev.off()
  }
}











jj <- data.frame(cancer_type = character(59), cancer_type_abb = character(59), median_clona_wgd = numeric(59), median_clona_non_wgd = numeric(59))

for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included$cancer_type_code[metadata_included$cancer_type == cancer_type])
  tmp_df <- global_timing_info_com[global_timing_info_com$cancer_type == cancer_type,]
  jj[i,1:2] <- c(cancer_type, cancer_type_abb)
  if (nrow(tmp_df[tmp_df$whole_genome_duplication & !is.na(tmp_df$total_clonal) & tmp_df$subclonal != 0,]) >= 5 & nrow(tmp_df[!(tmp_df$whole_genome_duplication) & !is.na(tmp_df$total_clonal) & tmp_df$subclonal != 0,]) >= 5){
    jj[i,3:4] <- c(median(tmp_df$total_clonal[tmp_df$whole_genome_duplication & !is.na(tmp_df$total_clonal) & tmp_df$subclonal != 0]/tmp_df$tmb[tmp_df$whole_genome_duplication & !is.na(tmp_df$total_clonal) & tmp_df$subclonal != 0], na.rm = T), median(tmp_df$total_clonal[!(tmp_df$whole_genome_duplication) & !is.na(tmp_df$total_clonal) & tmp_df$subclonal != 0]/tmp_df$tmb[!(tmp_df$whole_genome_duplication) & !is.na(tmp_df$total_clonal) & tmp_df$subclonal != 0], na.rm = T))
  } else {
    jj[i,3:4] <- rep(NA, times = 2)
  }
}


jj$cancer_type <- factor(jj$cancer_type)
jj$cancer_type_abb <- factor(jj$cancer_type_abb)

plot_median <- jj[!(is.na(jj$median_clona_non_wgd)) & !(is.na(jj$median_clona_wgd)),] %>% ggplot(aes(x = 100*median_clona_wgd, y = 100*median_clona_non_wgd, color = cancer_type)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb), max.overlaps = 20) +
  xlim(c(50,100)) +
  xlab("wgd median %") +
  ylim(c(50,100)) +
  ylab("non-wgd median %") +
  ggtitle("Percentage of Clonal mutations (mutationTimerR)") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw()



global_timing_info_com_cp <- global_timing_info_com
for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  tmp_df <- global_timing_info_com_cp[global_timing_info_com_cp$cancer_type == cancer_type,]
  if (nrow(tmp_df[tmp_df$whole_genome_duplication & !is.na(tmp_df$total_clonal) & tmp_df$subclonal != 0,]) >= 5 & nrow(tmp_df[!(tmp_df$whole_genome_duplication) & !is.na(tmp_df$total_clonal) & tmp_df$subclonal != 0,]) >= 5){
    print(cancer_type)
    print(nrow(tmp_df[tmp_df$whole_genome_duplication & !is.na(tmp_df$total_clonal) & tmp_df$subclonal_tmb != 0,]))
    print(nrow(tmp_df[!(tmp_df$whole_genome_duplication) & !is.na(tmp_df$total_clonal) & tmp_df$subclonal_tmb != 0,]))
  } else {
    global_timing_info_com_cp <- global_timing_info_com_cp[global_timing_info_com_cp$cancer_type != cancer_type,]
  }
}
global_timing_info_com_cp$cancer_type <- factor(global_timing_info_com_cp$cancer_type)

global_timing_info_com_cp <- global_timing_info_com_cp[!is.na(global_timing_info_com_cp$subclonal),]



plloott_box <- global_timing_info_com_cp[global_timing_info_com_cp$subclonal != 0,] %>% ggplot(aes(x = whole_genome_duplication, y = total_clonal/tmb, color = cancer_type)) + facet_wrap(~ cancer_type) +
  geom_boxplot() +
  # ggrepel::geom_text_repel(aes(label = cancer_type_abb)) +
  # xlim(c(50,100)) +
  xlab("WGD") +
  stat_compare_means(comparisons = list(c('FALSE', 'TRUE')), method = "wilcox.test", label = "p.signif") +
  ylab("Percentage of Clonal mutations (HMF Pipeline)") +
  ylim(c(0,1.5)) +
  ggtitle("Distribution of clonality percentage per cancer type (mutationTimerR)") 


com_plot <- ggarrange(plot_mean, plot_median, plloott_box, ncol=3 , common.legend = TRUE, legend="bottom")


for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig5-clonality-wgd-vs-nonwgd-per-cancer-mutationTimerR.png", height = 920, width = 1380)
    print(com_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig5-clonality-wgd-vs-nonwgd-per-cancer-mutationTimerR.pdf", height = 14, width = 21)
    print(com_plot)
    dev.off()
  }
}








# ================================================================================================================================================
# ================================================================================================================================================
# Somatic point mutations timing prim vs metastatic

global_timing_info_com_pruned <- global_timing_info_com[!is.na(global_timing_info_com$early_clonal) & !is.na(global_timing_info_com$late_clonal) & global_timing_info_com$early_clonal > 50 & global_timing_info_com$late_clonal > 50,]
nrow(global_timing_info_com_pruned)
for (cancer_type in unique(global_timing_info_com_pruned$cancer_type)) {
  tmp_df <- global_timing_info_com_pruned[global_timing_info_com_pruned$cancer_type == cancer_type,]
  if (nrow(tmp_df[tmp_df$is_metastatic,]) < 5 | nrow(tmp_df[!(tmp_df$is_metastatic),]) < 5) {
    global_timing_info_com_pruned <- global_timing_info_com_pruned[global_timing_info_com_pruned$cancer_type != cancer_type,]
  }
}

global_timing_info_com_pruned$early_to_late_clonal <- global_timing_info_com_pruned$early_clonal/global_timing_info_com_pruned$late_clonal

global_timing_info_com_pruned$cancer_type <- factor(global_timing_info_com_pruned$cancer_type, level = unique(global_timing_info_com_pruned$cancer_type))




# ************************************************************************************************************************************************
# fig 6 mutation timing per cancer type between prim and metas

global_mutation_timing_plot1 <- global_timing_info_com_pruned %>% ggplot(aes(x = cancer_type, y = log2(early_to_late_clonal), color = is_metastatic)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90))



for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig6-global-mutation-timing-prim-vs-met.png", height = 920, width = 1380)
    print(global_mutation_timing_plot1)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig6-global-mutation-timing-prim-vs-met.pdf", height = 14, width = 21)
    print(global_mutation_timing_plot1)
    dev.off()
  }
}


facet_name <- c(
  `FALSE` = "Primary",
  `TRUE` = "Metastatic"
)


global_mutation_timing_plot2 <- global_timing_info_com_pruned %>% ggplot(aes(x = cancer_type, y = log2(early_to_late_clonal))) +
  facet_wrap(~ is_metastatic, nrow = 2, labeller = as_labeller(facet_name)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Cancer Type") +
  ylab("Early-to-Late Clonal ratio (log2)")


# for (i in 1:2) {
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig6-global-mutation-timing-prim-and-metas.png", height = 920, width = 1380)
#     print(global_mutation_timing_plot2)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig6-global-mutation-timing-prim-and-metas.pdf", height = 14, width = 21)
#     print(global_mutation_timing_plot2)
#     dev.off()
#   }
# }
# 






global_mutation_timing_plot_per_cancer_type <- global_timing_info_com_pruned %>% ggplot(aes(x = is_metastatic, y = log2(early_to_late_clonal), color = is_metastatic)) +
  facet_wrap(~ cancer_type, scale = "free_y") +
  geom_boxplot() + 
  # ylim(c(-5,12)) +
  stat_compare_means(comparisons = list(c('FALSE', 'TRUE')), method = "wilcox.test", label = "p.value", size = 3) +
  ylab("Early-to-Late Clonal ratio (logg2)") +
  xlab("Metastastasis Status") +
  ggtitle("Mutational Timing Comparison")



# for (i in 1:2) {
#   if (i == 1) {
#     png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig6-global-mutation-timing-prim-vs-met-per-cancer.png", height = 920, width = 1380)
#     print(global_mutation_timing_plot_per_cancer_type)
#     dev.off()
#   }
#   if (i == 2) {
#     pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig6-global-mutation-timing-prim-vs-met-per-cancer.pdf", height = 14, width = 21)
#     print(global_mutation_timing_plot_per_cancer_type)
#     dev.off()
#   }
# }




com_plot <- ggarrange(global_mutation_timing_plot2, global_mutation_timing_plot_per_cancer_type, ncol = 2)


for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig6-global-mutation-timing-prim-vs-met-combined.png", height = 920, width = 1380)
    print(com_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig6-global-mutation-timing-prim-vs-met-combined.pdf", height = 14, width = 21)
    print(com_plot)
    dev.off()
  }
}



# ================================================================================================================================================
# ================================================================================================================================================
# WGD occurrence and timing vs tmb



head(wgd_timing_df)
str(wgd_timing_df)
table(wgd_timing_df$isWGD, wgd_timing_df$whole_genome_duplication, useNA = "always")
table(wgd_timing_df$isWGD, wgd_timing_df$molecular_timing, useNA = "always")
table(wgd_timing_df$whole_genome_duplication, wgd_timing_df$molecular_timing, useNA = "always")






# ************************************************************************************************************************************************
# fig 7 tmb per cancer type between wgd and non-wgd

ll <- data.frame(cancer_type = character(59), cancer_type_abb = character(59), cancer_type_abb_text = character(59), median_tmb_wgd = numeric(59), median_tmb_non_wgd = numeric(59))
cancer_types <- unique(metadata_included$cancer_type)


for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included$cancer_type_code[metadata_included$cancer_type == cancer_type])
  tmp_df <- wgd_timing_df[wgd_timing_df$cancer_type == cancer_type & !(wgd_timing_df$is_hypermutated),]
  ll[i,1:3] <- c(cancer_type, cancer_type_abb, cancer_type_abb)
  if (nrow(tmp_df[tmp_df$whole_genome_duplication,]) >= 5 & nrow(tmp_df[!(tmp_df$whole_genome_duplication),]) >= 5){
    res <- wilcox.test(tmp_df$tmb[tmp_df$whole_genome_duplication], tmp_df$tmb[!(tmp_df$whole_genome_duplication)])
    if (res$p.value > 0.05){
      ll[i,3] <- NA
    }
    ll[i,4:5] <- c(median(tmp_df$tmb[tmp_df$whole_genome_duplication], na.rm = T), median(tmp_df$tmb[!(tmp_df$whole_genome_duplication)], na.rm = T))
  } else {
    ll[i,4:5] <- rep(NA, times = 2)
  }
}

ll <- ll[!is.na(ll$median_tmb_wgd),]
ll$cancer_type <- paste0(ll$cancer_type, " (", ll$cancer_type_abb, ")")
ll$cancer_type <- factor(ll$cancer_type)
ll$cancer_type_abb <- factor(ll$cancer_type_abb)

plot_median <- ll[!(is.na(ll$median_tmb_non_wgd)) & !(is.na(ll$median_tmb_wgd)),] %>% ggplot(aes(x = log2(median_tmb_wgd), y = log2(median_tmb_non_wgd), color = cancer_type)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb_text), max.overlaps = 20) +
  xlim(c(10,18)) +
  xlab("TMB median (WGD+)") +
  ylim(c(10,18)) +
  ylab("TMB median (WGD-)") +
  ggtitle("TMB median") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw()



# creating the boxplot


wgd_timing_df_cp <- wgd_timing_df

for (cancer_type in cancer_types){
  tmp_df <- wgd_timing_df[wgd_timing_df$cancer_type == cancer_type & !(wgd_timing_df$is_hypermutated),]
  if (nrow(tmp_df[tmp_df$whole_genome_duplication,]) >= 5 & nrow(tmp_df[!tmp_df$whole_genome_duplication,]) >= 5){
    
  } else {
    wgd_timing_df_cp <- wgd_timing_df_cp[wgd_timing_df_cp$cancer_type != cancer_type,]
  }
}

wgd_timing_df_cp$cancer_type <- factor(wgd_timing_df_cp$cancer_type)


wgd_tmb_box_plot <- wgd_timing_df_cp[!(wgd_timing_df_cp$is_hypermutated),] %>% ggplot(aes(x = whole_genome_duplication, y = log2(tmb), color = whole_genome_duplication)) + facet_wrap(~ cancer_type) +
  geom_boxplot() +
  xlab("WGD") +
  stat_compare_means(comparisons = list(c('FALSE', 'TRUE')), method = "wilcox.test", label = "p.signif") +
  ylab("TMB (log2)") +
  ylim(c(5, 25)) 
  # ggtitle("Distribution of clonality percentage per cancer type (mutationTimerR)") 




com_plot <- ggarrange(plot_mean, plot_median, wgd_tmb_box_plot, ncol=3, common.legend = TRUE, legend="bottom")
com_plot <- annotate_figure(
  com_plot,
  bottom = "Hypermutated samples are excluded!")

for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig7-tmb-wgd-vs-nonwgd-per-cancer.png", height = 920, width = 1840)
    print(com_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig7-tmb-wgd-vs-nonwgd-per-cancer.pdf", height = 14, width = 28)
    print(com_plot)
    dev.off()
  }
}

# ************************************************************************************************************************************************
# fig 7.1 tmb per cancer type between wgd and non-wgd considering the metastasis status


zz <- data.frame(cancer_type = character(59), cancer_type_abb = character(59), cancer_type_abb_text_prim = character(59), cancer_type_abb_text_metas = character(59), median_tmb_wgd_prim = numeric(59), median_tmb_non_wgd_prim = numeric(59), median_tmb_wgd_metas = numeric(59), median_tmb_non_wgd_metas = numeric(59))
cancer_types <- unique(metadata_included$cancer_type)


for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included$cancer_type_code[metadata_included$cancer_type == cancer_type])
  tmp_df <- wgd_timing_df[wgd_timing_df$cancer_type == cancer_type & !(wgd_timing_df$is_hypermutated),]
  zz[i,1:4] <- c(cancer_type, cancer_type_abb, cancer_type_abb, cancer_type_abb)
  if (nrow(tmp_df[tmp_df$whole_genome_duplication & tmp_df$is_metastatic,]) >= 5 & nrow(tmp_df[tmp_df$whole_genome_duplication & !(tmp_df$is_metastatic),]) >= 5 & nrow(tmp_df[!(tmp_df$whole_genome_duplication) & tmp_df$is_metastatic,]) >= 5 & nrow(tmp_df[!(tmp_df$whole_genome_duplication) & !(tmp_df$is_metastatic),]) >= 5){
    res <- wilcox.test(tmp_df$tmb[tmp_df$whole_genome_duplication & !(tmp_df$is_metastatic)], tmp_df$tmb[!(tmp_df$whole_genome_duplication) & !(tmp_df$is_metastatic)])
    if (res$p.value > 0.05){
      zz[i,3] <- NA
    }
    
    res <- wilcox.test(tmp_df$tmb[tmp_df$whole_genome_duplication & tmp_df$is_metastatic], tmp_df$tmb[!(tmp_df$whole_genome_duplication) & tmp_df$is_metastatic])
    if (res$p.value > 0.05){
      zz[i,4] <- NA
    }
    
    zz[i,5:8] <- c(median(tmp_df$tmb[tmp_df$whole_genome_duplication & !(tmp_df$is_metastatic)], na.rm = T), median(tmp_df$tmb[!(tmp_df$whole_genome_duplication) & !(tmp_df$is_metastatic)], na.rm = T), 
                   median(tmp_df$tmb[tmp_df$whole_genome_duplication & tmp_df$is_metastatic], na.rm = T), median(tmp_df$tmb[!(tmp_df$whole_genome_duplication) & tmp_df$is_metastatic], na.rm = T))
  } else {
    zz[i,5:8] <- rep(NA, times = 4)
  }
}

zz <- zz[!is.na(zz$median_tmb_wgd_prim),]
zz$cancer_type <- paste0(zz$cancer_type, " (", zz$cancer_type_abb, ")")
zz$cancer_type <- factor(zz$cancer_type)
zz$cancer_type_abb <- factor(zz$cancer_type_abb)



plot_median_prim <- zz %>% ggplot(aes(x = log2(median_tmb_wgd_prim), y = log2(median_tmb_non_wgd_prim), color = cancer_type)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb_text_prim), max.overlaps = 20) +
  xlim(c(10,18)) +
  xlab("TMB median (WGD+)") +
  ylim(c(10,18)) +
  ylab("TMB median (WGD-)") +
  ggtitle("TMB median _ Primary") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw()


plot_median_metas <- zz %>% ggplot(aes(x = log2(median_tmb_wgd_metas), y = log2(median_tmb_non_wgd_metas), color = cancer_type)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb_text_metas), max.overlaps = 20) +
  xlim(c(10,18)) +
  xlab("TMB median (WGD+)") +
  ylim(c(10,18)) +
  ylab("TMB median (WGD-)") +
  ggtitle("TMB median _ Metastatic") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw()





# creating the boxplot


wgd_timing_df_cp <- wgd_timing_df

for (cancer_type in cancer_types){
  tmp_df <- wgd_timing_df[wgd_timing_df$cancer_type == cancer_type & !(wgd_timing_df$is_hypermutated),]
  if (nrow(tmp_df[tmp_df$whole_genome_duplication & tmp_df$is_metastatic,]) >= 5 & nrow(tmp_df[tmp_df$whole_genome_duplication & !(tmp_df$is_metastatic),]) >= 5 & nrow(tmp_df[!(tmp_df$whole_genome_duplication) & tmp_df$is_metastatic,]) >= 5 & nrow(tmp_df[!(tmp_df$whole_genome_duplication) & !(tmp_df$is_metastatic),]) >= 5){
    
  } else {
    wgd_timing_df_cp <- wgd_timing_df_cp[wgd_timing_df_cp$cancer_type != cancer_type,]
  }
}

wgd_timing_df_cp$cancer_type <- factor(wgd_timing_df_cp$cancer_type)
wgd_timing_df_cp$wgd_and_cohort <- NA
wgd_timing_df_cp$wgd_and_cohort[which(wgd_timing_df_cp$whole_genome_duplication & wgd_timing_df_cp$is_metastatic)] <- "WGD_METAS"
wgd_timing_df_cp$wgd_and_cohort[which(wgd_timing_df_cp$whole_genome_duplication & !(wgd_timing_df_cp$is_metastatic))] <- "WGD_PRIM"
wgd_timing_df_cp$wgd_and_cohort[which(!(wgd_timing_df_cp$whole_genome_duplication) & wgd_timing_df_cp$is_metastatic)] <- "NonWGD_METAS"
wgd_timing_df_cp$wgd_and_cohort[which(!(wgd_timing_df_cp$whole_genome_duplication) & !(wgd_timing_df_cp$is_metastatic))] <- "NonWGD_PRIM"

wgd_timing_df_cp$wgd_and_cohort <- factor(wgd_timing_df_cp$wgd_and_cohort, levels = c("NonWGD_PRIM", "WGD_PRIM", "NonWGD_METAS", "WGD_METAS"))


options(scripen = 999)

wgd_tmb_box_plot <- wgd_timing_df_cp[!(wgd_timing_df_cp$is_hypermutated),] %>% ggplot(aes(x = wgd_and_cohort, y = log2(tmb), color = whole_genome_duplication, fill = is_metastatic)) + facet_wrap(~ cancer_type) +
  geom_boxplot() +
  scale_fill_manual(values = c("black", "white")) +
  scale_color_manual(values = c("red", "blue")) +
  xlab("WGD") +
  stat_compare_means(comparisons = list(c('FALSE', 'TRUE')), method = "wilcox.test", label = "p.signif") +
  ylab("TMB (log2)") +
  ylim(c(10, 16)) +
  stat_compare_means(method = "anova") +
  theme(axis.text.x = element_text(angle = 90))
# ggtitle("Distribution of clonality percentage per cancer type (mutationTimerR)") 

com_plot <- ggarrange(plot_mean_prim, plot_median_prim, plot_mean_metas, plot_median_metas, nrow = 2, ncol=2, common.legend = TRUE, legend="bottom")
com_plot <- ggarrange(com_plot, wgd_tmb_box_plot, nrow = 1, ncol=2)


for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig7.1-tmb-wgd-vs-nonwgd-per-cancer-per-cohort.png", height = 920, width = 1840)
    print(com_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig7.1-tmb-wgd-vs-nonwgd-per-cancer-per-cohort.pdf", height = 14, width = 28)
    print(com_plot)
    dev.off()
  }
}


plot_median_ratios <- zz %>% ggplot(aes(x = median_tmb_wgd_metas/median_tmb_non_wgd_metas, y = median_tmb_wgd_prim/median_tmb_non_wgd_prim, color = cancer_type)) + 
  geom_point() +
  ggrepel::geom_text_repel(data = zz[!is.na(zz$cancer_type_abb_text_prim) & !is.na(zz$cancer_type_abb_text_metas),], aes(label = cancer_type_abb), max.overlaps = 20) +
  xlim(c(0,5)) +
  xlab("Ratio of TMB median WGD(+)/WGD(-) in Metastatic Cohort") +
  ylim(c(0,5)) +
  ylab("Ratio of TMB median WGD(+)/WGD(-) in Primary Cohort") +
  ggtitle("Ratio of TMB median change") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw() +
  geom_hline(yintercept = 1, color = "red", lty = 2, alpha = 0.5) +
  geom_vline(xintercept = 1, color = "red", lty = 2, alpha = 0.5)

for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig7.2-tmb-wgd-vs-nonwgd-per-cancer-per-cohort.png")
    print(plot_median_ratios)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig7.2-tmb-wgd-vs-nonwgd-per-cancer-per-cohort.pdf")
    print(plot_median_ratios)
    dev.off()
  }
}






# ************************************************************************************************************************************************
# fig 8 tmb per cancer type and wgd timing



# removing hyperutated samples results in better models (higher pearson correlations and r-squared values)

wgd_timing_df_cpp <- wgd_timing_df
cancer_types <- unique(metadata_included$cancer_type)

for (cancer_type in cancer_types){
  tmp_df <- wgd_timing_df[wgd_timing_df$whole_genome_duplication & !is.na(wgd_timing_df$molecular_timing) & wgd_timing_df$cancer_type == cancer_type & !(wgd_timing_df$is_hypermutated),]
  if (nrow(tmp_df) >= 5){
    
  } else {
    wgd_timing_df_cpp <- wgd_timing_df_cpp[wgd_timing_df_cpp$cancer_type != cancer_type,]
  }
}

wgd_timing_df_cpp$cancer_type <- factor(wgd_timing_df_cpp$cancer_type)



options(scripen = 999)
cancer_types <- as.character(unique(wgd_timing_df_cpp$cancer_type))


corr_df2 <- data.frame(cancer_type = character(24), cancer_type_abb = character(24), pearson_correlation = numeric(24), coefficient = numeric(24), p_val = numeric(24), adj_r_squared = numeric(24))
for (i in 1:length(unique(wgd_timing_df_cpp$cancer_type))){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included$cancer_type_code[metadata_included$cancer_type == cancer_type])
  tmp_df <- wgd_timing_df_cpp[wgd_timing_df_cpp$whole_genome_duplication & wgd_timing_df_cpp$cancer_type == cancer_type & !is.na(wgd_timing_df_cpp$molecular_timing) & !(wgd_timing_df_cpp$is_hypermutated),]
  p_corr <- cor(tmp_df$molecular_timing,
                log2(tmp_df$tmb))
  model <- lm(log2(tmb) ~ molecular_timing, data = tmp_df)
  sum_model <- summary(model)
  corr_df2[i,1:2] <- c(cancer_type, cancer_type_abb)
  corr_df2[i,3:6] <- c(p_corr, sum_model$coefficients[2,1], sum_model$coefficients[2,4], sum_model$adj.r.squared)
}


mod_summary_sign <- corr_df2$p_val
names(mod_summary_sign) <- corr_df2$cancer_type
mod_summary_stars <- NA                             # Named vector with significance stars
mod_summary_stars[mod_summary_sign < 0.1] <- "."
mod_summary_stars[mod_summary_sign < 0.05] <- "*"
mod_summary_stars[mod_summary_sign < 0.01] <- "**"
mod_summary_stars[mod_summary_sign < 0.001] <- "***"
mod_summary_stars[is.na(mod_summary_stars)] <- "n.s."

names(mod_summary_stars) <- names(mod_summary_sign)
corr_df2$p_val_sig_code <- as.vector(mod_summary_stars)
corr_df2$pearson_correlation <- round(corr_df2$pearson_correlation, digits = 2)
corr_df2$coefficient <- round(corr_df2$coefficient, digits = 2)



wgd_timing_tmb_box_plot <- wgd_timing_df_cpp[wgd_timing_df_cpp$whole_genome_duplication & !is.na(wgd_timing_df_cpp$molecular_timing) & !(wgd_timing_df_cpp$is_hypermutated),] %>% ggplot(aes(x = molecular_timing, y = log2(tmb))) + facet_wrap(~ cancer_type) +
  geom_point(size = 0.75) +
  ggrepel::geom_text_repel(data = corr_df2, x= 0.85, y = 25, aes(label = p_val_sig_code), size = 6) +
  ggrepel::geom_text_repel(data = corr_df2, x= 0.15, y = 4, aes(label = paste("Pearson corr.:", pearson_correlation, ", Line slope:", coefficient)), size = 3) +
  xlab("WGD Molecular Time") +
  geom_smooth(method='lm', formula= y~x, color = "red", se=F) +
  ylab("TMB (log2)") +
  ylim(c(5, 25)) +
  labs(title = "",
       subtitle = "",
       caption = "Hypermutated samples are excluded!")


for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig8-tmb-timing-wgd-vs-nonwgd-per-cancer.png", height = 1380, width = 1380)
    print(wgd_timing_tmb_box_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig8-tmb-timing-wgd-vs-nonwgd-per-cancer.pdf", height = 21, width = 21)
    print(wgd_timing_tmb_box_plot)
    dev.off()
  }
}


# plotting the statistics of the model

corr_df2$cancer_type <- paste0(corr_df2$cancer_type, " (", corr_df2$cancer_type_abb, ")")
corr_df2$cancer_type <- factor(corr_df2$cancer_type)
corr_df2$cancer_type_abb[which(corr_df2$p_val_sig_code == "n.s.")] <- NA

scatt_plot <- corr_df2 %>% ggplot(aes(x= coefficient, y = abs(pearson_correlation), size = -log10(p_val), color = cancer_type)) +
  geom_point() +
  xlim(c(-3,3))+
  ggrepel::geom_text_repel(aes(label = cancer_type_abb), max.overlaps = 20) +
  labs(title = "Correlation between WGD timing and TMB",
       subtitle = "",
       caption = "Hypermutated samples are excluded!")


for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig8-tmb-timing-wgd-vs-nonwgd-scatter-plot.png", height = 460, width = 920)
    print(scatt_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig8-tmb-timing-wgd-vs-nonwgd-scatter-plot.pdf", height = 7, width = 14)
    print(scatt_plot)
    dev.off()
  }
}



# 
# 
# corr_df2[corr_df2$p_val < 0.05,]
# 
# i <- 30
# cancer_type <- cancer_types[i]
# tmp_df <- wgd_timing_df_cpp[wgd_timing_df_cpp$whole_genome_duplication & wgd_timing_df_cpp$cancer_type == cancer_type & !is.na(wgd_timing_df_cpp$molecular_timing),]
# 
# 
# wgd_timing_tmb_box_plot <- tmp_df %>% ggplot(aes(x = molecular_timing, y = log2(tmb))) +
#   geom_point(size = 0.5) +
#   xlab("WGD Molecular Time") +
#   geom_smooth(method='lm', formula= y~x, color = "red", se=F) +
#   ylab("TMB (log2)")
# # ylim(c(5, 25)) 




# ************************************************************************************************************************************************
# fig 8.1 tmb per cancer type and wgd timing




# removing hyperutated samples results in better models (higher pearson correlations and r-squared values)





wgd_timing_df_cpp <- wgd_timing_df
cancer_types <- unique(metadata_included$cancer_type)

for (cancer_type in cancer_types){
  tmp_df <- wgd_timing_df[wgd_timing_df$whole_genome_duplication & !is.na(wgd_timing_df$molecular_timing) & wgd_timing_df$cancer_type == cancer_type & !(wgd_timing_df$is_hypermutated),]
  if (nrow(tmp_df[tmp_df$is_metastatic,]) >= 5 & nrow(tmp_df[!(tmp_df$is_metastatic),]) >= 5){
    
  } else {
    wgd_timing_df_cpp <- wgd_timing_df_cpp[wgd_timing_df_cpp$cancer_type != cancer_type,]
  }
}

wgd_timing_df_cpp$cancer_type <- factor(wgd_timing_df_cpp$cancer_type)

options(scripen = 999)
cancer_types <- as.character(unique(wgd_timing_df_cpp$cancer_type))


corr_df_prim <- data.frame(cancer_type = character(19), cancer_type_abb = character(19), pearson_correlation = numeric(19), coefficient = numeric(19), p_val = numeric(19), adj_r_squared = numeric(19))
for (i in 1:length(unique(wgd_timing_df_cpp$cancer_type))){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included$cancer_type_code[metadata_included$cancer_type == cancer_type])
  tmp_df <- wgd_timing_df_cpp[wgd_timing_df_cpp$whole_genome_duplication & wgd_timing_df_cpp$cancer_type == cancer_type & !is.na(wgd_timing_df_cpp$molecular_timing) & !(wgd_timing_df_cpp$is_hypermutated),]
  tmp_df <- tmp_df[!(tmp_df$is_metastatic),]
  p_corr <- cor(tmp_df$molecular_timing,
                log2(tmp_df$tmb))
  model <- lm(log2(tmb) ~ molecular_timing, data = tmp_df)
  sum_model <- summary(model)
  corr_df_prim[i,1:2] <- c(cancer_type, cancer_type_abb)
  corr_df_prim[i,3:6] <- c(p_corr, sum_model$coefficients[2,1], sum_model$coefficients[2,4], sum_model$adj.r.squared)
}

mod_summary_sign <- corr_df_prim$p_val
names(mod_summary_sign) <- corr_df_prim$cancer_type
mod_summary_stars <- NA                             # Named vector with significance stars
mod_summary_stars[mod_summary_sign < 0.1] <- "."
mod_summary_stars[mod_summary_sign < 0.05] <- "*"
mod_summary_stars[mod_summary_sign < 0.01] <- "**"
mod_summary_stars[mod_summary_sign < 0.001] <- "***"
mod_summary_stars[is.na(mod_summary_stars)] <- "n.s."

names(mod_summary_stars) <- names(mod_summary_sign)
corr_df_prim$p_val_sig_code <- as.vector(mod_summary_stars)
corr_df_prim$pearson_correlation <- round(corr_df_prim$pearson_correlation, digits = 2)
corr_df_prim$coefficient <- round(corr_df_prim$coefficient, digits = 2)



wgd_timing_tmb_box_plot_prim <- wgd_timing_df_cpp[wgd_timing_df_cpp$whole_genome_duplication & !is.na(wgd_timing_df_cpp$molecular_timing) & !(wgd_timing_df_cpp$is_hypermutated) & !(wgd_timing_df_cpp$is_metastatic),] %>% ggplot(aes(x = molecular_timing, y = log2(tmb))) + facet_wrap(~ cancer_type) +
  geom_point(size = 0.75) +
  ggrepel::geom_text_repel(data = corr_df_prim, x= 0.85, y = 25, aes(label = p_val_sig_code), size = 6) +
  ggrepel::geom_text_repel(data = corr_df_prim, x= 0.15, y = 4, aes(label = paste("Pearson corr.:", pearson_correlation, ", Line slope:", coefficient)), size = 3) +
  xlab("WGD Molecular Time") +
  geom_smooth(method='lm', formula= y~x, color = "red", se=F) +
  ylab("TMB (log2)") +
  ylim(c(5, 25)) +
  labs(title = "Relationship of WGD timing and TMB",
       subtitle = "Primary Cohort",
       caption = "Hypermutated samples are excluded!")










wgd_timing_df_cpp <- wgd_timing_df
cancer_types <- unique(metadata_included$cancer_type)

for (cancer_type in cancer_types){
  tmp_df <- wgd_timing_df[wgd_timing_df$whole_genome_duplication & !is.na(wgd_timing_df$molecular_timing) & wgd_timing_df$cancer_type == cancer_type & !(wgd_timing_df$is_hypermutated),]
  if (nrow(tmp_df[tmp_df$is_metastatic,]) >= 5 & nrow(tmp_df[!(tmp_df$is_metastatic),]) >= 5){
    
  } else {
    wgd_timing_df_cpp <- wgd_timing_df_cpp[wgd_timing_df_cpp$cancer_type != cancer_type,]
  }
}

wgd_timing_df_cpp$cancer_type <- factor(wgd_timing_df_cpp$cancer_type)

options(scripen = 999)
cancer_types <- as.character(unique(wgd_timing_df_cpp$cancer_type))


corr_df_metas <- data.frame(cancer_type = character(19), cancer_type_abb = character(19), pearson_correlation = numeric(19), coefficient = numeric(19), p_val = numeric(19), adj_r_squared = numeric(19))
for (i in 1:length(unique(wgd_timing_df_cpp$cancer_type))){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included$cancer_type_code[metadata_included$cancer_type == cancer_type])
  tmp_df <- wgd_timing_df_cpp[wgd_timing_df_cpp$whole_genome_duplication & wgd_timing_df_cpp$cancer_type == cancer_type & !is.na(wgd_timing_df_cpp$molecular_timing) & !(wgd_timing_df_cpp$is_hypermutated),]
  tmp_df <- tmp_df[tmp_df$is_metastatic,]
  p_corr <- cor(tmp_df$molecular_timing,
                log2(tmp_df$tmb))
  model <- lm(log2(tmb) ~ molecular_timing, data = tmp_df)
  sum_model <- summary(model)
  corr_df_metas[i,1:2] <- c(cancer_type, cancer_type_abb)
  corr_df_metas[i,3:6] <- c(p_corr, sum_model$coefficients[2,1], sum_model$coefficients[2,4], sum_model$adj.r.squared)
}

mod_summary_sign <- corr_df_metas$p_val
names(mod_summary_sign) <- corr_df_metas$cancer_type
mod_summary_stars <- NA                             # Named vector with significance stars
mod_summary_stars[mod_summary_sign < 0.1] <- "."
mod_summary_stars[mod_summary_sign < 0.05] <- "*"
mod_summary_stars[mod_summary_sign < 0.01] <- "**"
mod_summary_stars[mod_summary_sign < 0.001] <- "***"
mod_summary_stars[is.na(mod_summary_stars)] <- "n.s."

names(mod_summary_stars) <- names(mod_summary_sign)
corr_df_metas$p_val_sig_code <- as.vector(mod_summary_stars)
corr_df_metas$pearson_correlation <- round(corr_df_metas$pearson_correlation, digits = 2)
corr_df_metas$coefficient <- round(corr_df_metas$coefficient, digits = 2)



wgd_timing_tmb_box_plot_metas <- wgd_timing_df_cpp[wgd_timing_df_cpp$whole_genome_duplication & !is.na(wgd_timing_df_cpp$molecular_timing) & !(wgd_timing_df_cpp$is_hypermutated) & wgd_timing_df_cpp$is_metastatic,] %>% ggplot(aes(x = molecular_timing, y = log2(tmb))) + facet_wrap(~ cancer_type) +
  geom_point(size = 0.75) +
  ggrepel::geom_text_repel(data = corr_df_metas, x= 0.85, y = 25, aes(label = p_val_sig_code), size = 6) +
  ggrepel::geom_text_repel(data = corr_df_metas, x= 0.15, y = 4, aes(label = paste("Pearson corr.:", pearson_correlation, ", Line slope:", coefficient)), size = 3) +
  xlab("WGD Molecular Time") +
  geom_smooth(method='lm', formula= y~x, color = "red", se=F) +
  ylab("TMB (log2)") +
  ylim(c(5, 25)) +
  labs(title = "Relationship of WGD timing and TMB",
       subtitle = "Metastatic Cohort",
       caption = "Hypermutated samples are excluded!")






com_plot <- ggarrange(wgd_timing_tmb_box_plot_prim, wgd_timing_tmb_box_plot_metas, nrow = 2)




for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig8.1-tmb-timing-wgd-vs-nonwgd-scatter-plot-per-cohort.png", height = 920, width = 920)
    print(com_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig8.1-tmb-timing-wgd-vs-nonwgd-scatter-plot-per-cohort.pdf", height = 14, width = 14)
    print(com_plot)
    dev.off()
  }
}





corr_df_prim$cancer_type <- paste0(corr_df_prim$cancer_type, " (", corr_df_prim$cancer_type_abb, ")")
corr_df_prim$cancer_type <- factor(corr_df_prim$cancer_type)
corr_df_prim$cancer_type_abb[which(corr_df_prim$p_val_sig_code == "n.s.")] <- NA

scatt_plot_prim <- corr_df_prim %>% ggplot(aes(x= coefficient, y = abs(pearson_correlation), size = -log10(p_val), color = cancer_type)) +
  geom_point() +
  xlim(c(-4.5,4.5))+
  ggrepel::geom_text_repel(aes(label = cancer_type_abb), max.overlaps = 20) +
  labs(title = "Correlation between WGD timing and TMB",
       subtitle = "Primary Cohort",
       caption = "Hypermutated samples are excluded!") +
  geom_hline(yintercept = 0.2, color = "red", lty = 2, alpha = 0.5) +
  geom_vline(xintercept = 0, color = "red", lty = 2, alpha = 0.5) +
  xlab("WGD timing vs TMB (Line Slopes)") +
  ylab("Pearson Coefficient (ABS)")





corr_df_metas$cancer_type <- paste0(corr_df_metas$cancer_type, " (", corr_df_metas$cancer_type_abb, ")")
corr_df_metas$cancer_type <- factor(corr_df_metas$cancer_type)
corr_df_metas$cancer_type_abb[which(corr_df_metas$p_val_sig_code == "n.s.")] <- NA

scatt_plot_metas <- corr_df_metas %>% ggplot(aes(x= coefficient, y = abs(pearson_correlation), size = -log10(p_val), color = cancer_type)) +
  geom_point() +
  xlim(c(-4.5,4.5))+
  ggrepel::geom_text_repel(aes(label = cancer_type_abb), max.overlaps = 20) +
  labs(title = "Correlation between WGD timing and TMB",
       subtitle = "Metastatic Cohort",
       caption = "Hypermutated samples are excluded!") +
  geom_hline(yintercept = 0.2, color = "red", lty = 2, alpha = 0.5) +
  geom_vline(xintercept = 0, color = "red", lty = 2, alpha = 0.5) +
  xlab("WGD timing vs TMB (Line Slopes)") +
  ylab("Pearson Coefficient (ABS)")


com_plot <- ggarrange(scatt_plot_prim, scatt_plot_metas, ncol = 2, common.legend = T, legend = "bottom")



for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig8.2-tmb-timing-wgd-vs-nonwgd-scatter-plot-per-cohort.png", height = 690, width = 1265)
    print(com_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig8.2-tmb-timing-wgd-vs-nonwgd-scatter-plot-per-cohort.pdf", height = 10.5, width = 19.25)
    print(com_plot)
    dev.off()
  }
}


# ************************************************************************************************************************************************
# fig 9 Fraction of genome altered




# Fraction of genome altered

purple_purity_df <- read.csv(file = paste0(wd, "r-objects/purple-purity-fraction-changed.txt.gz"), header = T, stringsAsFactors = F, sep = "\t")


metadata_included <- merge(metadata_included, purple_purity_df, by = "sample_id")





metadata_included_cp <- metadata_included
cancer_types <- unique(metadata_included$cancer_type)

for (cancer_type in cancer_types){
  tmp_df <- metadata_included[metadata_included$cancer_type == cancer_type,]
  if (nrow(tmp_df[tmp_df$is_metastatic,]) >= 5 & nrow(tmp_df[!(tmp_df$is_metastatic),]) >= 5){
    
  } else {
    metadata_included_cp <- metadata_included_cp[metadata_included_cp$cancer_type != cancer_type,]
  }
}



metadata_included_cp$cancer_type <- paste0(metadata_included_cp$cancer_type, " (", metadata_included_cp$cancer_type_code , ")")




cancer_types <- unique(metadata_included_cp$cancer_type)
nr_row <- length(cancer_types)


ll <- data.frame(cancer_type = character(nr_row), cancer_type_abb = character(nr_row), cancer_type_abb_text = character(nr_row), mean_diploid_prim = numeric(nr_row), mean_diploid_metas = numeric(nr_row))


for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included_cp$cancer_type_code[metadata_included_cp$cancer_type == cancer_type])
  tmp_df <- metadata_included_cp[metadata_included_cp$cancer_type == cancer_type,]
  ll[i,1:3] <- c(cancer_type, cancer_type_abb, cancer_type_abb)
  
  res <- wilcox.test(tmp_df$diploidProportion[!(tmp_df$is_metastatic)], tmp_df$diploidProportion[tmp_df$is_metastatic])
  if (res$p.value > 0.05){
    ll[i,3] <- NA
  }
  ll[i,4:5] <- c(mean((1 - tmp_df$diploidProportion[!(tmp_df$is_metastatic)]), na.rm = T), mean((1- tmp_df$diploidProportion[tmp_df$is_metastatic]), na.rm = T))
}

# ll$cancer_type <- paste0(ll$cancer_type, " (", ll$cancer_type_abb, ")")
ll$cancer_type <- factor(ll$cancer_type)
ll$cancer_type_abb <- factor(ll$cancer_type_abb)

plot_mean <- ll %>% ggplot(aes(x = 100*mean_diploid_prim, y = 100*mean_diploid_metas, color = cancer_type)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb_text), max.overlaps = 20) +
  xlim(c(0,100)) +
  xlab("mean of Fraction of Genome Altered (%) _ Primary Cohort") +
  ylim(c(0,100)) +
  ylab("mean of Fraction of Genome Altered (%) _ Metastatic Cohort") +
  ggtitle("Altered Fraction mean") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw()



# creating the boxplot


for (cancer_type in cancer_types){
  tmp_df <- wgd_timing_df[wgd_timing_df$cancer_type == cancer_type & !(wgd_timing_df$is_hypermutated),]
  if (nrow(tmp_df[tmp_df$whole_genome_duplication,]) >= 5 & nrow(tmp_df[!tmp_df$whole_genome_duplication,]) >= 5){
    
  } else {
    wgd_timing_df_cp <- wgd_timing_df_cp[wgd_timing_df_cp$cancer_type != cancer_type,]
  }
}

wgd_timing_df_cp$cancer_type <- factor(wgd_timing_df_cp$cancer_type)


wgd_tmb_box_plot <- metadata_included_cp %>% ggplot(aes(x = is_metastatic, y = 100 * (1-diploidProportion), color = is_metastatic)) + facet_wrap(~ cancer_type_code) +
  geom_boxplot() +
  xlab("Cohort") +
  stat_compare_means(comparisons = list(c('FALSE', 'TRUE')), method = "wilcox.test", label = "p.signif") +
  ylab("Fraction Genome Altered (%)") +
  ylim(c(0, 120))
# ggtitle("Distribution of clonality percentage per cancer type (mutationTimerR)") 




com_plot <- ggarrange(plot_mean, plot_median, wgd_tmb_box_plot, ncol=3, common.legend = TRUE, legend="bottom")

com_plot <- annotate_figure(
  com_plot,
  bottom = "Hypermutated samples are excluded!")





for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig9-fraction-pf-genome-altered.png", height = 690, width = 1265)
    print(com_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig9-fraction-pf-genome-altered.pdf", height = 10.5, width = 19.25)
    print(com_plot)
    dev.off()
  }
}












metadata_included_cpp <- metadata_included
cancer_types <- unique(metadata_included$cancer_type)

for (cancer_type in cancer_types){
  tmp_df <- metadata_included[metadata_included$cancer_type == cancer_type,]
  if (nrow(tmp_df[tmp_df$is_metastatic & !(tmp_df$whole_genome_duplication),]) >= 5 & nrow(tmp_df[!(tmp_df$is_metastatic) & !(tmp_df$whole_genome_duplication),]) >= 5){
    
  } else {
    metadata_included_cpp <- metadata_included_cpp[metadata_included_cpp$cancer_type != cancer_type,]
  }
}



metadata_included_cpp$cancer_type <- paste0(metadata_included_cpp$cancer_type, " (", metadata_included_cpp$cancer_type_code , ")")




cancer_types <- unique(metadata_included_cpp$cancer_type)
nr_row <- length(cancer_types)


ll <- data.frame(cancer_type = character(nr_row), cancer_type_abb = character(nr_row), cancer_type_abb_text = character(nr_row), median_diploid_prim = numeric(nr_row), median_diploid_metas = numeric(nr_row))


for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included_cpp$cancer_type_code[metadata_included_cpp$cancer_type == cancer_type])
  tmp_df <- metadata_included_cpp[metadata_included_cpp$cancer_type == cancer_type & !(metadata_included_cpp$whole_genome_duplication),]
  ll[i,1:3] <- c(cancer_type, cancer_type_abb, cancer_type_abb)
  
  res <- wilcox.test(tmp_df$diploidProportion[!(tmp_df$is_metastatic)], tmp_df$diploidProportion[tmp_df$is_metastatic])
  if (res$p.value > 0.05){
    ll[i,3] <- NA
  }
  ll[i,4:5] <- c(median((1 - tmp_df$diploidProportion[!(tmp_df$is_metastatic)]), na.rm = T), median((1- tmp_df$diploidProportion[tmp_df$is_metastatic]), na.rm = T))
}

# ll$cancer_type <- paste0(ll$cancer_type, " (", ll$cancer_type_abb, ")")
ll$cancer_type <- factor(ll$cancer_type)
ll$cancer_type_abb <- factor(ll$cancer_type_abb)

plott_median <- ll %>% ggplot(aes(x = 100*median_diploid_prim, y = 100*median_diploid_metas, color = cancer_type)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb_text), max.overlaps = 20) +
  xlim(c(0,60)) +
  xlab("median of Fraction of Genome Altered (%) _ Primary Cohort") +
  ylim(c(0,60)) +
  ylab("median of Fraction of Genome Altered (%) _ Metastatic Cohort") +
  ggtitle("Altered Fraction median") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw()



# creating the boxplott



wgd_tmb_box_plott <- metadata_included_cpp[!(metadata_included_cpp$whole_genome_duplication),] %>% ggplot(aes(x = is_metastatic, y = 100 * (1-diploidProportion), color = is_metastatic)) + facet_wrap(~ cancer_type_code) +
  geom_boxplot() +
  xlab("Cohort") +
  stat_compare_means(comparisons = list(c('FALSE', 'TRUE')), method = "wilcox.test", label = "p.signif") +
  ylab("Fraction Genome Altered (%)") +
  ylim(c(0, 120))
# ggtitle("Distribution of clonality percentage per cancer type (mutationTimerR)") 




com_plott <- ggarrange(plott_mean, plott_median, wgd_tmb_box_plott, ncol=3, common.legend = TRUE, legend="bottom")

com_plott <- annotate_figure(
  com_plott,
  bottom = "WGD samples are excluded!")





for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig9-1-fraction-pf-genome-altered-without-wgd.png", height = 690, width = 1265)
    print(com_plott)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig9-1-fraction-pf-genome-altered-without-wgd.pdf", height = 10.5, width = 19.25)
    print(com_plott)
    dev.off()
  }
}







# ************************************************************************************************************************************************
# fig 10 ms timing


if (dir.exists("/hpc/cuppen/")){
  cosmic_sigs <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/cosmic-sigs.txt", sep = "\t", header = F)
} else {
  cosmic_sigs <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/cosmic-sigs.txt", sep = "\t", header = F)
}



if (dir.exists("/hpc/cuppen/")){
  ms_timing_df <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-df.txt.gz", sep = "\t", header = T, stringsAsFactors = F)
} else {
  ms_timing_df <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-df.txt.gz", sep = "\t", header = T, stringsAsFactors = F)
}

# df <- data.frame(sample_id = character(7049*51), cosmic_sig = character(7049*51), fold_change = numeric(7049*51))
# 
# 
# # length(metadata_included$sample_id)
# for (i in 1:length(metadata_included$sample_id)){
#   print(i)
#   sample_id <- metadata_included$sample_id[i]
#   i <- i + (50*(i-1))
#   
#   for (j in 1:length(cosmic_sigs$V1)){
#     
#     df[i+j-1, "sample_id"] <- sample_id
#     df[i+j-1, "cosmic_sig"] <- cosmic_sigs$V1[j]
#     late_count <-  ms_timing_df[ms_timing_df$sample_id == sample_id & ms_timing_df$timing == "clonal [late]",cosmic_sigs$V1[j]]
#     early_count <- ms_timing_df[ms_timing_df$sample_id == sample_id & ms_timing_df$timing == "clonal [early]",cosmic_sigs$V1[j]]
#     
#     if (!is.na(late_count) & !is.na(early_count) & late_count != 0 & early_count != 0){
#       df[i+j-1, "fold_change"] <- late_count/early_count
#     } else {
#       df[i+j-1, "fold_change"] <- NA
#     }
#   }
# }
# saveRDS(df, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/df.all.rds"))


df <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/df.all.rds"))


id8_sig <- df$sample_id[df$cosmic_sig == "ID18" & !is.na(df$fold_change)]

table(metadata_included$cancer_type[metadata_included$sample_id %in% id8_sig])


sbs88_sig <- df$sample_id[df$cosmic_sig == "SBS88" & !is.na(df$fold_change)]

table(metadata_included$cancer_type[metadata_included$sample_id %in% sbs88_sig])


summary(log2(df[!is.na(df$fold_change),"fold_change"]))


df$cosmic_sig <- factor(df$cosmic_sig)
# df[df$fold_change > 50 & !is.na(df$fold_change), "fold_change"] <- 50

gg_plot <- df[!is.na(df$fold_change),] %>% ggplot(aes(x = reorder(cosmic_sig,log2(fold_change),na.rm = TRUE), y = log2(fold_change))) +
  geom_boxplot() +
  geom_text(x = 40, y = -10, label = "Median of medians = -0.75", color = "red") +
  geom_text(x= 5, y = 7.5, label = "Late clonal", face = "bold") +
  geom_text(x= 5, y = -10, label = "Early clonal", face = "bold") +
  # scale_color_manual(guide = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = -0.7469, color = "red", lty = 2) +
  xlab("Cosmic Signature") +
  ylab("Fold Change (Log2)")




# Normalized

1/0.59589


df2 <- df

df2$fold_change[!is.na(df2$fold_change)] <- df2$fold_change[!is.na(df2$fold_change)] * 1.67862



gg_plot2 <- df2[!is.na(df2$fold_change),] %>% ggplot(aes(x = reorder(cosmic_sig,log2(fold_change),na.rm = TRUE), y = log2(fold_change))) +
  geom_boxplot() +
  geom_text(x= 5, y = 7.5, label = "Late clonal", face = "bold") +
  geom_text(x= 5, y = -10, label = "Early clonal", face = "bold") +
  # scale_color_manual(guide = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  xlab("Cosmic Signature") +
  ylab("Normalized Fold Change (Log2)")





for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig10-ms-timing.png", height = 480, width = 960)
    print(gg_plot)
    dev.off()
    
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig10-ms-timing-normalized.png", height = 480, width = 960)
    print(gg_plot2)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig10-ms-timing.pdf", height = 7, width = 14)
    print(gg_plot)
    dev.off()
    
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig10-ms-timing-normalized.pdf", height = 7, width = 14)
    print(gg_plot2)
    dev.off()
  }
}



### between metastatic and primary


df <- merge(df, metadata_included[,c("sample_id", "is_metastatic")], by = "sample_id")

head(df)


summary(log2(df[!is.na(df$fold_change),"fold_change"]))


dff <- df

cosmic_sigs <- as.character(unique(dff$cosmic_sig))

for (cosmic_sig in cosmic_sigs) {
  if (nrow(dff[dff$cosmic_sig == cosmic_sig & dff$is_metastatic & !is.na(dff$fold_change),]) < 3 | nrow(dff[dff$cosmic_sig == cosmic_sig & !(dff$is_metastatic) & !is.na(dff$fold_change),]) < 3){
    dff <- dff[dff$cosmic_sig != cosmic_sig,]
  }
}

dff$cosmic_sig <- factor(dff$cosmic_sig)

cosmic_sigs <- as.character(unique(dff$cosmic_sig))


test_df <- data.frame(cosmic_sig = character(48), p_value = numeric(48), effect_size = character(48))

for (i in 1:length(cosmic_sigs)) {
  cosmic_sig <- cosmic_sigs[i]
  print(cosmic_sig)
  test_df[i,1] <- cosmic_sig
  
  prim_fold_change <- dff[dff$cosmic_sig == cosmic_sig & !is.na(dff$fold_change) & !(dff$is_metastatic),"fold_change"]
  metas_fold_change <- dff[dff$cosmic_sig == cosmic_sig & !is.na(dff$fold_change) & dff$is_metastatic,"fold_change"]
  
  
  wilcox_res <- wilcox.test(prim_fold_change, metas_fold_change)
  test_df[i,2] <- wilcox_res$p.value
  print(wilcox_res)
  
  effsize_res <- cohen.d(prim_fold_change, metas_fold_change, hedges.correction=T)
  test_df[i,3] <- as.character(effsize_res$magnitude)
  print(effsize_res)
}



mod_summary_sign <- test_df$p_value
mod_summary_stars <- NA                             # Named vector with significance stars
mod_summary_stars[mod_summary_sign < 0.1] <- "."
mod_summary_stars[mod_summary_sign < 0.05] <- "*"
mod_summary_stars[mod_summary_sign < 0.01] <- "**"
mod_summary_stars[mod_summary_sign < 0.001] <- "***"
mod_summary_stars[is.na(mod_summary_stars)] <- "n.s."

test_df$p_val_sig <- mod_summary_stars

# test_df$p_val_sig[which(test_df$effect_size == "negligible")] <- NA
# test_df$p_val_sig[which(test_df$p_val_sig == "n.s.")] <- NA


dff$cosmic_sig <- factor(dff$cosmic_sig)

test_df$effect_size <- factor(test_df$effect_size, levels = c("large", "medium", "small", "negligible"))

# adding mean data
ms_df <- ddply(ms_timing_df[,-2], "sample_id", numcolwise(sum, na.rm = T))


ms_df[,2:52] <- round(ms_df[,2:52]/(rowSums(ms_df[,2:52])+0.000001), digits = 2)


sig_means <- data.frame(cosmic_sig = names(colMeans(ms_df[,2:52])), mean_rel_cont = as.numeric(colMeans(ms_df[,2:52])))

dff <- merge(dff, sig_means, by = "cosmic_sig")
dff$mean_rel_cont <- as.numeric(dff$mean_rel_cont)
dff$mean_rel_cont <- dff$mean_rel_cont*100

dff$mean_rel_cont_category <- NA
dff[which(dff$mean_rel_cont >= 5),"mean_rel_cont_category"] <- "High"
dff[which(dff$mean_rel_cont >= 1 & dff$mean_rel_cont < 5),"mean_rel_cont_category"] <- "Medium"
dff[which(dff$mean_rel_cont >= 0.1 & dff$mean_rel_cont < 1),"mean_rel_cont_category"] <- "Low"
dff[which(dff$mean_rel_cont < 0.1),"mean_rel_cont_category"] <- "Trivial"

dff$mean_rel_cont_category <- factor(dff$mean_rel_cont_category, levels = c("Trivial", "Low", "Medium", "High"))


gg_plott <- dff[!is.na(dff$fold_change),] %>% ggplot(aes(x = reorder(cosmic_sig,log2(fold_change),na.rm = TRUE), y = log2(fold_change), fill = is_metastatic)) +
  geom_boxplot(aes(alpha = mean_rel_cont_category)) +
  scale_alpha_discrete(labels = c("High (x>=5%)","Medium (1%<x<5%)","Low (0.1%<=x<1%)","Trivial (x<=<0.1%)"), name = "Mean Relative Contribution") +
  scale_fill_discrete(labels = c("Primary", "Metastatic"), name = "Cohort") +
  geom_text(x = 40, y = -10, label = "Median of medians = -0.75", color = "red") +
  geom_text(inherit.aes = F, data = test_df, aes(x = cosmic_sig, y = 9.5, label = p_val_sig, color = effect_size)) +
  scale_color_discrete(labels = c("Large", "Medium", "Small", "Negligible"), name = "Effect Size (Cohen's d)") +
  # stat_compare_means(aes(group = is_metastatic), method = "wilcox.test", label = "p.signif") +
  geom_text(x= 5, y = 7.5, label = "Late clonal", face = "bold") +
  geom_text(x= 5, y = -10, label = "Early clonal", face = "bold") +
  # scale_color_manual(guide = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = -0.7469, color = "red", lty = 2) +
  xlab("Cosmic Signature") +
  ylab("Fold Change (Log2)") +
  labs(caption = "Wilcox Test Used!", title = "Mutational processes during early and late clonal tumour evolution") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))



# Normalized

dff2 <- df

cosmic_sigs <- as.character(unique(df$cosmic_sig))

for (cosmic_sig in cosmic_sigs) {
  if (nrow(dff2[dff2$cosmic_sig == cosmic_sig & dff2$is_metastatic & !is.na(dff2$fold_change),]) < 3 | nrow(dff2[dff2$cosmic_sig == cosmic_sig & !(dff2$is_metastatic) & !is.na(dff2$fold_change),]) < 3){
    dff2 <- dff2[dff2$cosmic_sig != cosmic_sig,]
  }
}

dff2$cosmic_sig <- factor(dff2$cosmic_sig)
dff2$fold_change[!is.na(dff2$fold_change)] <- dff2$fold_change[!is.na(dff2$fold_change)] * 1.67862



dff2 <- merge(dff2, sig_means, by = "cosmic_sig")
dff2$mean_rel_cont <- as.numeric(dff2$mean_rel_cont)
dff2$mean_rel_cont <- dff2$mean_rel_cont*100

dff2$mean_rel_cont_category <- NA
dff2[which(dff2$mean_rel_cont >= 5),"mean_rel_cont_category"] <- "High"
dff2[which(dff2$mean_rel_cont >= 1 & dff2$mean_rel_cont < 5),"mean_rel_cont_category"] <- "Medium"
dff2[which(dff2$mean_rel_cont >= 0.1 & dff2$mean_rel_cont < 1),"mean_rel_cont_category"] <- "Low"
dff2[which(dff2$mean_rel_cont < 0.1),"mean_rel_cont_category"] <- "Trivial"

dff2$mean_rel_cont_category <- factor(dff2$mean_rel_cont_category, levels = c("Trivial", "Low", "Medium", "High"))




gg_plott2 <- dff2[!is.na(dff2$fold_change),] %>% ggplot(aes(x = reorder(cosmic_sig,log2(fold_change),na.rm = TRUE), y = log2(fold_change), fill = is_metastatic)) +
  geom_boxplot(aes(alpha = mean_rel_cont_category)) +
  scale_alpha_discrete(labels = c("High (x>=5%)","Medium (1%<x<5%)","Low (0.1%<=x<1%)","Trivial (x<=<0.1%)"), name = "Mean Relative Contribution") +
  scale_fill_discrete(labels = c("Primary", "Metastatic"), name = "Cohort") +
  geom_text(x= 5, y = 7.5, label = "Late clonal", face = "bold") +
  geom_text(x= 5, y = -10, label = "Early clonal", face = "bold") +
  geom_text(inherit.aes = F, data = test_df, aes(x = cosmic_sig, y = 9.5, label = p_val_sig, color = effect_size)) +
  scale_color_discrete(labels = c("Large", "Medium", "Small", "Negligible"), name = "Effect Size (Cohen's d)") +
  # scale_color_manual(guide = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  xlab("Cosmic Signature") +
  ylab("Normalized Fold Change (Log2)") +
  labs(caption = "Wilcox Test Used!", title = "Mutational processes during early and late clonal tumour evolution") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))




for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig10-ms-timing-per-cohort.png", height = 480, width = 960)
    print(gg_plott)
    dev.off()
    
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig10-ms-timing-per-cohort-normalized.png", height = 480, width = 960)
    print(gg_plott2)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig10-ms-timing-per-cohort.pdf", height = 7, width = 14)
    print(gg_plott)
    dev.off()
    
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig10-ms-timing-per-cohort-normalized.pdf", height = 7, width = 14)
    print(gg_plott2)
    dev.off()
  }
}





# ************************************************************************************************************************************************
# fig11 ms clonality purple


# df2 is when subclonal considered (subclonal + clonal [late]) and df3 when only subclonal category
# subclonal <-  ms_timing_df[ms_timing_df$sample_id == sample_id & ms_timing_df$timing == "clonal [late]",cosmic_sigs$V1[j]] + ms_timing_df[ms_timing_df$sample_id == sample_id & ms_timing_df$timing == "subclonal",cosmic_sigs$V1[j]]


# df3 <- data.frame(sample_id = character(7049*51), cosmic_sig = character(7049*51), fold_change = numeric(7049*51))
# 
# 
# # length(metadata_included$sample_id)
# for (i in 1:length(metadata_included$sample_id)){
#   print(i)
#   sample_id <- metadata_included$sample_id[i]
#   i <- i + (50*(i-1))
# 
#   for (j in 1:length(cosmic_sigs$V1)){
# 
#     df3[i+j-1, "sample_id"] <- sample_id
#     df3[i+j-1, "cosmic_sig"] <- cosmic_sigs$V1[j]
#     subclonal <- ms_timing_df[ms_timing_df$sample_id == sample_id & ms_timing_df$timing == "subclonal",cosmic_sigs$V1[j]]
#     clonal <- ms_timing_df[ms_timing_df$sample_id == sample_id & ms_timing_df$timing == "clonal [early]",cosmic_sigs$V1[j]] + ms_timing_df[ms_timing_df$sample_id == sample_id & ms_timing_df$timing == "clonal [NA]",cosmic_sigs$V1[j]] + ms_timing_df[ms_timing_df$sample_id == sample_id & ms_timing_df$timing == "clonal [late]",cosmic_sigs$V1[j]]
# 
#     if (!is.na(subclonal) & !is.na(clonal) & subclonal != 0 & clonal != 0){
#       df3[i+j-1, "fold_change"] <- subclonal/clonal
#     } else {
#       df3[i+j-1, "fold_change"] <- NA
#     }
#   }
# }
# saveRDS(df3, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/df3.all.rds"))


df3 <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/df3.all.rds"))

summary(log2(df3[!is.na(df3$fold_change),"fold_change"]))

df3 <- merge(df3, metadata_included[,c("sample_id", "is_metastatic")], by = "sample_id")

head(df3[!is.na(df3$fold_change),])





dff3 <- df3

cosmic_sigs <- as.character(unique(dff3$cosmic_sig))

for (cosmic_sig in cosmic_sigs) {
  if (nrow(dff3[dff3$cosmic_sig == cosmic_sig & dff3$is_metastatic & !is.na(dff3$fold_change),]) < 3 | nrow(dff3[dff3$cosmic_sig == cosmic_sig & !(dff3$is_metastatic) & !is.na(dff3$fold_change),]) < 3){
    dff3 <- dff3[dff3$cosmic_sig != cosmic_sig,]
  }
}

dff3$cosmic_sig <- factor(dff3$cosmic_sig)

cosmic_sigs <- as.character(unique(dff3$cosmic_sig))


# To normalize
summary(df3[!is.na(df3$fold_change),"fold_change"])
1/0.15954
dff3$fold_change[!is.na(dff3$fold_change)] <- dff3$fold_change[!is.na(dff3$fold_change)] * 6.268



row_no <- length(cosmic_sigs)
test_df <- data.frame(cosmic_sig = character(row_no), p_value = numeric(row_no), effect_size = character(row_no))

for (i in 1:length(cosmic_sigs)) {
  cosmic_sig <- cosmic_sigs[i]
  print(cosmic_sig)
  test_df[i,1] <- cosmic_sig
  
  prim_fold_change <- dff3[dff3$cosmic_sig == cosmic_sig & !is.na(dff3$fold_change) & !(dff3$is_metastatic),"fold_change"]
  metas_fold_change <- dff3[dff3$cosmic_sig == cosmic_sig & !is.na(dff3$fold_change) & dff3$is_metastatic,"fold_change"]
  
  
  wilcox_res <- wilcox.test(prim_fold_change, metas_fold_change)
  test_df[i,2] <- wilcox_res$p.value
  print(wilcox_res)
  
  effsize_res <- cohen.d(prim_fold_change, metas_fold_change, hedges.correction=T)
  test_df[i,3] <- as.character(effsize_res$magnitude)
  print(effsize_res)
}



mod_summary_sign <- test_df$p_value
mod_summary_stars <- NA                             # Named vector with significance stars
mod_summary_stars[mod_summary_sign < 0.1] <- "."
mod_summary_stars[mod_summary_sign < 0.05] <- "*"
mod_summary_stars[mod_summary_sign < 0.01] <- "**"
mod_summary_stars[mod_summary_sign < 0.001] <- "***"
mod_summary_stars[is.na(mod_summary_stars)] <- "n.s."

test_df$p_val_sig <- mod_summary_stars

# test_df$p_val_sig[which(test_df$effect_size == "negligible")] <- NA
# test_df$p_val_sig[which(test_df$p_val_sig == "n.s.")] <- NA


dff3$cosmic_sig <- factor(dff3$cosmic_sig)

test_df$effect_size <- factor(test_df$effect_size, levels = c("large", "medium", "small", "negligible"))

# adding mean data
ms_df <- ddply(ms_timing_df[,-2], "sample_id", numcolwise(sum, na.rm = T))


ms_df[,2:52] <- round(ms_df[,2:52]/(rowSums(ms_df[,2:52])+0.000001), digits = 2)


sig_means <- data.frame(cosmic_sig = names(colMeans(ms_df[,2:52])), mean_rel_cont = as.numeric(colMeans(ms_df[,2:52])))

dff3 <- merge(dff3, sig_means, by = "cosmic_sig")
dff3$mean_rel_cont <- as.numeric(dff3$mean_rel_cont)
dff3$mean_rel_cont <- dff3$mean_rel_cont*100

dff3$mean_rel_cont_category <- NA
dff3[which(dff3$mean_rel_cont >= 5),"mean_rel_cont_category"] <- "High"
dff3[which(dff3$mean_rel_cont >= 1 & dff3$mean_rel_cont < 5),"mean_rel_cont_category"] <- "Medium"
dff3[which(dff3$mean_rel_cont >= 0.1 & dff3$mean_rel_cont < 1),"mean_rel_cont_category"] <- "Low"
dff3[which(dff3$mean_rel_cont < 0.1),"mean_rel_cont_category"] <- "Trivial"

dff3$mean_rel_cont_category <- factor(dff3$mean_rel_cont_category, levels = c("Trivial", "Low", "Medium", "High"))


gg_plott <- dff3[!is.na(dff3$fold_change),] %>% ggplot(aes(x = reorder(cosmic_sig,log2(fold_change),na.rm = TRUE), y = log2(fold_change), fill = is_metastatic)) +
  geom_boxplot(aes(alpha = mean_rel_cont_category)) +
  scale_alpha_discrete(labels = c("High (x>=5%)","Medium (1%<x<5%)","Low (0.1%<=x<1%)","Trivial (x<=<0.1%)"), name = "Mean Relative Contribution") +
  scale_fill_discrete(labels = c("Primary", "Metastatic"), name = "Cohort") +
  geom_text(x = 30, y = -10, label = "Median of medians = -2.648", color = "red") +
  geom_text(inherit.aes = F, data = test_df, aes(x = cosmic_sig, y = 9.5, label = p_val_sig, color = effect_size)) +
  scale_color_discrete(labels = c("Large", "Medium", "Small", "Negligible"), name = "Effect Size (Cohen's d)") +
  # stat_compare_means(aes(group = is_metastatic), method = "wilcox.test", label = "p.signif") +
  geom_text(x= 5, y = 7.5, label = "Subclonal", face = "bold") +
  geom_text(x= 5, y = -7.5, label = "Clonal", face = "bold") +
  # scale_color_manual(guide = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = -2.648, color = "red", lty = 2) +
  xlab("Cosmic Signature") +
  ylab("Fold Change (Log2)") +
  ylim(c(-10,10)) +
  labs(caption = "Wilcox Test Used!", title = "Mutational processes during clonal and subclonal tumour evolution") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))



# for normalized plot


gg_plott2 <- dff3[!is.na(dff3$fold_change),] %>% ggplot(aes(x = reorder(cosmic_sig,log2(fold_change),na.rm = TRUE), y = log2(fold_change), fill = is_metastatic)) +
  geom_boxplot(aes(alpha = mean_rel_cont_category)) +
  scale_alpha_discrete(labels = c("High (x>=5%)","Medium (1%<x<5%)","Low (0.1%<=x<1%)","Trivial (x<=<0.1%)"), name = "Mean Relative Contribution") +
  scale_fill_discrete(labels = c("Primary", "Metastatic"), name = "Cohort") +
  geom_text(x= 5, y = 7.5, label = "Subclonal", face = "bold") +
  geom_text(x= 5, y = -7.5, label = "Clonal", face = "bold") +
  geom_text(inherit.aes = F, data = test_df, aes(x = cosmic_sig, y = 9.5, label = p_val_sig, color = effect_size)) +
  scale_color_discrete(labels = c("Large", "Medium", "Small", "Negligible"), name = "Effect Size (Cohen's d)") +
  # scale_color_manual(guide = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  xlab("Cosmic Signature") +
  ylab("Normalized Fold Change (Log2)") +
  # ylim(c(-10,10)) +
  labs(caption = "Wilcox Test Used!", title = "Mutational processes during clonal and subclonal tumour evolution") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))



for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig11-ms-clonality-per-cohort-mutatiotimer-literal.png", height = 480, width = 960)
    print(gg_plott)
    dev.off()
    
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig11-ms-clonality-per-cohort-normalized-mutatiotimer-literal.png", height = 480, width = 960)
    print(gg_plott2)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig11-ms-clonality-per-cohort-mutatiotimer-literal.pdf", height = 7, width = 14)
    print(gg_plott)
    dev.off()
    
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig11-ms-clonality-per-cohort-normalized-mutatiotimer-literal.pdf", height = 7, width = 14)
    print(gg_plott2)
    dev.off()
  }
}




## Fig11-1 making the same plot for binned purple clonality data (leaving out probable-clonal and probable-subclonal)



if (dir.exists("/hpc/cuppen/")){
  cosmic_sigs <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/cosmic-sigs.txt", sep = "\t", header = F)
} else {
  cosmic_sigs <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/cosmic-sigs.txt", sep = "\t", header = F)
}



if (dir.exists("/hpc/cuppen/")){
  ms_clonality_binned_purple <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-binned-purple-df.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
} else {
  ms_clonality_binned_purple <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-binned-purple-df.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
}


# dff <- data.frame(sample_id = character(7049*51), cosmic_sig = character(7049*51), fold_change = numeric(7049*51))
# 
# 
# # length(metadata_included$sample_id)
# for (i in 1:length(metadata_included$sample_id)){
#   print(i)
#   sample_id <- metadata_included$sample_id[i]
#   i <- i + (50*(i-1))
# 
#   for (j in 1:length(cosmic_sigs$V1)){
# 
#     dff[i+j-1, "sample_id"] <- sample_id
#     dff[i+j-1, "cosmic_sig"] <- cosmic_sigs$V1[j]
#     subclonal <- ms_clonality_binned_purple[ms_clonality_binned_purple$sample_id == sample_id & ms_clonality_binned_purple$timing == "subclonal",cosmic_sigs$V1[j]]
#     clonal <- ms_clonality_binned_purple[ms_clonality_binned_purple$sample_id == sample_id & ms_clonality_binned_purple$timing == "clonal",cosmic_sigs$V1[j]]
# 
#     if (!is.na(subclonal) & !is.na(clonal) & subclonal != 0 & clonal != 0){
#       dff[i+j-1, "fold_change"] <- subclonal/clonal
#     } else {
#       dff[i+j-1, "fold_change"] <- NA
#     }
#   }
# }
# saveRDS(dff, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/dff.all.ms-clonality-binned-purple-ratios.rds"))


dff <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/dff.all.ms-clonality-binned-purple-ratios.rds"))

dff <- merge(dff, metadata_included[,c("sample_id", "is_metastatic")], by = "sample_id")

head(dff[!is.na(dff$fold_change),])





dff3 <- dff

cosmic_sigs <- as.character(unique(dff3$cosmic_sig))

for (cosmic_sig in cosmic_sigs) {
  if (nrow(dff3[dff3$cosmic_sig == cosmic_sig & dff3$is_metastatic & !is.na(dff3$fold_change),]) < 3 | nrow(dff3[dff3$cosmic_sig == cosmic_sig & !(dff3$is_metastatic) & !is.na(dff3$fold_change),]) < 3){
    dff3 <- dff3[dff3$cosmic_sig != cosmic_sig,]
  }
}

dff3$cosmic_sig <- factor(dff3$cosmic_sig)

cosmic_sigs <- as.character(unique(dff3$cosmic_sig))


# To normalize
summary(log2(dff[!is.na(dff$fold_change),"fold_change"]))
summary(dff[!is.na(dff$fold_change),"fold_change"])
1/0.09413
dff3$fold_change[!is.na(dff3$fold_change)] <- dff3$fold_change[!is.na(dff3$fold_change)] * 10.62361



row_no <- length(cosmic_sigs)
test_df <- data.frame(cosmic_sig = character(row_no), p_value = numeric(row_no), effect_size = character(row_no))

for (i in 1:length(cosmic_sigs)) {
  cosmic_sig <- cosmic_sigs[i]
  print(cosmic_sig)
  test_df[i,1] <- cosmic_sig
  
  prim_fold_change <- dff3[dff3$cosmic_sig == cosmic_sig & !is.na(dff3$fold_change) & !(dff3$is_metastatic),"fold_change"]
  metas_fold_change <- dff3[dff3$cosmic_sig == cosmic_sig & !is.na(dff3$fold_change) & dff3$is_metastatic,"fold_change"]
  
  
  wilcox_res <- wilcox.test(prim_fold_change, metas_fold_change)
  test_df[i,2] <- wilcox_res$p.value
  print(wilcox_res)
  
  effsize_res <- cohen.d(prim_fold_change, metas_fold_change, hedges.correction=T)
  test_df[i,3] <- as.character(effsize_res$magnitude)
  print(effsize_res)
}



mod_summary_sign <- test_df$p_value
mod_summary_stars <- NA                             # Named vector with significance stars
mod_summary_stars[mod_summary_sign < 0.1] <- "."
mod_summary_stars[mod_summary_sign < 0.05] <- "*"
mod_summary_stars[mod_summary_sign < 0.01] <- "**"
mod_summary_stars[mod_summary_sign < 0.001] <- "***"
mod_summary_stars[is.na(mod_summary_stars)] <- "n.s."

test_df$p_val_sig <- mod_summary_stars

# test_df$p_val_sig[which(test_df$effect_size == "negligible")] <- NA
# test_df$p_val_sig[which(test_df$p_val_sig == "n.s.")] <- NA


dff3$cosmic_sig <- factor(dff3$cosmic_sig)

test_df$effect_size <- factor(test_df$effect_size, levels = c("large", "medium", "small", "negligible"))

# adding mean data
ms_df <- ddply(ms_clonality_binned_purple[,-2], "sample_id", numcolwise(sum, na.rm = T))


ms_df[,2:52] <- round(ms_df[,2:52]/(rowSums(ms_df[,2:52])+0.000001), digits = 2)


sig_means <- data.frame(cosmic_sig = names(colMeans(ms_df[,2:52])), mean_rel_cont = as.numeric(colMeans(ms_df[,2:52])))

dff3 <- merge(dff3, sig_means, by = "cosmic_sig")
dff3$mean_rel_cont <- as.numeric(dff3$mean_rel_cont)
dff3$mean_rel_cont <- dff3$mean_rel_cont*100

dff3$mean_rel_cont_category <- NA
dff3[which(dff3$mean_rel_cont >= 5),"mean_rel_cont_category"] <- "High"
dff3[which(dff3$mean_rel_cont >= 1 & dff3$mean_rel_cont < 5),"mean_rel_cont_category"] <- "Medium"
dff3[which(dff3$mean_rel_cont >= 0.1 & dff3$mean_rel_cont < 1),"mean_rel_cont_category"] <- "Low"
dff3[which(dff3$mean_rel_cont < 0.1),"mean_rel_cont_category"] <- "Trivial"

dff3$mean_rel_cont_category <- factor(dff3$mean_rel_cont_category, levels = c("Trivial", "Low", "Medium", "High"))


gg_plott <- dff3[!is.na(dff3$fold_change),] %>% ggplot(aes(x = reorder(cosmic_sig,log2(fold_change),na.rm = TRUE), y = log2(fold_change), fill = is_metastatic)) +
  geom_boxplot(aes(alpha = mean_rel_cont_category)) +
  scale_alpha_discrete(labels = c("High (x>=5%)","Medium (1%<x<5%)","Low (0.1%<=x<1%)","Trivial (x<=<0.1%)"), name = "Mean Relative Contribution") +
  scale_fill_discrete(labels = c("Primary", "Metastatic"), name = "Cohort") +
  geom_text(x = 30, y = -10, label = "Median of medians = -3.409", color = "red") +
  geom_text(inherit.aes = F, data = test_df, aes(x = cosmic_sig, y = 9.5, label = p_val_sig, color = effect_size)) +
  scale_color_discrete(labels = c("Large", "Medium", "Small", "Negligible"), name = "Effect Size (Cohen's d)") +
  # stat_compare_means(aes(group = is_metastatic), method = "wilcox.test", label = "p.signif") +
  geom_text(x= 5, y = 7.5, label = "Subclonal", face = "bold") +
  geom_text(x= 5, y = -7.5, label = "Clonal", face = "bold") +
  # scale_color_manual(guide = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = -2.648, color = "red", lty = 2) +
  xlab("Cosmic Signature") +
  ylab("Fold Change (Log2)") +
  ylim(c(-10,10)) +
  labs(caption = "Wilcox Test Used!", title = "Mutational processes during clonal and subclonal tumour evolution (PURPLE binned clonality info)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))



# for normalized plot


gg_plott2 <- dff3[!is.na(dff3$fold_change),] %>% ggplot(aes(x = reorder(cosmic_sig,log2(fold_change),na.rm = TRUE), y = log2(fold_change), fill = is_metastatic)) +
  geom_boxplot(aes(alpha = mean_rel_cont_category)) +
  scale_alpha_discrete(labels = c("High (x>=5%)","Medium (1%<x<5%)","Low (0.1%<=x<1%)","Trivial (x<=<0.1%)"), name = "Mean Relative Contribution") +
  scale_fill_discrete(labels = c("Primary", "Metastatic"), name = "Cohort") +
  geom_text(x= 5, y = 7.5, label = "Subclonal", face = "bold") +
  geom_text(x= 5, y = -10, label = "Clonal", face = "bold") +
  geom_text(inherit.aes = F, data = test_df, aes(x = cosmic_sig, y = 9.5, label = p_val_sig, color = effect_size)) +
  scale_color_discrete(labels = c("Large", "Medium", "Small", "Negligible"), name = "Effect Size (Cohen's d)") +
  # scale_color_manual(guide = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  xlab("Cosmic Signature") +
  ylab("Normalized Fold Change (Log2)") +
  # ylim(c(-10,10)) +
  labs(caption = "Wilcox Test Used!", title = "Mutational processes during clonal and subclonal tumour evolution (PURPLE binned clonality info)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))



for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig11-1-ms-clonality-per-cohort-purple-binned.png", height = 480, width = 960)
    print(gg_plott)
    dev.off()
    
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig11-1-ms-clonality-per-cohort-normalized-purple-binned.png", height = 480, width = 960)
    print(gg_plott2)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig11-1-ms-clonality-per-cohort-purple-binned.pdf", height = 7, width = 14)
    print(gg_plott)
    dev.off()
    
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig11-1-ms-clonality-per-cohort-normalized-purple-binned.pdf", height = 7, width = 14)
    print(gg_plott2)
    dev.off()
  }
}




# ************************************************************************************************************************************************
# fig 12 reproducing fig 2a of Gerstung et al.

library(tidyr)

global_timing_info_com_na_rm <- global_timing_info_com[!is.na(global_timing_info_com$subclonal),]

global_timing_info_com_subcl_0_rm <- global_timing_info_com_na_rm[global_timing_info_com_na_rm$subclonal != 0,]

suff_num_tiss <- names(table(global_timing_info_com_subcl_0_rm$tissue_type)[as.vector(table(global_timing_info_com_subcl_0_rm$tissue_type)) >= 10])

global_timing_info_com_subcl_0_rm <- global_timing_info_com_subcl_0_rm[global_timing_info_com_subcl_0_rm$tissue_type %in% suff_num_tiss,]



global_timing_info_com_subcl_0_rm[,c(3,4,5,6,7)] <- global_timing_info_com_subcl_0_rm[,c(3,4,5,6,7)]/rowSums(global_timing_info_com_subcl_0_rm[,c(3,4,5,7)])


all_order_sample <- vector()
global_timing_info_com_subcl_0_rm$early_plus_late <- rowSums(global_timing_info_com_subcl_0_rm[,c(3,4)])

for (tissue in unique(global_timing_info_com_subcl_0_rm$tissue_type)){
  tmp_df <- global_timing_info_com_subcl_0_rm[global_timing_info_com_subcl_0_rm$tissue_type == tissue,]
  tissue_order <- tmp_df[order(tmp_df$early_plus_late),"sample_id"]
  
  
  all_order_sample <- append(all_order_sample, tissue_order)
}


global_timing_info_com_subcl_0_rm$sample_id <- factor(global_timing_info_com_subcl_0_rm$sample_id, levels = all_order_sample)






global_timing_info_com_subcl_0_rm_tibble <-  gather(data = global_timing_info_com_subcl_0_rm, key = "clonality_category", value = "clonality_value", early_clonal, late_clonal, na_clonal, subclonal)


global_timing_info_com_subcl_0_rm_tibble$clonality_value <- 100*global_timing_info_com_subcl_0_rm_tibble$clonality_value
global_timing_info_com_subcl_0_rm_tibble$clonality_category <- factor(global_timing_info_com_subcl_0_rm_tibble$clonality_category)
global_timing_info_com_subcl_0_rm_tibble$clonality_category <- factor(global_timing_info_com_subcl_0_rm_tibble$clonality_category, levels = rev(levels(global_timing_info_com_subcl_0_rm_tibble$clonality_category)))
global_timing_info_com_subcl_0_rm_tibble$tissue_type <- factor(global_timing_info_com_subcl_0_rm_tibble$tissue_type)

global_timing_info_com_subcl_0_rm_tibble


plot_stacked_bar_plot <- global_timing_info_com_subcl_0_rm_tibble %>% ggplot(aes(x = sample_id, y = clonality_value, fill = clonality_category)) + facet_wrap(~ tissue_type, scale = "free_x") +
  geom_bar(position="fill", stat="identity", aes(color=clonality_category)) +
  scale_color_manual(values = c("yellow", "grey", "blue", "red")) +
  scale_fill_manual(values = c("yellow", "grey", "blue", "red")) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())




for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig12-timing-stacked-bar-plot.png", height = 960, width = 1360)
    print(plot_stacked_bar_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig12-timing-stacked-bar-plot.pdf", height = 14, width = 21)
    print(plot_stacked_bar_plot)
    dev.off()
  }
}




### Primary cohort





global_timing_info_com_na_rm <- global_timing_info_com[!is.na(global_timing_info_com$subclonal),]

global_timing_info_com_subcl_0_rm <- global_timing_info_com_na_rm[global_timing_info_com_na_rm$subclonal != 0,]

suff_num_tiss_prim <- names(table(global_timing_info_com_subcl_0_rm$tissue_type[!(global_timing_info_com_subcl_0_rm$is_metastatic)])[as.vector(table(global_timing_info_com_subcl_0_rm$tissue_type[!(global_timing_info_com_subcl_0_rm$is_metastatic)])) >= 10])
suff_num_tiss_metas <- names(table(global_timing_info_com_subcl_0_rm$tissue_type[global_timing_info_com_subcl_0_rm$is_metastatic])[as.vector(table(global_timing_info_com_subcl_0_rm$tissue_type[global_timing_info_com_subcl_0_rm$is_metastatic])) >= 10])

suff_num_tiss_shared <- intersect(suff_num_tiss_metas, suff_num_tiss_prim)

global_timing_info_com_subcl_0_rm <- global_timing_info_com_subcl_0_rm[global_timing_info_com_subcl_0_rm$tissue_type %in% suff_num_tiss_shared,]



global_timing_info_com_subcl_0_rm[,c(3,4,5,6,7)] <- global_timing_info_com_subcl_0_rm[,c(3,4,5,6,7)]/rowSums(global_timing_info_com_subcl_0_rm[,c(3,4,5,7)])

global_timing_info_com_subcl_0_rm_prim <- global_timing_info_com_subcl_0_rm[!(global_timing_info_com_subcl_0_rm$is_metastatic),]


all_order_sample <- vector()
global_timing_info_com_subcl_0_rm_prim$early_plus_late <- rowSums(global_timing_info_com_subcl_0_rm_prim[,c(3,4)])


for (tissue in unique(global_timing_info_com_subcl_0_rm_prim$tissue_type)){
  tmp_df <- global_timing_info_com_subcl_0_rm_prim[global_timing_info_com_subcl_0_rm_prim$tissue_type == tissue,]
  tissue_order <- tmp_df[order(tmp_df$early_plus_late),"sample_id"]
  
  
  all_order_sample <- append(all_order_sample, tissue_order)
}


global_timing_info_com_subcl_0_rm_prim$sample_id <- factor(global_timing_info_com_subcl_0_rm_prim$sample_id, levels = all_order_sample)






global_timing_info_com_subcl_0_rm_prim_tibble <-  gather(data = global_timing_info_com_subcl_0_rm_prim, key = "clonality_category", value = "clonality_value", early_clonal, late_clonal, na_clonal, subclonal)


global_timing_info_com_subcl_0_rm_prim_tibble$clonality_value <- 100*global_timing_info_com_subcl_0_rm_prim_tibble$clonality_value
global_timing_info_com_subcl_0_rm_prim_tibble$clonality_category <- factor(global_timing_info_com_subcl_0_rm_prim_tibble$clonality_category)
global_timing_info_com_subcl_0_rm_prim_tibble$clonality_category <- factor(global_timing_info_com_subcl_0_rm_prim_tibble$clonality_category, levels = rev(levels(global_timing_info_com_subcl_0_rm_prim_tibble$clonality_category)))
global_timing_info_com_subcl_0_rm_prim_tibble$tissue_type <- factor(global_timing_info_com_subcl_0_rm_prim_tibble$tissue_type)

global_timing_info_com_subcl_0_rm_prim_tibble


plot_stacked_bar_plot_prim <- global_timing_info_com_subcl_0_rm_prim_tibble %>% ggplot(aes(x = sample_id, y = clonality_value, fill = clonality_category)) + facet_wrap(~ tissue_type, scale = "free_x") +
  geom_bar(position="fill", stat="identity", aes(color=clonality_category)) +
  scale_color_manual(values = c("yellow", "grey", "blue", "red")) +
  scale_fill_manual(values = c("yellow", "grey", "blue", "red")) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  ggtitle("Timing and CLonality for Primary Cohort")






### Metastatic cohort



global_timing_info_com_na_rm <- global_timing_info_com[!is.na(global_timing_info_com$subclonal),]

global_timing_info_com_subcl_0_rm <- global_timing_info_com_na_rm[global_timing_info_com_na_rm$subclonal != 0,]

suff_num_tiss_prim <- names(table(global_timing_info_com_subcl_0_rm$tissue_type[!(global_timing_info_com_subcl_0_rm$is_metastatic)])[as.vector(table(global_timing_info_com_subcl_0_rm$tissue_type[!(global_timing_info_com_subcl_0_rm$is_metastatic)])) >= 10])
suff_num_tiss_metas <- names(table(global_timing_info_com_subcl_0_rm$tissue_type[global_timing_info_com_subcl_0_rm$is_metastatic])[as.vector(table(global_timing_info_com_subcl_0_rm$tissue_type[global_timing_info_com_subcl_0_rm$is_metastatic])) >= 10])

suff_num_tiss_shared <- intersect(suff_num_tiss_metas, suff_num_tiss_prim)

global_timing_info_com_subcl_0_rm <- global_timing_info_com_subcl_0_rm[global_timing_info_com_subcl_0_rm$tissue_type %in% suff_num_tiss_shared,]



global_timing_info_com_subcl_0_rm[,c(3,4,5,6,7)] <- global_timing_info_com_subcl_0_rm[,c(3,4,5,6,7)]/rowSums(global_timing_info_com_subcl_0_rm[,c(3,4,5,7)])

global_timing_info_com_subcl_0_rm_metas <- global_timing_info_com_subcl_0_rm[global_timing_info_com_subcl_0_rm$is_metastatic,]


all_order_sample <- vector()
global_timing_info_com_subcl_0_rm_metas$early_plus_late <- rowSums(global_timing_info_com_subcl_0_rm_metas[,c(3,4)])


for (tissue in unique(global_timing_info_com_subcl_0_rm_metas$tissue_type)){
  tmp_df <- global_timing_info_com_subcl_0_rm_metas[global_timing_info_com_subcl_0_rm_metas$tissue_type == tissue,]
  tissue_order <- tmp_df[order(tmp_df$early_plus_late),"sample_id"]
  
  
  all_order_sample <- append(all_order_sample, tissue_order)
}


global_timing_info_com_subcl_0_rm_metas$sample_id <- factor(global_timing_info_com_subcl_0_rm_metas$sample_id, levels = all_order_sample)






global_timing_info_com_subcl_0_rm_metas_tibble <-  gather(data = global_timing_info_com_subcl_0_rm_metas, key = "clonality_category", value = "clonality_value", early_clonal, late_clonal, na_clonal, subclonal)


global_timing_info_com_subcl_0_rm_metas_tibble$clonality_value <- 100*global_timing_info_com_subcl_0_rm_metas_tibble$clonality_value
global_timing_info_com_subcl_0_rm_metas_tibble$clonality_category <- factor(global_timing_info_com_subcl_0_rm_metas_tibble$clonality_category)
global_timing_info_com_subcl_0_rm_metas_tibble$clonality_category <- factor(global_timing_info_com_subcl_0_rm_metas_tibble$clonality_category, levels = rev(levels(global_timing_info_com_subcl_0_rm_metas_tibble$clonality_category)))
global_timing_info_com_subcl_0_rm_metas_tibble$tissue_type <- factor(global_timing_info_com_subcl_0_rm_metas_tibble$tissue_type)

global_timing_info_com_subcl_0_rm_metas_tibble


plot_stacked_bar_plot_metas <- global_timing_info_com_subcl_0_rm_metas_tibble %>% ggplot(aes(x = sample_id, y = clonality_value, fill = clonality_category)) + facet_wrap(~ tissue_type, scale = "free_x") +
  geom_bar(position="fill", stat="identity", aes(color=clonality_category)) +
  scale_color_manual(values = c("yellow", "grey", "blue", "red")) +
  scale_fill_manual(values = c("yellow", "grey", "blue", "red")) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  ggtitle("Timing and CLonality for Metastatic Cohort")





combined_plot <- ggarrange(plot_stacked_bar_plot_prim, plot_stacked_bar_plot_metas, ncol = 1, nrow = 2, common.legend = T, legend = "right")


for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig12-1-timing-stacked-bar-plot-per-cohort.png", height = 1360, width = 960)
    print(combined_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig12-1-timing-stacked-bar-plot-per-cohort.pdf", height = 21, width = 14)
    print(combined_plot)
    dev.off()
  }
}






# ************************************************************************************************************************************************
# using binned purple data

if (dir.exists("/hpc/cuppen/")){
  df_purple_clonality_binned <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/all-purple-clonality-binned.txt", sep = "\t", stringsAsFactors = F, header = FALSE)
} else {
  df_purple_clonality_binned <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/all-purple-clonality-binned.txt", sep = "\t", stringsAsFactors = F, header = FALSE)
}

colnames(df_purple_clonality_binned) <- c("sample_id", "clonal", "probably_clonal", "probably_subclonal", "subclonal")

df_purple_clonality_binned <- merge(df_purple_clonality_binned, metadata_included[,c("sample_id", "is_metastatic", "cancer_type", "cancer_type_code")], by = "sample_id")


df_purple_clonality_binned$tmb <- rowSums(df_purple_clonality_binned[,2:5])
df_purple_clonality_binned$clonal_summed <- rowSums(df_purple_clonality_binned[,2:3])

df_purple_clonality_binned$cancer_type <- factor(df_purple_clonality_binned$cancer_type)
df_purple_clonality_binned$cancer_type_code <- factor(df_purple_clonality_binned$cancer_type_code)


# Fig13
### Clonality from PURPLE with binned data

cancer_types <- unique(metadata_included$cancer_type)



ff <- data.frame(cancer_type = character(59), cancer_type_abb = character(59), median_clona_met = numeric(59), median_clona_prim = numeric(59))

for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included$cancer_type_code[metadata_included$cancer_type == cancer_type])
  tmp_df <- df_purple_clonality_binned[df_purple_clonality_binned$cancer_type == cancer_type,]
  
  ff[i,1:2] <- c(cancer_type, cancer_type_abb)
  if (nrow(tmp_df[tmp_df$is_metastatic,]) >= 5 & nrow(tmp_df[!(tmp_df$is_metastatic),]) >= 5){
    res <- wilcox.test(tmp_df$clonal_summed[tmp_df$is_metastatic]/tmp_df$tmb[tmp_df$is_metastatic], tmp_df$clonal_summed[!(tmp_df$is_metastatic)]/tmp_df$tmb[!(tmp_df$is_metastatic)])
    if (res$p.value > 0.05){
      ff[i,2] <- NA
    }
    ff[i,3:4] <- c(median(tmp_df$clonal_summed[tmp_df$is_metastatic]/tmp_df$tmb[tmp_df$is_metastatic], na.rm = T), median(tmp_df$clonal_summed[!(tmp_df$is_metastatic)]/tmp_df$tmb[!(tmp_df$is_metastatic)], na.rm = T))
  } else {
    ff[i,3:4] <- rep(NA, times = 2)
  }
}


ff$cancer_type <- factor(ff$cancer_type)
ff$cancer_type_abb <- factor(ff$cancer_type_abb)



plloott_median <- ff[!(is.na(ff$median_clona_met)) & !(is.na(ff$median_clona_prim)),] %>% ggplot(aes(x = 100*median_clona_prim, y = 100*median_clona_met, color = cancer_type)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb), max.overlaps = 20) +
  xlim(c(40,100)) +
  xlab("primary median %") +
  ylim(c(40,100)) +
  ylab("Metastatsis median %") +
  ggtitle("Percentage of Clonal mutations (PURPLE binned)") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw()








df_purple_clonality_binned_cp <- df_purple_clonality_binned


for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  tmp_df <- df_purple_clonality_binned_cp[df_purple_clonality_binned_cp$cancer_type == cancer_type,]
  if (nrow(tmp_df[tmp_df$is_metastatic,]) >= 5 & nrow(tmp_df[!(tmp_df$is_metastatic),]) >= 5){
    
  } else {
    df_purple_clonality_binned_cp <- df_purple_clonality_binned_cp[df_purple_clonality_binned_cp$cancer_type != cancer_type,]
  }
}

unique(df_purple_clonality_binned_cp$cancer_type)

df_purple_clonality_binned_cp$cancer_type <- factor(df_purple_clonality_binned_cp$cancer_type)





plloott_box <- df_purple_clonality_binned_cp %>% ggplot(aes(x = is_metastatic, y = clonal/tmb, color = cancer_type)) + facet_wrap(~ cancer_type) +
  geom_boxplot() +
  # ggrepel::geom_text_repel(aes(label = cancer_type_abb)) +
  # xlim(c(50,100)) +
  xlab("Metastatic") +
  stat_compare_means(comparisons = list(c('FALSE', 'TRUE')), method = "wilcox.test", label = "p.signif") +
  ylab("Percentage of Clonal mutations (PURPLE binned)") +
  ylim(c(0,1.5)) +
  ggtitle("Distribution of clonality percentage per cancer type (PURPLE Pipeline binned)")




com_plot <- ggarrange(plloott_mean, plloott_median, plloott_box,  ncol=3, common.legend = TRUE, legend="bottom")

for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig13-clonality-prim-vs-metas-per-cancer-purple-binned.png", height = 920, width = 1380)
    print(com_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig13-clonality-prim-vs-metas-per-cancer-purple-binned.pdf", height = 14, width = 21)
    print(com_plot)
    dev.off()
  }
}




# fig13-1
# I made a new summary of binned PURPLE clonality for which I set hard cutoffs instead of percentile cut-offs

if (dir.exists("/hpc/cuppen/")){
  df_purple_clonality_binned2 <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/all-purple-clonality-binned2.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
} else {
  df_purple_clonality_binned2 <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/all-purple-clonality-binned2.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
}


df_purple_clonality_binned2 <- merge(df_purple_clonality_binned2, metadata_included[,c("sample_id", "is_metastatic", "cancer_type", "cancer_type_code")], by = "sample_id")


df_purple_clonality_binned2$tmb <- rowSums(df_purple_clonality_binned2[,2:5])
df_purple_clonality_binned2$clonal_summed <- rowSums(df_purple_clonality_binned2[,2:3])

df_purple_clonality_binned2$cancer_type <- factor(df_purple_clonality_binned2$cancer_type)
df_purple_clonality_binned2$cancer_type_code <- factor(df_purple_clonality_binned2$cancer_type_code)




cancer_types <- unique(metadata_included$cancer_type)



vv <- data.frame(cancer_type = character(59), cancer_type_abb = character(59), median_clona_met = numeric(59), median_clona_prim = numeric(59))

for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included$cancer_type_code[metadata_included$cancer_type == cancer_type])
  tmp_df <- df_purple_clonality_binned2[df_purple_clonality_binned2$cancer_type == cancer_type,]
  
  vv[i,1:2] <- c(cancer_type, cancer_type_abb)
  if (nrow(tmp_df[tmp_df$is_metastatic,]) >= 5 & nrow(tmp_df[!(tmp_df$is_metastatic),]) >= 5){
    res <- wilcox.test(tmp_df$clonal_summed[tmp_df$is_metastatic]/tmp_df$tmb[tmp_df$is_metastatic], tmp_df$clonal_summed[!(tmp_df$is_metastatic)]/tmp_df$tmb[!(tmp_df$is_metastatic)])
    if (res$p.value > 0.05){
      vv[i,2] <- NA
    }
    vv[i,3:4] <- c(median(tmp_df$clonal_summed[tmp_df$is_metastatic]/tmp_df$tmb[tmp_df$is_metastatic], na.rm = T), median(tmp_df$clonal_summed[!(tmp_df$is_metastatic)]/tmp_df$tmb[!(tmp_df$is_metastatic)], na.rm = T))
  } else {
    vv[i,3:4] <- rep(NA, times = 2)
  }
}



vv$cancer_type <- factor(vv$cancer_type)
vv$cancer_type_abb <- factor(vv$cancer_type_abb)



plloott_median2 <- vv[!(is.na(vv$median_clona_met)) & !(is.na(vv$median_clona_prim)),] %>% ggplot(aes(x = 100*median_clona_prim, y = 100*median_clona_met, color = cancer_type)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb), max.overlaps = 20) +
  xlim(c(40,100)) +
  xlab("primary median %") +
  ylim(c(40,100)) +
  ylab("Metastatsis median %") +
  ggtitle("Percentage of Clonal mutations (PURPLE binned)") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw()



## Apparently it's better to work with percentiles than using hard cutoffs






# fig14
#### raw clonality fig (fig2a)


suff_num_tiss <- names(table(df_purple_clonality_binned$cancer_type)[as.vector(table(df_purple_clonality_binned$cancer_type[df_purple_clonality_binned$is_metastatic])) >= 10 & as.vector(table(df_purple_clonality_binned$cancer_type[!(df_purple_clonality_binned$is_metastatic)])) >= 10])

df_purple_clonality_binned_rm <- df_purple_clonality_binned[df_purple_clonality_binned$cancer_type %in% suff_num_tiss,]

df_purple_clonality_binned_rm$cancer_type <- factor(df_purple_clonality_binned_rm$cancer_type)


df_purple_clonality_binned_rm[,2:5] <- df_purple_clonality_binned_rm[,2:5]/df_purple_clonality_binned_rm[,8]


all_order_sample <- vector()
df_purple_clonality_binned_rm$clonal_order <- rowSums(df_purple_clonality_binned_rm[,c(2,3)])

for (tissue in unique(df_purple_clonality_binned_rm$cancer_type)){
  tmp_df <- df_purple_clonality_binned_rm[df_purple_clonality_binned_rm$cancer_type == tissue,]
  tissue_order <- tmp_df[order(tmp_df$clonal_order),"sample_id"]
  
  
  all_order_sample <- append(all_order_sample, tissue_order)
}


df_purple_clonality_binned_rm$sample_id <- factor(df_purple_clonality_binned_rm$sample_id, levels = all_order_sample)

str(df_purple_clonality_binned_rm)



df_purple_clonality_binned_rm_tibble <-  gather(data = df_purple_clonality_binned_rm, key = "clonality_category", value = "clonality_value", clonal, probably_clonal, probably_subclonal, subclonal)


df_purple_clonality_binned_rm_tibble$clonality_value <- 100*df_purple_clonality_binned_rm_tibble$clonality_value
df_purple_clonality_binned_rm_tibble$clonality_category <- factor(df_purple_clonality_binned_rm_tibble$clonality_category)
df_purple_clonality_binned_rm_tibble$clonality_category <- factor(df_purple_clonality_binned_rm_tibble$clonality_category, levels = rev(levels(df_purple_clonality_binned_rm_tibble$clonality_category)))
df_purple_clonality_binned_rm_tibble$cancer_type <- factor(df_purple_clonality_binned_rm_tibble$cancer_type)


plot_stacked_bar_plot_primary <- df_purple_clonality_binned_rm_tibble[!(df_purple_clonality_binned_rm_tibble$is_metastatic),] %>% ggplot(aes(x = sample_id, y = clonality_value, fill = clonality_category)) + facet_wrap(~ cancer_type, scale = "free_x") +
  geom_bar(position="fill", stat="identity", aes(color=clonality_category)) +
  scale_color_manual(values = c("yellow", "grey", "blue", "red")) +
  scale_fill_manual(values = c("yellow", "grey", "blue", "red")) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  ggtitle("Primary Cohort")

plot_stacked_bar_plot_metastatic <- df_purple_clonality_binned_rm_tibble[df_purple_clonality_binned_rm_tibble$is_metastatic,] %>% ggplot(aes(x = sample_id, y = clonality_value, fill = clonality_category)) + facet_wrap(~ cancer_type, scale = "free_x") +
  geom_bar(position="fill", stat="identity", aes(color=clonality_category)) +
  scale_color_manual(values = c("yellow", "grey", "blue", "red")) +
  scale_fill_manual(values = c("yellow", "grey", "blue", "red")) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  ggtitle("Metastatic Cohort")



  
  


combined_plot <- ggarrange(plot_stacked_bar_plot_primary, plot_stacked_bar_plot_metastatic, ncol = 1, nrow = 2, common.legend = T, legend = "right")


for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig14-timing-stacked-bar-plot-per-cohort-binned-purple.png", height = 1360, width = 960)
    print(combined_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig14-timing-stacked-bar-plot-per-cohort-binned-purple.pdf", height = 21, width = 14)
    print(combined_plot)
    dev.off()
  }
}




# ************************************************************************************************************************************************
# Exploring clonality per cohort (metastatic and primary and separately) and studying the effect of treatment




if (dir.exists("/hpc/cuppen/")){
  treatment_metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2 - pretreatment-13092021.tsv"), sep = "\t", header = T, stringsAsFactors = F)
} else {
  treatment_metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2 - pretreatment-13092021.tsv"), sep = "\t", header = T, stringsAsFactors = F)
}

treatment_metadata$treatment_none <- !apply(treatment_metadata[,9:14], MARGIN = 1, FUN = any)
metadata_included <- merge(metadata_included, treatment_metadata[,c(3,9:15)], by = "sample_id")

if (dir.exists("/hpc/cuppen/")){
  ms_clonality_binned_purple <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-binned-purple-df.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
} else {
  ms_clonality_binned_purple <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-binned-purple-df.txt.gz", sep = "\t", stringsAsFactors = F, header = T)
}

ms_clonality_binned_purple <- merge(ms_clonality_binned_purple, metadata_included[,c(1,25:31)], by = "sample_id") 


# Fig15: clonal percentage for each cohort per cancer type


df_purple_clonality_binned_cp <- df_purple_clonality_binned


df_purple_clonality_binned_cp$cancer_type <- paste0(df_purple_clonality_binned_cp$cancer_type, " (", df_purple_clonality_binned_cp$cancer_type_code , ")")


cancer_types <- unique(metadata_included$cancer_type)

for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  tmp_df <- df_purple_clonality_binned_cp[df_purple_clonality_binned_cp$cancer_type == cancer_type,]
  if (nrow(tmp_df[tmp_df$is_metastatic,]) >= 5 & nrow(tmp_df[!(tmp_df$is_metastatic),]) >= 5){
    
  } else {
    df_purple_clonality_binned_cp <- df_purple_clonality_binned_cp[df_purple_clonality_binned_cp$cancer_type != cancer_type,]
  }
}



plott_primary <- df_purple_clonality_binned_cp[!(df_purple_clonality_binned_cp$is_metastatic),] %>% ggplot(aes (x= reorder(cancer_type_code,log2(100*(clonal+probably_clonal)/tmb),na.rm = TRUE), y = 100*(clonal+probably_clonal)/tmb)) +
  geom_boxplot(aes(color = cancer_type)) +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Cancer Type") +
  ylab("Clonality Percentage") +
  ylim(c(30,100)) +
  # ylim(c(-10,10)) +
  labs(title = "Clonality per cancer type in primary cohort") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))




plott_metastatic <- df_purple_clonality_binned_cp[df_purple_clonality_binned_cp$is_metastatic,] %>% ggplot(aes (x= reorder(cancer_type_code,log2(100*(clonal+probably_clonal)/tmb),na.rm = TRUE), y = 100*(clonal+probably_clonal)/tmb)) +
  geom_boxplot(aes(color = cancer_type)) +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Cancer Type") +
  ylab("Clonality Percentage") +
  ylim(c(30,100)) +
  # ylim(c(-10,10)) +
  labs(title = "Clonality per cancer type in metastatic cohort") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))





combined_plot <- ggarrange(plott_primary, plott_metastatic, ncol = 1, nrow = 2, common.legend = T, legend = "right")

combined_plot <- annotate_figure(combined_plot, bottom = "Binned PURPLE clonality info used!")


for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig15-clonality-per-cohort-per-cancer-type-binned-purple.png", height = 1360, width = 1360)
    print(combined_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig15-clonality-per-cohort-per-cancer-type-binned-purple.pdf", height = 21, width = 21)
    print(combined_plot)
    dev.off()
  }
}



# ************************************************************************************************************************************************
# Fig16: clonal ratio of mutational signature for different treatment groups for each cohort separately



df_sig_cont <- df


df_sig_cont <- merge(df_sig_cont, df_purple_clonality_binned[,c(1,6:8)], by = "sample_id")
df_sig_cont <- merge(df_sig_cont, treatment_metadata[,c(3,9:15)], by = "sample_id")

df_sig_cont_tibb <- gather(df_sig_cont, key = "Treatment_type", value = "Used", 7:13)
df_sig_cont_tibb <- df_sig_cont_tibb[df_sig_cont_tibb$Used,]
df_sig_cont_tibb$Treatment_type <- factor(df_sig_cont_tibb$Treatment_type)




# remove signatures that are fewer than 3 times are active in a treatment group
df_sig_cont_tibb_cp <- df_sig_cont_tibb
nrow(df_sig_cont_tibb_cp)

for (i in 1:length(unique(df_sig_cont_tibb$cosmic_sig))){
  signature <- as.character(unique(df_sig_cont_tibb$cosmic_sig))[i]
  tmp <- df_sig_cont_tibb[df_sig_cont_tibb$cosmic_sig == signature  & !is.na(df_sig_cont_tibb$fold_change),]
  for (j in 1:length(unique(df_sig_cont_tibb$Treatment_type))){
    treatment <- as.character(unique(df_sig_cont_tibb$Treatment_type))[j]
    tmp2 <- tmp[tmp$Treatment_type == treatment,]
    if (nrow(tmp2) < 3){
      df_sig_cont_tibb_cp <- df_sig_cont_tibb_cp[(df_sig_cont_tibb_cp$cosmic_sig != signature | df_sig_cont_tibb_cp$Treatment_type != treatment),]
    }
  }
}


# Making a count dataframe


a <- length(unique(df_sig_cont_tibb_cp$cosmic_sig))
b <- length(unique(df_sig_cont_tibb_cp$Treatment_type))

count_df <- data.frame(cosmic_sig = character(a*b), Treatment_type = character(a*b), count = numeric(a*b))

first <- T
for (i in 1:a){
  signature <- as.character(unique(df_sig_cont_tibb_cp$cosmic_sig)[i])
  for (j in 1:b){
    if (first){
      t <- i + j -1
      first <- F
    } else {
      t <- (i-1)*b + j
    }
    treatment <- as.character(unique(df_sig_cont_tibb_cp$Treatment_type))[j]
    tmp <- df_sig_cont_tibb_cp[df_sig_cont_tibb_cp$cosmic_sig == signature & df_sig_cont_tibb_cp$Treatment_type == treatment & df_sig_cont_tibb_cp$is_metastatic & !is.na(df_sig_cont_tibb_cp$fold_change),]
    count_df[t,1:2] <- c(signature, treatment)
    count_df[t,3] <- nrow(tmp)
    
  }
}


str(count_df)
count_df$cosmic_sig <- factor(count_df$cosmic_sig)
count_df$Treatment_type <- factor(count_df$Treatment_type)
count_df$Treatment_type <- factor(count_df$Treatment_type, levels = rev(levels(df_sig_cont_tibb_cp$Treatment_type)))


count_df$count_label <- count_df$count

count_df <- count_df[which(count_df$count >= 3),]


# summary(df_sig_cont_tibb_cp$fold_change*1.6)
# 1/0.6
df_sig_cont_tibb_cp$fold_change <- df_sig_cont_tibb_cp$fold_change*1.6





#
gg_plot_sbs <- df_sig_cont_tibb_cp[df_sig_cont_tibb_cp$is_metastatic & !is.na(df_sig_cont_tibb_cp$fold_change) & str_detect(df_sig_cont_tibb_cp$cosmic_sig, "SBS"),] %>% ggplot(aes(x = reorder(cosmic_sig,log2(fold_change),na.rm = TRUE), y = log2(fold_change))) +
  geom_boxplot(aes(fill = Treatment_type)) +
  geom_text(inherit.aes = F, data = count_df[str_detect(count_df$cosmic_sig, "SBS"),], aes(x= cosmic_sig, y = 9, fill = Treatment_type, color = Treatment_type, label = count_label), angle = 90,  position = position_dodge(width = 1), size = 2) +
  geom_text(x= 5, y = 7.5, label = "Late clonal", face = "bold") +
  geom_text(x= 5, y = -10, label = "Early clonal", face = "bold") +
  guides(color = F) +
  # scale_color_manual(guide = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  xlab("Cosmic Signature") +
  ylab("Fold Change (Log2)") +
  ylim(c(-10,10)) +
  labs(title = "Mutational Signature Timing", caption = "SBS Signatures!") +  
  theme(plot.title = element_text(size = 15, hjust = 0.5))


gg_plot_not_sbs <- df_sig_cont_tibb_cp[df_sig_cont_tibb_cp$is_metastatic & !is.na(df_sig_cont_tibb_cp$fold_change) & !str_detect(df_sig_cont_tibb_cp$cosmic_sig, "SBS"),] %>% ggplot(aes(x = reorder(cosmic_sig,log2(fold_change),na.rm = TRUE), y = log2(fold_change))) +
  geom_boxplot(aes(fill = Treatment_type)) +
  geom_text(inherit.aes = F, data = count_df[!str_detect(count_df$cosmic_sig, "SBS"),], aes(x= cosmic_sig, y = 9, fill = Treatment_type, color = Treatment_type, label = count_label), angle = 90,  position = position_dodge(width = 1), size = 2) +
  geom_text(x= 5, y = 7.5, label = "Late clonal", face = "bold") +
  geom_text(x= 5, y = -10, label = "Sub clonal", face = "bold") +
  guides(color = F) +
  # scale_color_manual(guide = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  xlab("Cosmic Signature") +
  ylab("Fold Change (Log2)") +
  ylim(c(-10,10)) +
  labs(title = "Mutational Signature Timing", 
       caption = "DBS and ID Signatures!") +
  theme(plot.title = element_text(size = 15, hjust = 0.5))


combined_plot <- ggarrange(gg_plot_sbs, gg_plot_not_sbs, ncol = 1, nrow = 2, common.legend = T, legend = "right")

unique(df_sig_cont_tibb_cp$Treatment_type)
df_sig_cont_tibb_cp[df_sig_cont_tibb_cp$cosmic_sig == "SBS7a" & df_sig_cont_tibb_cp$Treatment_type == "treatment__Immunotherapy",]
count_df[count_df$cosmic_sig == "SBS7a" & count_df$Treatment_type == "treatment_none",]

for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig16-ms-timing-per-treatment-regimen.png", height = 920, width = 2300)
    print(combined_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig16-ms-timing-per-treatment-regimen.pdf", height = 14, width = 35)
    print(combined_plot)
    dev.off()
  }
}




# Fig17

# dff here is from the binned purple clonality info

df3_sig_clonality <- dff

df3_sig_clonality <- merge(df3_sig_clonality, df_purple_clonality_binned_cp[,c(1,7:8)], by = "sample_id")
df3_sig_clonality <- merge(df3_sig_clonality, treatment_metadata[,c(3,9:15)], by = "sample_id")



df3_sig_clonality_tibb <- gather(df3_sig_clonality, key = "Treatment_type", value = "Used", 7:13)
df3_sig_clonality_tibb <- df3_sig_clonality_tibb[df3_sig_clonality_tibb$Used,]
nrow(df3_sig_clonality_tibb)
df3_sig_clonality_tibb$Treatment_type <- factor(df3_sig_clonality_tibb$Treatment_type)
df3_sig_clonality_tibb$cosmic_sig <- factor(df3_sig_clonality_tibb$cosmic_sig)


df3_sig_clonality_tibb <- df3_sig_clonality_tibb[!is.na(df3_sig_clonality_tibb$fold_change),]







df3_sig_clonality_tibb_cp <- df3_sig_clonality_tibb
nrow(df3_sig_clonality_tibb_cp)
# 
for (i in 1:length(as.character(unique(df3_sig_clonality_tibb$cosmic_sig)))){
  signature <- as.character(unique(df3_sig_clonality_tibb$cosmic_sig))[i]
  tmp <- df3_sig_clonality_tibb[df3_sig_clonality_tibb$cosmic_sig == signature  & !is.na(df3_sig_clonality_tibb$fold_change),]
  for (j in 1:length(unique(df3_sig_clonality_tibb$Treatment_type))){
    treatment <- as.character(unique(df3_sig_clonality_tibb$Treatment_type))[j]
    tmp2 <- tmp[tmp$Treatment_type == treatment & tmp$is_metastatic,]
    if (nrow(tmp2) < 3){
      print(signature)
      print(treatment)
      print("========")
      df3_sig_clonality_tibb_cp <- df3_sig_clonality_tibb_cp[(df3_sig_clonality_tibb_cp$cosmic_sig != signature | df3_sig_clonality_tibb_cp$Treatment_type != treatment),]
    }
  }
}


signature = "SBS24"
treatment = "treatment__Other"
tmp <- df3_sig_clonality_tibb_cp[df3_sig_clonality_tibb_cp$cosmic_sig == signature  & !is.na(df3_sig_clonality_tibb_cp$fold_change),]
tmp2 <- tmp[tmp$Treatment_type == treatment & tmp$is_metastatic,]
nrow(tmp2) < 3
df3_sig_clonality_tibb_cp <- df3_sig_clonality_tibb_cp[(df3_sig_clonality_tibb_cp$cosmic_sig != signature | df3_sig_clonality_tibb_cp$Treatment_type != treatment),]
nrow(df3_sig_clonality_tibb_cp)


# Making a count dataframe

a <- length(unique(df3_sig_clonality_tibb_cp$cosmic_sig))
b <- length(unique(df3_sig_clonality_tibb_cp$Treatment_type))

count_df <- data.frame(cosmic_sig = character(a*b), Treatment_type = character(a*b), count = numeric(a*b))

first <- T
for (i in 1:a){
  signature <- as.character(unique(df3_sig_clonality_tibb_cp$cosmic_sig)[i])
  for (j in 1:b){
    if (first){
      t <- i + j -1
      first <- F
    } else {
      t <- (i-1)*b + j
    }
    treatment <- as.character(unique(df3_sig_clonality_tibb_cp$Treatment_type))[j]
    tmp <- df3_sig_clonality_tibb_cp[df3_sig_clonality_tibb_cp$cosmic_sig == signature & df3_sig_clonality_tibb_cp$Treatment_type == treatment & df3_sig_clonality_tibb_cp$is_metastatic & !is.na(df3_sig_clonality_tibb_cp$fold_change),]
    count_df[t,1:2] <- c(signature, treatment)
    count_df[t,3] <- nrow(tmp)
    
  }
}


str(count_df)
count_df$cosmic_sig <- factor(count_df$cosmic_sig)
count_df$Treatment_type <- factor(count_df$Treatment_type)

count_df$count_label <- count_df$count

count_df <- count_df[which(count_df$count >= 3),]





# summary(df3_sig_clonality_tibb_cp$fold_change)
# 1/0.2
# summary(df3_sig_clonality_tibb_cp$fold_change*5.9)
df3_sig_clonality_tibb_cp$fold_change <- df3_sig_clonality_tibb_cp$fold_change*5.9

# count_df$Treatment_type <- factor(count_df$Treatment_type, levels = rev(levels(df3_sig_clonality_tibb_cp$Treatment_type)))

#  
gg_plot_sbs <- df3_sig_clonality_tibb_cp[df3_sig_clonality_tibb_cp$is_metastatic & !is.na(df3_sig_clonality_tibb_cp$fold_change) & str_detect(df3_sig_clonality_tibb_cp$cosmic_sig, "SBS"),] %>% ggplot(aes(x = reorder(cosmic_sig,log2(fold_change),na.rm = TRUE), y = log2(fold_change), fill = Treatment_type)) +
  geom_boxplot() + # position = position_dodge(preserve = "single")
  geom_text(inherit.aes = F, data = count_df[str_detect(count_df$cosmic_sig, "SBS"),], aes(x= cosmic_sig, y = 6.5, fill = Treatment_type, color = Treatment_type, label = count_label), angle = 90,  position = position_dodge(width = 0.9), size = 3) +
  geom_text(x= 5, y = 7.5, label = "Subclonal", face = "bold") +
  geom_text(x= 5, y = -10, label = "Clonal", face = "bold") +
  guides(color = F) +
  # scale_color_manual(guide = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  xlab("Cosmic Signature") +
  ylab("Fold Change (Log2)") +
  ylim(c(-10,10)) +
  labs(title = "Mutational Signature Clonality", caption = "SBS signatures!") + # 
  theme(plot.title = element_text(size = 15, hjust = 0.5))



gg_plot_not_sbs <- df3_sig_clonality_tibb_cp[df3_sig_clonality_tibb_cp$is_metastatic & !is.na(df3_sig_clonality_tibb_cp$fold_change) & !str_detect(df3_sig_clonality_tibb_cp$cosmic_sig, "SBS"),] %>% ggplot(aes(x = reorder(cosmic_sig,log2(fold_change),na.rm = TRUE), y = log2(fold_change), fill = Treatment_type)) +
  geom_boxplot() + # position = position_dodge(preserve = "single")
  geom_text(inherit.aes = F, data = count_df[!str_detect(count_df$cosmic_sig, "SBS"),], aes(x= cosmic_sig, y = 6.5, fill = Treatment_type, color = Treatment_type, label = count_label), angle = 90,  position = position_dodge(width = 0.9), size = 3) +
  geom_text(x= 5, y = 7.5, label = "Subclonal", face = "bold") +
  geom_text(x= 5, y = -10, label = "Clonal", face = "bold") +
  guides(color = F) +
  # scale_color_manual(guide = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  xlab("Cosmic Signature") +
  ylab("Fold Change (Log2)") +
  ylim(c(-10,10)) +
  labs(title = "Mutational Signature Clonality", 
       caption = "DBS and ID signatures!") +
  theme(plot.title = element_text(size = 15, hjust = 0.5))

df3_sig_clonality_tibb[df3_sig_clonality_tibb$cosmic_sig == "SBS24",]
df3_sig_clonality_tibb_cp[df3_sig_clonality_tibb_cp$cosmic_sig == "SBS24",]
count_df[count_df$cosmic_sig == "SBS24",]


combined_plot <- ggarrange(gg_plot_sbs, gg_plot_not_sbs, ncol = 1, nrow = 2, common.legend = T, legend = "right")


for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig17-ms-clonality-per-treatment-regimen.png", height = 920, width = 2300)
    print(combined_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig17-ms-clonality-per-treatment-regimen.pdf", height = 14, width = 35)
    print(combined_plot)
    dev.off()
  }
}



# ************************************************************************************************************************************************
# Fig10-1: comparing the mutation timing of mutational signature between prim and metas in a dot plot per cancer type

library(gggibbous)
library(RColorBrewer)
library(gtools)


df <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/df.all.rds"))

df2 <- df


df2$fold_change[!is.na(df2$fold_change)] <- df2$fold_change[!is.na(df2$fold_change)] * 1.67862


df2 <- merge(df2, metadata_included[,c("sample_id", "cancer_type", "cancer_type_code" ,"is_metastatic")], by = "sample_id")
# 
sigs <- unique(df2$cosmic_sig)
cancerTypes <- unique(df2$cancer_type)
cancerTypesCodes <- unique(df2$cancer_type_code)
cohort <- unique(df2$is_metastatic)
nrow <- length(sigs)*length(cancerTypes)*length(cohort)



dff <- data.frame(cosmic_sig = character(nrow), cancer_types = character(nrow), cancer_types_code = character(nrow), is_metastatic = character(nrow), median_foldChange = numeric(nrow), mean_foldChange = numeric(nrow), number_of_samples = integer(nrow))


for (i in 1:length(sigs)){
  for (j in 1:length(cancerTypes)){
    jj <- j + ((length(cohort)-1)*(j-1))
    for (k in 1:length(cohort)){
      tmp_df <- df2[df2$cosmic_sig == as.character(sigs[i]) & df2$cancer_type == cancerTypes[j] & df2$is_metastatic == cohort[k] & !is.na(df2$fold_change),]
      dff[length(cancerTypes)*length(cohort)*(i-1)  + jj+k-1,1:4] <- c(as.character(sigs[i]), cancerTypes[j], cancerTypesCodes[j], cohort[k])
      if (nrow(tmp_df) >= 30){
        dff[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:7] <- c(median(tmp_df$fold_change, na.rm = T), mean(tmp_df$fold_change, na.rm = T), nrow(tmp_df))
      } else {
        dff[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:7] <- c(NA, NA, nrow(tmp_df))
      }
    }
  }
}

dff <- dff[!is.na(dff$median_foldChange),]

dff$is_metastatic <- as.logical(dff$is_metastatic)
dff_prim <- dff[!(dff$is_metastatic),]
dff_metas <- dff[dff$is_metastatic,]

for (i in 1:length(sigs)){
  for (j in 1:length(cancerTypes)){
    tmp <- dff[dff$cosmic_sig == sigs[i] & dff$cancer_types == cancerTypes[j],]
    if (nrow(tmp) == 2){

    } else {
      dff <- dff[dff$cosmic_sig != sigs[i] | dff$cancer_types != cancerTypes[j],]
    }
  }
}

saveRDS(dff, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/dff.fig10-1.rds"))



dff <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/dff.fig10-1.rds"))

dff2 <- dff


dff2$is_metastatic <- as.logical(dff2$is_metastatic)



# doing the statistical tests

df2 <- df2[!is.na(df2$fold_change),]


sigs <- as.character(unique(df2$cosmic_sig))
cancerTypes <- unique(df2$cancer_type)
cancerTypesCodes <- unique(df2$cancer_type_code)
nrow <- length(sigs)*length(cancerTypes)

for (i in 1:length(sigs)){
  for (j in 1:length(cancerTypes)){
    tmp <- df2[df2$cosmic_sig == sigs[i] & df2$cancer_type == cancerTypes[j],]
    tmp_prim <- tmp[!(tmp$is_metastatic),]
    tmp_metas <- tmp[tmp$is_metastatic,]
    
    if (nrow(tmp_prim) >= 10 & nrow(tmp_metas) >= 10){
      
    } else {
      df2 <- df2[df2$cosmic_sig != sigs[i] | df2$cancer_type != cancerTypes[j],]
    }
  }
}


sigs <- as.character(unique(df2$cosmic_sig))
cancerTypes <- unique(df2$cancer_type)
cancerTypesCodes <- unique(df2$cancer_type_code)
nrow <- length(sigs)*length(cancerTypes)


test_df <- data.frame(cosmic_sig = character(nrow), cancer_types = character(nrow), cancer_type_code = character(nrow), p_value = numeric(nrow), effect_size = character(nrow))



for (i in 1:length(sigs)){
  ii <- i + (j-1)*(i-1)
  for (j in 1:length(cancerTypes)){
    tmp <- df2[df2$cosmic_sig == sigs[i] & df2$cancer_type == cancerTypes[j],]
    tmp_prim <- tmp[!(tmp$is_metastatic),"fold_change"]
    tmp_metas <- tmp[tmp$is_metastatic,"fold_change"]
    
    test_df[ii+j-1,1:3] <- c(sigs[i], cancerTypes[j], cancerTypesCodes[j])
    if (nrow(tmp) > 0){
      wilcox_res <- wilcox.test(tmp_prim, tmp_metas)
      
      test_df[ii+j-1,4] <- wilcox_res$p.value
      
      effsize_res <- cohen.d(tmp_prim, tmp_metas, hedges.correction=T)
      test_df[ii+j-1,5] <- abs(effsize_res$estimate)
    } else {
      test_df[ii+j-1,4:5] <- rep(NA, 2)
      
    }
  }
}

test_df$cosmic_sig <- factor(test_df$cosmic_sig, levels = mixedsort(unique(test_df$cosmic_sig)))
dff2$cosmic_sig <- factor(dff2$cosmic_sig, levels = mixedsort(unique(dff2$cosmic_sig)))


dff2 <- merge(dff2, test_df, by = c("cosmic_sig", "cancer_types"))

dff2$effect_size <- as.numeric(dff2$effect_size)
dff2$signif_p <- dff2$p_value < 0.05

dff2$signif_p <- factor(dff2$signif_p, levels = c("FALSE", "TRUE"))

dot_plot <- dff2 %>%
  ggplot(aes(x=cosmic_sig, y = cancer_types_code, size = effect_size)) + 
  theme_bw() +
  gggibbous::geom_moon( 
    aes(ratio=0.5, right=is_metastatic, fill=log2(median_foldChange), linetype = signif_p), stroke = 0.3) +
  scale_fill_distiller(name = "log2(median of late clonal/early clonal)", palette='Spectral') +
  scale_size_continuous(name="Effect size (Cohen's d)", range=c(2,12)) +
  scale_linetype((name='Significance (< 0.05)')) +
  theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, size = 10), 
        axis.text.y=element_text(vjust=0.5, size = 10),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 30, hjust = 0.5),
        plot.caption = element_text(size = 10)) +
  labs(caption = "Left half: primary \n Right half: metastatic", title = "Mutational Signature Timing \n") +
  xlab("\nCosmic signature") +
  ylab ("Cancer types")





for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig10-1-1-ms-timing-per-cancer-dotplot-10-sample-citeria.png", height = 960, width = 1440)
    print(dot_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig10-1-1-ms-timing-per-cancer-dotplot-10-sample-citeria.pdf", height = 14, width = 21)
    print(dot_plot)
    dev.off()
  }
}



# plotting the two cohorts separately by (2.5% contribution threshold data)


if (dir.exists("/hpc/cuppen/")){
  cosmic_sigs <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/cosmic-sigs.txt", sep = "\t", header = F)
} else {
  cosmic_sigs <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/cosmic-sigs.txt", sep = "\t", header = F)
}



if (dir.exists("/hpc/cuppen/")){
  ms_timing_df <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-df-2.5-percent-threhold.txt.gz", sep = "\t", header = T, stringsAsFactors = F)
} else {
  ms_timing_df <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-df-2.5-percent-threshold.txt.gz", sep = "\t", header = T, stringsAsFactors = F)
}

df <- data.frame(sample_id = character(7049*51), cosmic_sig = character(7049*51), fold_change = numeric(7049*51))


# length(metadata_included$sample_id)
for (i in 1:length(metadata_included$sample_id)){
  print(i)
  sample_id <- metadata_included$sample_id[i]
  i <- i + (50*(i-1))

  for (j in 1:length(cosmic_sigs$V1)){

    df[i+j-1, "sample_id"] <- sample_id
    df[i+j-1, "cosmic_sig"] <- cosmic_sigs$V1[j]
    late_count <-  ms_timing_df[ms_timing_df$sample_id == sample_id & ms_timing_df$timing == "clonal [late]",cosmic_sigs$V1[j]]
    early_count <- ms_timing_df[ms_timing_df$sample_id == sample_id & ms_timing_df$timing == "clonal [early]",cosmic_sigs$V1[j]]

    if (!is.na(late_count) & !is.na(early_count) & late_count != 0 & early_count != 0){
      df[i+j-1, "fold_change"] <- late_count/early_count
    } else {
      df[i+j-1, "fold_change"] <- NA
    }
  }
}
saveRDS(df, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/df.all-2.5-percent-threshold.rds"))

df <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/df.all-2.5-percent-threshold.rds"))

df2 <-df

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
      tmp_df <- df2[df2$cosmic_sig == as.character(sigs[i]) & df2$cancer_type == cancerTypes[j] & df2$cohort == cohort[k] & !is.na(df2$fold_change),]
      dff[length(cancerTypes)*length(cohort)*(i-1)  + jj+k-1,1:4] <- c(as.character(sigs[i]), cancerTypes[j], cancerTypesCodes[j], cohort[k])
      if (nrow(tmp_df) >= 10){
        dff[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:7] <- c(median(tmp_df$fold_change, na.rm = T), mean(tmp_df$fold_change, na.rm = T), nrow(tmp_df))
      } else {
        dff[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:7] <- c(NA, NA, nrow(tmp_df))
      }
    }
  }
}

dff <- dff[!is.na(dff$median_foldChange),]



dff$cohort <- factor(dff$cohort)


# saveRDS(dff, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/dff.fig10-1.rds"))

dff2 <- dff


# dff2$is_metastatic <- as.logical(dff2$is_metastatic)

dff2_prim <- dff2[dff2$cohort == "PCAWG",]
dff2_metas <- dff2[dff2$cohort == "HMF",]




str(dff2_prim)
dot_plot_metas <- dff2_metas %>%
  ggplot(aes(x=cosmic_sig, y = cancer_types_code, fill=log2(median_foldChange)), size = 5) + 
  theme_bw() +
  gggibbous::geom_moon(ratio=1, stroke = 0.3) +
  scale_fill_distiller(name = "log2(median of late clonal/early clonal)", palette='Spectral') +
  # scale_size_continuous(name="Effect size (Cohen's d)", range=c(2,12)) +
  # scale_linetype((name='Significance (< 0.05)')) +
  theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, size = 10), 
        axis.text.y=element_text(vjust=0.5, size = 10),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 30, hjust = 0.5),
        plot.caption = element_text(size = 10)) +
  labs(caption = "10 sample + 2.5% threshold", title = "Mutational Signature Timing in HMF Tumors\n") +
  xlab("\nCosmic signature") +
  ylab ("Cancer types")

for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig10-3-ms-timing-per-cancer-dotplot-10-2.5-prim.png", height = 960, width = 1440)
    print(dot_plot_prim)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig10-3-ms-timing-per-cancer-dotplot-10-2.5-prim.pdf", height = 14, width = 21)
    print(dot_plot_prim)
    dev.off()
  }
}

for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig10-3-ms-timing-per-cancer-dotplot-10-2.5-metas.png", height = 960, width = 1440)
    print(dot_plot_metas)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig10-3-ms-timing-per-cancer-dotplot-10-2.5-metas.pdf", height = 14, width = 21)
    print(dot_plot_metas)
    dev.off()
  }
}



### fi10-2 using the mean fraction instead of median of ratios







if (dir.exists("/hpc/cuppen/")){
  cosmic_sigs <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/cosmic-sigs.txt", sep = "\t", header = F)
} else {
  cosmic_sigs <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/cosmic-sigs.txt", sep = "\t", header = F)
}



if (dir.exists("/hpc/cuppen/")){
  ms_timing_df <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-df.txt.gz", sep = "\t", header = T, stringsAsFactors = F)
} else {
  ms_timing_df <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-df.txt.gz", sep = "\t", header = T, stringsAsFactors = F)
}

# ddf <- data.frame(sample_id = character(7049*51), cosmic_sig = character(7049*51), fraction = numeric(7049*51))
# 
# 
# # length(metadata_included$sample_id)
# for (i in 1:length(metadata_included$sample_id)){
#   print(i)
#   sample_id <- metadata_included$sample_id[i]
#   i <- i + (50*(i-1))
# 
#   for (j in 1:length(cosmic_sigs$V1)){
# 
#     ddf[i+j-1, "sample_id"] <- sample_id
#     ddf[i+j-1, "cosmic_sig"] <- cosmic_sigs$V1[j]
#     late_count <-  ms_timing_df[ms_timing_df$sample_id == sample_id & ms_timing_df$timing == "clonal [late]",cosmic_sigs$V1[j]]
#     early_count <- ms_timing_df[ms_timing_df$sample_id == sample_id & ms_timing_df$timing == "clonal [early]",cosmic_sigs$V1[j]]
# 
#     if (!is.na(late_count) & !is.na(early_count) & late_count != 0 & early_count != 0){
#       ddf[i+j-1, "fraction"] <- late_count/(late_count + early_count)
#     } else {
#       ddf[i+j-1, "fold_change"] <- NA
#     }
#   }
# }
# 
# 
# saveRDS(ddf, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/ddf.all.rds"))

ddf <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/ddf.all.rds"))

# summary(ddf$fraction)


ddf2 <- ddf



ddf2 <- merge(ddf2, metadata_included[,c("sample_id", "cancer_type", "cancer_type_code" ,"is_metastatic")], by = "sample_id")


# 
# sigs <- unique(ddf2$cosmic_sig)
# cancerTypes <- unique(ddf2$cancer_type)
# cancerTypesCodes <- unique(ddf2$cancer_type_code)
# cohort <- unique(ddf2$is_metastatic)
# nrow <- length(sigs)*length(cancerTypes)*length(cohort)
# 
# 
# 
# ddff <- data.frame(cosmic_sig = character(nrow), cancer_types = character(nrow), cancer_types_code = character(nrow), is_metastatic = character(nrow), median_fraction = numeric(nrow), mean_fraction = numeric(nrow), number_of_samples = integer(nrow))
# 
# 
# for (i in 1:length(sigs)){
#   for (j in 1:length(cancerTypes)){
#     jj <- j + ((length(cohort)-1)*(j-1))
#     for (k in 1:length(cohort)){
#       tmp_df <- ddf2[ddf2$cosmic_sig == as.character(sigs[i]) & ddf2$cancer_type == cancerTypes[j] & ddf2$is_metastatic == cohort[k] & !is.na(ddf2$fraction),]
#       ddff[length(cancerTypes)*length(cohort)*(i-1)  + jj+k-1,1:4] <- c(as.character(sigs[i]), cancerTypes[j], cancerTypesCodes[j], cohort[k])
#       if (nrow(tmp_df) > 2){
#         ddff[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:7] <- c(median(tmp_df$fraction, na.rm = T), mean(tmp_df$fraction, na.rm = T), nrow(tmp_df))
#       } else {
#         ddff[length(cancerTypes)*length(cohort)*(i-1) + jj+k-1,5:7] <- c(NA, NA, nrow(tmp_df))
#       }
#     }
#   }
# }
# 
# ddff <- ddff[!is.na(ddff$median_fraction),]
# 
# 
# for (i in 1:length(sigs)){
#   for (j in 1:length(cancerTypes)){
#     tmp <- ddff[ddff$cosmic_sig == sigs[i] & ddff$cancer_types == cancerTypes[j],]
#     if (nrow(tmp) == 2){
# 
#     } else {
#       ddff <- ddff[ddff$cosmic_sig != sigs[i] | ddff$cancer_types != cancerTypes[j],]
#     }
#   }
# }
# 
# saveRDS(ddff, file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/ddff.fig10-2.rds"))

ddff <- readRDS(file = paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/r-objects/ms-timing-clonality-ratios/ddff.fig10-2.rds"))


ddff2 <- ddff


ddff2$is_metastatic <- as.logical(ddff2$is_metastatic)




# doing the statistical tests

ddf2 <- ddf2[!is.na(ddf2$fraction),]

sigs <- as.character(unique(ddf2$cosmic_sig))
cancerTypes <- unique(ddf2$cancer_type)
cancerTypesCodes <- unique(ddf2$cancer_type_code)
nrow <- length(sigs)*length(cancerTypes)

for (i in 1:length(sigs)){
  for (j in 1:length(cancerTypes)){
    tmp <- ddf2[ddf2$cosmic_sig == sigs[i] & ddf2$cancer_type == cancerTypes[j],]
    tmp_prim <- tmp[!(tmp$is_metastatic),]
    tmp_metas <- tmp[tmp$is_metastatic,]
    
    if (nrow(tmp_prim) > 2 & nrow(tmp_metas) > 2){
      
    } else {
      ddf2 <- ddf2[ddf2$cosmic_sig != sigs[i] | ddf2$cancer_type != cancerTypes[j],]
    }
  }
}


sigs <- as.character(unique(ddf2$cosmic_sig))
cancerTypes <- unique(ddf2$cancer_type)
cancerTypesCodes <- unique(ddf2$cancer_type_code)
nrow <- length(sigs)*length(cancerTypes)


test_df <- data.frame(cosmic_sig = character(nrow), cancer_types = character(nrow), cancer_type_code = character(nrow), p_value = numeric(nrow), effect_size = character(nrow))



for (i in 1:length(sigs)){
  ii <- i + (j-1)*(i-1)
  for (j in 1:length(cancerTypes)){
    tmp <- ddf2[ddf2$cosmic_sig == sigs[i] & ddf2$cancer_type == cancerTypes[j],]
    tmp_prim <- tmp[!(tmp$is_metastatic),"fraction"]
    tmp_metas <- tmp[tmp$is_metastatic,"fraction"]
    
    test_df[ii+j-1,1:3] <- c(sigs[i], cancerTypes[j], cancerTypesCodes[j])
    if (nrow(tmp) > 0){
      wilcox_res <- wilcox.test(tmp_prim, tmp_metas)
      
      test_df[ii+j-1,4] <- wilcox_res$p.value
      
      effsize_res <- cohen.d(tmp_prim, tmp_metas, hedges.correction=T)
      test_df[ii+j-1,5] <- abs(effsize_res$estimate)
    } else {
      test_df[ii+j-1,4:5] <- rep(NA, 2)
      
    }
  }
}

test_df$cosmic_sig <- factor(test_df$cosmic_sig, levels = mixedsort(unique(test_df$cosmic_sig)))


ddff2$cosmic_sig <- factor(ddff2$cosmic_sig, levels = mixedsort(unique(ddff2$cosmic_sig)))


ddff2 <- merge(ddff2, test_df, by = c("cosmic_sig", "cancer_types"))

ddff2$effect_size <- as.numeric(ddff2$effect_size)
ddff2$signif_p <- ddff2$p_value < 0.05

ddff2$signif_p <- factor(ddff2$signif_p, levels = c("FALSE", "TRUE"))

dot_plot <- ddff2 %>%
  ggplot(aes(x=cosmic_sig, y = cancer_types_code, size = effect_size)) + 
  theme_bw() +
  gggibbous::geom_moon( 
    aes(ratio=0.5, right=is_metastatic, fill=median_fraction, linetype = signif_p), stroke = 0.3) +
  scale_fill_distiller(name = "Median of late clonal fraction", palette='Spectral') +
  scale_size_continuous(name="Effect size (Cohen's d)", range=c(2,12)) +
  scale_linetype((name='Significance (< 0.05)')) +
  theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, size = 10), 
        axis.text.y=element_text(vjust=0.5, size = 10),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 30, hjust = 0.5),
        plot.caption = element_text(size = 10)) +
  labs(caption = "Left half: primary \n Right half: metastatic", title = "Mutational Signature Timing \n") +
  xlab("\nCosmic signature") +
  ylab ("Cancer types")




for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig10-2-ms-timing-median-fraction-per-cancer-dotplot.png", height = 960, width = 1440)
    print(dot_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig10-2-ms-timing-median-fraction-per-cancer-dotplot.pdf", height = 14, width = 21)
    print(dot_plot)
    dev.off()
  }
}





### fig18 clonality with updated softfilter pcawg





purple_timing <- read.csv(file = paste0(wd, "r-objects/all-updated-purple-timing.txt.gz"), stringsAsFactors = F, header = T, sep = "\t")

purple_timing <- merge(purple_timing, metadata_included[,c("sample_id", "cancer_type","cancer_type_code", "cohort")], by = "sample_id")

purple_timing$tmb <- rowSums(purple_timing[,6:7])



cancer_types <- unique(purple_timing$cancer_type)
cancer_type_abb <- unique(purple_timing$cancer_type_code)

ff <- data.frame(cancer_type = character(59), cancer_type_abb = character(59), cancer_type_complete = character(59), median_clona_met = numeric(59), median_clona_prim = numeric(59))

for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included$cancer_type_code[metadata_included$cancer_type == cancer_type])
  cancer_type_complete <- paste0(cancer_type, " (", cancer_type_abb , ")")
  tmp_df <- purple_timing[purple_timing$cancer_type == cancer_type & purple_timing$subclonal_tmb_80_cutoff != 0,]
  
  ff[i,1:3] <- c(cancer_type, cancer_type_abb, cancer_type_complete)
  # print(cancer_type)
  # print(head(tmp_df))
  # print(nrow(tmp_df[tmp_df$cohort == "PCAWG",]))
  # print(nrow(tmp_df[tmp_df$cohort == "HMF",]))
  if (nrow(tmp_df[tmp_df$cohort == "PCAWG",]) >= 5 & nrow(tmp_df[tmp_df$cohort == "HMF",]) >= 5){
    tmp_df_pcawg <- tmp_df$clonal_tmb_80_cutoff[tmp_df$cohort == "PCAWG"]/tmp_df$tmb[tmp_df$cohort == "PCAWG"]
    tmp_df_hmf <- tmp_df$clonal_tmb_80_cutoff[tmp_df$cohort == "HMF"]/tmp_df$tmb[tmp_df$cohort == "HMF"]
    res <- wilcox.test(tmp_df_pcawg, tmp_df_hmf)
    # print(res$p.value)
    if (res$p.value > 0.05){
      ff[i,2] <- NA
    }
    ff[i,4:5] <- c(median(tmp_df_hmf, na.rm = T), median(tmp_df_pcawg, na.rm = T))
  } else {
    ff[i,4:5] <- rep(NA, times = 2)
  }
}




ff$cancer_type <- factor(ff$cancer_type)
ff$cancer_type_abb <- factor(ff$cancer_type_abb)



plloott_median <- ff[!(is.na(ff$median_clona_met)) & !(is.na(ff$median_clona_prim)),] %>% ggplot(aes(x = 100*median_clona_prim, y = 100*median_clona_met, color = cancer_type_complete)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb), max.overlaps = 20) +
  xlim(c(50,100)) +
  xlab("primary median %") +
  ylim(c(50,100)) +
  ylab("Metastatsis median %") +
  ggtitle("Percentage of Clonal mutations (softfilter PCAWG + strictfilter HMF)") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw() +
  labs(caption = "SUBCL >= 0.8: subclonal \n SUBCL < 0.8: clonlal")





purple_timing_cp <- purple_timing


for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  tmp_df <- purple_timing_cp[purple_timing_cp$cancer_type == cancer_type & purple_timing_cp$subclonal_tmb_80_cutoff != 0,]
  if (nrow(tmp_df[tmp_df$cohort == "PCAWG",]) >= 5 & nrow(tmp_df[tmp_df$cohort == "HMF",]) >= 5){
    
  } else {
    purple_timing_cp <- purple_timing_cp[purple_timing_cp$cancer_type != cancer_type,]
  }
}


purple_timing_cp$cancer_type <- factor(purple_timing_cp$cancer_type)




plloott_box <- purple_timing_cp[purple_timing_cp$subclonal_tmb != 0,] %>% ggplot(aes(x = cohort, y = clonal_tmb_80_cutoff/tmb, color = cancer_type_code)) + facet_wrap(~ cancer_type_code) +
  geom_boxplot() +
  # ggrepel::geom_text_repel(aes(label = cancer_type_abb)) +
  # xlim(c(50,100)) +
  xlab("Cohort") +
  stat_compare_means(comparisons = list(c('PCAWG', 'HMF')), method = "wilcox.test", label = "p.signif") +
  ylab("Percentage of Clonal mutations (HMF Pipeline)") +
  ylim(c(0,1.5)) +
  ggtitle("Distribution of clonality percentage per cancer type (HMF Pipeline)") 



com_plot <- ggarrange(plloott_mean, plloott_median, plloott_box,  ncol=3, common.legend = TRUE, legend="bottom")

for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig18-clonality-softfilter-pcawg.png", height = 920, width = 1380)
    print(com_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig18-clonality-softfilter-pcawg.pdf", height = 14, width = 21)
    print(com_plot)
    dev.off()
  }
}


# using binned data



ff <- data.frame(cancer_type = character(59), cancer_type_abb = character(59), cancer_type_complete = character(59), mean_clona_met = numeric(59), mean_clona_prim = numeric(59))

for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included$cancer_type_code[metadata_included$cancer_type == cancer_type])
  cancer_type_complete <- paste0(cancer_type, " (", cancer_type_abb , ")")
  tmp_df <- purple_timing[purple_timing$cancer_type == cancer_type & purple_timing$subclonal_tmb_80_cutoff != 0,]
  
  ff[i,1:3] <- c(cancer_type, cancer_type_abb, cancer_type_complete)
  # print(cancer_type)
  # print(head(tmp_df))
  # print(nrow(tmp_df[tmp_df$cohort == "PCAWG",]))
  # print(nrow(tmp_df[tmp_df$cohort == "HMF",]))
  if (nrow(tmp_df[tmp_df$cohort == "PCAWG",]) >= 5 & nrow(tmp_df[tmp_df$cohort == "HMF",]) >= 5){
    tmp_df_pcawg <- tmp_df$clonal[tmp_df$cohort == "PCAWG"]/tmp_df$tmb[tmp_df$cohort == "PCAWG"]
    tmp_df_hmf <- tmp_df$clonal[tmp_df$cohort == "HMF"]/tmp_df$tmb[tmp_df$cohort == "HMF"]
    res <- wilcox.test(tmp_df_pcawg, tmp_df_hmf)
    # print(res$p.value)
    if (res$p.value > 0.05){
      ff[i,2] <- NA
    }
    ff[i,4:5] <- c(mean(tmp_df_hmf, na.rm = T), mean(tmp_df_pcawg, na.rm = T))
  } else {
    ff[i,4:5] <- rep(NA, times = 2)
  }
}




ff$cancer_type <- factor(ff$cancer_type)
ff$cancer_type_abb <- factor(ff$cancer_type_abb)



plloott_mean_binned1 <- ff[!(is.na(ff$mean_clona_met)) & !(is.na(ff$mean_clona_prim)),] %>% ggplot(aes(x = 100*mean_clona_prim, y = 100*mean_clona_met, color = cancer_type_complete)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb), max.overlaps = 20) +
  xlim(c(25,100)) +
  xlab("primary mean %") +
  ylim(c(25,100)) +
  ylab("Metastatsis mean %") +
  ggtitle("Percentage of Clonal mutations (softfilter PCAWG + strictfilter HMF)") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw() +
  labs(caption = "SUBCL > 0: subclonal \n SUBCL = 0: clonlal")






purple_timing_cp <- purple_timing


for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  tmp_df <- purple_timing_cp[purple_timing_cp$cancer_type == cancer_type & purple_timing_cp$subclonal_tmb_80_cutoff != 0,]
  if (nrow(tmp_df[tmp_df$cohort == "PCAWG",]) >= 5 & nrow(tmp_df[tmp_df$cohort == "HMF",]) >= 5){
    
  } else {
    purple_timing_cp <- purple_timing_cp[purple_timing_cp$cancer_type != cancer_type,]
  }
}


purple_timing_cp$cancer_type <- factor(purple_timing_cp$cancer_type)




plloott_box_binned1 <- purple_timing_cp[purple_timing_cp$subclonal_tmb_80_cutoff != 0,] %>% ggplot(aes(x = cohort, y = clonal/tmb, color = cancer_type_code)) + facet_wrap(~ cancer_type_code) +
  geom_boxplot() +
  # ggrepel::geom_text_repel(aes(label = cancer_type_abb)) +
  # xlim(c(50,100)) +
  xlab("Cohort") +
  stat_compare_means(comparisons = list(c('PCAWG', 'HMF')), method = "wilcox.test", label = "p.signif") +
  ylab("Percentage of Clonal mutations (HMF Pipeline)") +
  ylim(c(0,1.5)) +
  ggtitle("Distribution of clonality percentage per cancer type") 



com_plot <- ggarrange(plloott_mean_binned1, plloott_median_binned1, plloott_box_binned1,  ncol=3, common.legend = TRUE, legend="bottom")


for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig18-1-clonality-softfilter-pcawg-binned.png", height = 920, width = 1380)
    print(com_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig18-1-clonality-softfilter-pcawg-binned.pdf", height = 14, width = 21)
    print(com_plot)
    dev.off()
  }
}





# using binned data2



ff <- data.frame(cancer_type = character(59), cancer_type_abb = character(59), cancer_type_complete = character(59), median_clona_met = numeric(59), median_clona_prim = numeric(59))

for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  cancer_type_abb <- unique(metadata_included$cancer_type_code[metadata_included$cancer_type == cancer_type])
  cancer_type_complete <- paste0(cancer_type, " (", cancer_type_abb , ")")
  tmp_df <- purple_timing[purple_timing$cancer_type == cancer_type & purple_timing$subclonal_tmb_80_cutoff != 0,]
  
  ff[i,1:3] <- c(cancer_type, cancer_type_abb, cancer_type_complete)
  # print(cancer_type)
  # print(head(tmp_df))
  # print(nrow(tmp_df[tmp_df$cohort == "PCAWG",]))
  # print(nrow(tmp_df[tmp_df$cohort == "HMF",]))
  if (nrow(tmp_df[tmp_df$cohort == "PCAWG",]) >= 5 & nrow(tmp_df[tmp_df$cohort == "HMF",]) >= 5){
    tmp_df_pcawg <- (tmp_df$clonal[tmp_df$cohort == "PCAWG"] + tmp_df$probably_clonal[tmp_df$cohort == "PCAWG"])/tmp_df$tmb[tmp_df$cohort == "PCAWG"]
    tmp_df_hmf <- (tmp_df$clonal[tmp_df$cohort == "HMF"] + tmp_df$probably_clonal[tmp_df$cohort == "HMF"])/tmp_df$tmb[tmp_df$cohort == "HMF"]
    res <- wilcox.test(tmp_df_pcawg, tmp_df_hmf)
    # print(res$p.value)
    if (res$p.value > 0.05){
      ff[i,2] <- NA
    }
    ff[i,4:5] <- c(median(tmp_df_hmf, na.rm = T), median(tmp_df_pcawg, na.rm = T))
  } else {
    ff[i,4:5] <- rep(NA, times = 2)
  }
}



ff$cancer_type <- factor(ff$cancer_type)
ff$cancer_type_abb <- factor(ff$cancer_type_abb)



plloott_median_binned2 <- ff[!(is.na(ff$median_clona_met)) & !(is.na(ff$median_clona_prim)),] %>% ggplot(aes(x = 100*median_clona_prim, y = 100*median_clona_met, color = cancer_type_complete)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb), max.overlaps = 20) +
  xlim(c(50,100)) +
  xlab("primary median %") +
  ylim(c(50,100)) +
  ylab("Metastatsis median %") +
  ggtitle("Percentage of Clonal mutations (softfilter PCAWG + strictfilter HMF)") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw() +
  labs(caption = "SUBCL > quantile[50% of SUBCL != 0] & : subclonal \n SUBCL = 0 + SUBCL <= quantile[50% of SUBCL != 0]: clonlal")







purple_timing_cp <- purple_timing


for (i in 1:length(cancer_types)){
  cancer_type <- cancer_types[i]
  tmp_df <- purple_timing_cp[purple_timing_cp$cancer_type == cancer_type & purple_timing_cp$subclonal_tmb_80_cutoff != 0,]
  if (nrow(tmp_df[tmp_df$cohort == "PCAWG",]) >= 5 & nrow(tmp_df[tmp_df$cohort == "HMF",]) >= 5){
    
  } else {
    purple_timing_cp <- purple_timing_cp[purple_timing_cp$cancer_type != cancer_type,]
  }
}


purple_timing_cp$cancer_type <- factor(purple_timing_cp$cancer_type)




plloott_box_binned2 <- purple_timing_cp[purple_timing_cp$subclonal_tmb_80_cutoff != 0,] %>% ggplot(aes(x = cohort, y = (clonal+probably_clonal)/tmb, color = cancer_type_code)) + facet_wrap(~ cancer_type_code) +
  geom_boxplot() +
  # ggrepel::geom_text_repel(aes(label = cancer_type_abb)) +
  # xlim(c(50,100)) +
  xlab("Cohort") +
  stat_compare_means(comparisons = list(c('PCAWG', 'HMF')), method = "wilcox.test", label = "p.signif") +
  ylab("Percentage of Clonal mutations (HMF Pipeline)") +
  ylim(c(0,1.5)) +
  ggtitle("Distribution of clonality percentage per cancer type") 



com_plot <- ggarrange(plloott_mean_binned2, plloott_median_binned2, plloott_box_binned2,  ncol=3, common.legend = TRUE, legend="bottom")






for (i in 1:2) {
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig18-1-clonality-softfilter-pcawg-binned2.png", height = 920, width = 1380)
    print(com_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig18-1-clonality-softfilter-pcawg-binned2.pdf", height = 14, width = 21)
    print(com_plot)
    dev.off()
  }
}

