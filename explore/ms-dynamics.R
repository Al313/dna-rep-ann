### Purpose

# In this script I investigate the dynamics of mutational processes using mutational signature and mutational timing information.

### Loading required libraries

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(magrittr)
library(rcompanion)
library(epitools)

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
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/png/fig4-clonality-prim-vs-metas-per-cancer-mutationTimerR.png", height = 920, width = 1380)
    print(com_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/explore/figs-wgd/pdf/fig4-clonality-prim-vs-metas-per-cancer-mutationTimerR.pdf", height = 14, width = 21)
    print(com_plot)
    dev.off()
  }
}



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



