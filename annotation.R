

# Loading packages

library(dplyr)
library(magrittr)


# General assignments and definitions

`%notin%` <- Negate(`%in%`)



if (dir.exists("/hpc/cuppen/")){
  wd <- "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/"
} else {
  wd <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/"
}


# Reading in the metadata and the diplotype info

if (dir.exists("/hpc/cuppen/")){
  metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2-metadata_08072021.tsv"), sep = "\t", header = T, stringsAsFactors = F)
} else {
  metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2-metadata_08072021.tsv"), sep = "\t", header = T, stringsAsFactors = F)
}


metadata_included <- metadata[!(metadata$is_blacklisted),]





if (dir.exists("/hpc/cuppen/")){
  
  diplotypes <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/scripts/gene_ann_ber_ner/diplotypes_germPonFilt.txt.gz",
                         sep = "\t", header = T, stringsAsFactors = F)
  
} else {
  
  diplotypes <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/scripts/gene_ann_ber_ner/diplotypes_germPonFilt.txt.gz",
                         sep = "\t", header = T, stringsAsFactors = F)
  
}




# adding the pathway info

specific_pathways <- c("BER", "GG-NER", "NER", "NER", "NER", "NER", "NER", "NER", "NER", "TC-NER", "TC-NER", "GLY-BER", "MMR", "GLY-BER", "MMR", "MMR", "GLY-BER", "GLY-BER", "GLY-BER", "GLY-BER", "GLY-BER", "MMR", "GLY-BER", "POL_PCNA", "MMR", "POL_B", "POL_D1", "POL_E", "POL_K", "GLY-BER", "GLY-BER", "GLY-BER", "NER", "GG-NER", "NER/BER")
collective_pathways <- c("BER", "NER", "NER", "NER", "NER", "NER", "NER", "NER", "NER", "NER", "NER", "BER", "MMR", "BER", "MMR", "MMR", "BER", "BER", "BER", "BER", "BER", "MMR", "BER", "POL", "MMR", "POL", "POL", "POL", "POL", "BER", "BER", "BER", "NER", "NER", "NER/BER")

genes <- unique(diplotypes$hgnc_symbol)


gene_pathway <- data.frame(hgnc_symbol = genes, pathway = collective_pathways, sub_pathway = specific_pathways)

diplotypes %<>% left_join(gene_pathway, by = "hgnc_symbol")




# Reading in the signature contribution 

if (dir.exists("/hpc/cuppen/")){
  sig_cont_hmf <- readRDS(file = paste0(wd,"r-objects/list-sig-cont.rds"))
  sbs_mut_context_mat_hmf <- readRDS(file = paste0(wd, "r-objects/sbs_matrix-sig-context.rds"))
  dbs_mut_context_mat_hmf <- readRDS(file = paste0(wd, "r-objects/dbs_matrix-sig-context.rds"))
  id_mut_context_mat_hmf <- readRDS(file = paste0(wd, "r-objects/id_matrix-sig-context.rds"))
} else {
  sig_cont_hmf <- readRDS(file = paste0(wd,"r-objects/list-sig-cont.rds"))
  sbs_mut_context_mat_hmf <- readRDS(file = paste0(wd, "r-objects/sbs_matrix-sig-context.rds"))
  dbs_mut_context_mat_hmf <- readRDS(file = paste0(wd, "r-objects/dbs_matrix-sig-context.rds"))
  id_mut_context_mat_hmf <- readRDS(file = paste0(wd, "r-objects/id_matrix-sig-context.rds"))
}


sig_cont_hmf_combined <- append(sig_cont_hmf[[1]], sig_cont_hmf[[2]])


# 15 whitelisted samples were not annotated for their genes, probably due to missing germline info (.germline.vcf.gz in D0*R folder (purple is only for somatic calling)

metadata_included$sample_id[metadata_included$sample_id %notin% unique(diplotypes$sample)]

# 383 samples with diplotype info are filtered out and not included
unique(diplotypes$sample[diplotypes$sample %notin% metadata_included$sample_id])
unique(diplotypes$sample[diplotypes$sample %notin% metadata$sample_id])

diplotypes_included <- diplotypes[diplotypes$sample  %notin% unique(diplotypes$sample[diplotypes$sample %notin% metadata_included$sample_id]),]



#

table(diplotypes$diplotype_origin)
table(diplotypes$biall_type)

length(unique(diplotypes_included[diplotypes_included$pathway == "BER" & ((diplotypes_included$diplotype_origin %notin% c("germ-som", "som-som") & (diplotypes_included$a1.max_score >= 3 | diplotypes_included$a2.max_score >= 3)) | (diplotypes_included$diplotype_origin == "som-som" & (diplotypes_included$a1.max_score >= 3 | diplotypes_included$a2.max_score >= 3)) | (diplotypes_included$diplotype_origin == "germ-som" & (diplotypes_included$a1.max_score >= 4 | diplotypes_included$a2.max_score >= 3))), "sample"]))
length(unique(diplotypes_included[diplotypes_included$pathway == "BER" & diplotypes_included$diplotype_origin %notin% c("germ-som", "som-som") & (diplotypes_included$a1.max_score >= 3 | diplotypes_included$a2.max_score >= 3),"sample"]))
length(unique(diplotypes_included[diplotypes_included$pathway == "BER" & diplotypes_included$diplotype_origin == "germ-som"  & (diplotypes_included$a1.max_score >= 3 | diplotypes_included$a2.max_score >= 4),"sample"]))
length(unique(diplotypes_included[diplotypes_included$pathway == "BER" & diplotypes_included$diplotype_origin == "som-som" & (diplotypes_included$a1.max_score >= 4 | diplotypes_included$a2.max_score >= 4),"sample"]))

length(unique(diplotypes_included[(diplotypes_included$pathway == "BER" & diplotypes_included$diplotype_origin %notin% c("germ_som", "cnv_germ") & (diplotypes_included$a1.max_score >= 4 | diplotypes_included$a2.max_score >= 4)) | (diplotypes_included$pathway == "BER" & diplotypes_included$diplotype_origin %in% c("germ_som", "cnv_germ") & (diplotypes_included$a1.max_score >= 5 | diplotypes_included$a2.max_score >= 5)), "sample"]))

length(unique(diplotypes_included[(diplotypes_included$pathway == "BER" & diplotypes_included$diplotype_origin %notin% c("germ_som", "cnv_germ") & (diplotypes_included$a1.max_score >= 4 | diplotypes_included$a2.max_score >= 4)), "sample"]))

length(unique(diplotypes_included[(diplotypes_included$pathway == "BER" & diplotypes_included$diplotype_origin %in% c("germ_som", "cnv_germ") & (diplotypes_included$a1.max_score >= 5 | diplotypes_included$a2.max_score >= 5)), "sample"]))
table(diplotypes_included[(diplotypes_included$pathway == "BER" & diplotypes_included$diplotype_origin %in% c("germ_som", "cnv_germ") & (diplotypes_included$a1.max_score >= 5 | diplotypes_included$a2.max_score >= 5)), "biall_status"])
diplotypes_included[(diplotypes_included$pathway == "BER" & diplotypes_included$diplotype_origin %in% c("germ_som", "cnv_germ") & (diplotypes_included$a1.max_score >= 5 | diplotypes_included$a2.max_score >= 5)),]



diplotypes_included[diplotypes_included$pathway == "BER" & (diplotypes_included$a1.max_score >= 3 | diplotypes_included$a2.max_score >= 3) & diplotypes_included$biall_type != "none" & diplotypes_included$diplotype_origin == "germ_som",]



length(unique(diplotypes_included[diplotypes_included$diplotype_origin == "germ_som" & diplotypes_included$a1.max_score == 5 & diplotypes_included$pathway == "BER","sample"]))
length(unique(diplotypes_included[diplotypes_included$diplotype_origin == "cnv_germ" & diplotypes_included$a2.max_score == 5 & diplotypes_included$pathway == "BER",]))


for (i in unique(diplotypes_included[diplotypes_included$pathway == "BER" & (diplotypes_included$a1.max_score >= 4 | diplotypes_included$a2.max_score >= 4), "sample"])) {
  if (i %notin% metadata_included$sample_id) {
    print(i)
  }
  if (i %notin% names(sig_cont_hmf_combined)) {
    print(i)
  }
}


ber_df <- data.frame()


for (i in unique(diplotypes_included[diplotypes_included$pathway == "BER" & (diplotypes_included$a1.max_score >= 4 | diplotypes_included$a2.max_score >= 4), "sample"])) {
  ber_df[i,names(sig_cont_hmf_combined[[i]][["SBS"]])] <- as.vector(sig_cont_hmf_combined[[i]][["SBS"]])
}

nrow(ber_df)
ber_df_norm <- 100*(ber_df/rowSums(ber_df))


ber_df_norm_fil <- ber_df_norm[ber_df_norm$SBS30 > 10 | ber_df_norm$SBS36 > 10,]
nrow(ber_df_norm_fil)



### previous method


ber_df <- data.frame()

deficient_lenient_BER_diplotypes <- diplotypes_included[diplotypes_included$pathway == "BER" & (diplotypes_included$biall_status == "deep_deletion" | diplotypes_included$diplotype_origin == "cnv_germ" & (diplotypes_included$a1.max_score==5 & diplotypes_included$a2.max_score>=4) | diplotypes_included$diplotype_origin == "cnv_som" & (diplotypes_included$a1.max_score==5 & diplotypes_included$a2.max_score>=3) | diplotypes_included$diplotype_origin %in% c('germ_som','som_som') & (diplotypes_included$a1.max_score==5 & diplotypes_included$a2.max_score==5)), ]
nrow(deficient_lenient_BER_diplotypes)


for (i in unique(deficient_lenient_BER_diplotypes$sample)) {
  ber_df[i,names(sig_cont_hmf_combined[[i]][["SBS"]])] <- as.vector(sig_cont_hmf_combined[[i]][["SBS"]])
}

nrow(ber_df)
ber_df_norm <- 100*(ber_df/rowSums(ber_df))


ber_df_norm_fil <- ber_df_norm[ber_df_norm$SBS30 > 10 | ber_df_norm$SBS36 > 10,]
nrow(ber_df_norm_fil)


#



length(unique(diplotypes_included[diplotypes_included$hgnc_symbol == "POLE" & (diplotypes_included$a1.max_score >= 3 | diplotypes_included$a2.max_score >= 3), "sample"]))




for (i in unique(diplotypes_included[diplotypes_included$hgnc_symbol == "POLE" & (diplotypes_included$a1.max_score >= 3 | diplotypes_included$a2.max_score >= 3), "sample"])) {
  if (i %notin% metadata_included$sample_id) {
    print(i)
  }
  if (i %notin% names(sig_cont_hmf_combined)) {
    print(i)
  }
}


pole_df <- data.frame()



for (i in unique(diplotypes_included[diplotypes_included$hgnc_symbol == "POLE" & (diplotypes_included$a1.max_score >= 3 | diplotypes_included$a2.max_score >= 3), "sample"])) {
  pole_df[i,names(sig_cont_hmf_combined[[i]][["SBS"]])] <- as.vector(sig_cont_hmf_combined[[i]][["SBS"]])
}

nrow(pole_df)
pole_df_norm <- 100*(pole_df/rowSums(pole_df))

pole_df_norm_fil <- pole_df_norm[pole_df_norm$SBS10a > 10 | pole_df_norm$SBS10b > 10,]

nrow(pole_df_norm_fil)




all_sbs <- data.frame()



for (i in names(sig_cont_hmf_combined)) {
  all_sbs[i,names(sig_cont_hmf_combined[[i]][["SBS"]])] <- as.vector(sig_cont_hmf_combined[[i]][["SBS"]])
}

nrow(all_sbs)


all_sbs_norm <- 100*(all_sbs/rowSums(all_sbs, na.rm = T))
all_sbs_norm[all_sbs_norm$SBS10b > 10,]



any(rownames(all_sbs) == "NA")
any(rownames(all_sbs_norm) == "NA")
str(all_sbs_norm)
any(rownames(all_sbs_norm[all_sbs_norm$SBS10a > 10,]) == "NA")




all_sbs_sorted <- all_sbs[order(all_sbs$SBS30, decreasing = T),]




diplotypes[diplotypes$sample %in% unique(diplotypes[diplotypes$hgnc_symbol == "OGG1" & (diplotypes$a1.max_score >= 4 | diplotypes$a2.max_score >= 4),"sample"]) & diplotypes$a1.max_score < 4,][diplotypes[diplotypes$sample %in% unique(diplotypes[diplotypes$hgnc_symbol == "OGG1" & (diplotypes$a1.max_score >= 4 | diplotypes$a2.max_score >= 4),"sample"]) & diplotypes$a1.max_score < 4,"a2.max_score"] >= 4,]


diplotypes[diplotypes$pathway == "BER" & diplotypes$sample == "CPCT02030442T",]
