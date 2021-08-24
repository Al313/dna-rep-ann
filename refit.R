
# Loading libraries

library(MutationalPatterns)


# Loading mutSigExtractor package


if (dir.exists("/hpc/cuppen/")){
  
  base_dir <- "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/"
  devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mutSigExtractor/'))
  
} else {
  
  library(mutSigExtractor)
}



`%notin%` <- Negate(`%in%`)



if (dir.exists("/hpc/cuppen/")){
  wd <- "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/"
} else {
  wd <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/"
}


if (dir.exists("/hpc/cuppen/")){
  metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2-metadata_08072021.tsv"), sep = "\t", header = T, stringsAsFactors = F)
} else {
  metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2-metadata_08072021.tsv"), sep = "\t", header = T, stringsAsFactors = F)
}


metadata_included <- metadata[!(metadata$is_blacklisted),]



if (dir.exists("/hpc/cuppen/")){
  pcawg_meta <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/metadata-per-donor.Rds")
} else {
  pcawg_meta <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/metadata-per-donor.Rds")
}



if (dir.exists("/hpc/cuppen/")){
  hmf_meta <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj2/HMF/metadata.tsv", header = T,
                       sep = "\t", stringsAsFactors = F)
} else {
  hmf_meta <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj2/HMF/metadata.tsv", header = T,
                       sep = "\t", stringsAsFactors = F)
}




# all(metadata_included[metadata_included$cohort == "HMF", "sample_id"] %in% hmf_meta$sampleId)
# all(metadata_included[metadata_included$cohort == "PCAWG", "sample_id"] %in% pcawg_meta$icgc_donor_id)





sbs_signatures <- SBS_SIGNATURE_PROFILES_V3
dbs_signatures <- DBS_SIGNATURE_PROFILES
id_signatures <- INDEL_SIGNATURE_PROFILES

sig_cont_hmf <- list()
sig_cont_hmf[["HMF"]] <- list()
sig_cont_hmf[["PCAWG"]] <-list()


sbs_mut_context_mat_hmf <- matrix(nrow = nrow(metadata_included), ncol = nrow(sbs_signatures))
rownames(sbs_mut_context_mat_hmf) <- metadata_included$sample_id
colnames(sbs_mut_context_mat_hmf) <- rownames(sbs_signatures)

dbs_mut_context_mat_hmf <- matrix(nrow = nrow(metadata_included), ncol = nrow(dbs_signatures))
rownames(dbs_mut_context_mat_hmf) <- metadata_included$sample_id
colnames(dbs_mut_context_mat_hmf) <- rownames(dbs_signatures)

id_mut_context_mat_hmf <- matrix(nrow = nrow(metadata_included), ncol = nrow(id_signatures))
rownames(id_mut_context_mat_hmf) <- metadata_included$sample_id
colnames(id_mut_context_mat_hmf) <- rownames(id_signatures)



for (i in 1:nrow(metadata_included)){
  print(i)
  print(metadata_included$sample_id[i])
  
  if (metadata_included$cohort[i] == "HMF"){
    HMF <- T
    file_exists <- T
    sample_id <- metadata_included$sample_id[i]
  
  
  
    if (dir.exists("/hpc/cuppen/")){
      path <- paste0("/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/",
                                          hmf_meta$setName[hmf_meta$sampleId == sample_id], "/purple/", sample_id, ".purple.somatic.vcf.gz")
    } else {
  
  
    path <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/",
                   hmf_meta$setName[hmf_meta$sampleId == sample_id], "/purple/", sample_id, ".purple.somatic.vcf.gz")
    }
    
  } else if (metadata_included$cohort[i] == "PCAWG") {
    
    HMF <- F
    sample_id <- metadata_included$sample_id[i]
    
    if (dir.exists("/hpc/cuppen/")){
          path <- paste0("/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
                                              sample_id, "-from-jar/purple53/", sample_id, "T.purple.somatic.vcf.gz")
          if (file.exists(paste0("/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
                                 sample_id, "-from-jar/purple53/", sample_id, "T.purple.somatic.vcf.gz"))){
            file_exists <- T
          } else {
            file_exists <- F
          }
        } else {

        path <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
                       sample_id, "-from-jar/purple53/", sample_id, "T.purple.somatic.vcf.gz")
        if (file.exists(paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
                               sample_id, "-from-jar/purple53/", sample_id, "T.purple.somatic.vcf.gz"))){
          file_exists <- T
        } else {
          file_exists <- F
        }
        }
  }
  

  if (HMF | file_exists){
  
    snv_mut_matrix <- extractSigsSnv(path, output='contexts', vcf.filter='PASS', sample.name=sample_id)
    sbs_mut_context_mat_hmf[sample_id,] <- as.numeric(snv_mut_matrix)
  
  
    snv_mut_vector <- as.vector(snv_mut_matrix)
    names(snv_mut_vector) <- rownames(snv_mut_matrix)
    snv_sig_cont <- fitToSignatures(snv_mut_vector, sbs_signatures)
  
    dbs_mut_matrix <- extractSigsDbs(path, output='contexts', vcf.filter='PASS', sample.name=sample_id)
    dbs_mut_context_mat_hmf[sample_id,] <- as.numeric(dbs_mut_matrix)
  
    dbs_mut_vector <- as.vector(dbs_mut_matrix)
    names(dbs_mut_vector) <- rownames(dbs_mut_matrix)
    dbs_sig_cont <- fitToSignatures(dbs_mut_vector, dbs_signatures)
  
  
    # We will use the hmf set of ID mutation types which has 83 elements
    indel_mut_matrix <- extractSigsIndel(path, output='contexts', vcf.filter='PASS', method = "PCAWG", sample.name=sample_id)
    id_mut_context_mat_hmf[sample_id,] <- as.numeric(indel_mut_matrix)
  
    indel_mut_vector <- as.vector(indel_mut_matrix)
    names(indel_mut_vector) <- rownames(indel_mut_matrix)
    indel_sig_cont <- fitToSignatures(indel_mut_vector, id_signatures)
  
    if (HMF){
      sig_cont_hmf[[1]][[sample_id]] <- list(SBS = snv_sig_cont, DBS = dbs_sig_cont, ID = indel_sig_cont)
    } else {
      sig_cont_hmf[[2]][[sample_id]] <- list(SBS = snv_sig_cont, DBS = dbs_sig_cont, ID = indel_sig_cont)
    }
  } else {
    
    print(paste0("VCF file for ", sample_id, " does not exist!"))
    sbs_mut_context_mat_hmf[sample_id,] <- NA
    dbs_mut_context_mat_hmf[sample_id,] <- NA
    id_mut_context_mat_hmf[sample_id,] <- NA
    sig_cont_hmf[[2]][[sample_id]] <- NA
    
    }

}


if (dir.exists("/hpc/cuppen/")){
  saveRDS(sig_cont_hmf, file = paste0(wd,"r-objects/list-sig-cont.rds"))
  saveRDS(sbs_mut_context_mat_hmf, file = paste0(wd, "r-objects/sbs_matrix-sig-context.rds"))
  saveRDS(dbs_mut_context_mat_hmf, file = paste0(wd, "r-objects/dbs_matrix-sig-context.rds"))
  saveRDS(id_mut_context_mat_hmf, file = paste0(wd, "r-objects/id_matrix-sig-context.rds"))
} else {
  saveRDS(sig_cont_hmf, file = paste0(wd,"r-objects/list-sig-cont.rds"))
  saveRDS(sbs_mut_context_mat_hmf, file = paste0(wd, "r-objects/sbs_matrix-sig-context.rds"))
  saveRDS(dbs_mut_context_mat_hmf, file = paste0(wd, "r-objects/dbs_matrix-sig-context.rds"))
  saveRDS(id_mut_context_mat_hmf, file = paste0(wd, "r-objects/id_matrix-sig-context.rds"))
}



