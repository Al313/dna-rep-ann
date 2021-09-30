#' ########################################################### LOADING PACKAGES ###########################################################
#' 
#' library(pacman)
#' pacman::p_load(GenomicRanges, VariantAnnotation, reshape2, dplyr, BSgenome, BSgenome.Hsapiens.UCSC.hg19, devtools, nlme, foreach, RCurl,
#'                jsonlite, MutationalPatterns, seqminer, TxDb.Hsapiens.UCSC.hg19.knownGene, vcfR)
#' 
#' 
#' ########################################################### ASSIGNMENTS ####################################################################
#' 
#' # Assigning values
#' 
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' projectId <- "3c7d37894d2fa48ddc91ab162e23d6e3eda476e0e7f7a4d35ba261871ebf046b"
#' token <- "WrR6im6Xgf6aCHZw6GMPvGlKa4MyOUswHVnCTdfXYIgVGZbSbLzHnJPkDrjcKvs9UhFgOfihvF8aGs9nWkZ"
#' args <- commandArgs(trailingOnly=TRUE)
#' local <- "/home/ali313/Documents/studies/master/umc-project"
#' 
#' 
#' i <- 1
#' 
#' i <- as.numeric(args[1])
#' 
#' 
#' if (dir.exists("/hpc/cuppen/")){
#'   manifest <- read.csv(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/manifest_HMF_PCAWG.gene_ann.txt.gz", sep = "\t",
#'                        header = T, stringsAsFactors = F)
#' } else {
#'   manifest <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/manifest_HMF_PCAWG.gene_ann.txt.gz", sep = "\t",
#'                        header = T, stringsAsFactors = F)
#' }
#' 
#' sampleID <- manifest$sample[i]
#' print(sampleID)
#' 
#' if (dir.exists("/hpc/cuppen/")){
#'   path_to_vcf <- paste0(manifest[i, "dir"], manifest[i, "som_vcf"])
#' } else {
#'   path_to_vcf <- paste0(local, manifest[i, "dir"], manifest[i, "som_vcf"])
#' }
#' 
#' 
#' 
#' tadBedFileInput <- args[2]
#' binBedFileInput <- args[3]
#' repTimingFileInput <- args[4]
#' repTimingFileInputBoxtel <- args[5]
#' cpgInputFile <- args[6]
#' repOrientationInputFile <- args[7]
#' 
#' 
#' ########################################################### FUNCTIONS ########################################################################
#' 
#' # 1. Loading variants from database ------------------------------------------------------------------------------------------------------------------
#' # Connecting to and retrieving data from sparqling-genomics
#' 
#' SGSPARQL <- function (projectId, token, query, verbose = T){
#'   if (verbose){
#'     message("Connecting to server ...")
#'   }
#'   accumulator <- basicTextGatherer()
#'   accumulator$reset()
#'   
#'   curlPerform(url           = paste0("https://sg-n1.op.umcutrecht.nl/api/query/?project-id=", projectId),
#'               httpheader    = c("Accept"        = "application/json",
#'                                 "Authorization" = "Basic ZGxhYnVzZXI6SEIzdlMkV3VAeXV0YndSTg==",
#'                                 "Cookie"        = paste0("SGSession=", token),
#'                                 "Content-Type"  = "application/sparql-update"),
#'               customrequest = "POST",
#'               postfields    = query,
#'               writefunction = accumulator$update)
#'   
#'   jsonData    <- accumulator$value()
#'   jsonData
#'   data        <- fromJSON(jsonData)
#'   
#'   accumulator$reset()
#'   return (data)
#' }
#' 
#' SPARQL_query_SNVs_sampleID <- function(mutType, sampleID,projectId,token, verbose = T){
#'   # Check mut type
#'   if (verbose){
#'     message("Loading variants ...")
#'   }
#'   if (mutType != "SNV")
#'     stop("this is a query to retrieve SNVs")
#'   query <- sprintf(
#'     "
#'              PREFIX faldo:       <http://biohackathon.org/resource/faldo#>
#'              PREFIX sg:          <https://sparqling-genomics.org/0.99.12/>
#'              PREFIX vcf2rdf:     <sg://0.99.12/vcf2rdf/>
#'              PREFIX variant:     <sg://0.99.12/vcf2rdf/variant/>
#'              PREFIX info:        <sg://0.99.12/vcf2rdf/info/>
#'              PREFIX seq:         <sg://0.99.12/vcf2rdf/sequence/>
#'              
#'              SELECT DISTINCT STRAFTER(STR(?chromosome), 'nuccore/') AS ?chromosome
#'                              ?position
#'                              STRAFTER(STR(?REF), STR(seq:)) AS ?REF
#'                              STRAFTER(STR(?ALT), STR(seq:)) AS ?ALT
#'                              ?TNC
#'                              ?PURPLE_AF
#'                              ?KT
#'                              ?MH
#'                              ?REPS
#'                              ?REPC
#'              
#'              WHERE {
#'                GRAPH <hartwig://purple-somatic-vcf> {
#'               
#'                  ?origin  sg:filename             ?filename .
#'              
#'                  ?sample  sg:foundIn              ?origin ;
#'                           rdfs:label              ?sampleName .
#'              
#'                  ?variant sg:originatedFrom       ?origin ;
#'                           faldo:reference         ?chromosome ;
#'                           faldo:position          ?position ;
#'                           variant:REF             ?REF ;
#'                           variant:ALT             ?ALT ;
#'                           variant:FILTER          vcf2rdf:PASS ;
#'                           info:TNC                ?TNC ;
#'                           info:PURPLE_AF          ?PURPLE_AF .
#'                  OPTIONAL {
#'                   GRAPH <hartwig://purple-somatic-vcf> {
#'                     ?variant info:KT       ?opt_KT .
#'                     BIND(COALESCE(?opt_KT, '-') AS ?KT)
#'                   }
#'                  }
#'                  
#'                  OPTIONAL {
#'                   GRAPH <hartwig://purple-somatic-vcf> {
#'                     ?variant info:MH       ?opt_MH .
#'                     BIND(COALESCE(?opt_MH, '-') AS ?MH)
#'                   }
#'                 }
#'                 
#'                 OPTIONAL {
#'                   GRAPH <hartwig://purple-somatic-vcf> {
#'                     ?variant info:REP_S       ?opt_REPS .
#'                     BIND(COALESCE(?opt_REPS, '-') AS ?REPS)
#'                   }
#'                 }
#'                 
#'                  OPTIONAL {
#'                   GRAPH <hartwig://purple-somatic-vcf> {
#'                     ?variant info:REP_C       ?opt_REPC .
#'                     BIND(COALESCE(?opt_REPC, '-') AS ?REPC)
#'                   }
#'                  }
#'                 
#'                }
#'              
#'                FILTER (?sampleName = '%s'^^xsd:string)
#'                FILTER (STRLEN(STRAFTER(STR(?ALT), STR(seq:))) = 1)
#'                FILTER (STRLEN(STRAFTER(STR(?REF), STR(seq:))) = 1)
#'              }
#'   ",sampleID)
#'   results <- SGSPARQL (projectId, token, query)
#'   results <- unique(results) #remove duplicated mutations because of GT type in blood sample
#'   return(results)
#' }
#' 
#' 
#' # 2. Mutational signatures and their likelihoods ----------------------------------------------------------------------------------------------
#' 
#' fit_to_signatures_goldenratio = function(mut_matrix, signatures, type, cutoff, method = "golden-ratio-search", ...){
#'   # Check mutation type argument
#'   if (missing(type)) { type_default = TRUE }
#'   else { type_default = FALSE }
#'   type = check_mutation_type(type)
#'   
#'   # If signature object is a matrix, then look at "mut_matrix" for mutation type
#'   if (class(signatures) == "matrix")
#'   { 
#'     if (class(mut_matrix) == "matrix") 
#'     { 
#'       signatures = list("snv"=signatures) 
#'       mut_matrix = list("snv"=mut_matrix)
#'     }
#'     else 
#'     {
#'       signatures_list = list()
#'       for (m in names(mut_matrix))
#'       {
#'         if (all(rownames(mut_matrix[[m]]) %in% rownames(signatures))) 
#'         { 
#'           signatures_list[[m]] = signatures
#'           signatures = signatures_list
#'           break 
#'         }
#'       }
#'     }
#'     
#'     # If count matrix object is a matrix, then look at "signatures" for mutation type
#'   } else if (class(signatures) == "list")
#'   {
#'     if (class(mut_matrix) == "matrix") 
#'     {
#'       mut_list = list()
#'       for (m in names(signatures))
#'       {
#'         if (all(rownames(signatures[[m]]) %in% rownames(mut_matrix))) 
#'         { 
#'           mut_list[[m]] = mut_matrix
#'           mut_matrix = mut_list
#'           break 
#'         }
#'       }
#'       
#'       type = names(mut_matrix)
#'     }
#'   }
#'   
#'   if (class(mut_matrix) != "list")
#'     stop(paste("No list is given for 'mut_matrix' and mutation type",
#'                "could not be found in signature list"))
#'   
#'   # Get the mutation types asked for
#'   if (!type_default & !(all(type %in% names(mut_matrix))))
#'     stop("One or more mutation types are not found in count matrices")
#'   mut_matrix = mut_matrix[intersect(type, names(mut_matrix))]
#'   signatures = signatures[intersect(type, names(signatures))]
#'   
#'   # Extra arguments for whichSignatures()
#'   dots = list(...)
#'   
#'   # Solve the least squares error problem
#'   if (method == "least-squares")
#'   {
#'     if (missing(cutoff)) { cutoff = 0 }
#'     res = least_squares_error_fitting(mut_matrix, signatures, cutoff)
#'     
#'     # Solve the golden ratio search problem
#'   } else if (method == "golden-ratio-search")
#'   {
#'     # If signature.cutoff is not given, but cutoff is, then use 
#'     # signature.cutoff = cutoff in whichSignatures()
#'     if (!("signature.cutoff" %in% names(dots)) & !missing(cutoff)) 
#'     { res = golden_ratio_search_fitting(mut_matrix, signatures, signature.cutoff = cutoff, ...) }
#'     else 
#'     { res = golden_ratio_search_fitting(mut_matrix, signatures, ...) }
#'   }
#'   
#'   return(res)
#' }
#' 
#' check_mutation_type <- function(type){
#'   # Default type is "snv"
#'   if (missing(type)) {type = c("snv")}
#'   
#'   # Translate type to lower case
#'   type = unname(sapply(type, function(m) tolower(m)))
#'   
#'   if (any(type == "all")) { type = c("snv","dbs","indel") }
#'   else if (!all(type %in% c("snv","dbs","indel"))) 
#'   {
#'     stop("One or more mutation types given are unknown")
#'   }
#'   
#'   return(type)
#' }
#' 
#' golden_ratio_search_fitting <- function(mut_matrix, signatures, ...){
#'   # Check if mut_matrix and signatures are both lists
#'   if (class(mut_matrix) != "list" | isEmpty(names(mut_matrix)) | any(names(mut_matrix) == ""))
#'     stop("'mut_matrix' is not a list or some elements are not named")
#'   if (class(signatures) != "list" | isEmpty(names(signatures)) | any(names(signatures) == ""))
#'     stop("'signatures' is not a list or some elements are not named")
#'   
#'   mut_matrix_transposed = list()
#'   signatures_transposed = list()
#'   contribution = list()
#'   reconstructed = list()
#'   unknown = list()
#'   
#'   # For each mutation type in mut_matrix
#'   for (m in names(mut_matrix))
#'   {
#'     # Transpose the mutation matrix and signature matrix
#'     mut_matrix_transposed[[m]] = as.data.frame(t(mut_matrix[[m]]))
#'     signatures_transposed[[m]] = as.data.frame(t(signatures[[m]]))
#'     
#'     # For each sample get results from the golden ratio search
#'     result = sapply(rownames(mut_matrix_transposed[[m]]), function(n)
#'       whichSignatures(mut_matrix_transposed[[m]], sample.id = n, 
#'                       signatures_transposed[[m]], 
#'                       contexts.needed = TRUE, ...))
#'     
#'     # Write results of contribution, reconstructed and unknown as lists 
#'     # of mutation types
#'     contribution[[m]] = as.matrix(do.call(cbind, 
#'                                           lapply(1:ncol(result), function(i)
#'                                             data.frame(unlist(result[1,i])))))
#'     colnames(contribution[[m]]) = colnames(result)
#'     
#'     reconstructed[[m]] = as.matrix(do.call(cbind, 
#'                                            lapply(1:ncol(result), function(i) 
#'                                              data.frame(unlist(result[3,i])))))
#'     colnames(reconstructed[[m]]) = colnames(result)
#'     rownames(reconstructed[[m]]) = colnames(result[3][[1]])
#'     
#'     unknown[[m]] = as.matrix(do.call(cbind, lapply(1:ncol(result), function(i)
#'       data.frame(unlist(result[5,i])))))
#'     colnames(unknown[[m]]) = colnames(result)
#'   }
#'   
#'   if (length(contribution) == 1)
#'   {
#'     contribution = contribution[[1]]
#'     reconstructed = reconstructed[[1]]
#'     unknown = unknown[[1]]
#'   }
#'   
#'   res = list("contribution"=contribution, 
#'              "reconstructed"=reconstructed, 
#'              "unknown"=unknown)
#'   
#'   return(res)
#' }
#' 
#' tri_context_generator <-function(){
#'   tri_context <- vector()
#'   for (i in c("C>A","C>G","C>T","T>A","T>C","T>G")){
#'     first <- c(rep("A", 4), rep("C", 4), rep("G", 4), rep("T", 4))
#'     last <- rep(c("A", "C", "G","T"), 2)
#'     mut <- sprintf("[%s]", i)
#'     tri_context <- c(tri_context, paste0(first, mut, last))
#'   }
#'   return(tri_context)
#' }
#' 
#' mutation_likelihood_score <- function(mutation_df=df_SNVs, signatures=COSMIC_signatures,sampleID,
#'                                       colname_chrom,colname_position,colname_REF,colname_ALT,colname_TNC, verbose = T){
#'   options(warn = -1)
#'   if (verbose){
#'     message("Getting mutation signatures and their likelihood scores ...")
#'   }
#'   
#'   #create unique ID for each mut
#'   mutation_df <- mutation_df %>% tidyr::unite(mutID,c(colname_chrom,colname_position,colname_REF,colname_ALT), remove = FALSE) 
#'   
#'   df_input <- NULL
#'   df_input <- mutation_df %>% dplyr::select(colname_chrom,
#'                                             colname_position,
#'                                             colname_REF,
#'                                             colname_ALT,
#'                                             colname_TNC)
#'   df_input <- df_input %>% tidyr::unite(mutID,c(colname_chrom,colname_position,colname_REF,colname_ALT), remove = FALSE)
#'   colnames(df_input) <- c("mutID","chrom","position","ref","alt","trinucleotidecontext")
#'   
#'   df_input$ref_rev <- df_input$ref
#'   df_input$alt_rev <- df_input$alt
#'   df_input$trinucleotidecontext_rev <- df_input$trinucleotidecontext
#'   ref_A <- which(df_input$ref_rev == "A")
#'   ref_G <- which(df_input$ref_rev == "G")
#'   df_input$ref_rev[c(ref_A,ref_G)] <- as.character(reverseComplement(DNAStringSet(df_input$ref_rev[c(ref_A,ref_G)])))
#'   df_input$alt_rev[c(ref_A,ref_G)] <- as.character(reverseComplement(DNAStringSet(df_input$alt_rev[c(ref_A,ref_G)])))
#'   df_input$trinucleotidecontext_rev[c(ref_A,ref_G)] <- as.character(reverseComplement(DNAStringSet(df_input$trinucleotidecontext_rev[c(ref_A,ref_G)])))
#'   
#'   
#'   five_ref <- substr(df_input$trinucleotidecontext_rev, 1, 1)
#'   three_ref <- substr(df_input$trinucleotidecontext_rev, 3, 3)
#'   mutations <- sprintf("%s[%s>%s]%s", five_ref, df_input$ref_rev, df_input$alt_rev, three_ref)
#'   mut_mat <- matrix(table(factor(mutations, levels = tri_context_generator())))
#'   rownames(mut_mat) <- tri_context_generator()
#'   colnames(mut_mat) <- sampleID
#'   #plot_96_profile(mut_mat)
#'   
#'   rowsum <- colSums(signatures)
#'   relative_signatures <- signatures / rowsum
#'   lsq_contribution <- NULL
#'   #plot_96_profile(relative_signatures)
#'   
#'   
#'   ###option1###
#'   ###n_samples = dim(mut_mat)[2]
#'   ###n_signatures = dim(signatures)[2]
#'   ###lsq_contribution = matrix(NA, nrow=n_signatures, ncol=n_samples)
#'   ###lsq_reconstructed = matrix(NA, nrow=96, ncol=n_samples)
#'   ###for (i in 1:ncol(mut_mat)){
#'   ###  y = mut_mat[,i]
#'   ###  lsq = lsqnonneg(as.matrix(signatures), as.numeric(y))
#'   ###  lsq_contribution[,i] = lsq$x
#'   ###  lsq_reconstructed[,i] = as.matrix(signatures) %*% as.matrix(lsq$x)
#'   ###}
#'   
#'   ###option2###
#'   #' Find the linear combination of mutation signatures that most closely
#'   #' reconstructs the mutation matrix by using the golden ratio search
#'   #' algorithm implemented in the \link{deconstructSigs} package (Rosenthal et al. 2016).
#'   ###goldon <- fit_to_signatures_goldenratio(mut_matrix=list("snv" = as.matrix(mut_mat)),signatures= list("snv"=signatures) ,type="snv",method = "golden-ratio-search")
#'   ###goldon$unknown
#'   ###goldon <- goldon$contribution
#'   ###lsq_contribution <- goldon*nrow(df_input)
#'   
#'   ###option3###
#'   ##MutPatters
#'   strict_refit <- fit_to_signatures_strict(mut_mat, as.matrix(signatures), max_delta = 0.004)
#'   lsq_contribution <- strict_refit$fit_res$contribution
#'   
#'   
#'   
#'   
#'   # Add row and col names
#'   sample_names = colnames(mut_mat)
#'   signature_names = colnames(signatures)
#'   # mut_type_names = rownames(contribution_profile_192)
#'   #signature_names = sprintf("HMF_%s",LETTERS[1:16])
#'   mut_type_names = rownames(signatures)
#'   
#'   colnames(lsq_contribution) = sample_names
#'   rownames(lsq_contribution) = signature_names
#'   
#'   
#'   res = list(lsq_contribution)
#'   names(res) = c("contribution")
#'   
#'   absolute_contribution_sample_pattern <- foreach(i = 1:nrow(res$contribution), .combine = 'cbind') %do% {
#'     x <- relative_signatures[,i] * res$contribution[i,]
#'     x <- as.numeric(x)
#'     return(x)
#'   }
#'   
#'   relative_contribution_sample_pattern <- absolute_contribution_sample_pattern / rowSums(absolute_contribution_sample_pattern)
#'   colnames(relative_contribution_sample_pattern) <- signature_names
#'   rownames(relative_contribution_sample_pattern) <- rownames(mut_mat)
#'   #Check
#'   #rowSums(relative_contribution_sample_pattern)
#'   df_input_final <- foreach(x = rownames(relative_contribution_sample_pattern), .combine = "rbind") %do% {
#'     if (x %in% mutations){
#'       rows <- which(mutations == x)
#'       return(cbind(df_input[rows,], t(relative_contribution_sample_pattern[x,])))
#'     } else {
#'       return(NULL)
#'     }
#'   }
#'   
#'   df_input_final <- df_input_final %>% 
#'     mutate(position = as.character(position)) %>% 
#'     mutate(across(is.numeric, ~trimws(format(round(.,2), nsmall=2)))) %>%
#'     mutate(position = as.integer(position)) %>% 
#'     arrange(chrom,position)
#'   
#'   
#'   df_input_final$mut_sign <- colnames(df_input_final[,grepl("SBS",names(df_input_final))])[max.col(df_input_final[,grepl("SBS",names(df_input_final))],ties.method="first")]
#'   df_input_final$mut_sign_score <- apply(df_input_final[,grepl("SBS",names(df_input_final))], 1, max)
#'   
#'   mut_sign_contribution <- as.data.frame(t(table(df_input_final$mut_sign)))[,2:3]
#'   colnames(mut_sign_contribution) <- c("SBS","mut_sign")
#'   #  print("mutational contribution")
#'   #  print(mut_sign_contribution)
#'   #lsq_contribution %>% as.data.frame()%>% tibble::rownames_to_column(var = "SBS") %>% left_join(.,mut_sign_df) %>% as.matrix() %>% t()%>%as.data.frame()
#'   #write.table(df_input_final, sprintf("/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/analysis/Bastiaan/signatures_to_genes/df_inputs/breast_colon_all_muts/%s_contribution.txt", sampleid),row.names = F)
#'   df_out <- left_join(mutation_df, dplyr::select(df_input_final,mutID,mut_sign,mut_sign_score), by = "mutID")
#'   df_out$mutID <- NULL
#'   return(df_out)
#'   
#' }
#' 
#' get_known_signatures <- function (muttype = c("snv", "dbs", "indel", "tsb_snv"),
#'                                   source = c("COSMIC", "SIGNAL", "SPARSE"),
#'                                   sig_type = c("reference", "exposure", "tissue"),
#'                                   incl_poss_artifacts = FALSE,
#'                                   tissue_type = c(
#'                                     NA, "Biliary", "Bladder", "Bone",
#'                                     "Breast", "Cervix", "CNS",
#'                                     "Colorectal", "Esophagus", "Head",
#'                                     "Kidney", "Liver", "Lung",
#'                                     "Lymphoid", "Myeloid", "Ovary",
#'                                     "Pancreas", "Prostate", "Skin",
#'                                     "Stomach", "Thyroid", "Uterus"
#'                                   )) {
#'   # Validate arguments
#'   muttype <- match.arg(muttype)
#'   source <- match.arg(source)
#'   sig_type <- match.arg(sig_type)
#'   tissue_type <- match.arg(tissue_type)
#'   
#'   if (!.is_na(tissue_type) & sig_type != "tissue") {
#'     stop("tissue_type can only be used with `sig_type == 'tissue'`",
#'          call. = FALSE
#'     )
#'   }
#'   
#'   # Determine signature file name
#'   basename_sig <- paste0(muttype, "_", source, "_", sig_type, ".txt")
#'   fname_sig <- file.path("extdata", "signatures", basename_sig)
#'   fname_sig <- system.file(fname_sig, package = "MutationalPatterns")
#'   
#'   # Give error if file doesn't exist.
#'   if (!file.exists(fname_sig)) {
#'     stop(paste0(
#'       "The signature file: ", fname_sig, " does not exist.\n",
#'       "Look at the documentation of 'get_known_signatures()' for",
#'       " all the possible combinations of arguments."
#'     ),
#'     call. = FALSE
#'     )
#'   }
#'   
#'   # Read in signature file
#'   signatures <- read.table(fname_sig, sep = "\t", header = TRUE)
#'   
#'   
#'   # Remove meta columns
#'   if (muttype == "snv") {
#'     meta_cols <- c(1, 2)
#'   } else if (muttype == "tsb_snv") {
#'     meta_cols <- c(1, 2, 3)
#'   } else {
#'     meta_cols <- 1
#'   }
#'   signatures <- as.matrix(signatures[, -meta_cols, drop = FALSE])
#'   
#'   # Remove possible artifacts
#'   if (!incl_poss_artifacts) {
#'     if (source == "SIGNAL" & sig_type == "reference") {
#'       good_cols <- grep("Ref.Sig.N[0-9]{0-2}",
#'                         colnames(signatures),
#'                         invert = TRUE
#'       )
#'       signatures <- signatures[, good_cols, drop = FALSE]
#'     }
#'     
#'     if (source == "COSMIC" & muttype == "snv") {
#'       bad_sigs <- paste0("SBS", c(27, 43, seq(45, 60)))
#'       good_cols <- !colnames(signatures) %in% bad_sigs
#'       signatures <- signatures[, good_cols, drop = FALSE]
#'     }
#'   }
#'   
#'   # Select signatures of the specified tissue type
#'   if (!.is_na(tissue_type)) {
#'     tissue_cols <- grep(paste0("^", tissue_type, "_"), colnames(signatures))
#'     signatures <- signatures[, tissue_cols, drop = FALSE]
#'   }
#'   
#'   return(signatures)
#' }
#' 
#' .is_na <- function(x) {
#'   purrr::is_scalar_vector(x) && is.na(x)
#' }
#' 
#' 
#' # 3. Chromatin status --------------------------------------------------------------------------------------------------------------------------
#' 
#' # getTadChromatinStat <- function(vcf, bed, verbose = T){
#' #   if (verbose){
#' #     message("Getting chromatin status (from tads) ...")
#' #   }
#' # #  if (verbose) {message("Process starting ...")}
#' # #  vcf <- read.csv(file = vcfFile, sep = "\t", stringsAsFactors = F)
#' #   tmp_vcf <- vcf[,1:4]
#' #   coords <- paste("chr", tmp_vcf$chromosome, ":", tmp_vcf$position, "-", tmp_vcf$position, sep = "")
#' #   tabixOut <- c()
#' #   for (i in coords){
#' #     container <- seqminer::tabix.read(tabixFile = bed, i)
#' #     if (length(container) == 0){
#' #       tabixOut <- append(tabixOut, "NA\tNA\tNA\tNA")
#' #     } else{
#' #       tabixOut <- append(tabixOut, container)
#' #     }
#' #   }
#' #   tabixOut <- lapply(strsplit(tabixOut,'\t'), as.vector)
#' #   tabixOut <- as.data.frame(do.call(rbind,tabixOut), stringsAsFactors=F)
#' #   colnames(tabixOut) <- c("CHROM", "POS1", "POS2", "CHROMSTAT")
#' #   tmp_vcf <- cbind(tmp_vcf, tabixOut$CHROMSTAT)
#' #   colnames(tmp_vcf)[5] <- "CHROMATIN_STAT_TAD"
#' #   vcf <- cbind(vcf, tmp_vcf[5])
#' #   return(vcf)
#' # }
#' 
#' 
#' # getBinChromatinStat <- function(vcf, bed, verbose = T){
#' #   if (verbose){
#' #     message("Getting chromatin status (from bins) ...")
#' #   }
#' # #  if (verbose) {message("Process starting ...")}
#' # #  vcf <- read.csv(file = vcfFile, sep = "\t", stringsAsFactors = F)
#' #   tmp_vcf <- vcf[,1:4]
#' #   coords <- paste("chr", tmp_vcf$chromosome, ":", tmp_vcf$position, "-", tmp_vcf$position, sep = "")
#' #   tabixOut <- c()
#' #   for (i in coords){
#' #     container <- seqminer::tabix.read(tabixFile = bed, i)
#' #     if (length(container) == 0){
#' #       tabixOut <- append(tabixOut, "NA\tNA\tNA\tNA")
#' #     } else{
#' #       tabixOut <- append(tabixOut, container)
#' #     }
#' #   }
#' #   tabixOut <- lapply(strsplit(tabixOut,'\t'), as.vector)
#' #   tabixOut <- as.data.frame(do.call(rbind,tabixOut), stringsAsFactors=F)
#' #   colnames(tabixOut) <- c("CHROM", "POS1", "POS2", "CHROMSTAT")
#' #   tmp_vcf <- cbind(tmp_vcf, tabixOut$CHROMSTAT)
#' #   colnames(tmp_vcf)[5] <- "CHROMATIN_STAT_BIN"
#' #   vcf <- cbind(vcf, tmp_vcf[5])
#' #   return(vcf)
#' # }
#' 
#' 
#' 
#' # 3.1 Chromatin status optimized (without for loops) ----------------------------------------------------------------------------------------
#' 
#' getTadChromatinStatOpt <- function(vcf, bed, verbose = T){
#'   if (verbose){
#'     message("Getting chromatin status (from tads) ...")
#'     print("Getting chromatin status (from tads) ...")
#'   }
#'   #  if (verbose) {message("Process starting ...")}
#'   #  vcf <- read.csv(file = vcfFile, sep = "\t", stringsAsFactors = F)
#'   coords <- paste("chr", vcf$chromosome, ":", vcf$position, "-", vcf$position, sep = "")
#'   vcf$CHROMATIN_STAT_TAD <- NA
#'   tabixOut <- seqminer::tabix.read(tabixFile = bed, coords)
#'   tabixOut <- lapply(strsplit(tabixOut,'\t'), as.vector)
#'   tabixOut <- as.data.frame(do.call(rbind,tabixOut), stringsAsFactors=F)
#'   colnames(tabixOut) <- c("chrom", "chromStart", "chromEnd", "CHROMATIN_STAT_TAD")
#'   GenomeInfoDb::seqlevelsStyle(tabixOut$chrom)<- 'NCBI'
#'   tabixOut[,2] <- as.numeric(tabixOut[,2])
#'   tabixOut[,3] <- as.numeric(tabixOut[,3])
#'   indeces <- isOverlappingChromPos(vcf$chromosome, vcf$position, vcf$position, 
#'                                    tabixOut$chrom, tabixOut$chromStart, tabixOut$chromEnd)
#'   indeces[sapply(indeces, is.integer0)] <- NA
#'   names(indeces) <- paste0(vcf$chromosome, ":", vcf$position, "_", vcf$REF, "/", vcf$ALT)
#'   vcf[names(indeces[!is.na(indeces)]),11] <- tabixOut[,4]
#'   return(vcf)
#' }
#' 
#' getBinChromatinStatOpt <- function(vcf, bed, verbose = T){
#'   if (verbose){
#'     message("Getting chromatin status (from bins) ...")
#'     print("Getting chromatin status (from bins) ...")
#'   }
#'   #  if (verbose) {message("Process starting ...")}
#'   #  vcf <- read.csv(file = vcfFile, sep = "\t", stringsAsFactors = F)
#'   coords <- paste("chr", vcf$chromosome, ":", vcf$position, "-", vcf$position, sep = "")
#'   vcf$CHROMATIN_STAT_BIN <- NA
#'   tabixOut <- seqminer::tabix.read(tabixFile = bed, coords)
#'   tabixOut <- lapply(strsplit(tabixOut,'\t'), as.vector)
#'   tabixOut <- as.data.frame(do.call(rbind,tabixOut), stringsAsFactors=F)
#'   colnames(tabixOut) <- c("chrom", "chromStart", "chromEnd", "CHROMATIN_STAT_BIN")
#'   GenomeInfoDb::seqlevelsStyle(tabixOut$chrom)<- 'NCBI'
#'   # For some reason if you convert the below columns into numeric you will get error for some hmf vcf files!
#'   tabixOut[,2] <- as.integer(tabixOut[,2])
#'   tabixOut[,3] <- as.integer(tabixOut[,3])
#'   indeces <- isOverlappingChromPos(vcf$chromosome, vcf$position, vcf$position, 
#'                                    tabixOut$chrom, tabixOut$chromStart, tabixOut$chromEnd)
#'   indeces[sapply(indeces, is.integer0)] <- NA
#'   names(indeces) <- paste0(vcf$chromosome, ":", vcf$position, "_", vcf$REF, "/", vcf$ALT)
#'   vcf[names(indeces[!is.na(indeces)]),11] <- tabixOut[,4]
#'   return(vcf)
#' }
#' 
#' 
#' # 4. Replication timing ---------------------------------------------------------------------------------------------------------------------
#' 
#' # getRepliSeq <- function(vcf, bed, verbose = T){
#' #   if (verbose){
#' #     message("Getting replication timing ...")
#' #   }
#' #   #if (verbose) {message("Process starting ...")}
#' #   #vcf <- read.csv(file = vcfFile, sep = "\t", stringsAsFactors = F)
#' #   vcf <- vcf[,1:4]
#' #   coords <- paste("chr", vcf$chromosome, ":", vcf$position, "-", vcf$position, sep = "")
#' #   tabixOut <- c()
#' #   for (i in coords){
#' #     container <- seqminer::tabix.read(tabixFile = bed, i)
#' #     if (length(container) == 0){
#' #       tabixOut <- append(tabixOut, "NA\tNA\tNA\tNA\tNA")
#' #     } else{
#' #       tabixOut <- append(tabixOut, container)
#' #     }
#' #   }
#' #   tabixOut <- lapply(strsplit(tabixOut,'\t'), as.vector)
#' #   tabixOut <- as.data.frame(do.call(rbind,tabixOut), stringsAsFactors=F)
#' #   colnames(tabixOut) <- c("CHROM", "POS1", "POS2", "rep_timing_mean", "rep_timing_bin")
#' #   #vcf <- read.csv(file = vcfFile, sep = "\t", stringsAsFactors = F)
#' #   vcf <- cbind(vcf, tabixOut$rep_timing_bin)
#' #   colnames(vcf)[5] <- "rep_timing_bin"
#' #   #vcf$REPLISEQ <- as.numeric(vcf$REPLISEQ)
#' #   return(vcf)
#' # }
#' 
#' # 4.1 Replication timing optimized without for loops ------------------------------------------------------------------------------------------
#' 
#' 
#' 
#' getRepliSeqOpt <- function(vcf, bed, verbose = T){
#'   if (verbose){
#'     message("Getting replication timing ...")
#'     print("Getting replication timing ...")
#'   }
#'   #if (verbose) {message("Process starting ...")}
#'   #vcf <- read.csv(file = vcfFile, sep = "\t", stringsAsFactors = F)
#'   coords <- paste("chr", vcf$chromosome, ":", vcf$position, "-", vcf$position, sep = "")
#'   vcf$rep_timing_mean <- NA
#'   vcf$rep_timing_bin <- NA
#'   tabixOut <- seqminer::tabix.read(tabixFile = bed, coords)
#'   tabixOut <- lapply(strsplit(tabixOut,'\t'), as.vector)
#'   tabixOut <- as.data.frame(do.call(rbind,tabixOut), stringsAsFactors=F)
#'   colnames(tabixOut) <- c("chrom", "chromStart", "chromEnd", "rep_timing_mean", "rep_timing_bin")
#'   #vcf <- read.csv(file = vcfFile, sep = "\t", stringsAsFactors = F)
#'   GenomeInfoDb::seqlevelsStyle(tabixOut$chrom)<- 'NCBI'
#'   tabixOut[,2] <- as.numeric(tabixOut[,2])
#'   tabixOut[,3] <- as.numeric(tabixOut[,3])
#'   indeces <- isOverlappingChromPos(vcf$chromosome, vcf$position, vcf$position, 
#'                                    tabixOut$chrom, tabixOut$chromStart, tabixOut$chromEnd)
#'   indeces[sapply(indeces, is.integer0)] <- NA
#'   #vcf$REPLISEQ <- as.numeric(vcf$REPLISEQ)
#'   names(indeces) <- paste0(vcf$chromosome, ":", vcf$position, "_", vcf$REF, "/", vcf$ALT)
#'   vcf[names(indeces[!is.na(indeces)]),11:12] <- tabixOut[,4:5]
#'   return(vcf)
#' }
#' 
#' 
#' # 5. CpG islands ---------------------------------------------------------------------------------------------------------------------
#' 
#' # getCpgAnn <- function (vcf, bed, verbose = T){
#' #   if (verbose){
#' #     message("Getting CpG island status ...")
#' #   }
#' #   vcf$cpg <- NA
#' #   vcf$percpg <- NA
#' #   vcf$pergc <- NA
#' #   coords <- paste("chr", vcf$chromosome, ":", vcf$position, "-", vcf$position, sep = "")
#' #   tabixOut <- seqminer::tabix.read(tabixFile = bed, coords)
#' #   tabixOut <- lapply(strsplit(tabixOut,'\t'), as.vector)
#' #   tabixOut <- as.data.frame(do.call(rbind,tabixOut), stringsAsFactors=F)
#' #   tabixOut <- unique(tabixOut)
#' #   colnames(tabixOut) <- c("chrom", "chromStart", "chromEnd", "name", "length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp")
#' #   tabixOut$chromStart <- as.integer(tabixOut$chromStart)
#' #   tabixOut$chromEnd <- as.integer(tabixOut$chromEnd)
#' #   for (i in 1:nrow(vcf)){
#' #     for (j in 1:nrow(tabixOut)){
#' #       if (paste0("chr",vcf$chromosome[i]) == tabixOut$chrom[j]){
#' #         if (vcf$position[i] >= tabixOut$chromStart[j] & vcf$position[i] <= tabixOut$chromEnd[j]){
#' #           vcf[i,"cpg"] <- T
#' #           vcf[i,"percpg"] <- tabixOut$perCpg[j]
#' #           vcf[i,"pergc"] <- tabixOut$perGc[j]
#' #         }
#' #       }
#' #     }
#' #     if (is.na(vcf$cpg[i])){
#' #       vcf$cpg[i] <- F
#' #     }
#' #   }
#' #   vcf$percpg <- as.numeric(vcf$percpg)
#' #   vcf$pergc<- as.numeric(vcf$pergc)
#' #   return(vcf)
#' # }
#' 
#' 
#' 
#' # 5.1 CpG islands (without for loop) ------------------------------
#' 
#' isOverlapping <- function(start1, end1, start2, end2, verbose=F){
#'   # start1=genome_bins$start_linear
#'   # end1=genome_bins$end_linear
#'   #
#'   # start2=linearizeChromPos(cnv$chromosome, cnv$start)
#'   # end2=linearizeChromPos(cnv$chromosome, cnv$end)
#'   
#'   if( length(start1)!=length(end1) ){
#'     stop('start1 and end1 must be the same length')
#'   }
#'   if( length(start2)!=length(end2) ){
#'     stop('start2 and end2 must be the same length')
#'   }
#'   
#'   if(verbose){ pb <- txtProgressBar(max=length(start1), style=3) }
#'   lapply(1:length(start1), function(i){
#'     #i=1
#'     # out <- which(pmax(start1[i],start2) <= pmin(end1[i], end2))
#'     # if(length(out)!=0){ out } else { NA }
#'     if(verbose){ setTxtProgressBar(pb, i) }
#'     which(pmax(start1[i],start2) <= pmin(end1[i], end2))
#'   })
#' }
#' 
#' isOverlappingChromPos <- function(
#'   chrom1=NULL, start1=NULL, end1=NULL, chrom2=NULL, start2=NULL, end2=NULL,
#'   df1=NULL, df2=NULL
#' ){
#'   ## Debug
#'   # chrom1=bed$chrom
#'   # start1=bed$start
#'   # end1=bed$end
#'   #
#'   # chrom2=cnv$chrom
#'   # start2=cnv$start
#'   # end2=cnv$end
#'   #
#'   # df1 <- data.frame(chrom=chrom1, start=start1, end=end1, stringsAsFactors=F)
#'   # df2 <- data.frame(chrom=chrom2, start=start2, end=end2, stringsAsFactors=F)
#'   # colnames(df1)[1:3] <- colnames(df2)[1:3] <- c('chrom','start','end')
#'   
#'   if(!is.null(df1)){
#'     chrom1 <- df1[,1]
#'     start1 <- df1[,2]
#'     end1 <- df1[,3]
#'   }
#'   
#'   if(!is.null(df2)){
#'     chrom2 <- df2[,1]
#'     start2 <- df2[,2]
#'     end2 <- df2[,3]
#'   }
#'   
#'   GenomeInfoDb::seqlevelsStyle(chrom1)<- 'NCBI'
#'   GenomeInfoDb::seqlevelsStyle(chrom2)<- 'NCBI'
#'   
#'   chrom_lookup <- as.factor(unique(c(chrom1,chrom2)))
#'   
#'   chrom1 <- factor(chrom1, levels=chrom_lookup)
#'   chrom2 <- factor(chrom2, levels=chrom_lookup)
#'   
#'   pad_width <- nchar(max(c(start1, end1, start2, end2)))
#'   
#'   chromPosAsNumeric <- function(chrom, pos){
#'     as.numeric(paste0(
#'       as.integer(chrom),
#'       formatC(pos, width=pad_width, format='d', flag='0')
#'     ))
#'   }
#'   
#'   start1 <- chromPosAsNumeric(chrom1, start1)
#'   end1 <- chromPosAsNumeric(chrom1, end1)
#'   
#'   start2 <- chromPosAsNumeric(chrom2, start2)
#'   end2 <- chromPosAsNumeric(chrom2, end2)
#'   
#'   isOverlapping(start1, end1, start2, end2)
#' }
#' 
#' is.integer0 <- function(x)
#' {
#'   is.integer(x) && length(x) == 0L
#' }
#' 
#' getCpgAnnOpt <- function (vcf, bed, verbose = T){
#'   if (verbose){
#'     message("Getting CpG island status ...")
#'     print("Getting CpG island status ...")
#'   }
#'   vcf$cpg <- FALSE
#'   vcf$percpg <- NA
#'   vcf$pergc <- NA
#'   coords <- paste("chr", vcf$chromosome, ":", vcf$position, "-", vcf$position, sep = "")
#'   tabixOut <- seqminer::tabix.read(tabixFile = cpgInputFile, coords)
#'   if (length(tabixOut) != 0){   # Without this if statement we would get error for samples with low calls that wouldreturn nothing and the length of the tabiixOut is 0
#'   tabixOut <- lapply(strsplit(tabixOut,'\t'), as.vector)
#'   tabixOut <- as.data.frame(do.call(rbind,tabixOut), stringsAsFactors=F)
#'   colnames(tabixOut) <- c("chrom", "chromStart", "chromEnd", "name", "length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp")
#'   tabixOut$chromStart <- as.integer(tabixOut$chromStart)
#'   tabixOut$chromEnd <- as.integer(tabixOut$chromEnd)
#'   #tabixOut <- unique(tabixOut)
#'   GenomeInfoDb::seqlevelsStyle(tabixOut$chrom)<- 'NCBI'
#'   
#'   indeces <- isOverlappingChromPos(vcf$chromosome, vcf$position, vcf$position, 
#'                                    tabixOut$chrom, tabixOut$chromStart, tabixOut$chromEnd)
#'   indeces[sapply(indeces, is.integer0)] <- NA
#'   names(indeces) <- paste0(vcf$chromosome, ":", vcf$position, "_", vcf$REF, "/", vcf$ALT)
#'   vcf[names(indeces[!is.na(indeces)]),12:13] <- tabixOut[,8:9]
#'   vcf[,12] <- as.numeric(vcf[,12])
#'   vcf[,13] <- as.numeric(vcf[,13])
#'   vcf[names(indeces[!is.na(indeces)]),11] <- TRUE
#'   }
#'   return(vcf)
#' }
#' 
#' 
#' # 6. Trp orientation ---------------------------------------------------------------------------------------------------------------------
#' 
#' getTrpStrandAnn <- function(vcf_gr, verbose = T){
#'   if (verbose){
#'     message("Getting transcription strand orientation ...")
#'     print("Getting transcription strand orientation ...")
#'   }
#'   # # Reading in the data
#'   # vcf <- readVcf(file = path_to_vcf, genome = "hg19")
#'   # vcf_gr <- rowRanges(vcf)
#'   # 
#'   # # This is how you subset GRanges based on its metadata
#'   # vcf_gr <- vcf_gr[elementMetadata(vcf_gr)[,"FILTER"] == "PASS"]
#'   # 
#'   # Getting the gene body and transcription information
#'   suppressMessages(genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
#'   genes_hg19 <- keepStandardChromosomes(genes_hg19, pruning.mode = "coarse")
#'   genes_hg19 <- dropSeqlevels(genes_hg19, "chrM", pruning.mode = "coarse") 
#'   # Checking if the seqlevels match
#'   if (suppressWarnings(all(seqlevels(genes_hg19) == seqlevels(vcf_gr))) == FALSE){
#'     # If not matching then drop the "MT" chromosome level
#'     newStyle <- mapSeqlevels(seqlevels(vcf_gr), "UCSC") 
#'     vcf_gr <- renameSeqlevels(vcf_gr, newStyle)
#'     vcf_gr <- dropSeqlevels(vcf_gr, "chrM", pruning.mode = "coarse") 
#'   }
#'   
#'   # Annotating the variants with trp strand info
#'   trp_str_ann <- mut_strand(vcf_gr, genes_hg19, mode = "transcription")
#'   
#'   # Adding the trp str annotation to the GRanges object
#'   values(vcf_gr) <- cbind(values(vcf_gr), DataFrame(trp_str_ann))
#'   
#'   # This is because mut_strand function returns the strand info for "C" and "T" bases!
#'   
#'   # index_for_transcribed <- which((vcf_gr$REF == "A" | vcf_gr$REF == "G") & vcf_gr$trp_str_ann == "transcribed")
#'   # index_for_untranscribed <- which((vcf_gr$REF == "A" | vcf_gr$REF == "G") & vcf_gr$trp_str_ann == "untranscribed")
#'   # 
#'   # 
#'   # vcf_gr[index_for_transcribed,]$trp_str_ann <- "untranscribed"
#'   # vcf_gr[index_for_untranscribed,]$trp_str_ann <- "transcribed"
#'   
#'   # Returning the annotated GRnages
#'   return(vcf_gr)
#' }
#' 
#' 
#' # 7. Rep orientation ----------------------------------------------------------------------------------------------------------------------
#' 
#' getRepliStrandAnn <- function(vcf_gr, bed, verbose = T){
#'   if (verbose){
#'     message("Getting replication strand orientation ...")
#'     print("Getting replication strand orientation ...")
#'   }
#'   # Here I read in the actual vcf file and do not get the data from sparqling-gnomics dataset
#'   # vcf <- readVcf(file = vcfFile, genome = "hg19")
#'   # vcf_gr <- rowRanges(vcf)
#'   # This is how you subset GRanges based on its metadata
#'   # vcf_gr <- vcf_gr[elementMetadata(vcf_gr)[,"FILTER"] == "PASS"]
#'   # Making the seqlevels compatible
#'   newStyle <- mapSeqlevels(seqlevels(vcf_gr), "UCSC") 
#'   vcf_gr <- renameSeqlevels(vcf_gr, newStyle)
#'   # To get the granges with the same length
#'   vcf_gr_cp <- vcf_gr
#'   dropped_ranges <- vcf_gr_cp[seqnames(vcf_gr_cp) %in% c("chrY", "chrX")]
#'   if (length(dropped_ranges) != 0){
#'   dropped_ranges$rep_str_ann <- "-"}
#'   # We have to remove "chrX" and "chrY" as there is no information available for them in the replication direction file
#'   vcf_gr <- dropSeqlevels(vcf_gr, c("chrM", "chrX", "chrY"), pruning.mode = "coarse")
#'   # getting the replication strand file
#'   repli_strand_granges <- readRDS(file = bed)
#'   # Running mut_sig
#'   rep_str_ann <- mut_strand(vcf_gr, repli_strand_granges, mode = "replication")
#'   #  adding the annotation to the GRnages object
#'   values(vcf_gr) <- cbind(values(vcf_gr), DataFrame(rep_str_ann))
#'   # To get the granges with the same length
#'   grl <- GRangesList(vcf_gr, dropped_ranges)
#'   vcf_gr <- unlist(grl)
#'   # returning the annotated vcf
#'   return(vcf_gr)
#' }
#' 
#' 
#' # 8. Functional genomic units ---------------------------------------------------------------------------------------------------------------
#' 
#' getFuncUnitAnn <- function (vcf_gr, verbose = T) {
#'   if (verbose){
#'     message("Getting function units ...")
#'     print("Getting function units ...")
#'   }
#'   # Loading the txdb information
#'   txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'   if (length(seqlevels(txdb)) >= 25){
#'     # Removing extra chromosomes
#'     keepStandardChromosomes(txdb)
#'     newStyle <- mapSeqlevels(seqlevels(txdb), "NCBI")
#'     # Using "NCBI" seqnames and not "UCSC" ("chr1, ...)
#'     txdb <- renameSeqlevels(txdb, newStyle)
#'     # Dropping "MT" level
#'     dropSeqlevels(txdb, "MT")
#'   }
#'   # Getting the variants in collapsedVcf object
#'   #vcf <- readVcf(vcfFile, "hg19")
#'   # getting the GRnages object
#'   #vcf_gr <- rowRanges(vcf)
#'   # Remove variants that haven't passed PON filter
#'   #vcf_gr <- vcf_gr[elementMetadata(vcf_gr)[,"FILTER"] == "PASS"]
#'   # Dropping "MT" level
#'   # vcf_gr <- dropSeqlevels(vcf_gr, "MT", pruning.mode = "coarse")
#'   
#'   # Are seqlevels of "query" and "subject" compatible?
#'   intersect(seqlevels(vcf_gr), seqlevels(txdb))
#'   
#'   # Annotating with regard to functions
#'   loc_all <- suppressMessages(locateVariants(vcf_gr, txdb, AllVariants(), ignore.strand = T))
#'   
#'   # Merging ranges
#'   loc_all <- unique(loc_all)
#'   # Some variants are located in two different regions and with the following lines I remove them altogether
#'   duplicates <- unique(names(loc_all)[!isUnique(names(loc_all))])
#'   loc_all <- DataFrame(loc_all)
#'   
#'   loc_all_final <- loc_all[(!rownames(loc_all) %in% duplicates),]
#'   
#'   # That's the dumb way to do it which would tae up to 6 mins for 10,000 calls
#'   # vcf_gr$genomic_func <- NA
#'   # for (i in 1:length(vcf_gr)){
#'   #   if (names(vcf_gr)[i] %in% rownames(loc_all_final)){
#'   #     vcf_gr[i]$genomic_func <- as.character(loc_all_final[names(vcf_gr)[i], "LOCATION"])
#'   #   } else {
#'   #     vcf_gr[i]$genomic_func <- "NA"
#'   #   }
#'   # }
#'   
#'   # The smart way of doing that
#'   df_split <- split(vcf_gr, names(vcf_gr) %in% rownames(loc_all_final))
#'   df_split$'TRUE'$genomic_func <- loc_all_final$LOCATION
#'   if (length(df_split) != 1){  # For samples with very low number of variants the list might only have one item
#'   df_split$'FALSE'$genomic_func <- "NA"}
#'   annotated_vcf_gr <- unlist(df_split)
#'   names(annotated_vcf_gr) <- paste0(seqnames(annotated_vcf_gr), ":", start(ranges(annotated_vcf_gr)), "_",
#'                                     annotated_vcf_gr$ref, "/", annotated_vcf_gr$alt)
#'   annotated_vcf_gr <- annotated_vcf_gr[names(vcf_gr),]
#'   
#'   # return(vcf_gr)
#'   return(annotated_vcf_gr)
#' }
#' 
#' ########################################################### ANNOTATING #####################################################################
#' 
#' # 1. Loading variant calls from database ------------------------------------------------------------------------------------------------------
#' # TNC | PURPLE_AF | KT | MH | REPS | REPC
#' 
#' # Getting the data from the sparqling-genomics
#' # df_SNVs <- SPARQL_query_SNVs_sampleID("SNV",sampleID,projectId,token)
#' 
#' # Getting th data from the vcf files
#' 
#' #
#' # path_to_vcf <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/150720_HMFregCPCT_HMFxx5_HMFxx6_CPCT02020171/purple/CPCT02020171T.purple.somatic.vcf.gz"
#' # sampleID <- "CPCT02020171T"
#' 
#' vcf <- read.vcfR(path_to_vcf)
#' df_SNVs <- cbind(as.data.frame(getFIX(vcf)), INFO2df(vcf))
#' df_SNVs <- df_SNVs[df_SNVs$FILTER == "PASS",]
#' df_SNVs<- df_SNVs[,c("CHROM", "POS", "REF", "ALT", "TNC", "PURPLE_AF", "KT", "MH", "REP_C", "REP_S")]
#' df_SNVs$POS <- as.integer(df_SNVs$POS)
#' chrOrder <-c((1:22),"X","Y")
#' df_SNVs <- df_SNVs[order(factor(df_SNVs$CHROM, levels = chrOrder, ordered = T), df_SNVs$POS),]
#' 
#' rownames(df_SNVs) <- paste0(df_SNVs$CHROM, ":", df_SNVs$POS, "_", df_SNVs$REF, "/", df_SNVs$ALT)
#' df_SNVs <- df_SNVs[isUnique(paste0(df_SNVs$CHROM,":", df_SNVs$POS)),]
#' colnames(df_SNVs)[c(1,2)] <- c("chromosome", "position")
#' 
#' # Remove indels!
#' df_SNVs <- df_SNVs[nchar(df_SNVs$ALT) == 1 & nchar(df_SNVs$REF) == 1,]
#' 
#' # A problem that came to my attention is that there are some variants call that point to exactly the same position in the genome but their PURPLE_AF position is different
#' # That can cause problem (specifically when running isOverlapping function). THerefore with the command below I get rid of them. For example in sample
#' # "CPCT02070491T" there were two rows for the same position that differ in ALT and PURPLE_AF fields (chrX:98875822)
#' 
#' #nrow(df_SNVs)
#' df_SNVs <- df_SNVs[isUnique(paste0(df_SNVs$chromosome,":",df_SNVs$position)),]
#' 
#' # First I let's sort the data frame based on the chromosome levels 
#' # (this step is needed for making GRanges and by doing sorting here the the order of rows will be the same for all the annotation)
#' chrOrder <-c((1:22),"X","Y")
#' df_SNVs <- df_SNVs[order(factor(df_SNVs$chromosome, levels = chrOrder, ordered = T), df_SNVs$position),]
#' 
#' 
#' # setting the row names
#' rownames(df_SNVs) <- paste0(df_SNVs$chromosome, ":", df_SNVs$position, "_", df_SNVs$REF, "/", df_SNVs$ALT)
#' 
#' # message(head(df_SNVs))
#' # 2. Mutational Signature and their likelihood --------------------------------------------------------------------------------------------------
#' # mut_sign | mut_sign_score
#' 
#' #load cosmic signatures
#' COSMIC_signatures = as.data.frame(get_known_signatures(muttype = "snv", source = "COSMIC", sig_type = "reference"))
#' rownames(COSMIC_signatures) <- tri_context_generator()
#' df_SNVs_mut_likelihood <- suppressMessages(mutation_likelihood_score(df_SNVs,
#'                                                                      signatures=COSMIC_signatures, 
#'                                                                      sampleID, 
#'                                                                      colname_chrom="chromosome",
#'                                                                      colname_position="position",
#'                                                                      colname_REF="REF",
#'                                                                      colname_ALT="ALT",
#'                                                                      colname_TNC="TNC"))
#' 
#' 
#' # 3. Chromatin status --------------------------------------------------------------------------------------------------------------------------
#' # CHROMATIN_STAT_TAD | CHROMATIN_STAT_BIN
#' 
#' # tadBedFileInput <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotation-beds/chromatin-status-tads.bed.bgz"
#' # df_SNVs_chrom_stat_tad <- getTadChromatinStat(df_SNVs, tadBedFileInput)              # 16 secs
#' df_SNVs_chrom_stat_tad_opt <- getTadChromatinStatOpt(df_SNVs, tadBedFileInput)         # 5 secs
#' # all(df_SNVs_chrom_stat_tad[df_SNVs_chrom_stat_tad[,11] != "NA",11] == df_SNVs_chrom_stat_tad_opt[!is.na(df_SNVs_chrom_stat_tad_opt[,11]),11])
#' 
#' # binBedFileInput <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotation-beds/chromatin-status-bins.bed.bgz"
#' # df_SNVs_chrom_stat_bin <- getBinChromatinStat(df_SNVs, binBedFileInput)          # 15 secs
#' df_SNVs_chrom_stat_bin_opt <- getBinChromatinStatOpt(df_SNVs, binBedFileInput)   # 5 secs
#' # all(df_SNVs_chrom_stat_bin[df_SNVs_chrom_stat_bin[,11] != "NA",11] == df_SNVs_chrom_stat_bin_opt[!is.na(df_SNVs_chrom_stat_bin_opt[,11]),11])
#' 
#' 
#' # 4. Replication timing ---------------------------------------------------------------------------------------------------------------------
#' 
#' # repTimingFileInput <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotation-beds/repliseq_encode_mean_binned.bed.bgz"
#' # df_SNVs_rep_timing <- getRepliSeq(df_SNVs, repTimingFileInput)          # 50 secs
#' df_SNVs_rep_timing_opt <- getRepliSeqOpt(df_SNVs, repTimingFileInput)   # 17 secs
#' # table(df_SNVs_rep_timing$rep_timing_bin)
#' # table(df_SNVs_rep_timing_opt$rep_timing_bin)
#' # sum(is.na(df_SNVs_rep_timing_opt$rep_timing_bin))
#' 
#' # Boxtel paper data
#' # repTimingFileInputBoxtel <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotation-beds/boxtel_paper_all_RepliSeq_median.bed.bgz"
#' df_SNVs_rep_timing_opt_boxtel <- getRepliSeqOpt(df_SNVs, repTimingFileInputBoxtel)
#' 
#' # 5. CpG islands ---------------------------------------------------------------------------------------------------------------------------
#' # cpg | percpg | pergc
#' 
#' # cpgInputFile <- "/home/ali313/Documents/studies/master/umc-project/data/cpg-island/dataset/cpgIslandExtUnmasked.bed.bgz"
#' # df_SNVs_cpg_stat <- getCpgAnn(df_SNVs, cpgInputFile)              # 7 secs
#' #### df_SNVs_cpg_stat_opt <- getCpgAnnOpt(df_SNVs, cpgInputFile)         # 4 secs   
#' 
#' #all(df_SNVs_cpg_stat_opt[df_SNVs_cpg_stat_opt$cpg,11:13] == df_SNVs_cpg_stat[df_SNVs_cpg_stat$cpg,11:13])
#' 
#' 
#' # @ Making GRnages object -----------------------------------------------------------------------------------------------------------------
#' 
#' df_SNVs_gr <- GRanges(seqnames = df_SNVs$chromosome,ranges = IRanges(start = df_SNVs$position,end = df_SNVs$position),
#'                       ref = df_SNVs$REF, alt = df_SNVs$ALT) 
#' names(df_SNVs_gr) <- paste0(df_SNVs$chromosome, ":", df_SNVs$position, "_", df_SNVs$REF, "/", df_SNVs$ALT)
#' 
#' 
#' # 6. Trp orientation ---------------------------------------------------------------------------------------------------------------------
#' # trp_str_ann
#' 
#' df_SNVs_trp_str_gr <- getTrpStrandAnn(df_SNVs_gr)
#' 
#' 
#' # 7. Rep orientation ----------------------------------------------------------------------------------------------------------------------
#' # rep_str_ann
#' 
#' # repOrientationInputFile <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotation-beds/replication-direction-tableTerritories_Haradhvala_territories.rds.gz"
#' df_SNVs_rep_str_gr <- getRepliStrandAnn(df_SNVs_gr, repOrientationInputFile)
#' 
#' 
#' # 8. Functional genomic units -------------------------------------------------------------------------------------------------------------
#' # genomic_func
#' 
#' # the dumb way that would take 6 min for annotating 10,000 calls (look at the code, the commented out parts are the dumb way of doing it)
#' # That's the smart way of doing that which takes 15 secs
#' print("yeah")
#' df_SNVs_func_unit_gr <- getFuncUnitAnn(df_SNVs_gr)
#' print("yeah?")
#' 
#' # all(df_SNVs_func_unit_gr_smart == df_SNVs_func_unit_gr)
#' 
#' 
#' 
#' 
#' # @ getting all annotations together -----------------------------------------------------------------------------------------------------
#' # checking the number of rows and their order is the same:
#' 
#' # nrow(df_SNVs)
#' # nrow(df_SNVs_mut_likelihood)
#' # nrow(df_SNVs_chrom_stat_tad)
#' # nrow(df_SNVs_chrom_stat_bin)
#' # nrow(df_SNVs_rep_timing)
#' # nrow(df_SNVs_cpg_stat_opt)
#' # length(df_SNVs_trp_str_gr)
#' # length(df_SNVs_rep_str_gr)
#' # length(df_SNVs_func_unit_gr)
#' # all(rownames(df_SNVs_chrom_stat_tad) == rownames(df_SNVs))
#' # all(names(df_SNVs_func_unit_gr) == rownames(df_SNVs))
#' 
#' 
#' # final_ann_df <- cbind(df_SNVs_mut_likelihood, df_SNVs_chrom_stat_tad_opt$CHROMATIN_STAT_TAD, df_SNVs_chrom_stat_bin_opt$CHROMATIN_STAT_BIN,
#' #                       df_SNVs_rep_timing_opt$rep_timing_bin, df_SNVs_rep_timing_opt_boxtel$rep_timing_bin,
#' #                       df_SNVs_cpg_stat_opt[,11:13], mcols(df_SNVs_trp_str_gr)$trp_str_ann, mcols(df_SNVs_rep_str_gr)$rep_str_ann,
#' #                       mcols(df_SNVs_func_unit_gr)$genomic_func)
#' 
#' final_ann_df <- cbind(df_SNVs_mut_likelihood, df_SNVs_chrom_stat_tad_opt$CHROMATIN_STAT_TAD, df_SNVs_chrom_stat_bin_opt$CHROMATIN_STAT_BIN,
#'                       df_SNVs_rep_timing_opt$rep_timing_bin, df_SNVs_rep_timing_opt_boxtel$rep_timing_bin,
#'                       mcols(df_SNVs_trp_str_gr)$trp_str_ann, mcols(df_SNVs_rep_str_gr)$rep_str_ann,
#'                       mcols(df_SNVs_func_unit_gr)$genomic_func)
#' 
#' 
#' 
#' # colnames(final_ann_df)[c(1,2,13:22)] <- c("CHROM", "POS", "chromatin_status_tads", "chromatin_status_bins", "rep_timing", "rep_timing_boxtel", "cpg", "percpg",
#' #                                    "pergc", "trp_str_ann", "rep_str_ann", "genomic_func")
#' 
#' print("now?")
#' colnames(final_ann_df)[c(1,2,13:19)] <- c("CHROM", "POS", "chromatin_status_tads", "chromatin_status_bins", "rep_timing", "rep_timing_boxtel", 
#'                                           "trp_str_ann", "rep_str_ann", "genomic_func")
#' print("What about now?")
#' # This is a temporary solution for fixing the transcription strand orientation
#' index_for_transcribed <- which((final_ann_df$REF == "A" | final_ann_df$REF == "G") & final_ann_df$trp_str_ann == "transcribed")
#' index_for_untranscribed <- which((final_ann_df$REF == "A" | final_ann_df$REF == "G") & final_ann_df$trp_str_ann == "untranscribed")
#' 
#' final_ann_df[index_for_transcribed,]$trp_str_ann <- "untranscribed"
#' final_ann_df[index_for_untranscribed,]$trp_str_ann <- "transcribed"
#' 
#' 
#' # This is a temporary solution for fixing the replication strand orientation
#' 
#' index_for_right <- which((final_ann_df$REF == "A" | final_ann_df$REF == "G") & final_ann_df$rep_str_ann == "right")
#' index_for_left <- which((final_ann_df$REF == "A" | final_ann_df$REF == "G") & final_ann_df$rep_str_ann == "left")
#' 
#' final_ann_df[index_for_right,]$rep_str_ann <- "left"
#' final_ann_df[index_for_left,]$rep_str_ann <- "right"
#' 
#' 
#' # write.table(final_ann_df, "/home/ali313/Desktop/test.txt", sep = "\t", quote = F, row.names = F)
#' 
#' 
#' 
#' if (dir.exists("/hpc/cuppen/")){
#'   write.table(final_ann_df, file = gzfile(paste0("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/", sampleID, ".txt.gz")), quote = F, row.names = F, sep = "\t")
#' } else {
#'   write.table(final_ann_df, file = gzfile(paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/", sampleID, ".txt.gz")), quote = F, row.names = F, sep = "\t")
#' }
#' 



## Adding new ms annotations to the files

library(stringr)
library(dplyr)
local <- "/home/ali313/Documents/studies/master/umc-project"





if (dir.exists("/hpc/cuppen/")){
  wd <- "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/"
} else {
  wd <- paste0(local,"/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/")
}


if (dir.exists("/hpc/cuppen/")){
  metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2-metadata_25082021.tsv"), sep = "\t", header = T, stringsAsFactors = F)
} else {
  metadata <- read.csv(file = paste0(wd, "external-files/cancer_types_HMF_PCAWG_2-metadata_25082021.tsv"), sep = "\t", header = T, stringsAsFactors = F)
}

metadata_included <- metadata[!(metadata$is_blacklisted),]

metadata_included$tmb <- rowSums(metadata_included[,22:23])



# args[1]:args[2] # It takes longer than expected when run it on one node. Submit multiple jobs next time

for (i in 1:nrow(metadata_included)){
  print(i)
  sampleID <- metadata_included$sample_id[i]
  
  if (dir.exists("/hpc/cuppen/")){
    path_to_ann <- paste0("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/", sampleID, ".txt.gz")
  } else {
    path_to_ann <- paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/", sampleID, ".txt.gz")
  }
  
  
  if (dir.exists("/hpc/cuppen/")){
    path_to_ms <- paste0("/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/muts_assigned/", sampleID, ".txt.gz")
  } else {
    path_to_ms <- paste0(local, "/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/04_sigProfiler/sig_contrib/muts_assigned/", sampleID, ".txt.gz")
  }
  
  if (file.exists(path_to_ann) & file.exists(path_to_ms)){
    
    final_ann_df <- read.csv(file = path_to_ann, stringsAsFactors = F, header = T, sep = "\t")
    final_ann_df$CHROM <- as.character(final_ann_df$CHROM)
    
    ms_df <- read.csv(file = path_to_ms, stringsAsFactors = F, header = T, sep = "\t")
    ms_df$chrom <- as.character(ms_df$chrom)
    
    colnames(final_ann_df)[11:12] <- c("mut_sign_internship", "mut_sign_score_internship")
    
    colnames(ms_df) <- c("CHROM", "POS", "REF", "ALT", "CONTEXT", "mut_sign_updated", "mut_sign_score_updated")
    ms_df$CHROM <- sapply(str_split(ms_df$CHROM, pattern = "r"), "[[", 2)
    
    shared <- dplyr::inner_join(final_ann_df, ms_df[,c(1:4, 6:7)])
    
    
    if (dir.exists("/hpc/cuppen/")){
      write.table(shared, file = gzfile(paste0("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/updated/", sampleID, ".txt.gz")), quote = F, row.names = F, sep = "\t")
    } else {
      write.table(shared, file = gzfile(paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation/annotated-vcfs/updated/", sampleID, ".txt.gz")), quote = F, row.names = F, sep = "\t")
    }

  }
  
}
