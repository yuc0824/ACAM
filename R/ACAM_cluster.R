#' This is the description of the function ACAM_cluster.
#' @title ACAM_cluster
#' 
#' @description The function ACAM_cluster conducts the five clustering methods and output the four results of them which have the largest pairwise ARI.
#' 
#' @details The function ACAM_cluster conducts the five clustering methods and output the four results of them which have the largest pairwise ARI.
#' 
#' @param DF The input dataset. Make sure that the dataset is lognormalized. Note that the rows are cells, and the columns are marker genes.
#' @param threshold The threshold that any dataset with the sample size larger than this threshold should be used an accelerated clustering procedure. Default is 4000.
#' @param SEED The random seed.
#' @return The clustering results \it{cluster_results}.
#' @export 



# library(SAMEclustering)
# library(SingleCellExperiment)
# library(SC3)
# library(cidr)
# library(ADPclust)
# library(SIMLR)



ACAM_cluster <- function(DF, threshold = 4000, SEED = 123){
  if(nrow(DF) < threshold){
    clusters <- SAMEclustering::individual_clustering(inputTags = t(DF), mt_filter = TRUE,  mt.cutoff = 0.1, 
                                       percent_dropout = 0, SC3 = TRUE, gene_filter = FALSE, svm_num_cells = 5000, CIDR = TRUE, nPC.cidr = NULL, 
                                       Seurat = TRUE, nGene_filter = FALSE, low.genes = 0, high.genes = 0, nPC.seurat = 10, resolution = 0.7, 
                                       tSNE = TRUE, dimensions = 3, perplexity = 30, tsne_min_cells = 200, tsne_min_perplexity = 10, var_genes = NULL, 
                                       SIMLR = FALSE, diverse = FALSE, SEED = SEED)
  }else{
    clusters <- self_clustering(inputTags = t(DF), k_fixed = 15, SEED = SEED)
  }
  cluster_results <- Louvain(clusters)
  return(cluster_results)
}


self_clustering <- function(inputTags, percent_dropout = 10, tsne_min_cells = 200, dimensions = 3,
                            tsne_min_perplexity = 10, perplexity = 30, k_fixed = 0, SEED = 123){
  cluster_number <- NULL
  clusters <- NULL
  inputTags = as.matrix(inputTags)
  #SC3
  
  sc3OUTPUT <- sc3_SELF(inputTags = inputTags,percent_dropout = percent_dropout,k_fixed = k_fixed, SEED = SEED)
  clusters <- rbind(clusters, matrix(c(sc3OUTPUT), 
                                                   nrow = 1, byrow = TRUE))
  if(max(c(sc3OUTPUT)) > 0){
    cluster_number <- c(cluster_number, max(c(sc3OUTPUT)))
  }
  rm(sc3OUTPUT)
  #CIDR
  
  cidrOUTPUT <- cidr_SELF(inputTags = inputTags,percent_dropout = percent_dropout, SEED = SEED)
  clusters <- rbind(clusters, matrix(c(cidrOUTPUT@clusters), 
                                                   nrow = 1, byrow = TRUE))
  if(cidrOUTPUT@nCluster > 0){
    cluster_number <- c(cluster_number, cidrOUTPUT@nCluster)
  }
  rm(cidrOUTPUT)
  #Seurat
  
  seurat_output <- seurat_SELF(inputTags = inputTags,percent_dropout = percent_dropout)
  clusters <- rbind(clusters, matrix(c(seurat_output), 
                                                   nrow = 1, byrow = TRUE))
  if(max(!is.na(seurat_output)) > 0){
    cluster_number <- c(cluster_number, max(!is.na(seurat_output)))
  }
  rm(seurat_output)
  #tsne
  
  if (length(inputTags[1, ]) < tsne_min_cells) {
    perplexity = tsne_min_perplexity
  }
  tsne_kmeansOUTPUT <- tSNE_kmeans_SELF(inputTags = as.matrix(inputTags),
                                        percent_dropout = percent_dropout, k.min = 2, k.max = max(cluster_number),SEED = SEED)
  clusters <- rbind(clusters, matrix(c(tsne_kmeansOUTPUT$cluster), 
                                                   nrow = 1, byrow = TRUE))
  if(max(as.numeric(tsne_kmeansOUTPUT$cluster)) > 0){
    cluster_number <- c(cluster_number, max(as.numeric(tsne_kmeansOUTPUT$cluster)))
  }
  rm(tsne_kmeansOUTPUT)
  #SIMLR
  
  simlrOUTPUT <- SIMLR_SELF(inputTags = inputTags,percent_dropout = percent_dropout, 
                            k.min = 2, k.max = max(cluster_number), k_fixed = k_fixed, SEED = SEED)
  clusters <- rbind(clusters, simlrOUTPUT$y$cluster)
  rm(simlrOUTPUT)
  
  #Final
  message("Starting Selection...")
  rownames(clusters) <- c("SC3", "CIDR", 
                                 "Seurat", "tSNE+kmeans", "SIMLR")
  ARI = matrix(0, 5, 5)
  rownames(ARI) <- c("SC3", "CIDR", "Seurat", 
                     "tSNE+kmeans", "SIMLR")
  colnames(ARI) <- c("SC3", "CIDR", "Seurat", 
                     "tSNE+kmeans", "SIMLR")
  for (i in c("SC3", "CIDR", "Seurat", 
              "tSNE+kmeans", "SIMLR")) {
    for (j in c("SC3", "CIDR", "Seurat", 
                "tSNE+kmeans", "SIMLR")) {
      ARI[i, j] <- adjustedRandIndex(unlist(clusters[i, 
      ]), unlist(clusters[j, ]))
  }
  m1 <- which.min(apply(ARI, 1, var))
  clusters <- clusters[-m1, ]
  return(clusters)
  }
}


sc3_SELF <- function(inputTags, percent_dropout = 10, svm_num_cells = 1000, k_fixed, SEED = 123) 
{
  message("Performing SC3 clustering...")
  exp_cell_exprs <- NULL
  sc3OUTPUT <- NULL
  if (is.null(percent_dropout)) {
    inputTags <- inputTags
  }else {
    dropouts <- rowSums(inputTags == 0)/ncol(inputTags) * 100
    inputTags <- inputTags[-c(which(dropouts <= percent_dropout), 
                              which(dropouts >= 100 - percent_dropout)), ]
  }
  exp_cell_exprs <- SingleCellExperiment(assays = list(counts = inputTags))
  normcounts(exp_cell_exprs) <- t(t(inputTags)/colSums(inputTags)) * 1e+06
  logcounts(exp_cell_exprs) <- log2(normcounts(exp_cell_exprs) + 1)
  rowData(exp_cell_exprs)$feature_symbol <- rownames(exp_cell_exprs)
  exp_cell_exprs <- exp_cell_exprs[!duplicated(rowData(exp_cell_exprs)$feature_symbol), ]
  gc()
  
  if(k_fixed == 0){
    exp_cell_exprs <- sc3_estimate_k(exp_cell_exprs)
    gc()
    optimal_K <- metadata(exp_cell_exprs)$sc3$k_estimation
  }else{
    optimal_K <- k_fixed
  }
  if (ncol(inputTags) < svm_num_cells) {
    gc()
    exp_cell_exprs <- sc3(exp_cell_exprs, ks = optimal_K, biology = FALSE, gene_filter = F, n_cores = 1, kmeans_iter_max = 5000,rand_seed = SEED)
  }else if (ncol(inputTags) >= svm_num_cells) {
    gc()
    exp_cell_exprs <- sc3(exp_cell_exprs, ks = optimal_K, 
                          biology = FALSE, gene_filter = F, svm_max = svm_num_cells, 
                          svm_num_cells = svm_num_cells, n_cores = 1, kmeans_iter_max = 5000, rand_seed = SEED)
    gc()
    exp_cell_exprs <- sc3_run_svm(exp_cell_exprs, ks = optimal_K)
  }
  gc()
  p_Data <- colData(exp_cell_exprs)
  col_name <- paste("sc3_", optimal_K, "_clusters", 
                    sep = "")
  sc3OUTPUT <- p_Data[, grep(col_name, colnames(p_Data))]
  return(sc3OUTPUT)
}

cidr_SELF <- function(inputTags, percent_dropout = 10, nPC.cidr = NULL, SEED = 123) 
{
  message("Performing CIDR clustering...")
  set.seed(SEED)
  cidrOUTPUT <- NULL
  if (is.null(percent_dropout)) {
    inputTags_cidr <- inputTags
  }else {
    dropouts <- rowSums(inputTags == 0)/ncol(inputTags) * 100
    inputTags_cidr <- inputTags[-c(which(dropouts <= percent_dropout), 
                                   which(dropouts >= 100 - percent_dropout)), ]
  }
  cidrOUTPUT <- scDataConstructor(inputTags_cidr, tagType = "raw")
  gc()
  cidrOUTPUT <- determineDropoutCandidates(cidrOUTPUT)
  gc()
  cidrOUTPUT <- wThreshold(cidrOUTPUT)
  gc()
  cidrOUTPUT <- scDissim(cidrOUTPUT)
  gc()
  cidrOUTPUT <- scPCA(cidrOUTPUT, plotPC = F)
  gc()
  cidrOUTPUT <- nPC(cidrOUTPUT)
  gc()
  if (!is.null(nPC.cidr)) {
    cidrOUTPUT@nPC <- nPC.cidr
  }else {
    nPC.cidr <- cidrOUTPUT@nPC
  }
  cidrOUTPUT <- scCluster(cidrOUTPUT, nPC = nPC.cidr)
  return(cidrOUTPUT)
}

seurat_SELF <- function (inputTags, percent_dropout = 10, nGene_filter = FALSE, low.genes = 200, high.genes = 8000, 
                         nPC.seurat = NULL, resolution = 0.7, SEED = 123) 
{
  message("Performing Seurat clustering...")
  seuratOUTPUT <- NULL
  if (is.null(percent_dropout)) {
    inputTags <- inputTags
  }else {
    dropouts <- rowSums(inputTags == 0)/ncol(inputTags) * 100
    inputTags <- inputTags[-c(which(dropouts <= percent_dropout), 
                              which(dropouts >= 100 - percent_dropout)), ]
  }
  gc()
  seuratOUTPUT <- CreateSeuratObject(counts = inputTags, min.cells = 0, 
                                     min.features = 0, project = "single-cell clustering")
  if (nGene_filter == TRUE) {
    seuratOUTPUT <- subset(object = seuratOUTPUT, subset = nFeature_RNA > 
                             low.genes & nFeature_RNA < high.genes)
  }
  seuratOUTPUT = NormalizeData(object = seuratOUTPUT, normalization.method = "LogNormalize", 
                               scale.factor = 10000)
  gc()
  seuratOUTPUT = FindVariableFeatures(object = seuratOUTPUT, 
                                      selection.method = "vst", nfeatures = 2000)
  gc()
  all.genes <- rownames(seuratOUTPUT)
  seuratOUTPUT <- ScaleData(object = seuratOUTPUT, features = all.genes)
  if (nPC.seurat <= 20 || is.null(nPC.seurat)) {
    gc()
    seuratOUTPUT <- RunPCA(object = seuratOUTPUT, features = VariableFeatures(object = seuratOUTPUT), 
                           npcs = 20, seed.use = SEED, verbose = FALSE)
    gc()
    seuratOUTPUT <- FindNeighbors(seuratOUTPUT, dims = 1:20, 
                                  verbose = FALSE)
  }else {
    gc()
    seuratOUTPUT <- RunPCA(object = seuratOUTPUT, features = VariableFeatures(object = seuratOUTPUT), 
                           npcs = nPC.seurat, seed.use = SEED, verbose = FALSE)
    gc()
    seuratOUTPUT <- FindNeighbors(seuratOUTPUT, dims = 1:20, 
                                  verbose = FALSE)
  }
  gc()
  seuratOUTPUT <- FindClusters(object = seuratOUTPUT, resolution = resolution, 
                               verbose = FALSE)
  if (length(seuratOUTPUT@active.ident) < ncol(inputTags)) {
    seurat_output <- matrix(NA, ncol = ncol(inputTags), byrow = T)
    colnames(seurat_output) <- colnames(inputTags)
    seurat_retained <- t(as.matrix(as.numeric(seuratOUTPUT@active.ident)))
    colnames(seurat_retained) <- names(seuratOUTPUT@active.ident)
    for (i in 1:ncol(seurat_retained)) {
      seurat_output[1, colnames(seurat_retained)[i]] <- seurat_retained[1, 
                                                                        colnames(seurat_retained)[i]]
    }
  }else {
    seurat_output <- t(as.matrix(as.numeric(seuratOUTPUT@active.ident)))
  }
  return(seurat_output)
}

tSNE_kmeans_SELF <- function(inputTags, percent_dropout = 10, dimensions = 3, perplexity = 30, tsne_min_perplexity = 10,
                             k.min = 2, k.max = max(cluster_number), var_genes = NULL, SEED = 123) 
{
  message("Performing tsne-kmeans clustering...")
  input_lcpm <- NULL
  tsne_input <- NULL
  tsne_output <- NULL
  tsne_kmeansOUTPUT <- NULL
  adpOUTPUT <- NULL
  if (is.null(percent_dropout)) {
    inputTags_tsne <- inputTags
  }else {
    dropouts <- rowSums(inputTags == 0)/ncol(inputTags) * 100
    inputTags_tsne <- inputTags[-c(which(dropouts <= percent_dropout), 
                                   which(dropouts >= 100 - percent_dropout)), ]
  }
  tsne_input <- inputTags_tsne
  if (is.null(var_genes)) {
    set.seed(SEED)
    gc()
    tsne_output <- Rtsne(t(tsne_input), dims = dimensions, 
                         perplexity = perplexity, check_duplicates = FALSE)
  }else {
    se_genes = rep(NA, nrow(tsne_input))
    for (i in 1:nrow(tsne_input)) {
      se_genes[i] = sqrt(var(tsne_input[i, ])/length(tsne_input[i, ]))
    }
    decreasing_rank = order(se_genes, decreasing = TRUE)
    set.seed(SEED)
    gc()
    tsne_output <- Rtsne(t(tsne_input[decreasing_rank[1:var_genes], 
    ]), dims = dimensions, perplexity = perplexity)
  }
  gc()
  adpOUTPUT <- adpclust(tsne_output$Y, htype = "amise", 
                        centroids = "auto", nclust = k.min:k.max)
  gc()
  tsne_kmeansOUTPUT <- kmeans(tsne_output$Y, tsne_output$Y[adpOUTPUT$centers, ], adpOUTPUT$nclust, iter.max = 5000)
  gc()
  return(tsne_kmeansOUTPUT)
}


SIMLR_SELF <- function(inputTags, percent_dropout = 10, k.min = 2, k.max = max(cluster_number), k_fixed, SEED = 123) 
{
  message("Performing SIMLR clustering...")
  set.seed(SEED)
  if (is.null(percent_dropout)) {
    inputTags_simlr <- inputTags
  }else {
    dropouts <- rowSums(inputTags == 0)/ncol(inputTags) * 
      100
    inputTags_simlr <- inputTags[-c(which(dropouts <= percent_dropout), 
                                    which(dropouts >= 100 - percent_dropout)), ]
  }
  k_range <- k.min:k.max
  if(k_fixed == 0){
    best_k <- SIMLR_Estimate_Number_of_Clusters(inputTags_simlr, NUMC = k_range, cores.ratio = 1)
    k <- which.min(best_k$K1) + k_range[1] - 1
  }else{
    k = k_fixed
  }
  gc()
  if (dim(inputTags_simlr)[2] < 1000) {
    simlrOUTPUT <- SIMLR(inputTags_simlr, c = k.max, 
                         cores.ratio = 0)
  }else {
    simlrOUTPUT <- SIMLR_Large_Scale(log10(inputTags_simlr + 1), c = k)
  }
  return(simlrOUTPUT)
}


