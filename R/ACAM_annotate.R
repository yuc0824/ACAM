#' This is the description of the function ACAM_annotation.
#' @title ACAM_annotation
#'
#' @description The function ACAM_annotation annotates cells by the method ACAM.
#'
#' @details The function ACAM_annotation annotates cells by the method ACAM.
#'
#' @param DF The input dataset. Make sure that the dataset is lognormalized. Note that the rows are cells, and the columns are marker genes.
#' @param preprocess Whether the dataset need to be preprocessed. If true, DF will be lognormalized and zero expressed genes will be removed. Default is FALSE.
#' @param cluster_results The clustering results obtained from the function \code{ACAM_cluster}
#' @param gene.markers The species- and tissue-specific marker genes obtained from CellMatch database.
#' @param min_num The mininum number which any clusters have a larger size than this number will be considered as representative clusters. Default is 10.
#' @param oversample Whether we randomly select the cells until the size of the target group is the same as that of the adversarial group.
#' @param pca Whether PCA should be used in the dimension reduction.
#' @param pca.rank The dimension of which the dataset will be reduced into by the method PCA, if \code{pca == TRUE}.
#' @param umap.rank The dimension of which the dataset will be reduced into by the method umap.
#' @param umap.seed The random seed when conducting UMAP.
#' @param k.neighbors The number of neighbors considered in the final step kNN.
#' @param knn.seed The random seed when conducting kNN.
#' @return The annotation results \code{annotation_results}.
#' @export





# input data
ACAM_annotation <-function(DF,
                           preprocess = FALSE,
                           cluster_results,
                           gene.markers,
                           min_num = 10,
                           oversample = TRUE,
                           pca = FALSE,
                           pca.rank = 50,
                           umap.rank = 2,
                           umap.seed = 1,
                           k.neighbors = 1,
                           knn.seed = 1
                           )
{
  if(preprocess == TRUE){
    DF <- Seurat::LogNormalize(data = DF, scale.factor = 10000, verbose = F)
    zero_expressed <- which(apply(DF,2,sum) == 0)
    if(length(zero_expressed) >= 1){
      DF <- DF[,-zero_expressed]
    }
  }

  #representative clusters
  comb_vector <- rep(0,length(unique(cluster_results)))
  for(i in 1:length((cluster_results))){
    comb_vector[cluster_results[i]] <- comb_vector[cluster_results[i]] + 1
  }

  Y_min_in <- which(comb_vector >= min_num)
  Y_min_out <- which(comb_vector < min_num)
  Ycomb <- cluster_results
  Ycomb[which(Ycomb %in% Y_min_out)] <- 0
  Ycomb_in <- which(Ycomb != 0)
  Ycomb_out <- which(Ycomb == 0)


# eXtreme Gradient Boosting (XGBoost)
importance_names <- toupper(colnames(DF))
label <- NULL
Yclass.DF <- NULL
DFclass <- NULL
M <- NULL
X <- NULL
IM <- NULL
type <- NULL
info <- NULL
gain <- NULL
Gain <- NULL

for(i in 1:length(Y_min_in)){
  Yclass.DF[[i]] <- rep(0, length(Ycomb))
  Yclass.DF[[i]][which(Ycomb == Y_min_in[i])] <- 1

  Yclass.DF[[i]] <- as.numeric(Yclass.DF[[i]])
  #partial <- ceiling(length(Ycomb) / length(which(Ycomb == Y_min_in[i])) )
  DFclass <- cbind(Yclass.DF[[i]], DF)
 
#   #magnify
#   if(partial > 2){
#     for(k in 1:(partial- 2)){
#       DFclass <- rbind(DFclass, DFclass[which(Ycomb == Y_min_in[i]),])
#     }
#   }
  
  # oversample
  if(oversample == TRUE){
    sample_index <- sample(which(Ycomb == Y_min_in[i]), size = (length(Ycomb) - length(which(Ycomb == Y_min_in[i]))), replace = T)
    DFclass <- rbind(DFclass, DFclass[sample_index, ])
  }
  
  this_mean <- apply(DF[which(Ycomb == Y_min_in[i]), ],2,mean)
  other_mean <- apply(DF[which(Ycomb != Y_min_in[i]),],2,mean)
  if(length(which(this_mean > other_mean)) > 0){
    DFclass <- cbind(DFclass[,1],DFclass[,-1][,which(this_mean > other_mean)])

    #xgboost
    M[[i]] <- xgboost::xgb.DMatrix(data = as.matrix(DFclass[,-1]),label = DFclass[,1])
    X[[i]] <- xgboost::xgboost(M[[i]],
                           max.depth = 1,
                           eta = 0.5,
                           nround = 50,
                           objective = 'binary:hinge',
                           eval_metric = "auc")
    IM[[i]] <- xgboost::xgb.importance(importance_names[which(this_mean > other_mean)], model = X[[i]])
    IM[[i]][,2] <- as.vector(IM[[i]][,2])

    info[[i]] <- 0
    for(k in 1:nrow(IM[[i]])){
      info[[i]] <- c(info[[i]], which(as.data.frame(gene.markers)[,1]  == as.data.frame(IM[[i]])[k,1]))
    }
    info[[i]] <- info[[i]][-1]

    type[[i]] <- gene.markers[info[[i]],2]

    gain <- NULL
    mid_gene <- NULL
    for(j in 1:length(unique(type[[i]]))){
      mid_gene <- which(gene.markers[,2] == unique(type[[i]])[j])
      mid_index <- which(as.data.frame(IM[[i]])[,1] %in% as.vector(gene.markers[mid_gene,1]))
      gain[j] <- sum(IM[[i]][mid_index,2])

    }
    Gain[[i]] <- data.frame(unique(type[[i]]), gain)
    Gain[[i]] <- Gain[[i]][order(Gain[[i]][,2],decreasing = T),]
    if(nrow(Gain[[i]]) == 1){
      label[i] <- Gain[[i]][1,1]
    }else if(Gain[[i]][1,2] >= Gain[[i]][2,2]){
      label[i] <- Gain[[i]][1,1]
    }else{
      label[i] <- 'Unknown'
    }
    rm(gain)
  }else{
    label[i] <- 'Unknown'
  }
}

# kNN
if(pca == TRUE){
  DFpca <- stats::prcomp(DF, rank = pca.rank)$x
  set.seed(umap.seed)
  umap_DF <- uwot::umap(DFpca, n_components = umap.rank)
}else{
  set.seed(umap.seed)
  umap_DF <- uwot::umap(DF, n_components = umap.rank)
}


train.DF <- as.data.frame(cbind(Ycomb,umap_DF)[which(Ycomb !=0), ])
test.DF <- as.data.frame(cbind(Ycomb,umap_DF)[which(Ycomb ==0), ])
set.seed(knn.seed)
knn <- kknn::kknn(Ycomb ~.,
                   train = train.DF,
                   test = test.DF,
                   k = k.neighbors,
                   kernel = 'rectangular',
)

Ynew <- Ycomb
Ynew[which(Ycomb == 0)] <- knn$fitted.values
Yfinal <- rep(0,length(Ycomb))
for(i in 1:length(Y_min_in)){
  Yfinal[which(Ynew == Y_min_in[i])] <- label[i]
}
return(Yfinal)
}
