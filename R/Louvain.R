#' This is the description of the function Louvain.
#' @title Louvain
#'
#' @description The function Louvain assembles the result of five clustering methods used in the function \it{ACAM_cluster}.
#'
#' @details You can use this function to calculate x+1,then return the value of x+1.
#'
#' @param X The input dataset.
#' @return The clustering results \code{cluster_results}.



Louvain <- function(X){
  Method1.G <- matrix(0,nrow = ncol(X), ncol = ncol(X))
  Method2.G <- matrix(0,nrow = ncol(X), ncol = ncol(X))
  Method3.G <- matrix(0,nrow = ncol(X), ncol = ncol(X))
  Method4.G <- matrix(0,nrow = ncol(X), ncol = ncol(X))

  for(i in 1:ncol(X)){
    for(j in 1:i){
      if (X[1,i] == X[1,j]){
        Method1.G[i,j] <- 1
        Method1.G[j,i] <- 1
      }
      if (X[2,i] == X[2,j]){
        Method2.G[i,j] <- 1
        Method2.G[j,i] <- 1
      }
      if (X[3,i] == X[3,j]){
        Method3.G[i,j] <- 1
        Method3.G[j,i] <- 1
      }
      if (X[4,i] == X[4,j]){
        Method4.G[i,j] <- 1
        Method4.G[j,i] <- 1
      }
    }}
  total.G <- Method1.G + Method2.G + Method3.G + Method4.G
  diag(total.G) <- 0
  total.G <- total.G/4
  total.G[total.G < 1] <- 0

  CM_graph <- graph_from_adjacency_matrix(total.G, mode = 'undirected')
  CM_cluster <- cluster_louvain(CM_graph, weights = NULL)
  Y.result <- CM_cluster$membership
  return(Y.result)
}
