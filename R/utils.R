#' Calculate the mode or most frequent element (internal)
#'
#' @param v input vector 
#' @return value of mode
#' @keywords internal
#' 
Mode <- function(v) {
  uniqv <- unique(v)
  return ( uniqv[which.max(tabulate(match(v, uniqv)))] )
}


#' Calculate balancedness (internal)
#'
#' @param netobj graph network object
#' @param p number of assets (nodes)
#' @param q number of clusters
#' @return balancedness score
#' @import igraph
#' @keywords internal
#' 
Balancedness <- function(netobj, p, q) {
  clust_size <- igraph::components(netobj)$csize
  k <- length(clust_size)
  if (k<q)
  {
    clust_size <- c(rep(0,q-k), clust_size)
  }
  return( sum(abs(clust_size - p/q))/q )
}



#' Calculate normalized (0,1) balancedness (internal)
#'
#' @param netobj graph network object
#' @param p number of assets (nodes)
#' @param q number of clusters
#' @return normalized balancedness score
#' @import igraph
#' @keywords internal
#' 
Balancedness_norm <- function(netobj, p, q) {
  clust_size <- igraph::components(netobj)$csize
  k <- length(clust_size)
  if (k<q)
  {
    clust_size <- c(rep(0,q-k), clust_size)
  }
  num <- sum(abs(clust_size - p/q))
  clust_size <- c(rep(1,(q-1)), p-(q-1))
  dennum <- sum(abs(clust_size - p/q))
  return( 1 - num/dennum )
}


#' Calculate GINI (internal)
#'
#' @param netobj graph network object
#' @return GINI score
#' @keywords internal
#' 
GINI <- function(netobj) {
  clust_size <- igraph::components(netobj)$csize
  k <- length(clust_size)
  clust_size_avg <- mean(clust_size)
  num <- 0 
  for (i in 1:k)
  {
    for (j in 1:k) {
      num <- num + abs(clust_size[i]-clust_size[j]) }
  }
  return(num/ (2*k^2*clust_size_avg))
}



#' @title Evaluate Clustering Performance
#'
#' Clustering evaluation using multiple metrics
#'
#' @param net graph network object 
#' @param true_labels True cluster labels (optional)
#' @param p number of assets (nodes)
#' @param q number of clusters
#' @return Evaluation metrics
#' @export
#' 
evaluate_clustering <- function(net, true_labels, p, q) {
  labels_pred <- rep(0, p)
  memberships <- igraph::components(net)$membership
  for (j in 1:q){
    idx <- memberships %in% c(j)
    labels_pred[idx] <- Mode(true_labels[idx])
  }
  
  mask <- labels_pred != true_labels
  purity <- 1- sum(mask)/length(mask)
  
  
  labels_pred <- memberships
  labels_pred_sorted <- labels_pred 
  perms <- permn(c(1:q))
  acc_max <- 0
  for (k in 1:length(perms)){
    perm <- perms[[k]]
    for (j in 1:q){
      idx <- memberships %in% j
      labels_pred_sorted[idx] <- perm[j]
      
    }
    mask <- labels_pred_sorted != true_labels
    acc <- 1- sum(mask)/length(mask)
    if (acc>= acc_max) {
      acc_max <- acc
      ind_max <- k
    }
  }
  perm <- perms[[ind_max]]
  for (j in 1:q){
    idx <- labels_pred %in% j
    labels_pred_sorted[idx] <- perm[j]
  }
  
  NMI <- randnet::NMI(labels_pred_sorted, true_labels)
  
  ARI <- mclust::adjustedRandIndex(labels_pred_sorted, true_labels)
  mask <- labels_pred_sorted != true_labels
  accuracy <- 1- sum(mask)/length(mask)
  
  
  
  
  mod_gt <- modularity(net, true_labels)
  
  balanced_norm <- Balancedness_norm(net, p, q)
  GINI_metric <- GINI(net)
  
  metrics <- list( memberships = memberships,  labels_pred = labels_pred_sorted,  
                   accuracy = accuracy, purity = purity, 
                   mod = mod_gt,
                   balanced = balanced_norm, GINI = GINI_metric,
                   NMI = NMI, ARI = ARI)
  return( metrics )
}


# 
# spectral_clustering <- function(X, # matrix of data points
#                                 nn = 10, # the k nearest neighbors to consider
#                                 n_eig = 2) # m number of eignenvectors to keep
# {
#   mutual_knn_graph <- function(X, nn = 10)
#   {
#     D <- as.matrix( dist(X) ) # matrix of euclidean distances between data points in X
#     
#     # intialize the knn matrix
#     knn_mat <- matrix(0,
#                       nrow = nrow(X),
#                       ncol = nrow(X))
#     
#     # find the 10 nearest neighbors for each point
#     for (i in 1: nrow(X)) {
#       neighbor_index <- order(D[i,])[2:(nn + 1)]
#       knn_mat[i,][neighbor_index] <- 1 
#     }
#     
#     # Now we note that i,j are neighbors iff K[i,j] = 1 or K[j,i] = 1 
#     knn_mat <- knn_mat + t(knn_mat) # find mutual knn
#     
#     knn_mat[ knn_mat == 2 ] = 1
#     
#     return(knn_mat)
#   }
#   
#   graph_laplacian <- function(W, normalized = TRUE)
#   {
#     stopifnot(nrow(W) == ncol(W)) 
#     
#     g = colSums(W) # degrees of vertices
#     n = nrow(W)
#     
#     if(normalized)
#     {
#       D_half = diag(1 / sqrt(g) )
#       return( diag(n) - D_half %*% W %*% D_half )
#     }
#     else
#     {
#       return( diag(g) - W )
#     }
#   }
#   
#   W = mutual_knn_graph(X) # 1. matrix of similarities
#   L = graph_laplacian(W) # 2. compute graph laplacian
#   ei = eigen(L, symmetric = TRUE) # 3. Compute the eigenvectors and values of L
#   n = nrow(L)
#   return(list(Eigv = ei$vectors[,(n - n_eig):(n - 1)], W=W, L= L)) # return the eigenvectors of the n_eig smallest eigenvalues
#   
# }
