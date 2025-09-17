library(spectralGraphTopology)
library(igraph)
library(Matrix)
#' @title Plot a clustering graph from Laplacian matrix
#'
#' Plots a clustering graph based on a given Laplacian matrix and the ground-truth labels for the clusters. 
#' @param Laplacian r by r Laplacian matrix.
#' @param true_labels the numeric vector of true labels (of length r)
#' @param names the vector of names of the elements, e.g., stocks (of length r)
#' @param heavy_type a string which selects the statistical distribution of the data    .
#'        Valid values are "gaussian" or "student".
#' @param nu the degrees of freedom of the Student-t distribution.
#'        Must be a real number greater than 2.
#' @param sigma_e hyperparameter that controls graph weight sparsity and time-consistency
#' @param gamma hyperparameter that controls the sparsity of VAR coefficients in the variations of the weights
#' @param w0 initial vector of graph weights. Either a vector of length p(p-1)/2 or
#'        a string indicating the method to compute an initial value.
#' @param a0 initial value of the VAR coefficient
#' @param eta hyperparameter that controls the effect of the additional regularization to obtain a
#'        k-component graph
#' @param update_eta whether to update eta during the optimization.
#' @param d the nodes' degrees. Either a vector or a single value.
#' @param rho ADMM hyperparameter.
#' @param update_rho whether or not to update rho during the optimization.
#' @param maxiter maximum number of iterations.
#' @param reltol relative tolerance as a convergence criteria.
#' @param verbose whether or not to show a progress bar during the iterations.
#' @return The network object:
#' \item{\code{graph_net}}{the graph network object that can be used for plotting}
#' \item{\code{metric}}{the results of the clustering measures}
#' @import spectralGraphTopology
#' @import igraph
#' @import Matrix
#' @export

  plot_graph <- function(Laplacian, true_labels, names) {
  # number of nodes
  r <- nrow(Laplacian)
  

  w <- spectralGraphTopology:::Linv(Laplacian)
  adjacency <- spectralGraphTopology:::A(w)
  
  # number of components
  q = r - rankMatrix(Laplacian)
  
  graph_net <- graph_from_adjacency_matrix(adjacency, mode = "undirected", weighted = TRUE)
  
  
  # where do our predictions differ from true labels?
  metric <- evaluate_clustering(graph_net, true_labels, r, q)
  labels_pred_adj <-  metric$labels_pred_adj
  
  
  # if q>10 add some colors to the list
  colors <- c("#55efc4", "#ff7675", "#0984e3", "#a29bfe", "#B33771", "#48dbfb", "#FDA7DF", "#C4E538", "#0184a3", "#26dcfb")
  colors <- colors[1:q]
  
  
  
  # ground truth coloring 
  V(graph_net)$color <- c(colors[true_labels])
  V(graph_net)$type <- c(rep(FALSE, r))
  V(graph_net)$cluster <- c(true_labels)
  E(graph_net)$color <- apply(
    as.data.frame(get.edgelist(graph_net)), 1,
    function(x) {
      ifelse(V(graph_net)$cluster[x[1]] == V(graph_net)$cluster[x[2]],
             colors[V(graph_net)$cluster[x[1]]], "grey"
      )
    }
  )
  
  mask <- labels_pred_adj != true_labels
  node_labels <- names
  node_labels[!mask] <- NA
  label_colors <- rep("black",r)
  # label_colors[mask] <- "red"
  
  # plot network
  plot(graph_net,
       vertex.size = c(rep(4, r)),
       vertex.label = c(node_labels),
       vertex.label.color = label_colors,
       vertex.label.cex = 0.8, vertex.label.dist = 0.5,
       vertex.frame.color = c(colors[true_labels]),
       layout = layout_nicely(graph_net),
       vertex.label.family = "Helvetica", vertex.label.color = "black",
       vertex.shape = c(rep("circle", r)),
       edge.width = 3 * E(graph_net)$weight
  )
  
  
  # implied clusters
  true_labels_res <- true_labels
  true_labels <- igraph::components(graph_net)$membership
  V(graph_net)$color <- c(colors[true_labels])
  V(graph_net)$type <- c(rep(FALSE, r))
  V(graph_net)$cluster <- c(true_labels)
  E(graph_net)$color <- apply(
    as.data.frame(get.edgelist(graph_net)), 1,
    function(x) {
      ifelse(V(graph_net)$cluster[x[1]] == V(graph_net)$cluster[x[2]],
             colors[V(graph_net)$cluster[x[1]]], "grey"
      )
    }
  )
  
  
  node_labels <- rep(NA, r)
  # plot network
  plot(graph_net,
       vertex.size = c(rep(4, r)),
       vertex.label = c(node_labels),
       vertex.label.cex = 0.8, vertex.label.dist = 0.5,
       vertex.frame.color = c(colors[true_labels]),
       layout = layout_nicely(graph_net),
       vertex.label.family = "Helvetica", vertex.label.color = "black",
       vertex.shape = c(rep("circle", r)),
       edge.width = 3 * E(graph_net)$weight
  )
  
  
  
  # print("metric:\n")
  # print(metric$accuracy_adj)
  # print(metric$purity)
  # print(metric$mod_gt)
  # print(metric$ARI)
  
  
  results <- list(graph_net = graph_net, metric = metric)
  
  return(results)
}
