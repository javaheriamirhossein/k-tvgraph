#' @title Plot a clustering graph from Laplacian matrix
#'
#' Plots a clustering graph based on a given Laplacian matrix and the ground-truth labels for the clusters.
#' @param Laplacian r by r Laplacian matrix.
#' @param true_labels the numeric vector of true labels (of length r)
#' @param names the vector of names of the elements, e.g., stocks (of length r)
#' @return The network object:
#' \item{\code{graph_net}}{the graph network object that can be used for plotting}
#' \item{\code{metric}}{the results of the clustering measures}
#' @export

  plot_graph <- function(Laplacian, true_labels, names) {
  # number of nodes
  r <- nrow(Laplacian)


  w <- spectralGraphTopology:::Linv(Laplacian)
  adjacency <- spectralGraphTopology:::A(w)

  # number of components
  q = r - Matrix::rankMatrix(Laplacian)

  graph_net <- igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected", weighted = TRUE)


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
