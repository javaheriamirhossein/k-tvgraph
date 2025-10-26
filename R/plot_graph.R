#' @title Plot a clustering graph from adjacency matrix
#'
#' Plots a clustering graph based on a given adjacency matrix and the ground-truth labels for the clusters.
#' @param adjacency p by p adjacency matrix.
#' @param true_labels the numeric vector of true labels (of length p)
#' @param node_names the vector of names of the nodes, e.g., stocks (of length p)
#' @param implied_clusters whether to color the nodes based on implied clusters or not
#' @param Coords input graph node coordinates
#' @return The network object:
#' \item{\code{graph_net}}{the graph network object that can be used for plotting}
#' \item{\code{metric}}{the results of the clustering measures}
#' \item{\code{Coords}}{output graph node coordinates}
#' @export

plot_graph <- function(adjacency, true_labels, node_names = NULL, 
                       implied_clusters = FALSE,
                       verbose = FALSE,
                       Coords=NULL) {
  
  # number of nodes
  p <- nrow(adjacency)
  
  
  w <- spectralGraphTopology:::Ainv(adjacency)
  Laplacian <- spectralGraphTopology:::L(w)
  
  
  # number of components
  q = p - Matrix::rankMatrix(Laplacian)
  
  graph_net <- igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected", weighted = TRUE)
  
  
  # where do our predictions differ from true labels?
  clustering_res <- evaluate_clustering(graph_net, true_labels, p, q)
  labels_pred <-  clustering_res$labels_pred
  
  
  # if q>10 add some colors to the list
  colors <- c("#55efc4", "#ff7675", "#0984e3", "#a29bfe", "#B33771", "#48dbfb", "#FDA7DF", "#C4E538", "#0184a3", "#f9dcfb")
  colors <- colors[1:q]
  
  
  
  # ground truth coloring
  V(graph_net)$color <- c(colors[true_labels])
  V(graph_net)$type <- c(rep(FALSE, p))
  V(graph_net)$cluster <- c(true_labels)
  E(graph_net)$color <- apply(
    as.data.frame(get.edgelist(graph_net)), 1,
    function(x) {
      ifelse(V(graph_net)$cluster[x[1]] == V(graph_net)$cluster[x[2]],
             colors[V(graph_net)$cluster[x[1]]], "grey"
      )
    }
  )
  
  node_labels <- rep(NA, p)
  if (!is.null(node_names)){
    node_labels <- node_names
  }
  
  mask <- labels_pred != true_labels
  node_labels[!mask] <- NA
  
  if (!is.null(Coords)){
    layout <- Coords
  }
  else {
    layout <- layout_nicely(graph_net)
    Coords <- layout
  }
  
  label_colors <- rep("black",p)
  label_colors[mask] <- "red"
  
  
  
  if (!implied_clusters) {
    
    # plot network
    plot(graph_net,
         vertex.size = c(rep(4, p)),
         vertex.label = c(node_labels),
         vertex.label.color = label_colors,
         vertex.label.cex = 0.8, vertex.label.dist = 0.5,
         vertex.frame.color = c(colors[true_labels]),
         layout = layout,
         vertex.label.family = "Helvetica", vertex.label.color = "black",
         vertex.shape = c(rep("circle", p)),
         edge.width = 8 * E(graph_net)$weight
    )
  }
  
  else {
    # implied clusters
    true_labels <- igraph::components(graph_net)$membership
    
    true_labels_reordered <- rep(0, p)
    for (j in 1:q){
      idx <- true_labels %in% c(j)
      true_labels_reordered[idx] <- q-j+1
    }
    true_labels <- true_labels_reordered
    
    V(graph_net)$color <- c(colors[true_labels])
    V(graph_net)$type <- c(rep(FALSE, p))
    V(graph_net)$cluster <- c(true_labels)
    E(graph_net)$color <- apply(
      as.data.frame(get.edgelist(graph_net)), 1,
      function(x) {
        ifelse(V(graph_net)$cluster[x[1]] == V(graph_net)$cluster[x[2]],
               colors[V(graph_net)$cluster[x[1]]], "grey"
        )
      }
    )
    
    
    # plot network
    plot(graph_net,
         vertex.size = c(rep(4, p)),
         vertex.label = c(node_labels),
         vertex.label.cex = 0.8, vertex.label.dist = 0.5,
         vertex.frame.color = c(colors[true_labels]),
         layout = layout,
         vertex.label.family = "Helvetica", vertex.label.color = "black",
         vertex.shape = c(rep("circle", p)),
         edge.width = 8 * E(graph_net)$weight
    )
    
  }
  
  
  metric <- list(accuracy = clustering_res$accuracy, 
                 purity = clustering_res$purity, 
                 modularity = clustering_res$mod, 
                 balancedness = clustering_res$balanced, 
                 GINI = clustering_res$GINI,
                 ARI = clustering_res$ARI,
                 NMI = clustering_res$NMI)
  
  if (verbose) {
    cat("accuracy:", metric$accuracy,"\n")
    cat("purity:", metric$purity,"\n")
    cat("modularity:", metric$modularity,"\n")
    cat("balancedness:", metric$balancedness,"\n")
    cat("GINI:", metric$GINI,"\n")
    cat("ARI:", metric$ARI,"\n")
    cat("NMI:", metric$NMI,"\n")
  }
  
  results <- list(graph_net = graph_net, metric = metric, Coords = Coords)
  
  return(results)
}