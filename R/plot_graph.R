#' @title Plot a clustering graph from Laplacian matrix
#'
#' Plots a clustering graph based on a given Laplacian matrix and the ground-truth labels for the clusters.
#' @param Laplacian p by p Laplacian matrix.
#' @param true_labels the numeric vector of true labels (of length p)
#' @param node_names the vector of names of the nodes, e.g., stocks (of length p)
#' @param implied_clusters whether to color the nodes based on implied clusters or not
#' @param Coords input graph node coordinates
#' @return The network object:
#' \item{\code{graph_net}}{the graph network object that can be used for plotting}
#' \item{\code{metric}}{the results of the clustering measures}
#' \item{\code{Coords}}{output graph node coordinates}
#' @export

  plot_graph <- function(Laplacian, true_labels, node_names = NULL, 
                         implied_clusters = FALSE,
                         verbose = FALSE,
                         Coords=NULL) {
  
  # number of nodes
  p <- nrow(Laplacian)


  w <- spectralGraphTopology:::Linv(Laplacian)
  adjacency <- spectralGraphTopology:::A(w)

  # number of components
  q = p - Matrix::rankMatrix(Laplacian)

  graph_net <- igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected", weighted = TRUE)


  # where do our predictions differ from true labels?
  metric <- evaluate_clustering(graph_net, true_labels, p, q)
  labels_pred_adj <-  metric$labels_pred_adj


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
  
  mask <- labels_pred_adj != true_labels
  node_labels[!mask] <- NA
  
  if (!is.null(Coords)){
    layout <- Coords
  }
  else {
    layout <- layout_nicely(net_Fin)
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
         edge.width = 3 * E(graph_net)$weight
         )
  }
  
  else {
    # implied clusters
    true_labels <- igraph::components(graph_net)$membership
    
    true_labels_reordered <- rep(0, r)
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
         edge.width = 3 * E(graph_net)$weight
        )

  }

  
  if (verbose) {
    print("metric:\n")
    print(metric$accuracy_adj)
    print(metric$purity)
    print(metric$mod_gt)
    # print(metric$balanced_norm)
    # print(metric$GINI)
    print(metric$ARI)
   }

  results <- list(graph_net = graph_net, metric = metric, Coords = Coords)

  return(results)
}
