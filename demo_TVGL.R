library(fitHeavyTail)
library(xts)
library(quantmod)
library(igraph)
library(readr)
library(spectralGraphTopology)
library(combinat)
library(ggplot2)
library(reshape2)
# library(ktvgraph)


set.seed(42)


pdf(file = "Plots_TVGL_results.pdf")


# number of stocks
r <- 100

# number of sectors
q <- 8


# load SP500 stock prices_test into an xts table
stock_prices_orig <- readRDS("examples/stocks/sp500-data-2016-2020.rds")
stock_prices <- stock_prices_orig[1:1001,1:r]

winLen <-  200
Nday <- nrow(stock_prices)

Nwin <- Nday%/% winLen





# total nodes in the graph
colnames(stock_prices)[1:r]
#>   [1] "A"    "AAL"  "ABBV" "ABC"  "ABMD" "ABT"  "ADM"  "AEE"  "AEP"  "AES"
#>  [11] "AFL"  "AIG"  "AIV"  "AIZ"  "AJG"  "ALB"  "ALGN" "ALK"  "ALL"  "ALLE"
#>  [21] "ALXN" "AMCR" "AME"  "AMGN" "AMP"  "AMT"  "ANTM" "AON"  "AOS"  "APA"
#>  [31] "APD"  "ARE"  "ATO"  "AVB"  "AVY"  "AWK"  "AXP"  "BA"   "BAC"  "BAX"
#>  [41] "BDX"  "BEN"  "BIIB" "BIO"  "BK"   "BKR"  "BLK"  "BLL"  "BMY"  "BSX"
#>  [51] "BXP"  "C"    "CAG"  "CAH"  "CAT"  "CB"   "CBOE" "CBRE" "CCI"  "CE"
#>  [61] "CERN" "CF"   "CFG"  "CHD"  "CHRW" "CI"   "CINF" "CL"   "CLX"  "CMA"
#>  [71] "CME"  "CMI"  "CMS"  "CNC"  "CNP"  "COF"  "COG"  "COO"  "COP"  "COST"
#>  [81] "COTY" "CPB"  "CPRT" "CSX"  "CTAS" "CVS"  "CVX"  "D"    "DAL"  "DD"
#>  [91] "DE"   "DFS"  "DGX"  "DHR"  "DLR"  "DOV"  "DRE"  "DTE"  "DUK"  "DVA"


# compute log-returns
log_returns <- diff(log(stock_prices), na.pad = FALSE)





# # build network
SP500 <- read_csv("examples/stocks/SP500-sectors.csv")
stock_sectors <- SP500$GICS.Sector[SP500$Symbol %in% colnames(stock_prices)[1:r]]
stock_sectors_index <- as.numeric(as.factor(stock_sectors))



#----------------------------
## Online TV graph learning (proposed)
forget_fac <- 0.1
data_frame <- log_returns[1:winLen,]
S_cov <- cor(scale(data_frame))
w <- spectralGraphTopology:::w_init('naive', MASS::ginv(S_cov))
w0 <- w
A0 <- A(w)
A0 <- A0 / rowSums(A0)
w0 <- spectralGraphTopology:::Ainv(A0)
w0 = w0/sum(w0)

w_lagged <- w0

w_lagged <- 0

graphs_list <- vector("list", Nwin)


accuracy_adj_vec <- rep(0,Nwin)
purity_vec <- rep(0,Nwin)
mod_gt_vec <- rep(0,Nwin)
balanced_norm_vec <- rep(0,Nwin)
ARI_vec <- rep(0,Nwin)
GINI_vec <- rep(0,Nwin)
rank_mat <- rep(0,Nwin)


for (i in 1:Nwin){
  data_frame <- log_returns[((i-1)*winLen+1):(i*winLen),]
  nu <- fit_mvt(data_frame, nu = "MLE-diag-resampled")$nu
  graphs_list[[i]] <- learn_kcomp_heavytail_TV_graph_online(scale(data_frame), k = q, heavy_type = "student",
                                                     nu = nu,
                                                     sigma_e = exp(3),
                                                     w_lagged = w_lagged,
                                                     rho = 1,
                                                     d = 1,
                                                     w0 = w0,
                                                     update_eta = TRUE,
                                                     maxiter = 40,
                                                     verbose = TRUE)





  # ------------------------

  w <-  spectralGraphTopology:::Ainv(graphs_list[[i]]$adjacency)
  w_lagged <- w
  w0 <- w



  net_Fin <- graph_from_adjacency_matrix(graphs_list[[i]]$adjacency, mode = "undirected", weighted = TRUE)


  # where do predictions differ from GICS?
  metric <- evaluate_clustering(net_Fin, stock_sectors_index, r, q)



  accuracy_adj_vec[i] <- metric$accuracy_adj
  purity_vec[i] <- metric$purity
  mod_gt_vec[i] <- metric$mod_gt
  balanced_norm_vec[i] <- metric$balanced_norm
  GINI_vec[i] <- metric$GINI
  ARI_vec[i] <- metric$ARI
  rank_mat[i] <- rankMatrix(graphs_list[[i]]$laplacian)[1]

}


adjacency <- graphs_list[[Nwin]]$adjacency



for (j in 1:Nwin) {
  # adjacency <- graphs_list[[j]]$adjacency
  # net_Fin <- graph_from_adjacency_matrix(adjacency, mode = "undirected", weighted = TRUE)
 
  net_result <- plot_graph(graphs_list[[j]]$laplacian, stock_sectors_index, colnames(stock_prices))

  metric <- net_result$metric
  
  
#   # where do our predictions differ from GICS?
#   metric <- evaluate_clustering(net_Fin, stock_sectors_index, r, q)
#   
#   
#   labels_pred_adj <-  metric$labels_pred_adj
# 
# 
#   
#   colors <- c("#55efc4", "#ff7675", "#0984e3", "#a29bfe", "#B33771", "#48dbfb", "#FDA7DF", "#C4E538")
# 
#   plot_graph(graphs_list[[1]]$laplacian, stock_sectors_index)
# 
# 
#   # ground truth coloring
#   V(net_Fin)$color <- c(colors[stock_sectors_index])
#   V(net_Fin)$type <- c(rep(FALSE, r))
#   V(net_Fin)$cluster <- c(stock_sectors_index)
#   E(net_Fin)$color <- apply(
#     as.data.frame(get.edgelist(net_Fin)), 1,
#     function(x) {
#       ifelse(V(net_Fin)$cluster[x[1]] == V(net_Fin)$cluster[x[2]],
#              colors[V(net_Fin)$cluster[x[1]]], "grey"
#       )
#     }
#   )
# 
#   mask <- labels_pred_adj != stock_sectors_index
#   node_labels <- colnames(stock_prices)[1:r]
#   node_labels[!mask] <- NA
#   label_colors <- rep("black",r)
# 
# 
#   # plot network
#   plot(net_Fin,
#        vertex.size = c(rep(4, r)),
#        vertex.label = c(node_labels),
#        vertex.label.color = label_colors,
#        vertex.label.cex = 0.8, vertex.label.dist = 0.5,
#        vertex.frame.color = c(colors[stock_sectors_index]),
#        layout = layout_nicely(net_Fin),
#        vertex.label.family = "Helvetica", vertex.label.color = "black",
#        vertex.shape = c(rep("circle", r)),
#        edge.width = 3 * E(net_Fin)$weight
#   )
# 
# 
#   # implied clusters
#   stock_sectors_index_res <- stock_sectors_index
#   stock_sectors_index <- igraph::components(net_Fin)$membership
#   V(net_Fin)$color <- c(colors[stock_sectors_index])
#   V(net_Fin)$type <- c(rep(FALSE, r))
#   V(net_Fin)$cluster <- c(stock_sectors_index)
#   E(net_Fin)$color <- apply(
#     as.data.frame(get.edgelist(net_Fin)), 1,
#     function(x) {
#       ifelse(V(net_Fin)$cluster[x[1]] == V(net_Fin)$cluster[x[2]],
#              colors[V(net_Fin)$cluster[x[1]]], "grey"
#       )
#     }
#   )
# 
# 
#   node_labels <- colnames(stock_prices)[1:r]
#   # plot network
#   plot(net_Fin,
#        vertex.size = c(rep(4, r)),
#        vertex.label = c(node_labels),
#        vertex.label.cex = 0.8, vertex.label.dist = 0.5,
#        vertex.frame.color = c(colors[stock_sectors_index]),
#        layout = layout_nicely(net_Fin),
#        vertex.label.family = "Helvetica", vertex.label.color = "black",
#        vertex.shape = c(rep("circle", r)),
#        edge.width = 3 * E(net_Fin)$weight
#   )
# 
# 
#   stock_sectors_index <- stock_sectors_index_res
# 
}


# print the final values of the clustering measures
cat("ACC: ", metric$accuracy_adj, "\n")
cat("PUR: ", metric$purity, "\n")
cat("MOD: ", metric$mod_gt, "\n")
cat("ARI: ", metric$ARI, "\n")


# Plot the clustering measures vs frames number
metrics <- matrix(rep(0,4*Nwin), Nwin, 4)
metrics[,1] <- accuracy_adj_vec
metrics[,2] <- purity_vec
metrics[,3] <- mod_gt_vec
metrics[,4] <- ARI_vec

matplot(metrics, type = "b",pch=2,col = 1:4, ylab = "Metrics", xlab = "Frame")
names <- c("Accuracy", "Purity", "Modularity", "ARI")
legend("bottomleft", inset=0.01, legend=names, col=c(1:4),pch=15:18,
       bg= ("white"), horiz=F)



dev.off()
