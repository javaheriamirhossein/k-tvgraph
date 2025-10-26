library(fitHeavyTail)
library(xts)
library(quantmod)
library(igraph)
library(readr)
library(spectralGraphTopology)
library(combinat)
library(tvgraph)


set.seed(42)


pdf(file = "Plots_TVGL_results.pdf")


# number of stocks
p <- 100

# number of sectors (clusters)
q <- 8


# load SP500 stock prices_test into an xts table
stock_prices_orig <- readRDS("examples/stocks/sp500-data-2016-2020.rds")
stock_prices <- stock_prices_orig[1:1001,1:p]

# frame (window) length
winLen <-  200

Nday <- nrow(stock_prices)

Nwin <- Nday%/% winLen





# total nodes in the graph
colnames(stock_prices)[1:p]
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





# Load true labels
SP500 <- read_csv("examples/stocks/SP500-sectors.csv")
stock_sectors <- SP500$GICS.Sector[SP500$Symbol %in% colnames(stock_prices)[1:p]]
stock_sectors_index <- as.numeric(as.factor(stock_sectors))



#----------------------------
## Online TV graph learning (proposed)
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


accuracy_vec <- rep(0,Nwin)
purity_vec <- rep(0,Nwin)
modularity_vec <- rep(0,Nwin)
balanced_vec <- rep(0,Nwin)
ARI_vec <- rep(0,Nwin)
GINI_vec <- rep(0,Nwin)
rank_mat <- rep(0,Nwin)


for (i in 1:Nwin){
  data_frame <- log_returns[((i-1)*winLen+1):(i*winLen),]
  nu <- fit_mvt(data_frame, nu = "MLE-diag-resampled")$nu
  graphs_list[[i]] <- learn_kcomp_heavytail_TV_graph_online(scale(data_frame), k = q, heavy_type = "student",
                                                     nu = nu,
                                                     sigma_e = exp(0.1),
                                                     w_lagged = w_lagged,
                                                     rho = 3,
                                                     d = 1,
                                                     w0 = w0,
                                                     update_eta = TRUE,
                                                     maxiter = 40,
                                                     verbose = TRUE)





  # ------------------------

  w <-  spectralGraphTopology:::Ainv(graphs_list[[i]]$adjacency)
  w_lagged <- w
  w0 <- w



  graph_net <- graph_from_adjacency_matrix(graphs_list[[i]]$adjacency, mode = "undirected", weighted = TRUE)


  # where do predictions differ from GICS?
  metric <- evaluate_clustering(graph_net, stock_sectors_index, p, q)



  accuracy_vec[i] <- metric$accuracy
  purity_vec[i] <- metric$purity
  modularity_vec[i] <- metric$mod
  balanced_vec[i] <- metric$balanced
  GINI_vec[i] <- metric$GINI
  ARI_vec[i] <- metric$ARI
  rank_mat[i] <- Matrix::rankMatrix(graphs_list[[i]]$laplacian)[1]

}




Coords <- NULL
Ngraph <- length(graphs_list)


implied_clusters <- TRUE
for (i in 1:Ngraph){
  id <- Ngraph - i + 1
  gplt <- plot_graph(graphs_list[[id]]$adjacency, stock_sectors_index, colnames(stock_prices), Coords = Coords, 
                                 implied_clusters = implied_clusters )
  
  title(main =  paste("frame",id) )
  Coords <- gplt$Coords
}

# print the final values of the clustering measures
cat("ACC: ", metric$accuracy, "\n")
cat("PUR: ", metric$purity, "\n")
cat("MOD: ", metric$mod, "\n")
cat("ARI: ", metric$ARI, "\n")


# Plot the clustering measures vs frames number
metrics <- matrix(rep(0,4*Nwin), Nwin, 4)
metrics[,1] <- accuracy_vec
metrics[,2] <- purity_vec
metrics[,3] <- modularity_vec
metrics[,4] <- ARI_vec

matplot(metrics, type = "b",pch=2,col = 1:4, ylab = "Metrics", xlab = "Frame")
names <- c("Accuracy", "Purity", "Modularity", "ARI")
legend("bottomleft", inset=0.01, legend=names, col=c(1:4),pch=15:18,
       bg= ("white"), horiz=F)



dev.off()


