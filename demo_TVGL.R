library(fitHeavyTail)
library(xts)
library(quantmod)
library(igraph)
library(bipartite)
library(readr)
library(spectralGraphTopology)
library(fingraph)
library(combinat)
library(CVXR)
library(ggplot2)
library(reshape2)
library(PerformanceAnalytics)


set.seed(42)


pdf(file = "Plots_TVGL_results.pdf")


# load SP500 stock prices_test into an xts table
stock_prices_orig <- readRDS("examples/stocks/sp500-data-2016-2020.rds")
stock_prices <- stock_prices_orig[1:1001,1:100]
# number of sectors
q <- 8


winLen <-  200
Nday <- nrow(stock_prices) 

# number of stocks
r <- ncol(stock_prices) 
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

market_return <- Ad(getSymbols("^GSPC",
                           from = index(stock_prices[1]), to = index(stock_prices[nrow(stock_prices)]),
                           auto.assign = FALSE,
                           verbose = FALSE
))
market_log_returns <- diff(log(market_return), na.pad = FALSE)

# fit an univariate Student-t distribution to the log-returns of the market
# to obtain an estimate of the degrees of freedom (nu)
nu <- fit_mvt(market_log_returns, nu = "MLE-diag-resampled")$nu



# # build network
SP500 <- read_csv("examples/stocks/SP500-sectors.csv")
stock_sectors <- SP500$GICS.Sector[SP500$Symbol %in% colnames(stock_prices)[1:r]]
stock_sectors_return <- as.numeric(as.factor(stock_sectors))



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

components_graph <- vector("list", Nwin)

accuracy_vec <- rep(0,Nwin)
accuracy_adj_vec <- rep(0,Nwin)
purity_vec <- rep(0,Nwin)
mod_g_vec <- rep(0,Nwin)
mod_gt_vec <- rep(0,Nwin)
balanced_vec <- rep(0,Nwin)
balanced_norm_vec <- rep(0,Nwin)
GINI_vec <- rep(0,Nwin)
rank_mat <- rep(0,Nwin)

nu <- fit_mvt(scale(log_returns), nu="MLE-diag-resample")$nu

for (i in 1:Nwin){
  data_frame <- log_returns[((i-1)*winLen+1):(i*winLen),]
  graphs_list[[i]] <- learn_kcomp_heavytail_TV_graph_online(scale(data_frame), k = q, heavy_type = "student",
                                                     nu = nu,
                                                     sigma_e = 1,
                                                     w_lagged = w_lagged,
                                                     rho = 1,
                                                     d = 1,
                                                     w0 = w0,
                                                     update_eta = TRUE,
                                                     maxiter = 100,
                                                     verbose = TRUE)


  # # Fingraph --------------------------
  # graphs_list[[i]] <- learn_kcomp_heavytail_graph(scale(data_frame), k = q, heavy_type = "student",
  #                                                           nu = nu,
  #                                                           beta = 1e-8,
  #                                                           rho = 1,
  #                                                           d = 1,
  #                                                           w0 = w0,
  #                                                           update_beta = TRUE,
  #                                                           maxiter = 100,
  #                                                           verbose = TRUE)


  # # SGLA ----------------------------
  # S_cov <- cov(scale(data_frame))
  # graphs_list[[i]] <- learn_k_component_graph(S_cov, k = q, beta =10, verbose = TRUE,
  #                                                                                  fix_beta = FALSE, abstol = 1e-3)

#
#   # CLR -----------------------------
#   graphs_list[[i]] <- cluster_k_component_graph(t(scale(log_returns)), k = q, m = 10, lmd = 1, eigtol = 1e-09,
#                                                edgetol = 1e-06, maxiter = 1000)



  # ------------------------

  w <-  spectralGraphTopology:::Ainv(graphs_list[[i]]$adjacency)
  w_lagged <- w
  w0 <- w



  net_Fin <- graph_from_adjacency_matrix(graphs_list[[i]]$adjacency, mode = "undirected", weighted = TRUE)


  # where do predictions differ from GICS?
  metric <- evaluate_clustering(net_Fin, stock_sectors_return, r, q)


  components_graph[[i]] <- metric$labels_pred

  accuracy_adj_vec[i] <- metric$accuracy_adj
  purity_vec[i] <- metric$purity
  mod_g_vec[i] <- metric$mod_g
  mod_gt_vec[i] <- metric$mod_gt
  balanced_vec[i] <- metric$balanced
  balanced_norm_vec[i] <- metric$balanced_norm
  GINI_vec[i] <- metric$GINI

  rank_mat[i] <- rankMatrix(graphs_list[[i]]$laplacian)[1]

}


adjacency <- graphs_list[[Nwin]]$adjacency


#-------------------------------------
# Portfolio design

w_portf_list <- vector("list", Nwin)
for (j in 1:Nwin){

  data_frame <- log_returns[((j-1)*winLen+1):(j*winLen),]
  res_Student_t <- fit_mvt(scale(data_frame))
  mean_vec <- res_Student_t$mu
  Sigma_cov <- res_Student_t$cov
  gamma_2 <- 0.1
  u <- Variable(r)
  ret <- t(mean_vec) %*% u
  risk <- quad_form(u, Sigma_cov)
  smooth <- quad_form(u, graphs_list[[j]]$laplacian)
  obj <- ret  -gamma_2 * smooth
  constr <- list(u >= 0, sum(u) == 1)
  prob <- Problem(Maximize(obj), constr)
  result <- solve(prob)
  w_portf <- result$getValue(u)
  w_portf[w_portf<1e-3] <- 0
  w_portf_list[[j]]<- t(w_portf)
}


for (j in 1:Nwin) {
  adjacency <- graphs_list[[j]]$adjacency
  net_Fin <- graph_from_adjacency_matrix(adjacency, mode = "undirected", weighted = TRUE)
  
  
  # where do our predictions differ from GICS?
  metric <- evaluate_clustering(net_Fin, stock_sectors_return, r, q) 
  labels_pred_adj <-  metric$labels_pred_adj
  
  
  
  colors <- c("#55efc4", "#ff7675", "#0984e3", "#a29bfe", "#B33771", "#48dbfb", "#FDA7DF", "#C4E538")
  
  
  
  
  # ground truth coloring 
  V(net_Fin)$color <- c(colors[stock_sectors_return])
  V(net_Fin)$type <- c(rep(FALSE, r))
  V(net_Fin)$cluster <- c(stock_sectors_return)
  E(net_Fin)$color <- apply(
    as.data.frame(get.edgelist(net_Fin)), 1,
    function(x) {
      ifelse(V(net_Fin)$cluster[x[1]] == V(net_Fin)$cluster[x[2]],
             colors[V(net_Fin)$cluster[x[1]]], "grey"
      )
    }
  )
  
  mask <- labels_pred_adj != stock_sectors_return
  node_labels <- colnames(stock_prices)[1:r]
  node_labels[!mask] <- NA
  label_colors <- rep("black",r)

  
  # plot network
  plot(net_Fin,
       vertex.size = c(rep(4, r)),
       vertex.label = c(node_labels),
       vertex.label.color = label_colors,
       vertex.label.cex = 0.8, vertex.label.dist = 0.5,
       vertex.frame.color = c(colors[stock_sectors_return]),
       layout = layout_nicely(net_Fin),
       vertex.label.family = "Helvetica", vertex.label.color = "black",
       vertex.shape = c(rep("circle", r)),
       edge.width = 3 * E(net_Fin)$weight
  )
  
  
  # implied clusters
  stock_sectors_return_res <- stock_sectors_return
  stock_sectors_return <- igraph::components(net_Fin)$membership
  V(net_Fin)$color <- c(colors[stock_sectors_return])
  V(net_Fin)$type <- c(rep(FALSE, r))
  V(net_Fin)$cluster <- c(stock_sectors_return)
  E(net_Fin)$color <- apply(
    as.data.frame(get.edgelist(net_Fin)), 1,
    function(x) {
      ifelse(V(net_Fin)$cluster[x[1]] == V(net_Fin)$cluster[x[2]],
             colors[V(net_Fin)$cluster[x[1]]], "grey"
      )
    }
  )
  
  
  node_labels <- colnames(stock_prices)[1:r]
  # plot network
  plot(net_Fin,
       vertex.size = c(rep(4, r)),
       vertex.label = c(node_labels),
       vertex.label.cex = 0.8, vertex.label.dist = 0.5,
       vertex.frame.color = c(colors[stock_sectors_return]),
       layout = layout_nicely(net_Fin),
       vertex.label.family = "Helvetica", vertex.label.color = "black",
       vertex.shape = c(rep("circle", r)),
       edge.width = 3 * E(net_Fin)$weight
  )
  
  
  stock_sectors_return <- stock_sectors_return_res
  
}


print(metric$accuracy)
print(metric$accuracy_adj)
print(metric$purity)
print(metric$mod_g)
print(metric$mod_gt)
print(metric$balanced)
print(metric$balanced_norm)
print(metric$GINI)


metrics <- matrix(rep(0,4*Nwin), Nwin, 4)
metrics[,1] <- accuracy_adj_vec
metrics[,2] <- purity_vec
metrics[,3] <- mod_gt_vec
metrics[,4] <- balanced_norm_vec

matplot(metrics, type = "b",pch=2,col = 1:4, ylab = "Metrics")
names <- c("Accuracy", "Purity", "Modularity", "Balancedness")
legend("bottomleft", inset=0.01, legend=names, col=c(1:4),pch=15:18,
       bg= ("white"), horiz=F)


Nwin <- Nwin-1
Nday <- winLen*Nwin
prices_test <- stock_prices [(winLen+1):(winLen+Nday+1),1:r]
returns_test <- (prices_test/lag(prices_test) - 1)[-1]

mine_returns <- matrix(rep(0,4*Nday),4,Nday)

for (i in 1:Nwin) {
 for (j in ((i-1)*winLen+1):((i)*winLen) ) {
   if (j<=Nday) { 
     if (i<=1) {k1 <- i}
     mine_returns[1, j] <-  sum(returns_test[j, ] * w_portf_list[[k1]] )
     if (i<=2) {k2 <- i}
     mine_returns[2, j] <-  sum(returns_test[j, ] * w_portf_list[[k2]] )
     if (i<=3) {k3 <- i}
     mine_returns[3, j] <-  sum(returns_test[j, ] * w_portf_list[[k3]] )
     if (i<=4) {k4 <- i}
     mine_returns[4, j] <-  sum(returns_test[j, ] * w_portf_list[[k4]] )
     
   }
 }
}


all_returns <- matrix(rep(0,6*Nday),6,Nday)
for (i in 1:4) {
  all_returns[i,] <- cumsum(mine_returns[i,])
}

EWP_return <- rowSums(returns_test)/r
all_returns[5,] <- cumsum(EWP_return)

market_return <- market_log_returns[(winLen+1):(winLen+Nday)]
all_returns[6,] <- cumsum(market_return)

matplot(t(all_returns), type = "l",pch=1,col = 1:6, ylab = "Cum return", xlab= )
names <- c( "MGSP-1", "MGSP-2", "MGSP-3", "MGSP-4", "EWP", "Market")
legend("bottomleft", inset=0.01, legend=names, col=c(1:6),pch=15:(15+6-1),
       bg= ("white"), horiz=F)


all_returns <- all_returns + 1
df <- data.frame(index   = index(market_return),
                 MGSP_1 = all_returns[1,],
                 MGSP_2 = all_returns[2,],
                 MGSP_3 = all_returns[3,],
                 MGSP_4 = all_returns[4,],
                 EWP =  all_returns[5,],
                 Market = all_returns[6,]
)
#
names <- c( "MGSP_1", "MGSP_2", "MGSP_3", "MGSP_4","EWP", "Market")

molten_df <- melt(df, id.vars = "index", measure.vars = names)
ggplot(molten_df, aes(x = index, y = value, col = variable)) +
  geom_line() +
  ggtitle("Portfolio Returns") + labs(colour = "Portfolios") + ylab("Cumm. Return") + xlab("Date")

ggsave("ggplot_portfolios.pdf", width = 10, height = 5)





#-----------------------
# Evaluate the performance of the portfolios

# create empty portfolio matrix
ret_list <- vector(mode="list", length = Nwin-1)

for (i in 1:(Nwin-1)) {
  
}
rebal_indices <-  winLen*(1:(Nwin)-1) +1

for (i in 1:(Nwin)) {
  w_MVP_rolling <- returns_test
  w_MVP_rolling[] <- NA
  
  for (j in 1:(Nwin)) {
    if (j<(i+1)) {k1 <-j}
    w_MVP_rolling[rebal_indices[j], ] <- w_portf_list[[k1]]
  }
  w_MVP_rolling <- na.omit(w_MVP_rolling)
  tmp <- Return.portfolio(returns_test, weights = w_MVP_rolling, verbose = TRUE, rebalance_on = "days")
  ret_list[[i]] <- tmp$returns
}

w_EWP <- rep(1/r,r)
ret_EWP_static <- xts(returns_test %*% w_EWP, index(returns_test))
ret_EWP_static <- Return.portfolio(returns_test, weights = w_EWP, rebalance_on = "days")


#--------------
w_MVP_list <- vector("list", Nwin)
for (j in 1:Nwin){
  data_stocks <- log_returns[((j-1)*winLen+1):(j*winLen),]
  res_Student_t <- fit_mvt(scale(data_stocks))
  mean_vec <- res_Student_t$mu
  Sigma_cov <- res_Student_t$cov
  gamma_1 <- 0.1
  u <- Variable(r)
  ret <- t(mean_vec) %*% u
  risk <- quad_form(u, Sigma_cov)
  smooth <- quad_form(u, graphs_list[[j]]$laplacian)
  obj <- ret - gamma_1  * risk
  constr <- list(u >= 0, sum(u) == 1)
  prob <- Problem(Maximize(obj), constr)
  result <- solve(prob)
  w_portf <- result$getValue(u)
  w_portf[w_portf<1e-3] <- 0
  w_MVP_list[[j]]<- t(w_portf)
}


# create empty portfolio matrix
ret_list_MVP <- vector(mode="list", length = Nwin-1)
rebal_indices <-  winLen*(1:(Nwin)-1) +1
for (i in 1:(Nwin)) {
  w_MVP_rolling <- returns_test
  w_MVP_rolling[] <- NA
  for (j in 1:(Nwin)) {
    if (j<(i+1)) {k1 <-j}
    w_MVP_rolling[rebal_indices[j], ] <- w_MVP_list[[k1]]
  }
  w_MVP_rolling <- na.omit(w_MVP_rolling)
  tmp <- Return.portfolio(returns_test, weights = w_MVP_rolling, verbose = TRUE, rebalance_on = "days")
  ret_list_MVP[[i]] <- tmp$returns
}

ret_all <- cbind(ret_list[[1]],
                 ret_list[[2]],
                 ret_list[[3]],
                 ret_list[[4]],
                 ret_list_MVP[[4]],
                 ret_EWP_static)


colnames(ret_all) <- c("MGSP-1", "MGSP-2", "MGSP-3", "MGSP-4", "MVP Dynamic", "EWP")

{chart.CumReturns(ret_all, main = "Commulative return of portfolios",
                  geometric = FALSE, wealth.index = TRUE, legend.loc = "topleft", colorset = rich6equal)
  addEventLines(xts("", index(returns_test[winLen])), srt=90, pos=1, lwd = 2, col = "darkblue")
  addEventLines(xts("", index(returns_test[2*winLen])), srt=90, pos=1, lwd = 2, col = "darkred")
  addEventLines(xts("", index(returns_test[3*winLen])), srt=90, pos=1, lwd = 2, col = "darkgreen")
}


{chart.Drawdown(ret_all, main = "Drawdown of portfolios", 
                legend.loc = "bottomleft", colorset = rich10equal)
  addEventLines(xts("", index(returns_test[winLen])), srt=90, pos=1, lwd = 2, col = "darkblue")
  addEventLines(xts("", index(returns_test[2*winLen])), srt=90, pos=1, lwd = 2, col = "darkred")
  addEventLines(xts("", index(returns_test[3*winLen])), srt=90, pos=1, lwd = 2, col = "darkgreen")
}
print(round(t(table.AnnualizedReturns(ret_all)),2))




dev.off()
