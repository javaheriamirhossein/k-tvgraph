set.seed(42)
library(spectralGraphTopology)

test_that("test learn_kcomp_heavytail_TV_graph_online", {
  w1 <- c(1, 1, 1, 1, 1, 1)/3
  w2 <- c(1, 1, 1, 1, 1, 1)/3
  true_clust_num <- 2
  Laplacian <- block_diag(L(w1), L(w2))
  p <- ncol(Laplacian)
  nu <- 3
  X <- mvtnorm::rmvt(n = p * 500, delta = rep(0, p), sigma = ((nu-2)/nu) * MASS::ginv(Laplacian), df = nu)
  res <- learn_kcomp_heavytail_TV_graph_online(X, k = 2, rho = 1, heavy_type = "student", nu = nu, maxiter = 100)
  laplacian <- res$laplacian
  inferred_clust_num <- sum(spectralGraphTopology:::eigval_sym(laplacian) < 1e-10)
  expect_true(true_clust_num == inferred_clust_num)
})
