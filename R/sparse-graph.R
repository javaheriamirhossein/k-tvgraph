library(spectralGraphTopology)

#' @title Learns sparse Laplacian matrix of a connected graph
#'
#' Learns a connected graph via non-convex, sparse promoting regularization
#' functions such as MCP, SCAD, and re-weighted l1-norm.
#'
#' @param S a pxp sample covariance/correlation matrix, where p is the number
#'        of nodes of the graph
#' @param w0 initial estimate for the weight vector the graph.
#' @param alpha hyperparameter to control the level of sparsiness of the
#'        estimated graph
#' @param sparsity_type type of non-convex sparsity regularization. Available
#'        methods are: "mcp", "scad", "re-l1", and "none"
#' @param eps hyperparameter for the re-weighted l1-norm
#' @param eta learning rate
#' @param backtrack whether to update the learning rate using backtrack line
#'        search
#' @param maxiter maximum number of iterations
#' @param reltol relative tolerance on the Frobenius norm of the estimated
#'        Laplacian matrix as a stopping criteria
#' @param verbose whether or not to show a progress bar displaying the
#'        iterations
#' @return A list containing possibly the following elements:
#' \item{\code{laplacian}}{the estimated Laplacian Matrix}
#' \item{\code{adjacency}}{the estimated Adjacency Matrix}
#' \item{\code{maxiter}}{number of iterations taken to converge}
#' \item{\code{convergence}}{boolean flag to indicate whether or not the optimization converged}
#' \item{\code{elapsed_time}}{elapsed time recorded at every iteration}
#' @author Ze Vinicius, Jiaxi Ying, and Daniel Palomar
#' @export
#' @import spectralGraphTopology
learn_laplacian_pgd_connected <- function(X, w0 = NULL, alpha = 0, nu = 2, heavy_type = "student", sparsity_type = "none",
                                          eps = 1e-4, gamma = 2.001, q = 1, backtrack = TRUE,
                                          maxiter = 10000, reltol = 1e-5, verbose = TRUE) {
  # number of nodes
  S <- cov(scale(X))
  p <- nrow(S)
  n <- nrow(X)
  J <- matrix(1, p, p) / p
  
  Sinv <- MASS::ginv(S)
  

  
  
	# use safe initial learning rate
	eta <- 1 / (2*p)
  if (is.null(w0)) {
    w <- spectralGraphTopology:::Linv(Sinv)
    w[w < 0] <- 0
    w <- compute_initial_point(w, Sinv, eta)
  } else {
    w <- w0
    if(length(w) != .5 * p * (p-1)) stop(paste("dimension of the initial point must be ", .5 * p * (p-1)))
    if(any(w < 0)) stop("initial point must be nonnegative.")
    chol_status <- try(chol_factor <- chol(L(w) + J), silent = TRUE)
    chol_error <- ifelse(class(chol_status) == "try-error", TRUE, FALSE)
    if(chol_error[1]) stop("initial point provided must be a connected graph.")
  }
  Lw <- L(w)
  
  
  
  if (heavy_type == "student") {
    Sq <- vector(mode = "list", length = n)
    for (i in 1:n)
      Sq[[i]] <- t(X[i, ]) %*%  X[i, ] 
    
    S <- rep(0, p, p)
    for (q in 1:n)
      S <- S +   Sq[[q]] / ( sum(w * Lstar( Sq[[q]] )) + nu )
    
    S <- S * (p + nu)/n* nu
  } 
  # else if (heavy_type == "gaussian") {
  #   S <- cov(X)
  # }

  
  
  if (sparsity_type == "mcp") {
    H <- -(alpha + Lw / gamma) * (Lw >= -alpha*gamma)
    diag(H) <- rep(0, p)
    K <- S + H
  } else if (sparsity_type == "re-l1") {
    H <- alpha * (diag(p) - p * J)
    K <- S + H / ((-Lw + eps) ^ q)
  } else if (sparsity_type == "scad") {
    H <- -alpha * (Lw >= - alpha)
    H <- H + (-gamma * alpha - Lw) / (gamma - 1) * (Lw > -gamma*alpha) * (Lw < -alpha)
    diag(H) <- rep(0, p)
    K <- S + H
  }
  K <- S + H
  if (verbose)
    pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta",
                                     total = maxiter, clear = FALSE, width = 80)
  time_seq <- c(0)
  relerror_seq <- c()
  eta_seq <- c()
  start_time <- proc.time()[3]
  for (i in 1:maxiter) {
    w_prev <- w
    tryCatch(
      {gradient <- spectralGraphTopology:::Lstar(K - spectralGraphTopology:::inv_sympd(Lw + J))},
        error = function(err) {
          results <- list(laplacian = L(w_prev),
                          adjacency = A(w_prev),
                          maxiter = i,
                          convergence = FALSE,
                          elapsed_time = time_seq,
                          optimization_error = relerror_seq,
                          learning_rates = eta_seq)
          return(results)
      }
    )
    if (backtrack) {
      fun <- mle_pgd.obj(Lw = Lw, J = J, K = K)$obj_fun
      while(1) {
        wi <- w - eta * gradient
        wi[wi < 0] <- 0
        # compute the objective function at the updated value of w
        Lwi <- L(wi)
        aux <- mle_pgd.obj(Lw = Lwi, J = J, K = K)
        fun_t <- aux$obj_fun
        is_disconnected <- aux$is_disconnected
        # check whether the previous value of the objective function is
        # smaller than the current one
        if ((fun < fun_t - sum(gradient * (wi - w)) - (.5/eta)*norm(wi - w, '2')^2)
            | is_disconnected) {
          eta <- .5 * eta
        } else {
          eta <- 2 * eta
          break
        }
      }
    } else {
        wi <- w - eta * gradient
        wi[wi < 0] <- 0
    }
    if (verbose)
      pb$tick()
    relerror <- norm(L(wi) - Lw, 'F') / norm(Lw, 'F')
    has_converged <- (relerror < reltol) && (i > 1)
    relerror_seq <- c(relerror_seq, relerror)
    eta_seq <- c(eta_seq, eta)
    time_seq <- c(time_seq, proc.time()[3] - start_time)
    if (has_converged)
      break
    w <- wi
    Lw <- L(w)
    if (sparsity_type == "mcp") {
      H <- -(alpha + Lw / gamma) * (Lw >= -alpha*gamma)
      diag(H) <- rep(0, p)
      K <- S + H
    } else if (sparsity_type == "re-l1") {
      K <- S + H / ((-Lw + eps) ^ q)
    } else if (sparsity_type == "scad") {
      H <- -alpha * (Lw >= - alpha)
      H <- H + (-gamma * alpha - Lw) / (gamma - 1) * (Lw > -gamma*alpha) * (Lw < -alpha)
      diag(H) <- rep(0, p)
      K <- S + H
    }
  }
  results <- list(laplacian = L(wi),
                  adjacency = A(wi),
#                  obj_fun = nonconvex.obj(L(wi), J, S, alpha, gamma, sparsity_type),
                  maxiter = i,
                  convergence = has_converged,
                  elapsed_time = time_seq,
                  optimization_error = relerror_seq,
                  learning_rates = eta_seq)
  return(results)
}


mle_pgd.obj <- function(Lw, J, K) {
  chol_status <- try(chol_factor <- chol(Lw + J), silent = TRUE)
  chol_error <- ifelse(class(chol_status) == "try-error", TRUE, FALSE)
  if (chol_error[1]) {
    return(list(obj_fun = 1e16, is_disconnected = TRUE))
  } else {
     return(list(obj_fun = sum(Lw*K) - 2*sum(log(diag(chol_factor))), is_disconnected = FALSE))
  }
}


#nonconvex.obj <- function(Lw, J, S, alpha, gamma, sparsity_type) {
#  chol_factor <- chol(Lw + J)
#  if (sparsity_type == "mcp") {
#    mask <- (Lw >= -alpha*gamma)
#    H <- (-alpha * Lw  - .5 * Lw^2 / gamma) * mask + (.5 * gamma * alpha ^ 2) * (!mask)
#    diag(H) <- rep(0, nrow(H))
#  }
#  return(sum(Lw * S) - 2*sum(log(diag(chol_factor))) + sum(H))
#}

compute_initial_point <- function(w, Sinv, eta) {
  for (i in c(1:100)) {
    grad <- spectralGraphTopology:::Lstar(L(w) - Sinv)
    wi <- w - eta * grad
    wi[wi < 0] <- 0
    if (norm(w - wi, '2') / norm(w, '2') < 1e-4)
      break
    w <- wi
  }
  return(w + 1e-4)
}
