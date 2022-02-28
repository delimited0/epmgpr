get_polytope_constraints <- function(lb, ub, mu, L) {
  
  d <- length(lb)
  
  inf_idx <- c(is.infinite(lb), is.infinite(ub))
  
  A <- rbind(diag(d), -diag(d))[!inf_idx, ]
  A <- A %*% L
  
  b <- c(mu - lb, -mu + ub)[!inf_idx]
  
  return(list(A = A, b = b))
}

rtmvn <- function(n, mu, Sigma, lb, ub, method = "epess", initial = NULL,
                  N = 1, J = 1, verbose = FALSE) {
  
  if (is.null(initial)) {
    initial <- ifelse(is.finite(lb), lb + 1e-12, ifelse(is.finite(ub), ub - 1e-12, 0))
  }
  
  L <- t(chol(Sigma))
  
  constraints <- get_polytope_constraints(lb, ub, mu, L)
  
  rpolytmvn(n, mu, Sigma, constraints$A, b = constraints$b, lb, ub, 
            method, initial, L, N, J, verbose)
}


#'
#' if b provided, must have b = c(-lb, ub)
rpolytmvn <- function(n, mu, Sigma, A, b = NULL, lb, ub, method = "epess", 
                      initial = NULL, L = NULL, N, J, verbose = FALSE) {
  
  if (is.null(L)) {
    L <- t(chol(Sigma))
  }
  
  if (is.null(b)) {
    b <- c(-lb, ub)
  }
  
  d <- length(mu)
  moments <- epmgp::moments2(mu, Sigma, lb - mu, ub - mu, L)
  
  if (method == "epess") {
    std_samples <- t(sample_epess(n, moments$mu, chol(moments$Sigma), A, b, 
                                  J, N, initial, verbose))
  }
  
  samples <- std_samples %*% t(L) + matrix(rep(mu + moments$mu, n), 
                                           nrow = n, byrow = TRUE)
  samples[1, ] <- initial
  
  return(samples)
}