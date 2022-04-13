get_polytope_constraints <- function(lb, ub, A, mu, L) {
  
  d <- length(lb)
  
  inf_idx <- c(is.infinite(lb), is.infinite(ub))
  
  Amat <- rbind(A, -A)[!inf_idx, ]
  Amat <- Amat %*% L
  
  b <- c(-lb + A %*% mu, ub - A %*% mu)[!inf_idx]
  
  return(list(A = Amat, b = b))
}


#' @param N number of samples to draw uniformly per threshold level
#' @param J number of slices per iteration (recycling)
#' @export
rtmvn <- function(n, mu, Sigma, lb, ub, A = NULL, method = "epess", initial = NULL,
                  N = 1, J = 1, verbose = FALSE, burnin = 0) {
  
  L <- t(chol(Sigma))
  d <- length(mu)
  
  if (is.null(A)) {
    A = diag(d)
  }
  
  if (is.null(initial)) {
    A_inv <- MASS::ginv(A)
    initial <- A_inv %*% (lb + ub) / 2
  }
  
  # standardize the gaussian
  initial = solve(L, initial - mu)
  pc <- get_polytope_constraints(lb, ub, A, mu, L)
  moments <- moments2(rep(0, d), diag(d), lb - A %*% mu, ub - A %*% mu, A %*% L)
  
  total_samples = n + burnin
  std_samples = t( sample_epess(total_samples, moments$mu, chol(moments$Sigma), 
                             pc$A, pc$b,
                             J, N, initial, verbose) )
  
  # discard burn in samples
  std_samples = std_samples[(burnin+1):total_samples, ]
  
  # undo the approximate mean shift transformation
  std_samples = std_samples + matrix(rep(moments$mu, n), nrow = n, byrow = TRUE)
  
  # un-standardize samples
  samples = std_samples %*% t(L) + matrix(rep(mu, n), nrow = n, byrow = TRUE)
  
  return(samples)
  # rpolytmvn(n, mu, Sigma, constraints$A, b = constraints$b, lb, ub, 
  #           method, initial, L, N, J, verbose)
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
  moments <- moments2(mu, Sigma, lb - mu, ub - mu, L)
  
  if (method == "epess") {
    std_samples <- t(sample_epess(n, moments$mu, chol(moments$Sigma), A, b, 
                                  J, N, initial, verbose))
  }
  
  browser()
  samples <- std_samples %*% t(L) + 
    # matrix(rep(mu + moments$mu, n), nrow = n, byrow = TRUE)
    matrix(rep(mu, n), nrow = n, byrow = TRUE)
  
  # samples[1, ] <- initial
  
  return(samples)
}