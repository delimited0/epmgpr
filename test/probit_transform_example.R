library(epmgpr)

transform_check = readRDS('test/transform_check.RDS')
Xbeta = transform_check$Xbeta
Sigma = transform_check$Sigma
y = transform_check$y
A = transform_check$A
AXbeta = A %*% Xbeta
m = nrow(Xbeta)
trans_result = axisepmgp(
  m = rep(0, m), 
  K = A %*% tcrossprod(Sigma , A),
  lb = -AXbeta, 
  ub = rep(Inf, m), 
  100
)

A %*% trans_result$mu + Xbeta

original_result = epmgp(
  m = Xbeta, 
  K = Sigma, 
  C = t(A), 
  lb = rep(0, m), 
  ub = rep(Inf, m), 
  100
)

original_result$mu

original_result$Sigma


# transformed as input to polytope ep
trans_poly_result = epmgp(rep(0, m), A %*% tcrossprod(Sigma , A), diag(m), -AXbeta, rep(Inf, m), 100)

A %*% trans_poly_result$mu + Xbeta

A %*% trans_poly_result$Sigma %*% t(A)


cbind(
  original_result$mu,
  A %*% trans_poly_result$mu + Xbeta,
  A %*% trans_result$mu + Xbeta
)


# a convergent example -------------------------------------------------------
library(epmgpr)
working_check = readRDS('test/working_example.RDS')
Xbeta = working_check$Xbeta
Sigma = working_check$Sigma
y = working_check$y
A = working_check$A
AXbeta = A %*% Xbeta
m = nrow(Xbeta)

original_result = epmgp(
  m = Xbeta, 
  K = Sigma, 
  C = t(A), 
  lb = rep(0, m), 
  ub = rep(Inf, m), 
  100
)

# transformed as input to polytope ep
trans_poly_result = epmgp(rep(0, m), A %*% tcrossprod(Sigma , A), diag(m), -AXbeta, rep(Inf, m), 100)

# transformed in code 
trans_result = axisepmgp(
  m = rep(0, m), 
  K = A %*% tcrossprod(Sigma , A),
  lb = -AXbeta, 
  ub = rep(Inf, m), 
  100
)

cbind(
  original_result$mu,
  A %*% trans_poly_result$mu + Xbeta,
  A %*% trans_result$mu + Xbeta
)

cbind(
  A %*% (original_result$mu - Xbeta), 
  trans_poly_result$mu,
  trans_result$mu  
)




