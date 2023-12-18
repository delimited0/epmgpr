library(epmgpr)
d = 2
mu = rep(0, d)
Sigma = .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
lb = rep(0, d)
ub = rep(Inf, d)

orthant_lesson = didactic_moments(mu, Sigma, lb, ub)

orthant_lesson$tau_site_history
orthant_lesson$nu_site_history
orthant_lesson$tau_cavity_history
orthant_lesson$nu_cavity_history


estimate = axisepmgp(mu, Sigma, lb, ub, 200)
estimate$logZ
estimate$mu
estimate$Sigma

d = 2
mu = rep(0, d)
Sigma = .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
A = matrix(c(1, .2, .2, 1), nrow = 2)
lb = rep(0, d)
ub = rep(Inf, d)

lincon_lesson = didactic_moments(mu, Sigma, lb, ub, A)