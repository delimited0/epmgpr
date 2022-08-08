#' @export
didactic_moments = function(mu, Sigma, lb, ub, A = NULL, max_steps = 200) {
  
  if (is.null(A)) {
    result = didactic_axisepmgp(mu, Sigma, lb, ub, max_steps)
  }
  else {
    result = didactic_epmgp(mu, Sigma, A, lb, ub, max_steps)
  }
  
  iters = result$iters
  result$tau_site_history = result$tau_site_history[, 1:iters]
  result$nu_site_history = result$nu_site_history[, 1:iters]
  result$tau_cavity_history = result$tau_cavity_history[, 1:iters]
  result$nu_cavity_history = result$nu_cavity_history[, 1:iters]
  
  return(result)
}