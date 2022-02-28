#include "util.h"

const double EPS_CONVERGE = 1e-5;

// [[Rcpp::export]]
Rcpp::List axisepmgp(arma::vec m, arma::mat K, arma::vec lb, arma::vec ub) {
  
  arma::vec nu_site = arma::zeros(K.n_rows);
  arma::vec tau_site = arma::zeros(K.n_rows);
  
  // initialize q(x)
  arma::mat Sigma = K;
  arma::vec mu = (lb + ub) / 2;
  for (int i = 0; i < mu.n_elem; i++) {
    if (std::isinf(mu(i))) {
      mu(i) = copysign(1.0, mu(i)) * 100;
    }
  }
  
  // Rcpp::Rcout << "mu: " << mu << std::endl;
  
  arma::mat Kinvm = arma::solve(K, m);
  double logZ = arma::datum::inf;
  arma::vec mu_last = -arma::datum::inf * arma::ones(arma::size(mu));
  bool converged = false;
  int k = 1;
  
  // algorithm loop
  arma::vec tau_cavity;
  arma::vec nu_cavity;
  arma::mat L;
  arma::vec logz_hat;
  arma::vec sigma_hat;
  arma::vec mu_hat;
  
  while (!converged) {
    
    // Rcpp::Rcout << "Iteration " << k << " ==========" << std::endl;
    
    // make cavity distribution
    tau_cavity = 1 / arma::diagvec(Sigma) - tau_site;
    nu_cavity = mu / arma::diagvec(Sigma) - nu_site;
    
    // Rcpp::Rcout << "tau_cavity: " << tau_cavity << std::endl;
    // Rcpp::Rcout << "nu_cavity: " << nu_cavity << std::endl;
    
    // compute moments using truncated normals
    Rcpp::List moments = trunc_norm_moments(lb, ub, nu_cavity / tau_cavity, 1 / tau_cavity);
    arma::vec sigma_hat_out = moments["sigma_hat"];
    arma::vec logz_hat_out = moments["logz_hat"];
    arma::vec mu_hat_out = moments["mu_hat"];
    logz_hat = logz_hat_out;
    sigma_hat = sigma_hat_out;
    mu_hat = mu_hat_out;
    
    // Rcpp::Rcout << "sigma_hat: " << sigma_hat << std::endl;
    
    // update the site parameters
    arma::vec delta_tau_site = (1 / sigma_hat) - tau_cavity - tau_site;
    tau_site += delta_tau_site;
    nu_site = (mu_hat / sigma_hat) - nu_cavity;
    
    // Rcpp::Rcout << "tau_site: " << tau_site << std::endl;
    // Rcpp::Rcout << "nu_site: " << nu_site << std::endl;
    
    // enforce nonnegativity of tau_site
    if (arma::any(tau_site < 0)) {
      for (int i = 0; i < tau_site.n_elem; i++) {
        if (tau_site(i) > -1e-8) {
          tau_site(i) = 0.0;
        }
      }
    }
    
    // update q(x) Sigma and mu
    arma::mat S_site_half = arma::diagmat(arma::sqrt(tau_site));
    L = arma::chol(
      arma::eye(K.n_rows, K.n_cols) + S_site_half * K * S_site_half);
    arma::mat V = arma::solve(L.t(), S_site_half * K);
    Sigma = K - V.t() * V;
    mu = Sigma * (nu_site + Kinvm);
    
    // Rcpp::Rcout << "tau site: " << tau_site << std::endl;
    // Rcpp::Rcout << "tau cavity: " << tau_cavity << std::endl;
    // Rcpp::Rcout << "L: " << L << std::endl;
    
    // check convergence criteria
    if ((arma::norm(mu_last - mu)) < EPS_CONVERGE)
      converged = true;
    else
      mu_last = mu;
    k++;
  }
  
  if (logZ != -arma::datum::inf) {
    double lZ1 = 0.5 * arma::sum(arma::log(1 + tau_site / tau_cavity)) - 
      arma::sum(arma::log(arma::diagvec(L)));
    double lZ2 = 0.5 *arma::as_scalar(
      (nu_site - tau_site % m).t() * 
      (Sigma - arma::diagmat(1 / (tau_cavity + tau_site))) * 
      (nu_site - tau_site % m)
    );
    double lZ3 = 0.5 * arma::as_scalar(
      nu_cavity.t() * 
        arma::solve(arma::diagmat(tau_site) + arma::diagmat(tau_cavity), 
                    tau_site % nu_cavity / tau_cavity - 2 * nu_site)
    );
    double lZ4 = -0.5 * arma::as_scalar(
      (tau_cavity % m).t() * 
        arma::solve(arma::diagmat(tau_site) + arma::diagmat(tau_cavity), 
                    tau_site % m - 2 * nu_site)
    );

    // Rcpp::Rcout << "lz1: " << lZ1 << std::endl;
    // Rcpp::Rcout << "lz2: " << lZ2 << std::endl;
    // Rcpp::Rcout << "lz3: " << lZ3 << std::endl;
    // Rcpp::Rcout << "lz4: " << lZ4 << std::endl;
    // Rcpp::Rcout << "logzhat: " << logz_hat << std::endl;
    
    logZ = lZ1 + lZ2 + lZ3 + lZ4 + arma::sum(logz_hat);
  }
  
  Rcpp::List result = Rcpp::List::create(
    Rcpp::_["logZ"] = logZ,
    Rcpp::_["mu"] = mu,
    Rcpp::_["Sigma"] = Sigma
  );
  
  return result;
}