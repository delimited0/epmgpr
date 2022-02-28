#include "util.h"

const double EPS_CONVERGE = 1e-8;

// [[Rcpp::export]]
Rcpp::List epmgp(arma::vec m, arma::mat K, arma::mat C, arma::vec lb, arma::vec ub,
                 int max_steps) {
  int n = m.n_elem;
  int p = C.n_cols;
  
  arma::vec nu_site = arma::zeros(p);
  arma::vec tau_site = arma::zeros(p);
  
  arma::vec nu_cavity(p);
  arma::vec tau_cavity(p);
  arma::vec delta_tau_site(p);
  arma::vec delta_nu_site(p);
  arma::vec logz_hat(p);
  arma::vec mu_hat(p);
  arma::vec sigma_hat(p);
    
  // arma::vec frac_terms = 1 / arma::sum(arma::abs(C.t() * C), 0);
  arma::vec frac_terms = arma::ones(p);
  arma::vec damp_terms = arma::ones(p);
  double damp_redamp = 0.8;
  
  arma::mat Sigma = K;
  arma::mat mu = m;
  
  arma::mat L = arma::chol(K, "lower");
  
  double logz = arma::datum::inf;
  arma::vec mu_last = -arma::datum::inf * arma::ones(mu.n_elem);
  bool converged = false;
  int k = 1;
  bool restart_algorithm = false;
  
  while (!converged && k <= max_steps) {
    
    restart_algorithm = false;
    
    for (int j = 0; j < p; j++) {
      
      // make the cavity distribution
      arma::vec Cj = C.col(j);
      tau_cavity(j) = 
        1 / arma::as_scalar(Cj.t() * Sigma * Cj) - tau_site(j) / frac_terms(j);
      nu_cavity(j) = 
        arma::dot(Cj, mu) / (arma::as_scalar(Cj.t() * Sigma * Cj)) -
        nu_site(j) / frac_terms(j);
      
      // Rcpp::Rcout << "Made cavity distribution" << std::endl;
      // Rcpp::Rcout << "tau_cavity: " << tau_cavity << std::endl;
      // Rcpp::Rcout << "nu_cavity: " << nu_cavity << std::endl;

      if (tau_cavity(j) <= 0) {
        // problem negative cavity updates
        restart_algorithm = true;
        damp_terms(j) = damp_terms(j) * damp_redamp;
        break;
      }
      
      // compute moments using truncated normals
      Rcpp::List moments = trunc_norm_moments(
        lb.row(j), ub.row(j), 
        nu_cavity.row(j) / tau_cavity.row(j), 1. / tau_cavity.row(j));
      
      arma::vec sigma_hat_out = moments["sigma_hat"];
      arma::vec logz_hat_out = moments["logz_hat"];
      arma::vec mu_hat_out = moments["mu_hat"];
      
      sigma_hat(j) = sigma_hat_out(0);
      logz_hat(j) = logz_hat_out(0);
      mu_hat(j) = mu_hat_out(0);
      
      // Rcpp::Rcout << "Made moments" << std::endl;
      // Rcpp::Rcout << "sigma_hat: " << sigma_hat << std::endl;
      // Rcpp::Rcout << "logz_hat: " << logz_hat << std::endl;
      // Rcpp::Rcout << "mu_hat: " << mu_hat << std::endl;
      
      if (sigma_hat(j) == 0) {
        // the algorithm has found a 0 weight dimension, terminate
        converged = true;
        logz = -arma::datum::inf;
        mu = arma::datum::nan;
        Sigma = arma::datum::nan;
        break;
      }
      
      delta_tau_site(j) = 
        damp_terms(j) * (frac_terms(j) * (1 / sigma_hat(j) - tau_cavity(j)) - tau_site(j));
      delta_nu_site(j) = 
        damp_terms(j) * (frac_terms(j) * (mu_hat(j) / sigma_hat(j) - nu_cavity(j)) - nu_site(j));
      tau_site(j) += delta_tau_site(j);
      nu_site(j) += delta_nu_site(j);
    
      // Rcpp::Rcout << "Updated sites" << std::endl;
      // Rcpp::Rcout << "delta_tau_site: " << delta_tau_site << std::endl;
      // Rcpp::Rcout << "delta_nu_site: " << delta_nu_site << std::endl;
      // Rcpp::Rcout << "tau_site: " << tau_site << std::endl;
      // Rcpp::Rcout << "nu_site: " << nu_site << std::endl;
      
      if (tau_site(j) < 0) {
        // if result negative, either due to numerical precision or error
        if (tau_site(j) > -1e-6)
          tau_site(j) = 0;
      }
      
      // update q(x) (Sigma and mu)
      arma::vec sc = Sigma * Cj;
      Sigma -= 
        (delta_tau_site(j) / (1 + delta_tau_site(j) * arma::dot(Cj, sc))) *
        (sc * sc.t());
      mu += 
        ((delta_nu_site(j) - delta_tau_site(j) * arma::dot(Cj, mu)) / 
         (1 + delta_tau_site(j) * arma::dot(Cj, sc))) * sc;

      // Rcpp::Rcout << "Updated q(x)" << std::endl;
      // Rcpp::Rcout << "sc: " << sc << std::endl;
      // Rcpp::Rcout << "Sigma: " << Sigma << std::endl;
      // Rcpp::Rcout << "mu: " << mu << std::endl;
      
    }
    
    // Rcpp::Rcout << "mu_last: " << mu_last << std::endl;
    
    // check convergence criteria
    if (arma::norm(mu_last - mu, 2) < EPS_CONVERGE)
      converged = true;
    else
      mu_last = mu;
    
    k++;
    
    // if sites are oscillating, restart everything
    if (restart_algorithm) {
      k = 1;
      converged = false;
      damp_terms = damp_redamp * damp_terms;
      nu_site = arma::zeros(p, 1);
      tau_site = arma::zeros(p, 1);
      mu = m;
      Sigma = K;
    }
  }
  
  // compute logz, the EP MGP result, from q(x)
  if (logz != -arma::datum::inf) {
    arma::vec mu_cavity = nu_cavity / tau_cavity;
    arma::vec sigma_cavity = 1 / tau_cavity;
    arma::mat lz_det_mat = arma::eye(n, n);
    for (int j = 0; j < p; j++) {
      arma::vec lc = L.t() * C.col(j);
      lz_det_mat += tau_site(j) * lc * lc.t();
    }
    
    // Rcpp::Rcout << "lzdetmat: " << lz_det_mat << std::endl;
    
    double logdet_value;
    double sign;
    arma::log_det(logdet_value, sign, lz_det_mat);
    
    double lz0 = -.5 * logdet_value;
    double lz1 = 
      -.5 * std::pow(arma::norm(arma::solve(L, m)), 2) + 
       .5 * std::pow(arma::norm(arma::solve(L, mu)), 2) + 
       .5 * arma::as_scalar(arma::sum(tau_site % (arma::pow((C.t() * mu), 2))));
    double lz2 = 
      arma::sum(frac_terms % logz_hat) + 
      .5 * arma::sum(frac_terms % arma::log1p(tau_site % sigma_cavity / frac_terms));
    double lz3 = 
      .5 * arma::sum((arma::pow(mu_cavity, 2) % tau_site % frac_terms - 
                      2 * mu_cavity % nu_site % frac_terms -
                      arma::pow(nu_site, 2) % sigma_cavity) / 
                      (frac_terms + tau_site % sigma_cavity));
    
    // Rcpp::Rcout << "lz: " << lz0 << ", " << lz1 << ", " << lz2 << ", " << lz3 << std::endl;
    
    logz = lz0 + lz1 + lz2 + lz3;
  }
  
  Rcpp::List result = Rcpp::List::create(
    Rcpp::_["logZ"] = logz,
    Rcpp::_["mu"] = mu,
    Rcpp::_["Sigma"] = Sigma,
    Rcpp::_["logz_hat"] = logz_hat,
    Rcpp::_["nu_cavity"] = nu_cavity,
    Rcpp::_["tau_cavity"] = tau_cavity,
    Rcpp::_["nu_site"] = nu_site,
    Rcpp::_["tau_site"] = tau_site
  );
  
  return result;
}
