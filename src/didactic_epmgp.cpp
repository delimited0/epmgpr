#include "util.h"

const double EPS_CONVERGE = 1e-5;

// [[Rcpp::export]]
Rcpp::List didactic_axisepmgp(arma::vec m, arma::mat K, arma::vec lb, arma::vec ub,
                              int max_steps) {
  
  arma::mat nu_site_history = arma::zeros(K.n_rows, max_steps);
  arma::mat tau_site_history = arma::zeros(K.n_rows, max_steps);
  arma::mat nu_cavity_history = arma::zeros(K.n_rows, max_steps);
  arma::mat tau_cavity_history = arma::zeros(K.n_rows, max_steps);
    
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
  
  while (!converged && k <= max_steps) {
    
    // Rcpp::Rcout << "Iteration " << k << " ==========" << std::endl;
    
    // make cavity distribution
    tau_cavity = 1 / arma::diagvec(Sigma) - tau_site;
    nu_cavity = mu / arma::diagvec(Sigma) - nu_site;
    
    tau_cavity_history.col(k-1) = tau_cavity;
    nu_cavity_history.col(k-1) = nu_cavity;
    
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
    
    tau_site_history.col(k-1) = tau_site;
    nu_site_history.col(k-1) = nu_site;
    
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
    Rcpp::_["Sigma"] = Sigma,
    Rcpp::_["tau_site_history"] = tau_site_history,
    Rcpp::_["nu_site_history"] = nu_site_history,
    Rcpp::_["tau_cavity_history"] = tau_cavity_history,
    Rcpp::_["nu_cavity_history"] = nu_cavity_history,
    Rcpp::_["iters"] = k-1
  );
  
  return result;
}


// [[Rcpp::export]]
Rcpp::List didactic_epmgp(arma::vec m, arma::mat K, arma::mat C, arma::vec lb, arma::vec ub,
                 int max_steps) {
  int n = m.n_elem;
  int p = C.n_cols;
  
  arma::mat nu_site_history = arma::zeros(K.n_rows, max_steps);
  arma::mat tau_site_history = arma::zeros(K.n_rows, max_steps);
  arma::mat nu_cavity_history = arma::zeros(K.n_rows, max_steps);
  arma::mat tau_cavity_history = arma::zeros(K.n_rows, max_steps);
  
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
      
      tau_cavity_history.col(k-1) = tau_cavity;
      nu_cavity_history.col(k-1) = nu_cavity;
      
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
      
      tau_site_history.col(k-1) = tau_site;
      nu_site_history.col(k-1) = nu_site;
      
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
    Rcpp::_["tau_site_history"] = tau_site_history,
    Rcpp::_["nu_site_history"] = nu_site_history,
    Rcpp::_["tau_cavity_history"] = tau_cavity_history,
    Rcpp::_["nu_cavity_history"] = nu_cavity_history,
    Rcpp::_["iters"] = k-1
  );
  
  return result;
}