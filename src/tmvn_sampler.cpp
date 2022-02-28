#include "tmvn_sampler.h"

double TmvnSampler::lpdf_tmvn(const arma::vec & x) {
  
  arma::vec eval = F * x + g;
  
  double log_density;
  if (arma::any(eval < 0))
    log_density = -arma::datum::inf;
  else {
    double q = arma::dot(x, x);
    double c = this->dim * std::log(2. * arma::datum::pi);
    log_density = -(c + q) / 2.;
    
    // arma::vec xmu = x - this->ep_mean;
    // arma::rowvec Q = xmu.t() * this->ep_chol_inv;
    // double q = arma::as_scalar(Q * Q.t());
    // double c = this->dim * std::log(2. * arma::datum::pi) + 
      // 2 * arma::sum(arma::log(arma::diagvec(this->ep_chol)));
    // log_density = -(c + q) / 2.;
  }
  
  return log_density;
}

