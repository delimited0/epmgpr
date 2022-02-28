#ifndef TMVN_SAMPLER_H
#define TMVN_SAMPLER_H

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

class TmvnSampler {
public:
  
  TmvnSampler(arma::mat F, arma::vec g, arma::vec ep_mean, arma::mat ep_chol) : 
    dim(F.n_cols), n_constr(F.n_rows), F(F), g(g),  ep_mean(ep_mean), 
    ep_chol(ep_chol) {
    
    ep_chol_inv = arma::inv(arma::trimatu(ep_chol));
  }
  
  double lpdf_tmvn(const arma::vec & x);
  
  int dim;
  int n_constr;
  arma::mat F;
  arma::vec g;
  arma::vec ep_mean;
  arma::mat ep_chol;  // upper cholesky factor
  arma::mat ep_chol_inv;
};


#endif