#ifndef UNIFORM_EPESS_H
#define UNIFORM_EPESS_H

#include "tmvn_sampler.h"

class UniformEPESS : public TmvnSampler {
public:
  arma::vec curr_sample;
  double curr_log_like;
  arma::mat ep_cov_inv;
  int J;  // number of thresholds per ellipse
  int N;  // number of samples per angle proposal range
  
  UniformEPESS(arma::vec initial, arma::vec ep_mean, arma::mat ep_chol, 
               arma::mat F, arma::vec g, int J, int N) : 
    curr_sample(initial), J(J), N(N), TmvnSampler(F, g, ep_mean, ep_chol) {
    
    ep_cov_inv = ep_chol_inv * ep_chol_inv.t();
    curr_log_like = this->pseudo_llik(initial);
  }
  
  arma::vec wall_hitting(const arma::vec & nu);
  
  arma::mat sample();
  
  double pseudo_llik(const arma::vec & x);
};

class ThreshRegion {
public:
  double a0;
  double a1;
  double a2;
  double a3;
  double a4;
  
  ThreshRegion(double a0, double a1, double a2, double a3, double a4) :
    a0(a0), a1(a1), a2(a2), a3(a3), a4(a4) {
  }
  
  double eval(const double x);
};



#endif 
