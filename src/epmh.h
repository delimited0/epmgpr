#ifndef EPMH_H
#define EPMH_H

#include "tmvn_sampler.h"

class EPMH : public TmvnSampler {
public:
  arma::vec curr_sample;
  double curr_log_like;
  
  EPMH(arma::vec initial, arma::vec ep_mean, arma::mat ep_chol,
       arma::mat F, arma::vec g) :
    curr_sample(initial), TmvnSampler(F, g, ep_mean, ep_chol) {
    
    curr_log_like = this->lpdf_tmvn(initial);
  }
};

#endif