#include "epmh.h"

// [[Rcpp::export]]
arma::mat sample_epmh(int n_samples, arma::vec ep_mean, arma::mat ep_chol,
                      arma::mat F, arma::vec g, arma::vec initial) {
  EPMH epmh = EPMH(initial, ep_mean, ep_chol, F, g);
  arma::mat samples = arma::zeros(F.n_cols, n_samples);
  samples.col(0) = initial;
  
  arma::vec proposal;
  double prop_log_like;
  for (int i = 1; i < n_samples; i++) {
    
    // propose until a point is accepted
    while (true) {
      proposal = epmh.ep_mean + epmh.ep_chol.t() * arma::randn(epmh.dim);
      prop_log_like = epmh.lpdf_tmvn(proposal);
    
      if (prop_log_like - epmh.curr_log_like > std::log(arma::randu())) {
        
        epmh.curr_log_like = prop_log_like;
        epmh.curr_sample = proposal;
        
        break;
      }
    }
    
    samples.col(i) = epmh.curr_sample;
  }
  
  return samples;
}