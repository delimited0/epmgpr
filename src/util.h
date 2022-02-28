#ifndef UTIL_H
#define UTIL_H

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

double erfcx (double x);
Rcpp::List trunc_norm_moments(arma::vec lb_in, arma::vec ub_in, 
                              arma::vec mu_in, arma::vec sigma_in);

#endif