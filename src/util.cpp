#include "util.h"

// scaled complementary error function
// from https://stackoverflow.com/a/39777361
// [[Rcpp::export]]
double erfcx (double x) {
  double a, d, e, m, p, q, r, s, t;

  a = fmax (x, 0.0 - x); // NaN preserving absolute value computation

  /* Compute q = (a-4)/(a+4) accurately. [0,INF) -> [-1,1] */
  m = a - 4.0;
  p = a + 4.0;
  r = 1.0 / p;
  q = m * r;
  t = fma (q + 1.0, -4.0, a);
  e = fma (q, -a, t);
  q = fma (r, e, q);

  /* Approximate (1+2*a)*exp(a*a)*erfc(a) as p(q)+1 for q in [-1,1] */
  p =             0x1.edcad78fc8044p-31;  //  8.9820305531190140e-10
  p = fma (p, q,  0x1.b1548f14735d1p-30); //  1.5764464777959401e-09
  p = fma (p, q, -0x1.a1ad2e6c4a7a8p-27); // -1.2155985739342269e-08
  p = fma (p, q, -0x1.1985b48f08574p-26); // -1.6386753783877791e-08
  p = fma (p, q,  0x1.c6a8093ac4f83p-24); //  1.0585794011876720e-07
  p = fma (p, q,  0x1.31c2b2b44b731p-24); //  7.1190423171700940e-08
  p = fma (p, q, -0x1.b87373facb29fp-21); // -8.2040389712752056e-07
  p = fma (p, q,  0x1.3fef1358803b7p-22); //  2.9796165315625938e-07
  p = fma (p, q,  0x1.7eec072bb0be3p-18); //  5.7059822144459833e-06
  p = fma (p, q, -0x1.78a680a741c4ap-17); // -1.1225056665965572e-05
  p = fma (p, q, -0x1.9951f39295cf4p-16); // -2.4397380523258482e-05
  p = fma (p, q,  0x1.3be1255ce180bp-13); //  1.5062307184282616e-04
  p = fma (p, q, -0x1.a1df71176b791p-13); // -1.9925728768782324e-04
  p = fma (p, q, -0x1.8d4aaa0099bc8p-11); // -7.5777369791018515e-04
  p = fma (p, q,  0x1.49c673066c831p-8);  //  5.0319701025945277e-03
  p = fma (p, q, -0x1.0962386ea02b7p-6);  // -1.6197733983519948e-02
  p = fma (p, q,  0x1.3079edf465cc3p-5);  //  3.7167515521269866e-02
  p = fma (p, q, -0x1.0fb06dfedc4ccp-4);  // -6.6330365820039094e-02
  p = fma (p, q,  0x1.7fee004e266dfp-4);  //  9.3732834999538536e-02
  p = fma (p, q, -0x1.9ddb23c3e14d2p-4);  // -1.0103906603588378e-01
  p = fma (p, q,  0x1.16ecefcfa4865p-4);  //  6.8097054254651804e-02
  p = fma (p, q,  0x1.f7f5df66fc349p-7);  //  1.5379652102610957e-02
  p = fma (p, q, -0x1.1df1ad154a27fp-3);  // -1.3962111684056208e-01
  p = fma (p, q,  0x1.dd2c8b74febf6p-3);  //  2.3299511862555250e-01

  /* Divide (1+p) by (1+2*a) ==> exp(a*a)*erfc(a) */
  d = a + 0.5;
  r = 1.0 / d;
  r = r * 0.5;
  q = fma (p, r, r); // q = (p+1)/(1+2*a)
  t = q + q;
  e = (p - q) + fma (t, -a, 1.0); // residual: (p+1)-q*(1+2*a)
  r = fma (e, r, q);

  /* Handle argument of infinity */
  if (a > 0x1.fffffffffffffp1023) r = 0.0;

  /* Handle negative arguments: erfcx(x) = 2*exp(x*x) - erfcx(|x|) */
  if (x < 0.0) {
    s = x * x;
    d = fma (x, x, -s);
    e = exp (s);
    r = e - r;
    r = fma (e, d + d, r);
    r = r + e;
    if (e > 0x1.fffffffffffffp1023) r = e; // avoid creating NaN
  }
  return r;
}

// [[Rcpp::export]]
Rcpp::List trunc_norm_moments(arma::vec lb_in, arma::vec ub_in, 
                              arma::vec mu_in, arma::vec sigma_in) {
  int d = lb_in.n_elem;
  arma::vec logz_hat_out(d);
  arma::vec z_hat_out(d);
  arma::vec mu_hat_out(d);
  arma::vec sigma_hat_out(d);
  
  for (int i = 0; i < d; i++) {
    
    double lb = lb_in(i);
    double ub = ub_in(i);
    double mu = mu_in(i);
    double sigma = sigma_in(i);
    
    double logz_hat_other_tail;
    double logz_hat;
    double mean_const;
    double var_const;
    
    // establish bounds
    double a = (lb - mu) / std::sqrt(2 * sigma);
    double b = (ub - mu) / std::sqrt(2 * sigma);
    
    // Rcpp::Rcout << "a: " << a << std::endl;
    // Rcpp::Rcout << "b: " << b << std::endl;
    
    // stable calculation
    
    // problem case 1
    if (std::isinf(a) && std::isinf(b)) {
      // check the sign
      if (copysign(1.0, a) == copysign(1.0, b)) {
        // integrating from inf to inf, should be 0
        logz_hat_out(i) = -arma::datum::inf;
        z_hat_out(i) = 0.0;
        mu_hat_out(i) = a;
        sigma_hat_out(i) = 0.0;
        continue;
      }
      else {
        logz_hat_out(i) = 0.0;
        z_hat_out(i) = 1.0;
        mu_hat_out(i) = mu;
        sigma_hat_out(i) = sigma;
        continue;
      }
    }
    
    // problem case 2
    else if (a > b) {
      // bounds are wrong, return 0 by convention
      logz_hat_out(i) = -arma::datum::inf;
      z_hat_out(i) = 0;
      mu_hat_out(i) = mu;
      sigma_hat_out(i) = 0;
      continue;
    }
    
    // typical case 1
    else if (a == -arma::datum::inf) {
      // integrating up to b
      if (b > 26.0) {
        // will be very close to 1
        logz_hat_other_tail = std::log(0.5) + std::log(erfcx(b)) - std::pow(b, 2);
        logz_hat = std::log1p(-std::exp(logz_hat_other_tail));
      }
      else {
        // b is small enough 
        logz_hat = std::log(0.5) + std::log(erfcx(-b)) - std::pow(b, 2);
      }
      
      mean_const = -2.0 / erfcx(-b);
      var_const = (-2.0 / erfcx(-b)) * (ub + mu);
    }
    
    // typical case 2
    else if (b == arma::datum::inf) {
      // Rcpp::Rcout << "handling unbounded upper constraint" << std::endl;
      // Rcpp::Rcout << "a: " << a << std::endl;
      // integrate from a to inf
      if (a < -26.0) {
        // will be very close to 1
        logz_hat_other_tail = std::log(0.5) + std::log(erfcx(-a)) - std::pow(a, 2);
        logz_hat = std::log1p(-std::exp(logz_hat_other_tail));
      }
      else {
        // should be stable
        logz_hat = std::log(0.5) + std::log(erfcx(a)) - std::pow(a, 2);
      }
      
      mean_const = 2.0 / erfcx(a);
      var_const = (2.0 / erfcx(a)) * (lb + mu);
    }
    
    // typical case 3
    else {
      // range from a to b, need stable exponent calculation
      double exp_a2b2 = std::exp(std::pow(a, 2) - std::pow(b, 2));
      double exp_b2a2 = std::exp(std::pow(b, 2) - std::pow(a, 2));
      
      if (copysign(1.0, a) == copysign(1.0, b)) {
        // exploit symmetry in problem to make calculations stable for erfcx
        double maxab = std::max(std::abs(a), std::abs(b));
        double minab = std::min(std::abs(a), std::abs(b));
        
        logz_hat = 
          std::log(0.5) - std::pow(minab, 2) + 
          std::log( std::abs( 
              std::exp( -(std::pow(maxab, 2) - std::pow(minab, 2))) * 
                erfcx(maxab) - 
                erfcx(minab) ) );
        
        double erfcx_a = erfcx(std::abs(a));
        double erfcx_b = erfcx(std::abs(b));
        mean_const = 2 * copysign(1.0, a) * (
          1 / (( erfcx_a - exp_a2b2 * erfcx_b )) -
            1 / (( exp_b2a2 * erfcx_a - erfcx_b ))
        );
        var_const = 2 * copysign(1.0, a) * (
          (lb + mu) / (erfcx_a - exp_a2b2 * erfcx_b) - 
            (ub + mu) / (exp_b2a2 * erfcx_a - erfcx_b)
        );
      }
      
      else {
        // the signs are different, so b > a and b >= 0 and a <= 0
        // Rcpp::Rcout << "Second branch" << std::endl;
        if (std::abs(b) >= std::abs(a)) { 
          if (a >= -26.0) {
            // do things normally
            // Rcpp::Rcout << "do things normally branch" << std::endl;
            logz_hat = std::log(0.5) - std::pow(a, 2) + std::log( 
              erfcx(a) - std::exp(-(std::pow(b, 2) - std::pow(a, 2))) * erfcx(b)  
            );
            // Rcpp::Rcout << "a: " << a << std::endl;
            // Rcpp::Rcout << "b: " << b << std::endl;
            // Rcpp::Rcout << "log - a^2: " << std::log(0.5) - std::pow(a, 2) << std::endl;
            // Rcpp::Rcout << "thing: " << erfcx(a) - std::exp(-(std::pow(b, 2) - std::pow(a, 2))) * erfcx(b) <<
            //   std::endl;
            // Rcpp::Rcout << "log thing: " << std::log( 
            //     erfcx(a) - std::exp(-(std::pow(b, 2) - std::pow(a, 2))) * erfcx(b)  
            // ) << std::endl;
            // Rcpp::Rcout << "exp(-(b^2 - a^2)): " << std::exp(-(std::pow(b, 2) - std::pow(a, 2))) << std::endl;
            // Rcpp::Rcout << "logz_hat: " << logz_hat << std::endl;
            // Rcpp::Rcout << "erfcx(a): " << erfcx(a) << std::endl;
            // Rcpp::Rcout << "erfcx(b): " << erfcx(b) << std::endl;
            
            mean_const = 2 * (
              1 / (erfcx(a) - exp_a2b2 * erfcx(b)) -
                1 / (exp_b2a2 * erfcx(a) - erfcx(b))
            );
            var_const = 2 * (
              (lb + mu) / (erfcx(a) - exp_a2b2 * erfcx(b)) - 
                (ub + mu) / (exp_b2a2 * erfcx(a) - erfcx(b))
            );
            
            // Rcpp::Rcout << "mean_const: " << mean_const << std::endl;
            // Rcpp::Rcout << "var_const: " << var_const << std::endl;
          }
          
          else {
            // a is too small, so put in something close to 2 instead 
            // Rcpp::Rcout << "Third branch" << std::endl;
            logz_hat = std::log(0.5) - std::pow(b, 2) + std::log( 
              erfcx(-b) - std::exp(-(std::pow(a, 2) - std::pow(b, 2))) * erfcx(-a)
            );
            
            mean_const = 2 * (
              1 / (erfcx(a) - exp_a2b2 * erfcx(b)) -
                1 / (std::exp(std::pow(b, 2)) * 2 - erfcx(b))
            );
            var_const = 2 * (
              (lb + mu) / (erfcx(a) - exp_a2b2 * erfcx(b)) -
                (ub + mu) / (std::exp(std::pow(b, 2)) * 2 - erfcx(b))
            );
          }
        }
        
        else {
          // abs(a) is bigger than abs(b), so reverse the calculation
          if (b <= 26.0) {
            // do things normally but mirrored across 0
            logz_hat = std::log(0.5) - std::pow(b, 2) + std::log(
              erfcx(-b) - std::exp(-(std::pow(a, 2) - std::pow(b, 2))) * erfcx(-a)
            );
            
            mean_const = -2 * (
              1 / (erfcx(-a) - exp_a2b2 * erfcx(-b)) -
                1 / (exp_b2a2 * erfcx(-a) - erfcx(-b))
            );
            var_const = -2 * (
              (lb + mu) / (erfcx(-a) - exp_a2b2 * erfcx(-b)) -
                (ub + mu) / (exp_b2a2 * erfcx(-a) - erfcx(-b)) 
            );
          }
          else {
            // b is too big, put something close to 2 instead
            logz_hat = std::log(0.5) + std::log(
              2. - std::exp(-std::pow(a, 2)) * erfcx(-a) - std::exp(-std::pow(b, 2)) * erfcx(b)
            );
            
            mean_const = -2 * (
              1 / (erfcx(-a) - std::exp(std::pow(a, 2)) * 2) - 
                1 / (exp_b2a2 * erfcx(-a) - erfcx(-b))
            );
            var_const = -2 * (
              (lb + mu) / (erfcx(-a) - std::exp(std::pow(a, 2)) * 2) - 
                (ub + mu) / (exp_b2a2 * erfcx(-a) - erfcx(-b))
            );
          }
        }
      }
    }
    
    // Rcpp::Rcout << "Log z hat: " << logz_hat << std::endl;
    
    double z_hat = std::exp(logz_hat);
    double mu_hat = mu + mean_const * std::sqrt(sigma / (2 * arma::datum::pi));
    double sigma_hat = 
      sigma + var_const * sqrt(sigma / (2 * arma::datum::pi)) +
      std::pow(mu, 2) - std::pow(mu_hat, 2);      
    
    logz_hat_out(i) = logz_hat;
    z_hat_out(i) = z_hat;
    mu_hat_out(i) = mu_hat;
    sigma_hat_out(i) = sigma_hat;
  }
  
  Rcpp::List result = Rcpp::List::create(
    Rcpp::_["logz_hat"] = logz_hat_out,
    Rcpp::_["z_hat"] = z_hat_out,
    Rcpp::_["mu_hat"] = mu_hat_out,
    Rcpp::_["sigma_hat"] = sigma_hat_out
  );
  
  return result;
}