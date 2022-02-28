#include "uniform_epess.h"
#include <complex>

const double TOL = 1e-10;
const int MAX_REJECT = 100;

// [[Rcpp::export]]
arma::vec range_intersection(arma::vec first, arma::vec second) {
  arma::vec intersection = arma::zeros(second.n_elem + (first.n_elem - 2));
  
  int k = 0;
  while (!first.is_empty() && !second.is_empty()) {
    
    if (first(0) > second(0)) {
      arma::vec temp = second;
      second = first;
      first = temp;
    }
    
    if (first(1) < second(0)) {
      first.shed_rows(0, 1);
      continue;
    }
    else if (first(1) == second(0)) {
      intersection(k) = second(0);
      intersection(k+1) = second(0);
      k += 2;
      
      first.shed_rows(0, 1);
      continue;
    }
    else {
      if (first(1) == second(1)) {
        intersection(k) = second(0);
        intersection(k+1) = second(1);
        k += 2;
        
        first.shed_rows(0, 1);
        second.shed_rows(0, 1);
      }
      else if (first(1) < second(1)) {
        intersection(k) = second(0);
        intersection(k+1) = first(1);
        k += 2;
        
        first.shed_rows(0, 1);
      }
      else {
        intersection(k) = second(0);
        intersection(k+1) = second(1);
        k += 2;
        
        second.shed_rows(0, 1);
      }
    }
  }
  
  arma::vec result = intersection.subvec(0, k-1);
  return result;
}

double simulate(arma::vec slice_range) {
  // sample angle uniformly from angle slice range
  
  arma::vec diff = arma::zeros(slice_range.n_elem);
  
  for (int i = 0; i < slice_range.n_elem; i += 2) {
    diff(i) = slice_range(i+1) - slice_range(i);
  }
  
  // Rcpp::Rcout << "Diff: " << diff << std::endl;
  
  double total = arma::sum(diff);
  diff /= total;
  
  // double u = arma::randu();
  // int idx;
  // for (int i = 0; i < diff.n_elem; i++) {
  //   if (u < diff(i)) {
  //     idx = i;
  //     break;
  //   }
  // }
  Rcpp::IntegerVector choices = Rcpp::seq(0, diff.n_elem-1);
  int idx = Rcpp::sample(choices, 1, false,
                         Rcpp::NumericVector(diff.begin(), diff.end()))[0];
  
  return arma::randu() * (slice_range(idx+1) - slice_range(idx)) + slice_range(idx);
}

double ThreshRegion::eval(const double x) {
  double val = 
    a4 * std::pow(std::cos(x), 2.) + 
    a3 * std::cos(x) * std::sin(x) +
    a2 * std::sin(x) +
    a1 * std::cos(x) + 
    a0;
  return val;
}

arma::vec UniformEPESS::wall_hitting(const arma::vec & nu) {
  arma::vec a = nu;
  arma::vec b = this->curr_sample;
  
  arma::vec fa = this->F * a;
  arma::vec fb = this->F * b;
  
  arma::vec U = arma::sqrt(arma::square(fa) + arma::square(fb));
  arma::vec phi = arma::atan2(-fa, fb);
  
  arma::vec gplus = this->g + this->F * this->ep_mean;
  arma::uvec pn = arma::find(arma::abs(gplus / U) < 1);
  
  arma::vec angle_slice;
  if (pn.n_elem > 0) {
    
    arma::vec phn = phi(pn);
    
    // time when coordinates hit walls
    arma::vec t1 = -phn + arma::acos(-gplus(pn) / U(pn));  
    arma::vec t2 = -phn + 2 * arma::datum::pi - arma::acos(-gplus(pn) / U(pn));
    
    // t2 is always greater than t1
    arma::vec range = {0.0, t1(0), t2(0), 2. * arma::datum::pi};
    arma::vec new_range;
    if (t1.n_elem > 1) {
      for (int i = 1; i < t1.n_elem; i++) {
        new_range = {0., t1(i), t2(i), 2. * arma::datum::pi};
        range = range_intersection(range, new_range);
      }
    }
    
    angle_slice = range;    
  }
  else {
    angle_slice = {0.0, 2 * arma::datum::pi};
  }
  
  return angle_slice;
} 

arma::mat UniformEPESS::sample() {
  
  arma::mat samples(this->dim, this->N * this->J);

  arma::vec nu = this->ep_chol.t() * arma::randn(this->dim);
  
  // Rcpp::Rcout << "nu: " << nu << std::endl;
  
  arma::vec angle_slice = wall_hitting(nu);
  
  // Rcpp::Rcout << "Angle slice: " << angle_slice << std::endl;
  
  for (int j = 0; j < J; j++) {
    
    // random threshold
    double hh = std::log(arma::randu()) + this->curr_log_like;
    
    // griding threshold
    // double hh = std::log( (double) (j+1) / J) + this->curr_log_like;
    
    arma::mat mat = this->ep_cov_inv - arma::eye(this->dim, this->dim);
    
    double a0 = -2 * hh + arma::as_scalar(nu.t() * mat * nu);
    double a1 = -2 * arma::as_scalar(this->curr_sample.t() * this->ep_mean);
    double a2 = -2 * arma::as_scalar(nu.t() * this->ep_mean);
    double a3 = 2 * arma::as_scalar(this->curr_sample.t() * mat * nu);
    double a4 = arma::as_scalar(this->curr_sample.t() * mat * this->curr_sample) - (a0 + 2 * hh);
    ThreshRegion thresh_region = ThreshRegion(a0, a1, a2, a3, a4);
    
    double a = std::pow(a4, 2) + std::pow(a3, 2);
    double b = (2 * a1 * a4) + (2 * a2 * a3);
    double c = std::pow(a1, 2) + 2 * a0 * a4 - std::pow(a3, 2) + std::pow(a2, 2);
    double d = (2 * a0 * a1) - (2 * a2 * a3);
    double e = std::pow(a0, 2) - std::pow(a2, 2);
    
    double p = (8*a*c - 3*std::pow(b, 2)) / (8*std::pow(a, 2));
    double q = (std::pow(b, 3) - 4*a*b*c + 8*d*std::pow(a, 2)) / (8*std::pow(a, 3));
    
    double delta0 = std::pow(c, 2) - 3*b*d + 12*a*e;
    double delta1 = 2*std::pow(c, 3) - 9*b*c*d + 27*std::pow(b, 2)*e + 27*a*std::pow(d, 2) - 72*a*c*e;
    
    arma::cx_double temp = 
      std::pow( arma::cx_double( std::pow(delta1, 2.) - 4*std::pow(delta0, 3.) , 0.), .5);
    arma::cx_double Q = std::pow( (delta1 + temp) / 2. , 1./3.);
    arma::cx_double S = .5 * std::pow( (Q + (delta0 / Q)) / (3 * a) - (2./3.) * p , .5);
    
    arma::cx_double k1 = std::pow((q / S) - 2*p - 4.*std::pow(S, 2), .5);
    arma::cx_double k2 = std::pow(-(q / S) - 2*p - 4.*std::pow(S, 2), .5);
    
    arma::cx_double x1 = -(b / (4*a)) - S + .5*k1;
    arma::cx_double x2 = -(b / (4*a)) - S - .5*k1;
    arma::cx_double x3 = -(b / (4*a)) + S + .5*k2;
    arma::cx_double x4 = -(b / (4*a)) + S - .5*k2;
    
    arma::cx_vec complex_roots = {x1, x2, x3, x4};
    arma::uvec np = arma::find(arma::abs(arma::imag(complex_roots)) < TOL);
    
    // Rcpp::Rcout << "complex roots: " << complex_roots << std::endl;
    // Rcpp::Rcout << "np: " << np << std::endl;
    // Rcpp::Rcout << "complex roots (np): " << complex_roots(np) << std::endl;
    
    arma::vec roots = arma::sort( arma::real(complex_roots(np)) );
    
    // Rcpp::Rcout << "Roots: " << roots << std::endl;
    
    // finding angle ranges
    arma::vec exact_range;
    arma::vec range_1;
    if (roots.n_elem == 0) {
      // no real roots
      exact_range = angle_slice;
    }
    else {
      // finding the root range
      arma::uvec ab = arma::find((roots >= -1.) && (roots <= 1.));
      arma::vec real_roots = roots(ab);
      
      // Rcpp::Rcout << "real roots: " << real_roots << std::endl;
      
      if (real_roots.n_elem == 0) {
        // no real roots between [-1, 1]
        exact_range = angle_slice;
      }
      else {
      
        arma::vec theta(real_roots.n_elem + 2);
        theta(0) = 0.0;
        theta(theta.n_elem - 1) = 2. * arma::datum::pi;
        
        arma::vec actual_theta(real_roots.n_elem);
        for (int k = 0; k < real_roots.n_elem; k++) {
          arma::vec values = { std::acos(real_roots(k)), 2. * arma::datum::pi - std::acos(real_roots(k)) };
          arma::vec thresh_vec = { std::abs(thresh_region.eval(values(0))), 
                                   std::abs(thresh_region.eval(values(1))) };
          arma::uword p = arma::index_min(thresh_vec);
          actual_theta(k) = values(p);
        }
        theta.subvec(1, theta.n_elem-2) = arma::sort(actual_theta);
        
        // Rcpp::Rcout << "Theta: " << theta << std::endl;
        
        if (theta.n_elem == 2) {
          
          range_1 = {0.0, 2. * arma::datum::pi};
        }
        else {
          
          // Rcpp::Rcout << "In block 2" << std::endl;
          for (int i = 0; i < theta.n_elem - 1; i++) {
            double point = (theta(i) + theta(i+1)) / 2;
            if (thresh_region.eval(point) > 0) {
              range_1 = arma::join_vert(range_1, theta.subvec(i, i+1));
            }
          }
          
          // Rcpp::Rcout << "Finish block 2, range1: " << range_1 << std::endl;
        }
        
        exact_range = range_intersection(angle_slice, range_1);
      }
    }
     
    // ?? what was the point of previous code?
    // exact_range = angle_slice;
    
    // Rcpp::Rcout << "sample uniformly from range " << exact_range << std::endl;
    
    arma::vec proposal;
    double proposal_log_like;
    // Rcpp::Rcout << "sample matrix size: " << arma::size(samples) << std::endl;
    for (int i = j*this->N; i < (j+1)*N; i++) {
      
      int rejects = 0;
      
      // Rcpp::Rcout << "Proposing from range: " << exact_range << std::endl;
      // Rcpp::Rcout << "acceptance threshold: " << hh << std::endl;
      while (true) {
        
        double phi = simulate(exact_range);
        proposal = this->curr_sample * std::cos(phi) + nu * std::sin(phi);
        proposal_log_like = this->pseudo_llik(proposal);
        
        // Rcpp::Rcout << "proposal: " << proposal << std::endl;
        // Rcpp::Rcout << "proposal log like: " << proposal_log_like << std::endl;
        
        // arma::vec eval = F * proposal + g;
        // if (arma::any(eval < 0))
        //   Rcpp::Rcout << "out of bounds" << std::endl;
        
        if (proposal_log_like > hh) {
          // Rcpp::Rcout << "------" << std::endl;
          break;
        }
        
        rejects++;
      }
      
      // Rcpp::Rcout << "i: " << i << std::endl;
      // Rcpp::Rcout << "saving to matrix, proposal: " << proposal << std::endl;
      this->curr_log_like = proposal_log_like;
      this->curr_sample = proposal;
      samples.col(i) = this->curr_sample;
      // Rcpp::Rcout << "sample stored" << std::endl;
    }
  }
  
  return samples;
}

double UniformEPESS::pseudo_llik(const arma::vec & x) {
  
  arma::vec q = this->ep_chol_inv * x;
  double lpdf_prior =  
    -.5 * arma::as_scalar(
        arma::dot(q, q) + this->dim*std::log(2*arma::datum::pi) +
        2 * arma::sum(arma::log(arma::diagvec(this->ep_chol))) );
  
  return this->lpdf_tmvn(x + this->ep_mean) - lpdf_prior;
}

// [[Rcpp::export]]
arma::mat sample_epess(int n_samples, arma::vec ep_mean, arma::mat ep_chol, 
                       arma::mat F, arma::vec g, int J, int N, arma::vec initial,
                       bool verbose) {
  
  UniformEPESS epess = UniformEPESS(initial, ep_mean, ep_chol, F, g, J, N);
  arma::mat samples = arma::zeros(F.n_cols, n_samples);
  samples.col(0) = initial;
  
  arma::mat output;
  for (int i = 1; i < (n_samples - N * J + 1); i += N*J) {
    
    if (verbose) Rcpp::Rcout << "Iteration: " << i << std::endl;
    
    output = epess.sample();
    samples.cols(i, i + N*J - 1) = output;
    epess.curr_sample = output.col(N*J - 1);
  } 
  
  return samples;
}

// [[Rcpp::export]]
arma::vec test_wall_hit(arma::vec nu, arma::vec initial, 
                        arma::vec ep_mean, arma::mat ep_chol,
                        arma::mat F, arma::vec g, int J, int N) {
  UniformEPESS epess = UniformEPESS(initial, ep_mean, ep_chol, F, g, J, N);
  return epess.wall_hitting(nu);
}
