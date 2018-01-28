// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Scraped from here to get a way to set seed:
// http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
void set_seed(unsigned int seed) {

  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);

}

double sign(double x) {
  if (x > 0.0) {
    return 1.0;
  } else if (x < 0.0) {
    return -1.0;
  } else {
    return 0.0;
  }
}

//' Thresholding function for penalized regression with power penalty
//'
//' @name powThresh
//'
//' @description Gives the value of \eqn{\beta} that minimizes: \cr
//' \deqn{(\beta - z)^2/2 + \lambda|\beta|^q} \cr
//' for fixed \eqn{z}, \eqn{\lambda > 0} and \eqn{q > 0}. \cr \cr
//' For \eqn{q \le 1}, uses thresholding function given in Marjanovic and Solo (2014).
//' For \eqn{q > 1}, computed bisection.
//'
//' @usage \code{powThresh(z =  1, lambda = 1, q = 1)}
//'
//' @param \code{z} vector of values to find the minimizing beta for (\code{powThresh} is applied coordinatewise when \code{length(z)} > 1)
//' @param \code{lambda} scalar multiplier applied to penalty
//' @param \code{q} scalar power of penalty function
//'
//' @source For \eqn{q\le 1}, uses thresholding function given in Marjanovic and Solo (2014).
//' Marjanovic, Goran, and Victor Solo. "\eqn{l_q} Sparsity Penalized Linear Regression With Cyclic Descent." IEEE Transactions on Signal Processing 62.6 (2014): 1464-1475.
//'
// [[Rcpp::export]]
NumericVector powThresh(NumericVector z, double lambda, double q) {

  int p = z.size();
  NumericVector js(p, 0.0);
  double betaLambda = pow(2.0*lambda*(1.0 - q), 1.0/(2.0 - q));
  double hLambda = betaLambda + lambda*q*pow(betaLambda, q - 1.0);
  double betaBar0 = 0.0;
  double betaBar1 = 0.0;
  double l = 0.0;
  double u = 0.0;
  double m = 0.0;

  if (q < 1) {
    for (int j = 0; j < p; j++) {

      if (fabs(z[j]) < hLambda) {
        js[j] = 0.0;
      } else if (z[j] == 0.0) {
        js[j] = std::max(0.0, betaLambda);
      } else {

        betaBar0 = fabs(z[j]) + 1.0;
        betaBar1 = fabs(z[j]);
        while(fabs(betaBar0 - betaBar1) > pow(10.0, -14.0)) {
          betaBar0 = betaBar1;
          betaBar1 = fabs(z[j]) - lambda*q*pow(betaBar0, q - 1.0);
        }

        js[j] = sign(z[j])*betaBar1;

      }

    }
  } else {

    for (int j = 0; j < p; j++) {

      l = 0.0;
      u = fabs(z[j]);

      while (u - l > pow(10.0, -7.0)) {
        m = (l + (u - l)/2.0) + lambda*q*pow(l + (u - l)/2.0, q - 1.0);
        if (m > abs(z[j])) {
          u = l + (u - l)/2.0;
        } else if (fabs(m - abs(z[j])) < pow(10.0, -14.0)) {
          u = l = (l + (u - l)/2.0);
        } else {
          l = l + (u - l)/2.0;
        }
      }

      js[j] = sign(z[j])*(u + l)/2.0;
    }

  }

  return js;
}

