#ifndef BVCUTILITIES_H
#define BVCUTILITIES_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double rtnorm(double mu, double sigma);
arma::vec rtnorm1(int n, arma::vec mean , arma::vec sd ,arma::vec lower_bound,arma::vec upper_bound);
double rtnorm0(double mu, double sigma);
double rtnorm2(double a, bool lb, double mu, double sigma);
double rinvgaussian(double mu, double lambda);
double rinvGauss(double mu, double lambda);
double rinvGauss_1(double mu, double lambda);
arma::vec mvrnormCpp(const arma::vec& mu, const arma::mat& sigma);
arma::vec mvrnormCpp(const arma::vec& mu, const arma::mat& sigma, double tol);
#endif
