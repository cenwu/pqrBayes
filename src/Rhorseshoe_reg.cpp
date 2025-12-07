#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List RBhorse_reg(arma::mat xx, arma::vec y, arma::mat W, int maxSteps, arma::vec hatAlpha, arma::vec hatBeta, double hatSigma, double hatcSq, arma::vec hatV, arma::vec hatSg, arma::vec hatnu, arma::mat invSigAlpha0, double hatTauSq, double xi1, double xi2, double hatZeta, double a, double b, double c, double d, int progress)
{
  unsigned int n = xx.n_rows, s = xx.n_cols, q = W.n_cols;
  arma::mat gsAlpha(maxSteps, q),
  gsBeta(maxSteps, s),
  gsV(maxSteps, n),
  gsnu(maxSteps, s),
  gsKappa(maxSteps, s),
  gsSg(maxSteps, s);
  
  arma::vec 
  gsTauSq(maxSteps),
  gsSigma(maxSteps),
  gscSq(maxSteps),
  gsZeta(maxSteps),
  gsMeff(maxSteps),
  gsMSE(maxSteps);
  
  arma::mat varAlpha, tWWoV(q,q), temp;
  arma::vec res, RWoV(q), meanAlpha, muV,kappa(s);
  double lambV, xi1Sq = std::pow(xi1, 2), xi2Sq = std::pow(xi2, 2), XgXgoV, RXgoV, meanG, varG, ResSqoV, sj, meff;
  
  
  for (int k = 0; k < maxSteps; k++) {
    // Rcpp::Rcout << "alpha" << std::endl;
    res = y - xx * hatBeta - xi1*hatV;
    tWWoV = (W.each_col()/hatV).t() * W;
    RWoV = arma::sum(W.each_col()% (res/hatV), 0).t();
    varAlpha = arma::inv_sympd(tWWoV*hatSigma/xi2Sq + invSigAlpha0);
    meanAlpha = varAlpha * RWoV * hatSigma / xi2Sq;
    hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
    res -= W * hatAlpha;
    gsAlpha.row(k) = hatAlpha.t();
    
    // Rcpp::Rcout << "v" << std::endl;
    res += xi1*hatV;
    lambV = hatSigma*xi1Sq/xi2Sq + 2*hatSigma;
    muV = arma::sqrt((xi1Sq+2*xi2Sq) / arma::square(res));
    for(unsigned int i = 0; i<n; i++){
      bool flag = true;
      while(flag){
        hatV(i) = 1/rinvGauss(muV(i), lambV);
        if(hatV(i)<=0 || std::isinf(hatV(i)) || std::isnan(hatV(i))){
          if(progress != 0) Rcpp::Rcout << "hatV(i) <= 0 or nan or inf" << std::endl; 
          Rcpp::checkUserInterrupt();
        }else{
          flag = false;
        }
      }
    }
    res -= xi1*hatV;
    gsV.row(k) = hatV.t();
    
    
    // Rcpp::Rcout << "S" << std::endl;
    for(unsigned int j = 0; j<s; j++){
      double rateSg = 1/hatnu(j)+pow(hatBeta(j),2)/(2*hatTauSq);
      hatSg(j) = 1/R::rgamma(1, 1/rateSg);
    }
    gsSg.row(k) = hatSg.t();
    
    // Rcpp::Rcout << "beta" << std::endl;
    for(unsigned int j=0; j<s; j++){
      res += xx.col(j) * hatBeta(j);
      XgXgoV = arma::as_scalar((xx.col(j)/hatV).t() * xx.col(j));
      double Sg = hatcSq*hatSg(j)/(hatcSq+hatTauSq*hatSg(j));
      varG = 1/(XgXgoV*hatSigma/xi2Sq + 1/(hatTauSq*Sg));
      
      RXgoV = arma::sum(xx.col(j) % (res/hatV));
      meanG = varG * RXgoV * hatSigma / xi2Sq;
      hatBeta(j) = R::rnorm(meanG, sqrt(varG));
      res -= xx.col(j) * hatBeta(j);
    }
    gsBeta.row(k) = hatBeta.t();
    // csq|
    double shapec = (s+c)/2;
    double ratec = 0.5*(d+arma::accu(square(hatBeta)));
    hatcSq = 1/R::rgamma(shapec, 1/ratec);
    gscSq(k) = hatcSq;
    
    // Rcpp::Rcout << "meff" << std::endl;
    for (unsigned int j = 0; j < s; j++) {
      sj = arma::accu(arma::square(xx.col(j)) / (xi2Sq * (1/hatSigma)* hatV));
      kappa(j) = 1/ (1 + hatTauSq * hatSg(j) * sj);
    }
    meff = arma::accu(1- kappa);
    
    gsMeff(k) = meff;
    gsKappa.row(k) = kappa.t();
    
    // tausq|
    double shapeS = (s+1)/2;
    double rateS = 1/hatZeta + arma::accu(square(hatBeta)%(1/hatSg))/2;
    hatTauSq = 1/R::rgamma(shapeS, 1/rateS);
    gsTauSq(k) = hatTauSq;
    
    //nu
    for(unsigned int j = 0; j<s; j++){
      double ratenu = 1+1/hatSg(j);
      hatnu(j) = 1/R::rgamma(1, 1/ratenu);
    }
    gsnu.row(k) = hatnu.t();
    
    // zeta|
    double ratezeta = 1+1/hatTauSq;
    hatZeta = 1/R::rgamma(1, 1/ratezeta);
    gsZeta(k) = hatZeta;
    
    // Rcpp::Rcout << "tau" << std::endl;
    double shape = a + 3*n/2;
    ResSqoV = arma::accu(arma::square(res)/hatV);
    double rate = b + arma::accu(hatV) + ResSqoV/(2*xi2Sq);
    hatSigma = R::rgamma(shape, 1/rate);
    gsSigma(k) = hatSigma;
    
    
    gsMSE(k) = arma::mean(arma::abs(res));
    if(k % 100 == 0){
      Rcpp::checkUserInterrupt();
    }
    if(progress != 0 && k % progress == 0){
      Rcpp::Rcout << "\nIter." << k << "  MAD: " << gsMSE(k) << std::endl;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("GS.alpha") = gsAlpha,
                            Rcpp::Named("GS.beta") = gsBeta,
                            Rcpp::Named("GS.sigma") = gsSigma,
                            Rcpp::Named("GS.v") = gsV,
                            Rcpp::Named("GS.s") = gsSg,
                            Rcpp::Named("GS.kappa") = gsKappa,
                            Rcpp::Named("GS.meff") = gsMeff,
                            Rcpp::Named("GS.mad") = gsMSE);
}