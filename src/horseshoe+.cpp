#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List Bhorse_plus(arma::mat xx, arma::vec y, arma::mat W, int maxSteps, arma::vec hatBeta, arma::vec hatAlpha, arma::vec hatlamSq, arma::vec hatetaSq, arma::vec hatnu_1, arma::vec hatnu_2,arma::mat invSigAlpha0, double hatTauSq, double hatSigmaSq, double hatZeta, double a, double b,int progress)
{
  unsigned int n = xx.n_rows, s = xx.n_cols, clc = W.n_cols;
  arma::mat gsAlpha(maxSteps, clc),
  gsBeta(maxSteps, s),
  gsnu_1(maxSteps, s),
  gsnu_2(maxSteps, s),
  gsetaSq(maxSteps, s),
  gsKappa(maxSteps, s),
  gslamSq(maxSteps, s);
  
  arma::vec gsTauSq(maxSteps),
  gsSigmaSq(maxSteps),
  gsZeta(maxSteps),
  gsMeff(maxSteps),
  gsMSE(maxSteps);
  
  arma::mat tWW = W.t()*W, varAlpha;
  arma::vec res, meanAlpha,kappa(s);
  double tempS, varRs, meanRs, sj, meff;
  
  arma::vec tBrBrDiag = sum(square(xx), 0).t();
  
  for (int k = 0; k < maxSteps; k++) {
    // alpha|
    varAlpha = arma::inv(tWW/hatSigmaSq + invSigAlpha0);
    res = y - xx * hatBeta;
    meanAlpha = varAlpha * (W.t() * res/hatSigmaSq);
    hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
    res -= W * hatAlpha;
    gsAlpha.row(k) = hatAlpha.t();
    
    //Rcpp::Rcout << "beta" << std::endl;
    for(unsigned int j=0; j<s; j++){		
      tempS = 1/(tBrBrDiag(j) + 1/(hatTauSq*hatlamSq(j)));
      varRs = hatSigmaSq * tempS;
      res += xx.col(j) * hatBeta(j);
      meanRs = arma::as_scalar(tempS * xx.col(j).t() * res);
      hatBeta(j) = R::rnorm(meanRs, std::sqrt(varRs));
      res -= xx.col(j) * hatBeta(j);
    }
    gsBeta.row(k) = hatBeta.t();
    
    //Rcpp::Rcout << "sigma2" << std::endl;
    // sigma.sq|
    double shapeSig = a+(n+s)/2;
    double rateSig = b+0.5*(arma::accu(arma::square(res)) + 
                          arma::accu(square(hatBeta) % (1/(hatTauSq*hatlamSq))));
    hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
    gsSigmaSq(k) = hatSigmaSq;
    
    
    // lamdasq|
    //Rcpp::Rcout << "lambda2" << std::endl;
    for(unsigned int j = 0; j<s; j++){
      double ratelamSq = 1/hatnu_1(j)+pow(hatBeta(j),2)/(2*hatTauSq*hatSigmaSq);
      hatlamSq(j) = 1/R::rgamma(1, 1/ratelamSq);
    }
    gslamSq.row(k) = hatlamSq.t();
    
    // Rcpp::Rcout << "meff" << std::endl;
    for (unsigned int j = 0; j < s; j++) {
      sj = arma::accu(arma::square(xx.col(j)));
      kappa(j) = 1/ (1 + (1/ hatSigmaSq)* hatTauSq * sj*hatlamSq(j));
    }
    meff = arma::accu(1- kappa);
    gsMeff(k) = meff;
    gsKappa.row(k) = kappa.t();
    
   // Rcpp::Rcout << "tau2" << std::endl;
    double shapeS = (s+1)/2;
    double rateS = 1/hatZeta + arma::accu(square(hatBeta)%(1/hatlamSq))/(2*hatSigmaSq);
    hatTauSq = 1/R::rgamma(shapeS, 1/rateS);
    gsTauSq(k) = hatTauSq;
    
    //nu_1
    //Rcpp::Rcout << "nu1" << std::endl;
    for(unsigned int j = 0; j<s; j++){
      double ratenu_1 = 1/hatetaSq(j)+1/hatlamSq(j);
      hatnu_1(j) = 1/R::rgamma(1, 1/ratenu_1);
    }
    gsnu_1.row(k) = hatnu_1.t();
    
    //nu_2
    //Rcpp::Rcout << "nu2" << std::endl;
    for(unsigned int j = 0; j<s; j++){
      double ratenu_2 = 1+1/hatetaSq(j);
      hatnu_2(j) = 1/R::rgamma(1, 1/ratenu_2);
    }
    
    //hatetaSq
    //Rcpp::Rcout << "eta2" << std::endl;
    for(unsigned int j = 0; j<s; j++){
      double rateEta = 1/hatnu_1(j)+1/hatnu_2(j);
      hatetaSq(j) = 1/R::rgamma(1, 1/rateEta);
    }
    gsetaSq.row(k) = hatetaSq.t();
    
    // zeta|
    //Rcpp::Rcout << "zeta" << std::endl;
    double ratezeta = 1+1/hatTauSq;
    hatZeta = 1/R::rgamma(1, 1/ratezeta);
    gsZeta(k) = hatZeta;
    
    // MSE
    gsMSE(k) = arma::mean(arma::square(res));
    if(k % 100 == 0){
      Rcpp::checkUserInterrupt();
    }
    if(progress != 0 && k % progress == 0){
      Rcpp::Rcout << "\nIter." << k << "  mse: " << gsMSE(k) << std::endl;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("GS.alpha") = gsAlpha,
                            Rcpp::Named("GS.beta") = gsBeta,
                            Rcpp::Named("GS.TAUsq") = gsTauSq,
                            Rcpp::Named("GS.lambda.sq") = gslamSq,
                            Rcpp::Named("GS.sigma.sq") = gsSigmaSq,
                            Rcpp::Named("GS.kappa") = gsKappa,
                            Rcpp::Named("GS.meff") = gsMeff,
                            Rcpp::Named("GS.mse") = gsMSE);
}
