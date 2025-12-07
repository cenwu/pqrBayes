Robust_lin <- function(g, y, e, quant, iterations,prior, hyper,debugging){
 
  p = dim(g)[2]
  
  x = cbind(1,g)
  n = length(y); 
  
  if(!is.null(e)){
    q = dim(e)[2]
    invSigAlpha0 = diag(10^-3, (1+q))
    hatAlpha = rep(1,(1+q))
    w = cbind(1,e)
  }else{
    q = 0
    hatAlpha=rep(1,1)
    invSigAlpha0 = diag(10^-3, 1)
    w = matrix(1,n,1)
  }
  
  xi1 = (1-2*quant)/(quant*(1-quant))
  xi2 = sqrt(2/(quant*(1-quant)))
  hatBeta=rep(1,p)
  hatSg=rep(1,p)
  hatTau = 1
  hatV = rep(1,n) 
  hatEtaSq = 1
  hatPi = 0.5
  
  sh0_1 = ifelse(is.null(hyper$a0), 1, hyper$a0)
  sh0_0 = ifelse(is.null(hyper$b0), 1, hyper$b0)
  
  a = ifelse(is.null(hyper$c1), 1, hyper$c1)
  b = ifelse(is.null(hyper$c2), 1, hyper$c2)
  
  r = ifelse(is.null(hyper$d2), 1, hyper$d2)
  hatTauSq=1
  hatcSq = 1
  hatlamSq = rep(1,p)
  hatZeta = 1
  hatnu = rep(1,p)
  hatnu_1 = rep(1,p)
  hatnu_2 = rep(1,p)
  hatetaSq = rep(1,p)
  hatSigma = 1
  hatSigmaSq = 1
  progress=0
  c = 0.002
  d = 90
  progress = ifelse(debugging, 10^(floor(log10(iterations))-1), 0)
  if(prior == "SS"){
    fit=BRLSS(g, y, w, iterations, hatAlpha, hatBeta, hatTau, hatV, hatSg, invSigAlpha0, hatPi, hatEtaSq, xi1, xi2, r, a, b,sh0_1,sh0_0, progress)}
  else if (prior =="HS"){
    fit = RBhorse (g, y, w, iterations, hatAlpha, hatBeta, hatSigma, hatV, hatSg, hatnu, invSigAlpha0, hatTauSq, xi1, xi2, hatZeta, a, b, progress)
  }
  else if (prior =="HS+" ){
    fit = RBhorse_plus(g, y, w, iterations, hatAlpha, hatBeta, hatSigma, hatV, hatSg, hatetaSq, hatnu_1,hatnu_2, invSigAlpha0, hatTauSq, xi1, xi2, hatZeta, a, b, progress)
  }
  else if (prior =="RHS"){
    fit = RBhorse_reg (g, y, w, iterations, hatAlpha, hatBeta, hatSigma, hatcSq, hatV, hatSg, hatnu, invSigAlpha0, hatTauSq, xi1, xi2, hatZeta, a, b, c, d, progress)
  }
  else if(prior == "Laplace"){fit= BRL(g, y, w, iterations, hatAlpha, hatBeta, hatTau, hatV, hatSg, invSigAlpha0, hatEtaSq, xi1, xi2, r, a, b, progress)}
  else{stop("The specified prior is currently not supported.")}
  out = list(fit=fit,iterations=iterations)
  return(out)
}
