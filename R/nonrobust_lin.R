NonRobust_lin <- function(g, y, e,iterations,prior,debugging){

  p = dim(g)[2]
  
  x = cbind(1,g)
  n = length(y); 

  if(!is.null(e)){
    q = dim(e)[2]
    hatAlpha=rep(1,(1+q))
    invSigAlpha0 = diag(10^-3, (1+q))
    w = cbind(1,e)
  }else{
    q = 0
    hatAlpha=rep(1,1)
    invSigAlpha0 = diag(10^-3, 1)
    w = matrix(1,n,1)
  }
  hatBeta = rep(1,p)
  hatInvTauSq = rep(1,p)
  hatLambdaSq=1
  hatSigmaSq=1
  aStar=1 
  bStar=1
  alpha=1 
  gamma=1
  hatPi=0.5
  sh1=1
  sh0=1
  hatTauSq=1
  hatcSq = 1
  hatSg = rep(1,p)
  hatlamSq = rep(1,p)
  hatZeta = 1
  hatnu = rep(1,p)
  hatnu_1 = rep(1,p)
  hatnu_2 = rep(1,p)
  hatetaSq = rep(1,p)
  hatSigma = 1
  hatSigmaSq = 1
  c = 0.002
  d = 90
  progress=0
  
  progress = ifelse(debugging, 10^(floor(log10(iterations))-1), 0)
  if(prior == "SS"){
    fit=BLSS(g, y, w, iterations, hatAlpha, hatBeta, hatInvTauSq, invSigAlpha0, hatPi, hatLambdaSq, hatSigmaSq, aStar, bStar, alpha, gamma, sh1, sh0, progress)}
  else if (prior =="HS"){
    fit = Bhorse(g, y, w, iterations, hatBeta, hatAlpha, hatlamSq, hatnu,invSigAlpha0, hatTauSq, hatSigmaSq, hatZeta, aStar, bStar, progress)
  }
  else if (prior =="HS+" ){
    fit = Bhorse_plus(g, y, w, iterations, hatBeta, hatAlpha, hatlamSq, hatetaSq, hatnu_1, hatnu_2, invSigAlpha0, hatTauSq, hatSigmaSq, hatZeta, aStar, bStar, progress)
  }
  else if (prior =="RHS"){
    fit = Bhorse_reg(g, y, w, iterations, hatBeta, hatAlpha, hatlamSq, hatnu,invSigAlpha0, hatTauSq, hatcSq, hatSigmaSq, hatZeta, aStar, bStar, c, d, progress)
  }
  else if (prior == "Laplace"){fit=BL(g, y, w, iterations, hatBeta, hatAlpha, hatInvTauSq, invSigAlpha0, hatLambdaSq, hatSigmaSq, aStar, bStar, alpha, gamma, progress)}
  else{stop("The specified prior is currently not supported.")}
  out=list(fit=fit,iterations=iterations)
  return(out)
}
