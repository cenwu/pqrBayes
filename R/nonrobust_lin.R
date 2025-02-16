NonRobust_lin <- function(g, y, e,iterations,sparse,debugging){
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
  hatLambdaSqStar=1
  hatSigmaSq=1
  aStar=1 
  bStar=1
  alpha=1 
  gamma=1
  hatPi=0.5
  sh1=1
  sh0=1
  progress=0
  
  progress = ifelse(debugging, 10^(floor(log10(iterations))-1), 0)
  if(sparse){fit=BLSS(g, y, w, iterations, hatAlpha, hatBeta, hatInvTauSq, invSigAlpha0, hatPi, hatLambdaSq, hatSigmaSq, aStar, bStar, alpha, gamma, sh1, sh0, progress)}
  else{fit=BL(g, y, w, iterations, hatBeta, hatAlpha, hatInvTauSq, invSigAlpha0, hatLambdaSqStar, hatSigmaSq, aStar, bStar, alpha, gamma, progress)}
  out=list(fit=fit,iterations=iterations)
  return(out)
}
