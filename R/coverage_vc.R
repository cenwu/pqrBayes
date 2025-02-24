coverage_vc = function(obj,coefficient,u.grid){
  kn = obj$obj$kn
  degree = obj$obj$degree
  d=kn+degree+1
  u.star = seq(0, 1, length=kn+2)[-c(1,kn+2)]
  Knots.star = as.numeric(stats::quantile(u.grid, u.star))
  pi.star = splines::bs(u.grid, knots=Knots.star, intercept=TRUE, degree=degree)[,1:(d)]
  # 2.5% quantile of posterior samples
  
  c1_25.C=rep(0,d)
  for (i in 1:d) {
    c1_25.C[i]=stats::quantile(obj$coefficients$GS.alpha[,i],0.025)
  }
  
  c2_25.C=rep(0,dim(obj$coefficients$GS.beta)[2])
  for (i in 1:dim(obj$coefficients$GS.beta)[2]) {
    c2_25.C[i]=stats::quantile(obj$coefficients$GS.beta[,i],0.025)
  }
  
  coeffmatrix.C1_25=as.matrix(cbind(c1_25.C,matrix(c2_25.C,nrow = d)))
  
  
  # 97.5% quantile of posterior samples
  
  c1_975.C=rep(0,d)
  for (i in 1:d) {
    c1_975.C[i]=stats::quantile(obj$coefficients$GS.alpha[,i],0.975)
  }
  
  c2_975.C=rep(0,dim(obj$coefficients$GS.beta)[2])
  for (i in 1:dim(obj$coefficients$GS.beta)[2]) {
    c2_975.C[i]=stats::quantile(obj$coefficients$GS.beta[,i],0.975)
  }
  
  coeffmatrix.C2_975=as.matrix(cbind(c1_975.C,matrix(c2_975.C,nrow = d)))
  gamma.var1 = pi.star%*% coeffmatrix.C1_25[,1:ncol(coefficient)]
  gamma.var2 = pi.star%*% coeffmatrix.C2_975[,1:ncol(coefficient)]
  CI = matrix(,nrow(coefficient),ncol(coefficient))
  for(i in 1:ncol(coefficient)){
    CI[,i] = ifelse(coefficient[,i]>=gamma.var1[,i]&coefficient[,i]<=gamma.var2[,i],1,0)
  }
  C = apply(CI, 2, function(col) ifelse(all(col == 1), 1, 0))
  return(C)
}
