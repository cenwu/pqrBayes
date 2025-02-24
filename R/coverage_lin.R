coverage_lin = function(obj,coefficient,u.grid=NULL){
  C = c()
  for(j in 1:dim(obj$coefficients$GS.beta)[2])
  {
    t1=as.matrix(obj$coefficients$GS.beta[,j])
    c_j = as.matrix(stats::quantile(t1,c(0.025,0.975)))
    C[j] = ifelse(coefficient[j]>=c_j[1]&coefficient[j]<=c_j[2],1,0)
  }
  
  return(C)
}
