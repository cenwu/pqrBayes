coverage_lin = function(obj,coefficient,u.grid=NULL){
  iterations = obj$obj$iterations
  C = c()
  for(j in 1:dim(obj$coefficients$GS.beta)[2])
  {
    t1=as.matrix(obj$coefficients$GS.beta[,j])
    t1 = t1[seq(iterations/2+1, iterations,1),]
    c_j = as.matrix(stats::quantile(t1,c(0.025,0.975)))
    C[j] = ifelse(coefficient[j]>=c_j[1]&coefficient[j]<=c_j[2],1,0)
  }
  
  return(C)
}
