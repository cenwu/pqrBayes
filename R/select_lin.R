linselect <- function(obj,sparse){
  iterations = obj$obj$iterations
  if(sparse){
    method="Sparse"
    mpm <- function(x)
    {
      if (mean(x) >= 0.5) {1}
      else {0}
    }
    sg1 = obj$coefficients$GS.beta
    sg1[which(sg1!=0)]=1
    id=c()
    for(j in 1:ncol(sg1)){
      t1=as.matrix(sg1[,j])
      t1 = t1[seq(iterations/2+1, iterations,1),]
      q_t1 = mpm(t1)
      id = c(id,q_t1)
    }
    linselect=list(method=method,id=id)
  }else{
    method="Nonsparse"
    fun <- function(x)
    {
      pp = prod(x)
      if(sign(pp)==1) {1}
      else {0}
    }
    sg1 = obj$coefficients$GS.beta
    q_t1=c()
    for(j in 1:ncol(sg1)){
      t1=as.matrix(sg1[,j])
      t1 = t1[seq(iterations/2+1, iterations,1),]
      q_t1 = as.matrix(stats::quantile(t1,c(0.025,0.975)))
    }
    id = apply(q_t1, 2, fun)
    
    linselect=list(method=method,id=id)
  }
  
  #class(VCselect)="VCselect"
  return(linselect)
  
}