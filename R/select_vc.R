VCselect <- function(obj,sparse){
  if(sparse){
    method="Sparse"
    idgene=obj$obj$idgene
    iterations = obj$obj$iterations
    id=which(idgene>iterations/4)
    VCselect=list(method=method,id=id,idgene = idgene)
  }else{
    method="Nonsparse"
    kn = obj$obj$kn
    degree = obj$obj$degree
    d=kn+degree+1
    iterations = obj$obj$iterations
    # 2.5% quantile of posterior samples
    
    c2.C=rep(0,dim(obj$coefficients$GS.beta)[2])
    for (i in 1:dim(obj$coefficients$GS.beta)[2]) {
      c2.C[i]=stats::quantile(obj$coefficients$GS.beta[(iterations/2+1):iterations,i],0.025)
    }
    
    coeffmatrix.C1=matrix(c2.C,nrow = d)
    
    # 97.5% quantile of posterior samples
    
    c2.C=rep(0,dim(obj$coefficients$GS.beta)[2])
    for (i in 1:dim(obj$coefficients$GS.beta)[2]) {
      c2.C[i]=stats::quantile(obj$coefficients$GS.beta[(iterations/2+1):iterations,i],0.975)
    }
    
    coeffmatrix.C2=matrix(c2.C,nrow = d)
    
    idgene=rep(0,dim(coeffmatrix.C1)[2])
    for (i in 1:dim(coeffmatrix.C1)[2]) {
      for (j in 1:d) {
        if(coeffmatrix.C1[j,i]*coeffmatrix.C2[j,i]>0){idgene[i]=1}
      }
    }
    id=which(idgene==1)
    VCselect=list(method=method,id=id,idgene = idgene)
  }
  
  #class(VCselect)="VCselect"
  return(VCselect)
  
}