
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pqrBayes

> Bayesian Penalized Quantile Regression

<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/pqrBayes)](https://cran.r-project.org/package=pqrBayes)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/pqrBayes)](https://www.r-pkg.org:443/pkg/pqrBayes)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/last-month/pqrBayes)](https://www.r-pkg.org:443/pkg/pqrBayes)

<!-- badges: end -->


Bayesian regularized quantile regression utilizing sparse priors to 
    impose exact sparsity leads to efficient Bayesian shrinkage estimation, variable 
    selection and statistical inference. In this package, we have implemented robust
    Bayesian variable selection with spike-and-slab priors under high-dimensional
    linear regression models ([Fan et al. (2024)](https://doi.org/10.3390/e26090794) and 
    [Ren et al. (2023)](https://doi.org/10.1111/biom.13670), and regularized quantile varying
    coefficient models ([Zhou et al.(2023)](https://doi.org/10.1016/j.csda.2023.107808)). In particular, 
    valid robust Bayesian inferences under both models in the presence of heavy-tailed errors
    can be validated on finite samples. Additional models including robust Bayesian 
    group LASSO and robust Bayesian binary LASSO are also included. The Markov Chain Monte Carlo (MCMC) algorithms 
    of the proposed and alternative models are implemented in C++. 
   

## How to install

 - To install from github, run these two lines of code in R

<!-- end list -->

    install.packages("devtools")
    devtools::install_github("cenwu/pqrBayes")

- Released versions of pqrBayes are available on CRAN
    <!-- [(link)](https://cran.r-project.org/package=pqrBayes) --> ,
    and can be installed within R via

<!-- end list -->

    install.packages("pqrBayes")

## Example 1 (Robust Bayesian Inference for Sparse Linear Regression)

#### Data Generation for Linear Model 

    Data <- function(n,p,quant){
      sig1 = matrix(0,p,p)
      diag(sig1)=1
      for (i in 1: p)
      {
      for (j in 1: p)
      {
      sig1[i,j]=0.5^abs(i-j)
      }
     }
    xx = MASS::mvrnorm(n,rep(0,p),sig1)
    x = cbind(1,xx)
    error=rt(n,2) -quantile(rt(n,2),probs = quant) # can also be changed to normal error for non-robust setting
    beta = c(0,1,1.5,2,rep(0,p-3))
    betaa = beta[-1]
    y = x%*%beta+error
    dat = list(y=y, x=xx, beta=betaa)
    return(dat)
    }
#### 95% empirical coverage probabilities for linear regression coefficients

    n=100; p=500; rep=1000;
    quant = 0.5; # focus on median for Bayesian inference

    CI_RBLSS = CI_RBL = CI_BLSS = CI_BL= matrix(0,rep,p)

    for (h in 1:rep) {
    dat = Data(n,p,quant)
    y = dat$y
    g = dat$x
    coefficient = dat$beta

    # an intercept is automatically included by the package

    fit = pqrBayes(g, y, u=NULL, d=NULL,e=NULL,quant=quant, iterations=10000, burn.in = NULL, spline = NULL,robust = TRUE, sparse=TRUE, model = "linear", hyper=NULL,debugging=FALSE)

    coverage = coverage(fit,coefficient,u.grid=NULL, model = "linear")

    fit1 = pqrBayes(g, y, u=NULL,d=NULL, e=NULL,quant=quant, iterations=10000, burn.in = NULL, spline = NULL, robust = TRUE, sparse=FALSE, model = "linear", hyper=NULL,debugging=FALSE)

    coverage1 = coverage(fit1,coefficient,u.grid=NULL, model = "linear")
  
    fit2 = pqrBayes(g, y, u=NULL,d=NULL, e=NULL,quant=quant, iterations=10000, burn.in = NULL, spline = NULL, robust = FALSE, sparse=TRUE, model= "linear", hyper=NULL,debugging=FALSE)
   
    coverage2 = coverage(fit2,coefficient,u.grid=NULL, model = "linear")
  
    fit3 = pqrBayes(g, y, u=NULL,d=NULL, e=NULL,quant=quant, iterations=10000, burn.in = NULL, spline = NULL, robust = FALSE, sparse=FALSE, model = "linear", hyper=NULL,debugging=FALSE)
    
    coverage3 = coverage(fit3,coefficient,u.grid=NULL, model = "linear")
  
    CI_RBLSS[h,] = coverage
    CI_RBL[h,]   = coverage1
    CI_BLSS[h,]  = coverage2
    CI_BL[h,]    = coverage3
    cat("Replicate = ", h, "\n")
    
    }
    # the intercept has not been regularized
    cp_RBLSS =  colMeans(CI_RBLSS)[1:3] # 95% empirical coverage probabilities for coefficients under the robust linear model
    cp_BLSS  =  colMeans(CI_BLSS)[1:3]
    cp_RBL   =  colMeans(CI_RBL)[1:3]
    cp_BL    =  colMeans(CI_BL)[1:3]

## Example 2 (Robust Bayesian Inference for Sparse Varying Coefficients)
#### Data Generation for the Varying Coefficient Model 

    Data <- function(n,p,quant){
      sig1 = matrix(0,p,p)
      diag(sig1)=1
      for (i in 1: p)
      {
      for (j in 1: p)
      {
      sig1[i,j]=0.5^abs(i-j)
      }
     }
    x = MASS::mvrnorm(n,rep(0,p),sig1)
    x = cbind(1,x)
    error=rt(n,2) -quantile(rt(n,2),probs = quant)
    u = runif(n,0.01,0.99)
    gamma0 = 2+2*sin(u*2*pi)
    gamma2 = -6*u*(1-u)
    gamma1 = 2*exp(2*u-1)
    gamma3= -4*u^3
    y = gamma1*x[,2] + gamma2*x[,3]  + gamma3*x[,4] + gamma0 + error
    dat = list(y=y, u=u, x=x, gamma=cbind(gamma0,gamma1,gamma2,gamma3))
    return(dat)
    }
#### 95% empirical coverage probabilities for sparse varying coefficients

    n=250; p=100; # the actual dimension after basis expansion is 505
    rep=200;
    quant = 0.5; # focus on median for Bayesian inference

    CI_RBGLSS = CI_RBGL = CI_BGLSS = CI_BGL= c()

    for (h in 1:rep) {
    dat = Data(n,p,quant)
    y = dat$y
    u = dat$u
    x = dat$x
    g = x[,-1]
    kn=2
    degree=2
    u.grid = (1:200)*0.005
    gamma_0_grid = 2+2*sin(2*u.grid*pi)
    gamma_1_grid = 2*exp(2*u.grid-1)
    gamma_2_grid = -6*u.grid*(1-u.grid)
    gamma_3_grid = -4*u.grid^3
    coefficient = cbind(gamma_0_grid,gamma_1_grid,gamma_2_grid,gamma_3_grid)

    # a varying intercept is automatically included by the package

    fit = pqrBayes(g, y, u, d=NULL,e=NULL,quant=quant, iterations=10000, burn.in = NULL, spline = list(kn=2,degree=2), robust = TRUE, sparse=TRUE, model = "VC", hyper=NULL,debugging=FALSE)
    
    coverage = coverage(fit,coefficient,u.grid, model = "VC")

    fit1 = pqrBayes(g, y, u, d=NULL,e=NULL,quant=quant, iterations=10000, burn.in = NULL, spline = list(kn=2, degree=2), robust = TRUE, sparse=FALSE, model = "VC", hyper=NULL,debugging=FALSE)
    
    coverage1 = coverage(fit1,coefficient,u.grid, model = "VC")
  
    fit2 = pqrBayes(g, y, u, d = NULL,e=NULL,quant=quant, iterations=10000, burn.in = NULL, spline = list(kn=2, degree=2), robust = FALSE, sparse=TRUE, model = "VC", hyper=NULL,debugging=FALSE)
    
    coverage2 = coverage(fit2,coefficient,u.grid, model = "VC")
  
    fit3 = pqrBayes(g, y, u, d=NULL,e=NULL,quant=quant, iterations=10000, burn.in = NULL, spline = list(kn=2, degree=2), robust = FALSE, sparse=FALSE, model = "VC", hyper=NULL,debugging=FALSE)
   
    coverage3 = coverage(fit3,coefficient,u.grid,model = "VC")
  
    CI_RBGLSS = rbind(CI_RBGLSS,coverage)
    CI_RBGL   = rbind(CI_RBGL,coverage1)
    CI_BGLSS  = rbind(CI_BGLSS,coverage2)
    CI_BGL    = rbind(CI_BGL,coverage3)
    cat("Replicate = ", h, "\n")
    
    }
    # the varying intercept has not been regularized
    cp_RBGLSS =  colMeans(CI_RBGLSS) # 95% empirical coverage probabilities for the varying coefficients under the default setting
    cp_BGLSS  =  colMeans(CI_BGLSS)
    cp_RBGL   =  colMeans(CI_RBGL)
    cp_BGL    =  colMeans(CI_BGL)

## Example 3 (Bayesian Shrinkage Estimation for Robust Bayesian Group LASSO)

#### Data Generation for Group LASSO 

    Data <- function(n,p,quant){
      sig1 = matrix(0,p,p)
      diag(sig1)=1
      for (i in 1: p)
      {
      for (j in 1: p)
      {
      sig1[i,j]=0.5^abs(i-j)
      }
     }
    xx = MASS::mvrnorm(n,rep(0,p),sig1)
    x = cbind(1,xx)
    error=rt(n,2) -quantile(rt(n,2),probs = quant) # can also be changed to normal error for non-robust setting
    beta = c(0,1,1.5,2,0,0,0,0.5,0.55,0.6,rep(0,p-9))
    betaa = beta[-1]
    y = x%*%beta+error
    dat = list(y=y, x=xx, beta=betaa)
    return(dat)
    }
#### robust Bayesian shrinkage estimation under group LASSO

    n=100; p=300;
    quant = 0.5; # focus on median for Bayesian estimation
    dat = Data(n,p,quant)
    y = dat$y
    g = dat$x
    coefficient = dat$beta

    # an intercept is automatically included by the package
    # the intercept has not been regularized
    fit = pqrBayes(g, y, u=NULL,d=3, e=NULL,quant=quant, iterations=10000, burn.in = NULL, spline = NULL,robust = TRUE, sparse=TRUE, model = "group", hyper=NULL,debugging=FALSE)
    estimation_1 = estimation.pqrBayes(fit,coefficient,model="group")
    coeff_est_1 = estimation_1$coeff.est    
    mse_1 = estimation_1$error$MSE
    
    fit1 = pqrBayes(g, y, u=NULL,d=3, e=NULL,quant=quant, iterations=10000, burn.in = NULL, spline = NULL, robust = TRUE, sparse=FALSE, model = "group", hyper=NULL,debugging=FALSE)
    estimation_2 = estimation.pqrBayes(fit1,coefficient,model="group")
    coeff_est_2 = estimation_2$coeff.est    
    mse_2 = estimation_2$error$MSE 
    
    fit2 = pqrBayes(g, y, u=NULL,d=3, e=NULL,quant=quant, iterations=10000, burn.in = NULL, spline = NULL, robust = FALSE, sparse=TRUE, model= "group", hyper=NULL,debugging=FALSE)
    estimation_3 = estimation.pqrBayes(fit2,coefficient,model="group")
    coeff_est_3 = estimation_3$coeff.est    
    mse_3 = estimation_3$error$MSE    
    
    fit3 = pqrBayes(g, y, u=NULL,d=3, e=NULL,quant=quant, iterations=10000, burn.in = NULL, spline = NULL, robust = FALSE, sparse=FALSE, model = "group", hyper=NULL,debugging=FALSE)
    estimation_4 = estimation.pqrBayes(fit3,coefficient,model="group")
    coeff_est_4 = estimation_4$coeff.est    
    mse_4 = estimation_4$error$MSE    

## Methods

This package provides implementation for methods from
  
  - Fan, K., Subedi, S., Yang, G., Lu, X., Ren, J. and Wu, C. (2024). Is Seeing Believing? A Practitioner's Perspective on High-dimensional Statistical Inference in Cancer Genomics Studies. [Entropy, 26(9),794](https://doi.org/10.3390/e26090794)
  - Zhou, F., Ren, J., Ma, S. and Wu, C. (2023). The Bayesian Regularized Quantile Varying Coefficient Model.  [Computational Statistics & Data Analysis, 107808](https://doi.org/10.1016/j.csda.2023.107808)
  - Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y., and Wu, C. (2023). Robust Bayesian variable selection for gene–environment interactions. [Biometrics, 79(2), 684-694](https://doi.org/10.1111/biom.13670)
