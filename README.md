
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pqrBayes

> Bayesian Penalized Quantile Regression

<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/pqrBayes)](https://cran.r-project.org/package=pqrBayes)
[![Codecov test
coverage](https://codecov.io/gh/cenwu/pqrBayes/branch/master/graph/badge.svg)](https://app.codecov.io/gh/cenwu/pqrBayes?branch=master)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/pqrBayes)](https://www.r-pkg.org:443/pkg/pqrBayes)
[![R-CMD-check](https://github.com/cenwu/pqrBayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cenwu/pqrBayes/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->


The quantile varying coefficient model is robust to data heterogeneity, 
    outliers and heavy-tailed distributions in the response variable. In addition, 
    it can flexibly model dynamic patterns of regression coefficients through 
    nonparametric varying coefficient functions. In this package, we have implemented 
    the Gibbs samplers of the penalized Bayesian quantile varying coefficient model with 
    spike-and-slab priors [Zhou et al.(2023)]<doi:10.1016/j.csda.2023.107808> for efficient 
    Bayesian shrinkage estimation, variable selection and statistical inference. In particular,
    valid Bayesian inferences on sparse quantile varying coefficient functions can be validated 
    on finite samples. The Markov Chain Monte Carlo (MCMC) algorithms of the proposed
    and alternative models can be efficiently performed by using the package.   

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

## Example
#### Data Generation    
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
    x = mvrnorm(n,rep(0,p),sig1)
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
#### 95% coverage probability for the varying coefficients
    n=250; p=100; rep=200;
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

    fit = pqrBayes(g, y, u, e=NULL,quant=quant, iterations=10000, kn=2, degree=2, robust = TRUE, sparse=TRUE, hyper=NULL,debugging=FALSE)
    posterior=fit$posterior
    coverage = coverage(posterior,coefficient,u.grid,kn,degree)

    fit1 = pqrBayes(g, y, u, e=NULL,quant=quant, iterations=10000, kn=2, degree=2, robust = TRUE, sparse=FALSE, hyper=NULL,debugging=FALSE)
    posterior1=fit1$posterior
    coverage1 = coverage(posterior1,coefficient,u.grid,kn,degree)
  
    fit2 = pqrBayes(g, y, u, e=NULL,quant=quant, iterations=10000, kn=2, degree=2, robust = FALSE, sparse=TRUE, hyper=NULL,debugging=FALSE)
    posterior2=fit2$posterior
    coverage2 = coverage(posterior2,coefficient,u.grid,kn,degree)
  
    fit3 = pqrBayes(g, y, u, e=NULL,quant=quant, iterations=10000, kn=2, degree=2, robust = FALSE, sparse=FALSE, hyper=NULL,debugging=FALSE)
    posterior3=fit3$posterior
    coverage3 = coverage(posterior3,coefficient,u.grid,kn,degree)
  
    CI_RBGLSS = rbind(CI_RBGLSS,coverage)
    CI_RBGL   = rbind(CI_RBGL,coverage1)
    CI_BGLSS  = rbind(CI_BGLSS,coverage2)
    CI_BGL    = rbind(CI_BGL,coverage3)
    cat("iteration = ", h, "\n")
    
    }
    # the intercept gamma_0 has not been regularized
    cp_RBGLSS =  colMeans(CI_RBGLSS) # 95% coverage probabilities for the varying coeffcients under the default setting
    cp_BGLSS  =  colMeans(CI_BGLSS)
    cp_RBGL   =  colMeans(CI_RBGL)
    cp_BGL    =  colMeans(CI_BGL)

## News

### pqrBayes 1.0.4 \[2025-01-24\]

- Added functions to compute the coverage probability of varying coefficients.
- Added examples to obtain the inference results reported in the reference.

### pqrBayes 1.0.3 \[2024-12-20\]

- Fixed the issue of no output in the VCselect() function and added examples.
- Updated the list of pqrBayes output objects.
- Added non-robust methods and corresponding examples.
- Updated the documentation.

## Methods

This package provides implementation for methods proposed in

  - Zhou, F., Ren, J., Ma, S. and Wu, C. (2023). The Bayesian Regularized Quantile Varying Coefficient Model.  {\emph{Computational Statistics & Data Analysis}, 107808} \doi{10.1016/j.csda.2023.107808}
