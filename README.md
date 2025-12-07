
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


Bayesian regularized quantile regression utilizing two major classes of shrinkage priors 
    (the spike-and-slab priors and the horseshoe family of priors) leads to efficient Bayesian 
    shrinkage estimation, variable selection and valid statistical inference. In this package, we have implemented robust
    Bayesian variable selection with spike-and-slab priors under high-dimensional
    linear regression models ([Fan et al. (2024)](https://doi.org/10.3390/e26090794) and 
    [Ren et al. (2023)](https://doi.org/10.1111/biom.13670)), and regularized quantile varying
    coefficient models ([Zhou et al.(2023)](https://doi.org/10.1016/j.csda.2023.107808)). In   
    particular, valid robust Bayesian inferences under both models in the presence of heavy-tailed
    errors can be validated on finite samples. Additional models with spike-and-slab priors 
    include robust Bayesian group LASSO and robust binary Bayesian LASSO ([Fan and Wu (2025)](https://doi.org/10.1002/sta4.70078)). 
    Besides, robust sparse Bayesian regression with the horseshoe family of (horseshoe, horseshoe+ and 
    regularized horseshoe) priors has also been implemented and yielded valid inference results under heavy-tailed model errors
    ([Fan et al.(2025)](https://doi.org/10.48550/arXiv.2507.10975)). The Markov Chain Monte Carlo (MCMC) 
    algorithms of the proposed and alternative models are implemented in C++. 
   

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
    
    
## Example 1 (Robust Bayesian Inference for Sparse Linear Regression with Spike-and-Slab Priors)

#### Data Generation for Sparse Linear Model 

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
#### 95% empirical coverage probabilities for sparse linear regression coefficients

    n=100; p=500; rep=1000;
    quant = 0.5; # focus on median for Bayesian inference
    
    CI_RBLSS = CI_RBL = CI_BLSS = CI_BL= matrix(0,rep,p)

    for (h in 1:rep) {
    dat = Data(n,p,quant)
    y = dat$y
    g = dat$x
    coefficient = dat$beta

    # an intercept not subject to regularization is automatically included by the package
    
    # RBLSS: robust Bayesian LASSO with spike-and-slab priors (Ren et al., Biometrics, 2023)
    fit = pqrBayes(g, y ,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL ,robust = TRUE, prior = "SS", model = "linear", hyper=NULL,debugging=FALSE)
    coverage = coverage(fit,coefficient,u.grid=NULL, model = "linear")
    
    # RBL: Bayesian quantile LASSO (Li, Xi & Lin, Bayesian Analysis, 2010)
    fit1 = pqrBayes(g, y,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = TRUE, prior = "Laplace", model = "linear", hyper=NULL,debugging=FALSE)
    coverage1 = coverage(fit1,coefficient,u.grid=NULL, model = "linear")
    
    # BLSS: Bayesian LASSO with spike-and-slab priors (Ren et al., Biometrics, 2023)
    fit2 = pqrBayes(g, y ,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL ,robust = FALSE, prior = "SS", model = "linear", hyper=NULL,debugging=FALSE)
    coverage2 = coverage(fit2,coefficient,u.grid=NULL, model = "linear")
    
    # Bayesian LASSO  (Park and Casella, JASA, 2008)
    fit3 = pqrBayes(g, y,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = FALSE, prior = "Laplace", model = "linear", hyper=NULL,debugging=FALSE)
    coverage3 = coverage(fit3,coefficient,u.grid=NULL, model = "linear")
    
    CI_RBLSS[h,] = coverage
    CI_RBL[h,]    = coverage1
    CI_BLSS[h,]   = coverage2
    CI_BL[h,]     = coverage3
    
    cat("Replicate = ", h, "\n")
    
    }
    # the intercept has not been regularized
    cp_RBLSS =  colMeans(CI_RBLSS)[1:3] # 95% empirical coverage probabilities for coefficients under the robust linear model
    cp_RBL    =  colMeans(CI_RBL)[1:3]
    cp_BLSS   =  colMeans(CI_BLSS)[1:3] # 95% empirical coverage probabilities for coefficients under the non-robust linear model
    cp_BL     =  colMeans(CI_BL)[1:3]
    
#### Bayesian shrinkage estimaton via robust Bayesian LASSO with spike-and-slab priors

    n=100; p=500;
    quant = 0.5; # focus on median for Bayesian estimation
    dat = Data(n,p,quant)
    y = dat$y
    g = dat$x
    coefficient = dat$beta

    # an intercept not subject to regularization is automatically included by the package

    # RBLSS: robust Bayesian LASSO with spike-and-slab priors (Ren et al., Biometrics, 2023)

    fit = pqrBayes(g, y ,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL ,robust = TRUE, prior = "SS", model = "linear", hyper=NULL,debugging=FALSE)

    fit$coefficients$GS.beta # posterior samples for regression coefficients
    estimation_1 = estimation.pqrBayes(fit,coefficient,model="linear")
    coeff_est_1 = estimation_1$coeff.est    
    mse_1 = estimation_1$error$MSE
    
    # RBL: Bayesian quantile LASSO (Li, Xi & Lin, Bayesian Analysis, 2010)

    fit1 = pqrBayes(g, y,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = TRUE, prior = "Laplace", model = "linear", hyper=NULL,debugging=FALSE)

    fit1$coefficients$GS.beta # posterior samples for regression coefficients
    estimation_2 = estimation.pqrBayes(fit1,coefficient,model="linear")
    coeff_est_2 = estimation_2$coeff.est    
    mse_2 = estimation_2$error$MSE
    
    # BLSS: Bayesian LASSO with spike-and-slab priors (Ren et al., Biometrics, 2023)
    
    fit2 = pqrBayes(g, y ,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL ,robust = FALSE, prior = "SS", model = "linear", hyper=NULL,debugging=FALSE)
    
    estimation_3 = estimation.pqrBayes(fit2,coefficient,model="linear")
    coeff_est_3 = estimation_3$coeff.est
    mse_3 = estimation_3$error$MSE
    
    # Bayesian LASSO  (Park and Casella, JASA, 2008)
    
    fit3 = pqrBayes(g, y,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = FALSE, prior = "Laplace", model = "linear", hyper=NULL,debugging=FALSE)
    
    estimation_4 = estimation.pqrBayes(fit3,coefficient,model="linear")
    coeff_est_4 = estimation_4$coeff.est
    mse_4 = estimation_4$error$MSE
    
## Example 2 (Robust Bayesian Inference for Sparse Linear Regression with the Horseshoe Family of Priors)

#### Data Generation for Sparse Linear Model 

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
    
    CI_RBHS = CI_RBHS_plus = CI_RBRHS = CI_BHS= CI_BHS_plus = CI_BRHS =matrix(0,rep,p)

    for (h in 1:rep) {
    dat = Data(n,p,quant)
    y = dat$y
    g = dat$x
    coefficient = dat$beta

    # an intercept not subject to regularization is automatically included by the package
    
    
    # RBHS: robust Bayesian Regression with horseshoe priors
    fit1 = pqrBayes(g, y ,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = TRUE, prior = "HS", model = "linear", hyper=NULL,debugging=FALSE)
    coverage1 = coverage(fit1,coefficient,u.grid=NULL, model = "linear")
    
    # RBHS+: robust Bayesian Regression with horseshoe plus priors
    fit2 = pqrBayes(g, y ,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = TRUE, prior = "HS+", model = "linear", hyper=NULL,debugging=FALSE)
    coverage2 = coverage(fit2,coefficient,u.grid=NULL, model = "linear")
    
    # RBRHS: robust Bayesian Regression with regularized horseshoe priors
    fit3 = pqrBayes(g, y ,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = TRUE, prior = "RHS", model = "linear", hyper=NULL,debugging=FALSE)
    coverage3 = coverage(fit3,coefficient,u.grid=NULL, model = "linear")
    
    # BHS: Bayesian Regression with horseshoe priors
    fit4 = pqrBayes(g, y ,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = FALSE, prior = "HS", model = "linear", hyper=NULL,debugging=FALSE)
    coverage4=coverage(fit4,coefficient,u.grid=NULL, model = "linear")
    
    # BHS+: Bayesian Regression with horseshoe plus priors
    fit5 = pqrBayes(g, y ,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = FALSE, prior = "HS+", model = "linear", hyper=NULL,debugging=FALSE)
    coverage5 = coverage(fit5,coefficient,u.grid=NULL, model = "linear")
    
    # BRHS: Bayesian Regression with regularized horseshoe priors
    fit6 = pqrBayes(g, y ,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = FALSE, prior = "RHS", model = "linear", hyper=NULL,debugging=FALSE)
    coverage6 = coverage(fit6,coefficient,u.grid=NULL, model = "linear")
    
    
    CI_RBHS[h,]   = coverage1
    CI_RBHS_plus[h,]  = coverage2
    CI_RBRHS[h,]   = coverage3
    CI_BHS[h,]    = coverage4
    CI_BHS_plus[h,]   = coverage5
    CI_BRHS[h,]  = coverage6
    
    cat("Replicate = ", h, "\n")
    
    }
    # the intercept has not been regularized
  
    cp_RBHS  =  colMeans(CI_RBHS)[1:3]
    cp_RBHS_plus   =  colMeans(CI_RBHS_plus)[1:3]
    cp_RBRHS    =  colMeans(CI_RBRHS)[1:3]
    cp_BHS    =  colMeans(CI_BHS)[1:3]
    cp_BHS_plus   =  colMeans(CI_BHS_plus)[1:3]
    cp_BRHS   =  colMeans(CI_BRHS)[1:3]

#### robust Bayesian shrinkage with the horseshoe family of priors

    n=100; p=500;
    quant = 0.5; # focus on median for Bayesian estimation
    dat = Data(n,p,quant)
    y = dat$y
    g = dat$x
    coefficient = dat$beta

    # an intercept not subject to regularization is automatically included by the package

    # RBHS: robust Bayesian Regression with horseshoe priors

    fit1 = pqrBayes(g, y ,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = TRUE, prior = "HS", model = "linear", hyper=NULL,debugging=FALSE)
    
    fit1$coefficients$GS.beta # posterior samples for regression coefficients
    estimation_1 = estimation.pqrBayes(fit1,coefficient,model="linear")
    coeff_est_1 = estimation_1$coeff.est    
    mse_1 = estimation_1$error$MSE

    # RBHS+: robust Bayesian Regression with horseshoe plus priors

    fit2 = pqrBayes(g, y ,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = TRUE, prior = "HS+", model = "linear", hyper=NULL,debugging=FALSE)

    fit2$coefficients$GS.beta # posterior samples for regression coefficients
    estimation_2 = estimation.pqrBayes(fit2,coefficient,model="linear")
    coeff_est_2 = estimation_2$coeff.est    
    mse_2 = estimation_2$error$MSE

    # RBRHS: robust Bayesian Regression with regularized horseshoe priors

    fit3 = pqrBayes(g, y ,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = TRUE, prior = "RHS", model = "linear", hyper=NULL,debugging=FALSE)
 
    fit3$coefficients$GS.beta # posterior samples for regression coefficients
    estimation_3 = estimation.pqrBayes(fit3,coefficient,model="linear")
    coeff_est_3 = estimation_3$coeff.est    
    mse_3 = estimation_3$error$MSE
    
    # BHS: Bayesian Regression with horseshoe priors
    fit4 = pqrBayes(g, y ,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = FALSE, prior = "HS", model = "linear", hyper=NULL,debugging=FALSE)
 
    estimation_4 = estimation.pqrBayes(fit4,coefficient,model="linear")
    coeff_est_4 = estimation_4$coeff.est
    mse_4 = estimation_4$error$MSE
    
    # BHS+: Bayesian Regression with horseshoe plus priors
    fit5 = pqrBayes(g, y ,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = FALSE, prior = "HS+", model = "linear", hyper=NULL,debugging=FALSE)
 
    estimation_5 = estimation.pqrBayes(fit5,coefficient,model="linear")
    coeff_est_5 = estimation_5$coeff.est
    mse_5 = estimation_5$error$MSE
    
    # BRHS: Bayesian Regression with regularized horseshoe priors
    fit6 = pqrBayes(g, y ,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = FALSE, prior = "RHS", model = "linear", hyper=NULL,debugging=FALSE)

    estimation_6 = estimation.pqrBayes(fit6,coefficient,model="linear")
    coeff_est_6 = estimation_6$coeff.est
    mse_6 = estimation_6$error$MSE
    

## Example 3 (Robust Bayesian Inference for Sparse Varying Coefficients)
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

    CI_BQRVCSS = CI_BQRVC = CI_BVCSS = CI_BVC= c()

    for (h in 1:rep) {
    dat = Data(n,p,quant)
    y = dat$y
    u = dat$u
    x = dat$x
    g = x[,-1]
    u.grid = (1:200)*0.005
    gamma_0_grid = 2+2*sin(2*u.grid*pi)
    gamma_1_grid = 2*exp(2*u.grid-1)
    gamma_2_grid = -6*u.grid*(1-u.grid)
    gamma_3_grid = -4*u.grid^3
    coefficient = cbind(gamma_0_grid,gamma_1_grid,gamma_2_grid,gamma_3_grid)

    # a varying intercept not subject to regularization is automatically included by the package
    
    # BQRVCSS: Bayesian regularized quantile VC model with spike-and-slab priors (Zhou et al., CSDA, 2023)

    fit = pqrBayes(g, y,e=NULL, d = NULL, quant=quant, iterations=10000, burn.in = NULL, robust = TRUE, prior="SS", model = "VC", hyper=NULL,debugging=FALSE)
    coverage = coverage(fit,coefficient,u.grid, model = "VC")
    
    # BQRVC: Bayesian regularized quantile VC model (Zhou et al., CSDA, 2023)

    fit1 = pqrBayes(g, y,e=NULL, d =NULL,quant=quant, iterations=10000, burn.in = NULL, robust = TRUE, prior = "Laplace", model = "VC", hyper=NULL,debugging=FALSE)
    coverage1 = coverage(fit1,coefficient,u.grid, model = "VC")
    
    # BVCSS: Bayesian regularized VC model with spike-and-slab priors (Zhou et al., CSDA, 2023)
  
    fit2 = pqrBayes(g, y,e=NULL, d =NULL, quant=quant, iterations=10000, burn.in = NULL, robust = FALSE, prior = "SS", model = "VC", hyper=NULL,debugging=FALSE)
    coverage2 = coverage(fit2,coefficient,u.grid, model = "VC")
    
    # BVC: Bayesian regularized VC model (Zhou et al., CSDA, 2023)
  
    fit3 = pqrBayes(g, y,e=NULL, d =NULL, quant=quant, iterations=10000, burn.in = NULL, robust = FALSE, prior = "Laplace", model = "VC", hyper=NULL,debugging=FALSE)
    coverage3 = coverage(fit3,coefficient,u.grid,model = "VC")
    
    CI_BQRVCSS = rbind(CI_BQRVCSS,coverage)
    CI_BQRVC   = rbind(CI_BQRVC,coverage1)
    CI_BVCSS   = rbind(CI_BVCSS,coverage2)
    CI_BVC     = rbind(CI_BVC,coverage3)
    cat("Replicate = ", h, "\n")
    
    }
    # the varying intercept has not been regularized
    cp_BQRVCSS =  colMeans(CI_BQRVCSS) # 95% empirical coverage probabilities for the varying coefficients under the default setting
    cp_BQRVC   =  colMeans(CI_BQRVC)
    cp_BVCSS   =  colMeans(CI_BVCSS)
    cp_BVC     =  colMeans(CI_BVC)

## Example 4 (Bayesian Shrinkage Estimation for Robust Bayesian Group LASSO)

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
#### robust Bayesian shrinkage with spike and slab prior estimation under group LASSO

    n=100; p=300;
    quant = 0.5; # focus on median for Bayesian estimation
    dat = Data(n,p,quant)
    y = dat$y
    g = dat$x
    coefficient = dat$beta

    # an intercept not subject to regularization is automatically included by the package

    # RBGLSS, (Ren et al. Biometrics, 2023)
    fit = pqrBayes(g, y, e=NULL, d=3, quant=quant, iterations=10000, burn.in = NULL,robust = TRUE, prior = "SS", model = "group", hyper=NULL,debugging=FALSE)
    
    fit$coefficients$GS.beta # posterior samples for regression coefficients 
    estimation_1 = estimation.pqrBayes(fit,coefficient,model="group")
    coeff_est_1 = estimation_1$coeff.est    
    mse_1 = estimation_1$error$MSE
    
    # RBGL: Bayesian Quantile Group LASSO (Li, Xi, & Lin, Bayesian Analysis, 2010)
    fit1 = pqrBayes(g, y, e=NULL,d=3, quant=quant, iterations=10000, burn.in = NULL, robust = TRUE, prior = "Laplace", model = "group", hyper=NULL,debugging=FALSE)
    
    fit1$coefficients$GS.beta # posterior samples for regression coefficients
    estimation_2 = estimation.pqrBayes(fit1,coefficient,model="group")
    coeff_est_2 = estimation_2$coeff.est    
    mse_2 = estimation_2$error$MSE 
    
    # BGLSS: Bayesian group LASSO with spike-and-slab priors (Xu & Ghosh, Bayesian Analysis, 2015) 
    fit2 = pqrBayes(g, y, d=3, e=NULL,quant=quant, iterations=10000, burn.in = NULL, robust = FALSE, prior = "SS", model= "group", hyper=NULL,debugging=FALSE)
    estimation_3 = estimation.pqrBayes(fit2,coefficient,model="group")
    coeff_est_3 = estimation_3$coeff.est    
    mse_3 = estimation_3$error$MSE    
    
    # BGL: Bayesian group LASSO (Casella et al., Bayesian Analysis, 2010)
    fit3 = pqrBayes(g, y,d=3, e=NULL,quant=quant, iterations=10000, burn.in = NULL, robust = FALSE, prior = "Laplace", model = "group", hyper=NULL,debugging=FALSE)
    estimation_4 = estimation.pqrBayes(fit3,coefficient,model="group")
    coeff_est_4 = estimation_4$coeff.est    
    mse_4 = estimation_4$error$MSE    

## Methods

This package provides implementation for methods from
  
  - Fan, K., Subedi, S., Yang, G., Lu, X., Ren, J. and Wu, C. (2024). Is Seeing Believing? A Practitioner's Perspective on High-dimensional Statistical Inference in Cancer Genomics Studies. [Entropy, 26(9),794](https://doi.org/10.3390/e26090794)
  - Zhou, F., Ren, J., Ma, S. and Wu, C. (2023). The Bayesian Regularized Quantile Varying Coefficient Model.  [Computational Statistics & Data Analysis, 107808](https://doi.org/10.1016/j.csda.2023.107808)
  - Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y., and Wu, C. (2023). Robust Bayesian variable selection for gene–environment interactions. [Biometrics, 79(2), 684-694](https://doi.org/10.1111/biom.13670)
  - Fan, K. and Wu, C. (2025). A New Robust Binary Bayesian LASSO. [Stat, 14(3), e70078.](https://doi.org/10.1002/sta4.70078)
  - Fan, K., Srijana, S., Dissanayake, V. and Wu, C. (2025). Robust Bayesian high-dimensional variable selection and inference with the horseshoe family of priors. [arXiv:2507.10975](https://doi.org/10.48550/arXiv.2507.10975)
