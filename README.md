
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pqrBayes

> Bayesian Penalized Quantile Regression

The quantile varying coefficient model is robust to data heterogeneity, 
    outliers and heavy-tailed distributions in the response variable due to the check
    loss function in quantile regression. In addition, it can flexibly model the dynamic
    pattern of regression coefficients through nonparametric varying coefficient 
    functions. Although high dimensional quantile varying coefficient model has been 
    examined extensively in the frequentist framework, the corresponding Bayesian variable 
    selection methods have rarely been developed. In this package, we have implemented 
    the Gibbs samplers of the penalized Bayesian quantile varying coefficient model with 
    the spike-and-slab priors [Zhou et al.(2023)]<doi.org/10.1016/j.csda.2023.107808>. 
    The Markov Chain Monte Carlo (MCMC) algorithms of the proposed
    and alternative models can be efficiently performed by using the package.


## How to install

    - To install from Github, run these two lines of code in R

<!-- end list -->

    install.packages("devtools")
    devtools::install_github("cenwu/pqrBayes")

  - Released versions of pqrBayes are available on CRAN
    [(link)](https://cran.r-project.org/package=pqrBayes), and can be
    installed within R via

<!-- end list -->

    install.packages("pqrBayes")



## Methods

This package provides implementation for methods proposed in

  - Zhou, F., Ren, J., Ma, S. and Wu, C. (2023). The Bayesian Regularized Quantile Varying Coefficient Model.  Computational Statistics & Data Analysis, 107808 \doi{10.1016/j.csda.2023.107808}
