################################################################################
## Normal Gamma copula 

## Real data example


## Written by Rosy Oh

## Last updated 28 Mar 2024

################################################################################

# install.packages("dplyr")
# install.packages("VineCopula")
# install.packages("GB2")
# install.packages("metRology")
# install.packages("copula")
# install.packages("xlsx")
# install.packages("actuar")

library(VineCopula)
library(metRology) # for pt.scaled function 
library(copula)
library(MASS) # for glm.nb
library(actuar) # qztnbinom # truncated negative binomial
library(dplyr)
library(xlsx)


source("functions_NGcopula.r")


################################################################################
## real dataset 

source("data_NGcopula.r") # prepare dataset for NG copula

dim(train_dat) #[1] 817  11
round(table(train_dat$N)/nrow(train_dat),3)

mo1 = glm.nb(N~.-M-S-PolicyNum-exposure+offset(log(exposure)), 
             data=train_dat)
# summary(mo1)
mo2 = glm(M~.-N-S-PolicyNum-exposure, data=subset(train_dat, N>0), 
          family = Gamma(link='log'))
# summary(mo2)
d1 = dim(model.matrix(mo2))[2]; d1


#---------------------------------------------------
# estimation 

# par =(beta1, beta2, r, psi, k)
fit_indep = optim(par=c(coef(mo1), coef(mo2), log(mo1$theta), 
                        0.1, -0.1),
                  fn=log_lik_indep,
                  n=train_dat$N, m=train_dat$M, X=model.matrix(mo1), 
                  w=train_dat$exposure, 
                  method="BFGS", hessian=TRUE,
                  control=list(maxit=10000,fnscale=-1))
# fit_indep$par

#------------------------------------------------------
fit_IFMng  = optim(par=c(1), log_lik_IFMng_onlya ,
                   n=train_dat$N, m=train_dat$M, X=model.matrix(mo1), 
                   w=train_dat$exposure, 
                   B1 =fit_indep$par[1     :   d1 ],
                   B2 =fit_indep$par[(d1+1):(2*d1)],
                   r  =exp(fit_indep$par[d1*2+1]),
                   psi=exp(fit_indep$par[d1*2+2]),
                   k  =exp(fit_indep$par[d1*2+3]),
                   #hessian=TRUE, 
                   method="Brent", upper=10, lower=-10,
                   control=list(fnscale=-1))
# fit_IFMng$par


# par =(beta1, beta2, r, psi, k, a)
fit_jointng = optim(par=c(fit_indep$par, fit_IFMng$par), 
                     fn=log_lik_jointng_onlya,
                     n=train_dat$N, m=train_dat$M, X=model.matrix(mo1),
                     w=train_dat$exposure,  
                     hessian=TRUE, method="BFGS", 
                     control=list(maxit=10000,fnscale=-1))
# fit_jointng$convergence
# fit_jointng$par

#------------------------------------------------------
fit_IFMgs  = optim(par=0, fn=log_lik_IFMgs,
                    n=train_dat$N, m=train_dat$M, 
                    X=model.matrix(mo1), w=train_dat$exposure, 
                    B1 =fit_indep$par[1     :   d1 ],
                    B2 =fit_indep$par[(d1+1):(2*d1)],
                    r  =exp(fit_indep$par[d1*2+1]),
                    psi=exp(fit_indep$par[d1*2+2]),
                    k  =exp(fit_indep$par[d1*2+3]),
                    method="Brent", upper=10, lower=-10,
                   control=list(maxit=10000,fnscale=-1))
# fit_IFMgs 

fit_jointgs = optim(par=c(fit_indep$par, fit_IFMgs$par),
                    fn=log_lik_jointgs,
                    n=train_dat$N, m=train_dat$M, X=model.matrix(mo1),
                    w=train_dat$exposure,
                    hessian=TRUE, method="BFGS", 
                    control=list(maxit=10000,fnscale=-1))

# fit_jointgs$convergence
# fit_jointgs$par

#------------------------------------------------------
fit_IFMjoe  = optim(par=0, fn=log_lik_IFMjoe,
                    n=train_dat$N, m=train_dat$M, X=model.matrix(mo1),
                    w=train_dat$exposure, 
                    B1 =fit_indep$par[1     :   d1 ],
                    B2 =fit_indep$par[(d1+1):(2*d1)],
                    r  =exp(fit_indep$par[d1*2+1]),
                    psi=exp(fit_indep$par[d1*2+2]),
                    k  =exp(fit_indep$par[d1*2+3]),
                    method="Brent", upper=10, lower=-10,
                    control=list(maxit=10000,fnscale=-1))
# fit_IFMjoe

fit_jointjoe = optim(par=c(fit_indep$par, fit_IFMjoe$par),
                     fn=log_lik_jointjoe,
                     n=train_dat$N, m=train_dat$M, X=model.matrix(mo1),
                     w=train_dat$exposure, 
                     method="BFGS", hessian=TRUE,
                     control=list(maxit=10000,fnscale=-1))
# fit_jointjoe$convergence
# fit_jointjoe$par



#####################################################################
## validation 

dim(test_dat) # 716  11
head(test_dat)
test_dat$exposure = 1

# simulated dataset
mo1_test = glm.nb(N~.-M-S-PolicyNum-exposure, data=test_dat)
summary(mo1_test)


# total loss 

B=5000 # nrow(train_dat)#
indep_Sv=c();gs_Sv=c();joe_Sv=c();ng_Sv=c()
pb1 <- txtProgressBar(min = 1, max = B, style = 3)
for(i in 1:B){
  indep_simdat = simDat(Xmat=model.matrix(mo1_test), w=test_dat$exposure,
                        param=fit_indep$par, cop="indep")
  gs_simdat = simDat(Xmat=model.matrix(mo1_test), w=test_dat$exposure,
                     param=fit_jointgs$par, cop="gs")
  joe_simdat = simDat(Xmat=model.matrix(mo1_test), w=test_dat$exposure, 
                      param=fit_jointjoe$par, cop="joe")
  ng_simdat = simDat(Xmat=model.matrix(mo1_test), w=test_dat$exposure,
                     param=fit_jointng$par, cop="ng")
  
  indep_Sv[i] = sum(indep_simdat$S)
  gs_Sv[i] = sum(gs_simdat$S)
  joe_Sv[i] = sum(joe_simdat$S)
  ng_Sv[i] = sum(ng_simdat$S)
  setTxtProgressBar(pb1, i) 
}


