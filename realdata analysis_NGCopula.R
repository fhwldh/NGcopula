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


#------------------------------------------------------
## summarize result : Table 4


# Independent
res_indep = summary_tab(fit_indep, d1=d1, 
                        extra=NULL, cop="indep")
# Gaussian
res_gs = summary_tab(fit_jointgs, d1=d1, 
                     extra="rho", cop="gs")
# Joe
res_joe = summary_tab(fit_jointjoe, d1=d1, 
                      extra="rho", cop="joe")
# NG
res_ng =summary_tab(fit_jointng, d1=d1, 
                    extra="alpha", cop="ng")

# goodness of fit 
gof = cbind(res_indep$gof, res_gs$gof, res_joe$gof, res_ng$gof)
colnames(gof) = c("indep","gaussian","joe", "ng")
gof

# estimation result 
est_tab = cbind(rbind(res_indep$est_tab,rep(NA,2)), 
                res_gs$est_tab, 
                res_joe$est_tab, 
                res_ng$est_tab)
row.names(est_tab)[18] = "copula"
# est_tab


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


#------------------------------------------------------------
## result of validation : Table 5, Figure3, Figure4


# total loss of test dataset
log(sum(test_dat$S)) # 14.01836


# Table 5 
round(rbind(summary(log(indep_Sv)), summary(log(gs_Sv)),
            summary(log(joe_Sv)), summary(log(ng_Sv)) ),3)

round(rbind(sd(log(indep_Sv)),
            sd(log(gs_Sv)),
            sd(log(joe_Sv)),
            sd(log(ng_Sv)) ),3)    

    
#--------------------------------------------------
# Figure 3
# pdf("tail_logS_val.pdf", width=6, height=5)
par(mfrow=c(4,1), mar=c(4,4,2,2))
d_log_indep_Sv = density(log(indep_Sv), width=0.7)
idx1 = d_log_indep_Sv$x>14.5
d_log_gs_Sv = density(log(gs_Sv), width=0.7)
idx2 = d_log_gs_Sv$x>14.5
d_log_joe_Sv = density(log(joe_Sv), width=0.7)
idx3 = d_log_joe_Sv$x>14.5
d_log_ng_Sv = density(log(ng_Sv), width=0.7)
idx4 = d_log_ng_Sv$x>14.5

par(mfrow=c(1,1))
plot(d_log_indep_Sv$x[idx1], d_log_indep_Sv$y[idx1], 
     type="l", xlim=c(14.5, 21), ylim=c(0,0.3),
     ylab="Density", xlab="Total loss (log)", lwd=2)
lines(d_log_gs_Sv$x[idx2], d_log_gs_Sv$y[idx2], lty=2, col=2, lwd=2)
lines(d_log_joe_Sv$x[idx3], d_log_joe_Sv$y[idx3], lty=3, col=3, lwd=2)
lines(d_log_ng_Sv$x[idx4], d_log_ng_Sv$y[idx4], lty=4, col=4, lwd=2)

legend("topright", c("Independent","Gaussian","Joe", "NG"),
       lty=c(1:4), col=c(1:4), lwd=2)
# dev.off()




#--------------------------------------------------
# Figure 4
# pdf("hist_logS_val.pdf", width=8, height=8)
par(mfrow=c(4,1), mar=c(4,4,2,2))
hist(log(indep_Sv), breaks=50, freq=F, col="white",
     main="Independent model", xlab="", 
     ylim=c(0,1.4), xlim=c(12,21))
abline(v=log(sum(test_dat$S)), col="tomato", lwd=2)
abline(v=mean(log(indep_Sv)), col="navy", lty=2, lwd=2)
legend("topright", c("Actual total loss", "Mean estimated total loss"), 
       col=c("tomato","navy"), lty=c(1,2), lwd=2)
hist(log(gs_Sv), breaks=60, freq=F, col="white",
     main="Gaussian model", xlab="",  
     ylim=c(0,1.4), xlim=c(12,21))
abline(v=log(sum(test_dat$S)), col="tomato", lwd=2)
abline(v=mean(log(gs_Sv)), col="navy", lty=2, lwd=2)
hist(log(joe_Sv), breaks=50, freq=F, col="white",
     main="Joe copula model", xlab="",  
     ylim=c(0,1.4), xlim=c(12,21))
abline(v=log(sum(test_dat$S)), col="tomato", lwd=2)
abline(v=mean(log(joe_Sv)), col="navy", lty=2, lwd=2)
hist(log(ng_Sv), breaks=50,  freq=F, col="white",
     main="NG-copula model", ylim=c(0,1.4), xlim=c(12,21),
     xlab=expression(paste("Total loss (log)")))
abline(v=log(sum(test_dat$S)), col="tomato", lwd=2)
abline(v=mean(log(ng_Sv)), col="navy", lty=2, lwd=2)
# dev.off()
