################################################################################
## Functions: NG copula 

## Last updated : 28 Mar 2024
## Written by Rosy Oh 

################################################################################
#----------------------------------------------------------------------

library(GB2) # Generalized Beta Distribution
library(metRology)


#----------------------------------------------------------------------

sim_uv_ng = function(N=10000, param){
  
  a=1/param[1]; b=param[2]
  
  u = runif(N)
  tmp1 = qgamma(u, shape=a, rate=b)
  tmp2 = rnorm(N, 0, sqrt(1/tmp1))
  v = pt.scaled(tmp2, df=2*a, mean=0, sd=sqrt(b/a))
  
  ng_sample = cbind(u,v)
  return(ng_sample)
}


#--------------------------------------------
# zero-truncated negative binomial distribution 

# NB mean: r(1-p)/p
# pztnbinom(n, size=w*r, prob=r/(lambda1+r))
# assume n is a vector and X is matrix

f1_p = function(n, X, w, B1, r){
  lambda1 = w*exp(X%*%B1)
  result= dztnbinom(n, size=w*r, prob=r*w/(lambda1+r*w))
  return(result)
}


F1_p = function(n, X, w, B1, r){
  lambda1 = w*exp(X%*%B1)
  result= pztnbinom(n, size=w*r, prob=r*w/(lambda1+r*w))
  return(result)
}

F1inv_p = function(pp, X, w, B1, r){
  lambda1 = w*exp(X%*%B1)
  result= qztnbinom(pp, size=w*r, prob=r*w/(lambda1+r*w))
  return(result)
}


#--------------------------------------------
# GB2 : generalized beta distribution of the second kind 
# GB2(a=shape1, b=scale, p=shape2, q=shape3)

F2inv = function(pp, X, B2, w, psi, k) {
  result = qgb2(pp, 1, exp(X%*%B2), w*psi, w*k+1)
  return(c(result)) 
}

F2 = function(m, X, B2, w, psi, k) {
  result = pgb2(m, 1, exp(X%*%B2), w*psi, w*k+1)
  return(c(result))
}

f2 = function(m, X, B2, w, psi, k) {
  result <- dgb2(m, 1, exp(X%*%B2), w*psi, w*k+1)
  return(c(result) )
}


#--------------------------------------------
## log likelihood function 
 
# par =(beta1, beta2, r, psi, k, a, b)


log_lik_jointng_onlya <- function(n, m, X, w, param_vec, b=1){
  
  d1 = dim(X)[2]
  B1=param_vec[1:d1]
  B2=param_vec[(d1+1):(2*d1)]
  
  r  =exp(param_vec[(2*d1)+1])
  psi=exp(param_vec[(2*d1)+2])
  k  =exp(param_vec[(2*d1)+3])
  a  = 1/exp(param_vec[(2*d1)+4])
  
  lamb = w*exp(X%*%B1)
  p= 1-dnbinom(0, mu=lamb, size=w*r)
  idx = n>0
  tmp_h_ng = h_ng(n[idx], m[idx], X[idx, ], w[idx],  
                  B1, B2, r, a, b, psi,k)
  hng = ifelse(tmp_h_ng==0, 1e-26, tmp_h_ng)
  temp1 = sum(log(  p[ idx]*hng))
  
  temp2 = sum(log(1-p[!idx]))
  
  return((temp1+temp2))
}


log_lik_jointgs = function(param_vec, n, m, X, w){
  
  d1 = dim(X)[2]
  B1=param_vec[1:d1]
  B2=param_vec[(d1+1):(2*d1)]
  
  r  =exp( param_vec[(2*d1)+1])
  psi=exp( param_vec[(2*d1)+2])
  k  =exp( param_vec[(2*d1)+3])
  rho=atan(param_vec[(2*d1)+4])*2/pi
  
  lamb = w*exp(X%*%B1)
  p= 1-dnbinom(0, mu=lamb, size=w*r)
  # p= 1-(1+lamb/r/w)^(-r*w)
  idx = n>0
  tmp_h_gs = h_gs(n[idx], m[idx], X[idx, ], w[idx],  
                  B1, B2, r, rho, psi,k)
  hgs = ifelse(tmp_h_gs==0, 1e-26, tmp_h_gs) 

  temp1 = sum(log(  p[ idx]*hgs))
  temp2 = sum(log(1-p[!idx]))
  
  return(temp1+temp2)
}


log_lik_jointjoe <- function(param_vec, n, m, X, w){
  
  d1 = dim(X)[2]
  B1=param_vec[1:d1]
  B2=param_vec[(d1+1):(2*d1)]
  
  r  =exp( param_vec[(2*d1)+1])
  psi=exp( param_vec[(2*d1)+2])
  k  =exp( param_vec[(2*d1)+3])
  rho=exp( param_vec[(2*d1)+4])+1
  
  lamb = w*exp(X%*%B1)
  p= 1-dnbinom(0, mu=lamb, size=w*r)
  idx = n>0
  
  tmp_h_joe = h_joe(n[idx], m[idx], X[idx, ], w[idx],  
                  B1, B2, r, rho, psi,k)
  hjoe = ifelse(tmp_h_joe==0, 1e-26, tmp_h_joe) 
  temp1 = sum(log(  p[ idx]*hjoe))
  temp2 = sum(log(1-p[!idx]))
  
  return(temp1+temp2) }


#--------------------------------------------
h_ng = function(n, m, X, w, B1, B2, r, a, b, psi, k){

  # F2m = F2(m, X, B2, psi, k)
  F2m = F2(m, X, B2, w, psi, k)
  qt_F2 = qt.scaled(F2m, df=2*a, mean=0, sd=sqrt(b/a))
  
  temp1 = pgamma(qgamma(F1_p(n, X, w, B1, r), shape=a, rate=b), 
                 shape=a+0.5, rate=b+0.5*(qt_F2)^2)-  
    pgamma(qgamma(F1_p(n-1, X, w, B1, r), shape=a, rate=b), 
           shape=a+0.5, rate=b+0.5*(qt_F2)^2)
  temp2 = f2(m, X, B2, w, psi, k)
  
  result = temp1*temp2

  return(result)
}


# BiCopHfunc : Evaluate the conditional distribution function
# hfunc1: h_1(u2|u1) conditional distribution function
# hfunc2: h_2(u1|u2)
h_gs = function(n, m, X, w, B1, B2, r, rho, psi, k){
  
  F2m = F2(m, X, B2, w, psi, k)
  temp1 = BiCopHfunc(F1_p(n, X, w, B1, r),F2m, 
                     family = 1, par=rho)$hfunc2 -
    BiCopHfunc(F1_p(n-1, X, w, B1, r),F2m, 
               family = 1, par=rho)$hfunc2
  temp2 = f2(m, X, B2, w, psi, k)
  
  result = temp1*temp2
  return(result) 
  }

h_joe= function(n, m, X, w, B1, B2, r, rho, psi, k){
  
  F2m = F2(m, X, B2, w, psi, k)
  temp1 = BiCopHfunc(F1_p(n, X, w, B1, r),F2m, 
                     family = 6, par=rho)$hfunc2 -
    BiCopHfunc(F1_p(n-1, X, w, B1, r),F2m, 
               family = 6, par=rho)$hfunc2
  temp2 = f2(m, X, B2, w, psi, k)
  result = temp1*temp2

  return(result) 
}


#--------------------------------------------
h_indep<- function(n, m, X, w, B1, B2, r,  psi, k){
  temp1 = f1_p(n  , X, w, B1, r)
  temp2 = f2(m, X, B2, w, psi, k)
  
  result = temp1*temp2
  return(result)
}


log_lik_indep <- function(n, m, X, w, param_vec){
  
  d1 = dim(X)[2]
  B1=param_vec[1:d1]
  B2=param_vec[(d1+1):(2*d1)]
  
  r  =exp(param_vec[(2*d1)+1])
  psi=exp(param_vec[(2*d1)+2])
  k  =exp(param_vec[(2*d1)+3])
  
  lamb = w*exp(X%*%B1)
  # p= c(1-(1+lamb/r/w)^(-r*w)) # 1-P(X=0)
  p = 1-dnbinom(0, mu=lamb, size=w*r)
  idx = n>0
  temp1 = sum(log(  p[ idx]) + log(h_indep(n[idx], m[idx], X[idx, ], 
                                            w[idx],  B1, B2,
                                            r, psi,k))) 
  temp2 = sum(log(1-p[!idx]))
  return((temp1+temp2) )
  }

#--------------------------------------------
log_lik_IFMng <- function(n, m, X, w, B1, B2, r, psi, k, 
                          param_vec){ #, b=1
  
  a  = 1/exp(param_vec[1])
  b  = exp(param_vec[2])
  
  lamb = w*exp(X%*%B1)
  p = 1-dnbinom(0, mu=lamb, size=w*r)
  idx = n>0
  temp1 = sum(log(  p[ idx]) +
                log(h_ng(n[idx], m[idx], X[idx, ], w[idx],  
                         B1, B2,r, a,b, psi,k)))
  temp2 = sum(log(1-p[!idx]))
  return(temp1+temp2) }


log_lik_IFMng_onlya <- function(n, m, X, w, B1, B2, r, psi, k, 
                                param_vec, b=1){ 
  a  = 1/exp(param_vec[1])
 
  lamb = w*exp(X%*%B1)

  p = 1-dnbinom(0, mu=lamb, size=w*r)
  idx = n>0
  temp1 = sum(log(  p[ idx]) +log(h_ng(n[idx], m[idx], X[idx, ], w[idx],  
                                       B1, B2,
                                       r, a,b, psi,k)))
  temp2 = sum(log(1-p[!idx]))
  return(temp1+temp2)
}


log_lik_IFMgs = function(n, m, X, w, B1, B2, r, psi, k, 
                         param_vec){
  
  rho=2/pi*atan(param_vec)
  
  lamb = w*exp(X%*%B1)
  p=1-dnbinom(0, mu=lamb, size=w*r)
  idx = n>0
  temp1 = sum(log(  p[ idx]) +
                log(h_gs(n[idx], m[idx], X[idx, ], w[idx],  
                         B1, B2, r, rho, psi,k)))
  temp2 = sum(log(1-p[!idx]))
  
  return(temp1+temp2) }


log_lik_IFMjoe <- function(n, m, X, w, B1, B2, r, psi, k, param_vec){
  
  rho=exp(param_vec)+1
  
  lamb = w*exp(X%*%B1)
  p=1-dnbinom(0, mu=lamb, size=w*r)
  idx = n>0
  temp1 = sum(log(  p[ idx]) +
                log(h_joe(n[idx], m[idx], X[idx, ], w[idx],  
                          B1, B2, r, rho, psi,k)))
  temp2 = sum(log(1-p[!idx]))
  return(temp1+temp2) 
}


#----------------------------------------------------------------------
simDat = function(Xmat, w,  param, cop, b = 1){
  param = as.numeric(param)
  
  
  d1 = dim(Xmat)[2]
  
  B1 =param[1     :   d1 ]
  B2 =param[(d1+1):(2*d1)]
  r  = exp(param[d1*2+1])
  psi= exp(param[d1*2+2])
  k  = exp(param[d1*2+3])
  
  I=nrow(Xmat)
  lamb = w*exp(Xmat%*%B1)
  
  p = 1-dnbinom(0, mu=lamb, size=w*r)
  R = rbinom(I, size=1, prob=p) # R ~ Ber(p)
  idx_0 = which(R==0)
  idx_not0 = which(R==1); Iplus = length(idx_not0)
  
  if(cop=="indep"){
    uv= rCopula(Iplus, indepCopula(dim=2))
  }else  if(cop == "gs"){
    theta = atan(param[d1*2+4])*2/pi
    uv= rCopula(Iplus, normalCopula(theta, dim=2))
  }else if(cop=="joe"){
    theta = exp(param[d1*2+4])+1
    uv = rCopula(Iplus, joeCopula(theta, dim=2))
  }else{
    a = 1/exp(param[d1*2+4]) # alpha 
    u = runif(Iplus)
    tmp1 = qgamma(u, shape=a, rate=b)
    tmp2 = rnorm(length(tmp1), 0, sqrt(1/tmp1))
    v = pt.scaled(tmp2, df=2*a, mean=0, sd=sqrt(b/a))

    uv = cbind(u,v)
  }
  
  Xplus = Xmat[idx_not0,]
  wplus = w[idx_not0]
  # for n>0
  Nplus = F1inv_p(uv[,1], X=Xplus, w=wplus, B1=B1, r=r)
  # for m>0
  Mplus=c()
  for(i in 1:nrow(Xplus)){
    Mplus[i] = F2inv(uv[i,2], X=Xplus[i,], w=wplus[i], B2=B2, psi=psi, k=k)  
  }

  N=rep(0, I); M=rep(0, I)
  N[idx_not0]=Nplus; M[idx_not0]=Mplus
  simdat = data.frame(Xmat,
                      w=w, idx_not0=R,
                      N=N, M=M,S=N*M)%>%
    filter(!is.na(M))%>%filter(!is.na(N))
  
  return(simdat)
  
}


#---------------------------------------------------------------------

simDat_NGcopula = function(Xmat, param, b = 1 ){
  param = as.numeric(param)
  
  B1 = param[1:3];
  B2 = param[4:6];
  r = param[7]; psi=param[8]; k=param[9]
  a = 1/param[10] # alpha 
 
  I=nrow(Xmat)
  w=rep(1, I)
  
  lamb = w*exp(Xmat%*%B1)
  
  p = 1-dnbinom(0, mu=lamb, size=w*r)
  R = rbinom(I, size=1, prob=p) # R ~ Ber(p)
  idx_0 = which(R==0)
  idx_not0 = which(R==1); Iplus = length(idx_not0)
  
  
  u = runif(Iplus)
  tmp1 = qgamma(u, shape=a, rate=b)
  tmp2 = rnorm(length(tmp1), 0, sqrt(1/tmp1))
  v = pt.scaled(tmp2, df=2*a, mean=0, sd=sqrt(b/a))
  
  uv = cbind(u,v)
  
  Xplus = Xmat[idx_not0,]
  
  # for n>0
  Nplus = F1inv_p(uv[,1], X=Xplus, w=rep(1, Iplus), B1=B1, r=r)
  # for m>0
  Mplus=c()
  for(i in 1:nrow(Xplus)){
    Mplus[i] = F2inv(uv[i,2], X=Xplus[i,], w=w[i], B2=B2, psi=psi, k=k)  
  } 
  
 
  N=rep(0, I); M=rep(0, I)
  N[idx_not0]=Nplus; M[idx_not0]=Mplus
  simdat = data.frame(x1 = unname(Xmat[,2]), 
                       x2 = unname(Xmat[,3]), 
                       N=N, M=M, 
                       w=w, idx_not0=R)%>%
    filter(!is.na(M))%>%filter(!is.na(N))
  
  return(simdat)
  
}

#----------------------------------------------------------------------
Sim_NG_est = function(Xmat, param, B){
  # runif(1)
  set.seed(766) # 101, 742
  par.name = names(param)
  param = as.numeric(param)
  
  init.time=proc.time()
  est.our.B = matrix(NA, ncol=length(param), nrow=B)
  se.our.B = matrix(NA, ncol=length(param), nrow=B)
  pb1 <- txtProgressBar(min = 1, max = B, style = 3)
  n.error = 0
  for(bb in 1:B){###################################################### <- start iteration 
    ## generate data set 
    sim_dat = simDat_NGcopula(Xmat=Xmat, param=param)
    
    # par = c(beta1, beta2, r, psi, k)
    fit_indep = optim(par=c(c(-1,0.1,0.1),
                            c(log(mean(sim_dat$M)),0.1,0.1), 
                            0, 1, 1),#coef(mo1), coef(mo2), log(mo1$theta)
                      fn=log_lik_indep,
                      n=sim_dat$N, m=sim_dat$M, X=Xmat, 
                      w=sim_dat$w,  # exposure
                      method="BFGS", #hessian=TRUE,
                      control=list(maxit=10000,fnscale=-1))
    d1 = ncol(Xmat)
    
    sim_N = ifelse(sum(sim_dat$idx_not0)<1000, sum(sim_dat$idx_not0), 1000)
    tmp_idx = sample(1:sum(sim_dat$idx_not0), sim_N)
    tmp_dat = sim_dat%>%filter(N>0)%>%.[tmp_idx,]
    tmp_v=F2(tmp_dat$M, X=model.matrix(~tmp_dat$x1+tmp_dat$x2), w=sim_dat$w,
       B2=fit_indep$par[(d1+1):(2*d1)], 
       psi=exp(fit_indep$par[d1*2+2]), k=exp(fit_indep$par[d1*2+3]))
    cc = cor(tmp_dat[tmp_v<0.5,]$N, tmp_dat[tmp_v<0.5,]$M,
             method="kendall", use="complete.obs")
    rm(tmp_idx, tmp_dat, tmp_v)
    
    if(cc>0.95){
      init.theta=30
    }else if(cc<0.95 & cc>0.8){
      init.theta=10
    }else if(cc<0.8 & cc>0.5){
      init.theta=5
    }else if(cc<0.5 & cc>0.3){
      init.theta=1
    }else{
      init.theta=0.1
    }

    fit_IFMng  = optim(par=c(log(init.theta),1), fn=log_lik_IFMng, 
                       n=sim_dat$N, m=sim_dat$M, X=Xmat,
                       w=sim_dat$w,
                       B1 =fit_indep$par[1     :   d1 ],
                       B2 =fit_indep$par[(d1+1):(2*d1)],
                       r  =exp(fit_indep$par[d1*2+1]),
                       psi=exp(fit_indep$par[d1*2+2]),
                       k  =exp(fit_indep$par[d1*2+3]),
                       control=list(fnscale=-1))
    fit = try(optim(par=round(c(fit_indep$par, fit_IFMng$par[1]),1), 
                    fn=log_lik_jointng_onlya,
                    n=sim_dat$N, m=sim_dat$M, X=Xmat, 
                    w=sim_dat$w,  # exposure
                    hessian=TRUE, method="BFGS", 
                    control=list(maxit=10000,fnscale=-1)))
    
    while('try-error' %in% class( fit) ){
      n.error=n.error+1
      sim_dat = simDat_NGcopula(Xmat=Xmat, param=param)
      
      # par = c(beta1, beta2, r, psi, k)
      fit_indep = optim(par=c(c(-1,0.1,0.1),
                              c(log(mean(sim_dat$M)),0.1,0.1), 
                              0, 1, 1),
                        fn=log_lik_indep,
                        n=sim_dat$N, m=sim_dat$M, X=Xmat, 
                        w=sim_dat$w,  # exposure
                        method="BFGS", 
                        control=list(maxit=10000,fnscale=-1))
      d1 = ncol(Xmat)
      
      sim_N = ifelse(sum(sim_dat$idx_not0)<1000,
                     sum(sim_dat$idx_not0), 1000)
      tmp_idx = sample(1:sum(sim_dat$idx_not0), sim_N)
      tmp_dat = sim_dat%>%filter(N>0)%>%.[tmp_idx,]
      tmp_v=F2(tmp_dat$M, X=model.matrix(~tmp_dat$x1+tmp_dat$x2),w=sim_dat$w,
               B2=fit_indep$par[(d1+1):(2*d1)],
               psi=exp(fit_indep$par[d1*2+2]), k=exp(fit_indep$par[d1*2+3]))
      cc = cor(tmp_dat[tmp_v<0.5,]$N, tmp_dat[tmp_v<0.5,]$M,
               method="kendall", use="complete.obs")
      rm(tmp_idx, tmp_dat, tmp_v)
    
      if(cc>0.90){
        init.theta=30
      }else if(cc<0.90 & cc>0.8){
        init.theta=10
      }else if(cc<0.8 & cc>0.5){
        init.theta=5
      }else if(cc<0.5 & cc>0.3){
        init.theta=1
      }else{
        init.theta=0.1
      }

      fit_IFMng  = optim(par=c(log(init.theta),1), fn=log_lik_IFMng, #c(init.alpha, 1)
                         n=sim_dat$N, m=sim_dat$M, X=Xmat, 
                         w=sim_dat$w, 
                         B1 =fit_indep$par[1     :   d1 ],
                         B2 =fit_indep$par[(d1+1):(2*d1)],
                         r  =exp(fit_indep$par[d1*2+1]),
                         psi=exp(fit_indep$par[d1*2+2]),
                         k  =exp(fit_indep$par[d1*2+3]),
                         #hessian=TRUE, 
                         # method="Brent", upper=10, lower=-10,
                         control=list(fnscale=-1))
      
      fit = try(optim(par=round(c(fit_indep$par, fit_IFMng$par[1]),1), 
                      fn=log_lik_jointng_onlya,
                      n=sim_dat$N, m=sim_dat$M, X=Xmat, 
                      w=sim_dat$w,  # exposure, 
                      hessian=TRUE, method="BFGS", 
                      control=list(maxit=10000,fnscale=-1)))
    }
    
    d1 = ncol(Xmat)
    est = c(fit$par[1:(2*d1)], exp(fit$par[(2*d1+1):length(fit$par)]))
    trs = c(rep(1,2*d1), exp(fit$par[(2*d1+1):length(fit$par)]))
    se = sqrt(diag(solve(-fit$hessian)))*trs
    est_tab = round(data.frame(est=est, se=se),4)
    row.names(est_tab)=par.name
    
    est.our.B[bb,] = est_tab$est
    se.our.B[bb,]= est_tab$se
 
    setTxtProgressBar(pb1, bb) 
  } ## end iteration 
  time= proc.time()- init.time
 
  mean.est.our = colMeans(est.our.B)
  bias=mean.est.our-param
  mae=rowMeans(abs(t(est.our.B)-param ))
  mse=rowMeans((t(est.our.B)-param )^2)
  rb = bias/param*100
  measure= cbind(param, mean.est.our, bias, mae, mse, rb )
  
  result = list(est.B = est.our.B, measure=measure, 
                n.error=n.error, time= time[3], 
                se.B = se.our.B)
}


