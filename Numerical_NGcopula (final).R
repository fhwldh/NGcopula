################################################################################
## Normal-Gamma copula 

## Figure1, numerical study, Figure2 

## Written by Rosy Oh

## Last updated  28 March 2024
################################################################################

# install.packages("VineCopula")
# install.packages("GB2")
# install.packages("metRology")
# install.packages("copula")

library(VineCopula)
library(metRology) # for pt.scaled function 
library(copula)
library(MASS) # for glm.nb
library(actuar) # qztnbinom # truncated negative binomial
library(dplyr)
library(tidyr)


source("functions_NGcopula.r")

################################################################################
## Figure 1 


# pdf("scatter_NGcopula_theta.pdf", width=8, height=5)
par(mfrow=c(2,3), mar=c(4,4,2,2))
theta= c(0.1, 0.5, 1, 5, 10, 30) 
beta=1
tmp = matrix(NA, nrow=length(theta), ncol=3)
pb <- txtProgressBar(min = 1, max = length(theta), style = 3)
for(i in 1:length(theta)){
  uv = sim = sim_uv_ng(N=10000, param=c(theta[i],beta))
  plot(uv[,1], uv[,2], pch='.',
       xlab = expression(U), ylab = expression(V),
       main=paste(expression(theta),"=",theta[i]))
  setTxtProgressBar(pb, i) 
}
# dev.off()



par(mfrow=c(1,1))


################################################################################
## Numerical study : estimation 



#---------------------------------------------------
# simulated dataset 


# runif(1)
set.seed(367)


I=5000

x1 = rbinom(I, p=0.5, size=1);  x2 = rbinom(I, p=0.5, size=1)
X=model.matrix(~x1+x2)
class(X)
head(X)

## parameter (Table 1)

nsetting=12
param.set= data.frame(
  b10=rep(-1.5, nsetting), b11 = rep(0.5, nsetting),
  b12 = c(rep(1, 6),rep(2, 6)), b20 = rep(6.5, nsetting),
  b21 = rep(0.1, nsetting), b22 = rep(0.5, nsetting),
  r = rep(0.2, nsetting),   psi = rep(1.6, nsetting),
  k = rep(0.2, nsetting),  theta = rep(c(0.1, 0.5, 1, 5, 10, 30),2))


#---------------------------------------------------
## simulation 

B=100 # 500

# B1, B2, r, psi, k, alpha 
pb <- txtProgressBar(min = 1, max = nrow(param.set), style = 3)
sim.res = list()
for(i in 1:nrow(param.set)){
  param = param.set[i,]
  print(paste0("iter=", i))
  print(paste(colnames(param.set),":",param))
  sim.res[[i]] = Sim_NG_est(Xmat=X, param=param, B=B)
  print(sim.res[[i]]$measure)
}


sum(sapply(sim.res, "[[", "time"))/3600 #4.25 hr
sapply(sim.res, "[[", "n.error")
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 # theta 


#---------------------------------------------------
# result : Table2, Table3

RB = matrix(NA, nrow=nrow(param.set), ncol=ncol(param.set))
MSE = matrix(NA, nrow=nrow(param.set), ncol=ncol(param.set))
for(j in 1:nrow(param.set)){
  RB[j,] =  sim.res[[j]]$measure[,"rb"]
  MSE[j,] = sim.res[[j]]$measure[,"mse"]
}
round(RB,2) # relative bias
round(MSE,4) # mean squared error



#---------------------------------------------------
## simulation for total loss 

# other dataset 

set.seed(295)
I=5000
x1 = rbinom(I, p=0.5, size=1);  x2 = rbinom(I, p=0.5, size=1)
X=model.matrix(~x1+x2)

B=5000 
ng_Sv = list()
pb1 <- txtProgressBar(min = 1, max = nrow(param.set), style = 3)
for(j in 1:nrow(param.set)){
  ng_Sv_tmp=c()
  set.seed(5756)
  for(i in 1:B){
    param = param.set[j,]
    param[(ncol(X)*2+1):length(param)] = 
      log(param[(ncol(X)*2+1):length(param)] )
    
    ng_simdat = simDat(Xmat=X, w=rep(1,nrow(X)),
                       param=param, cop="ng")

    ng_Sv_tmp[i] = sum(ng_simdat$S)
   
  }
  ng_Sv[[j]] = ng_Sv_tmp
  setTxtProgressBar(pb1, j) 
}


#-----------------------------------------------------------------------
## Figure2 

# scenario 1 ~ scenario 6
# pdf("tail_logS_sim1.pdf", width=5, height=4)
par(mfrow=c(1,1), mar=c(4,4,2,2))

d_logSv = density(log(ng_Sv[[1]]), width=0.6)
idx1 = d_logSv$x>16 & d_logSv$x<18 
mycol= hcl.colors(6)
plot(d_logSv$x[idx1], d_logSv$y[idx1], type="l", col=mycol[1],
     xlim=c(16.5, 17.2), ylim=c(0,2),  lwd=2,
     ylab="Density", xlab="Total loss (log)")
for(i in 2:6){
  d_logSv2 = density(log(ng_Sv[[i]]), width=0.6)
  idx2 = d_logSv2$x>16 & d_logSv2$x<18 
  lines(d_logSv2$x[idx2], d_logSv2$y[idx2], lty=i, col=mycol[i], lwd=2)
}
# dev.off()


# scenario 7 ~ scenario 12
# pdf("tail_logS_sim2.pdf", width=5, height=4)
par(mfrow=c(1,1), mar=c(4,4,2,2))

d_logSv = density(log(ng_Sv[[7]]), width=0.6)
idx1 = d_logSv$x>17 & d_logSv$x<18 
mycol= hcl.colors(6)
plot(d_logSv$x[idx1], d_logSv$y[idx1], type="l", col=mycol[1],
     xlim=c(17.2, 18), ylim=c(0,2), lwd=2,
     ylab="Density", xlab="Total loss (log)")
for(i in 8:12){
  d_logSv2 = density(log(ng_Sv[[i]]), width=0.6)
  idx2 = d_logSv2$x>17 & d_logSv2$x<18 
  lines(d_logSv2$x[idx2], d_logSv2$y[idx2], lty=i, col=mycol[i-6], lwd=2)
}
legend("topright", 
       c(expression(paste(theta,"=0.1  ")),expression(paste(theta,"=0.5  ")),
         expression(paste(theta,"=1 ")),expression(paste(theta,"=5 ")),
         expression(paste(theta,"=10 ")),expression(paste(theta,"=30 ")) ),
       lty=c(1:6), col=hcl.colors(6), lwd=2)
# dev.off()




