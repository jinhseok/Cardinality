#Bayesian analysis of children's CP-knower status at time 2
#see mult_appendix for model details

library(MCMCpack)
set.seed(8285)
setwd("~/Desktop/Analysis/Cardinality")
data=read.csv(file="~/Desktop/Analysis/PS_Data/Y1_ABC_clean.csv")
#save only variables that are needed to X
a=data$Num_T2/100 #ANS acc at T2
t=data$MagBox_T2/100  #OTS acc at T2
c=data$CPknow_T2 #0 if non-CP, 1 if CP-knower
X=cbind(a,t,c)
X=na.omit(X)
a=X[,1]
t=X[,2]
c=X[,3]
N=length(c)
rm(data)


#MCMC sampler stuff
M=10000
alpha=1:M
beta=1:M
a.sdTune=.6
b.sdTune=.45
a.counter=0
b.counter=0

#exponential priors on alpha and beta
a.lambda=.5
b.lambda=.5

#start values
alpha[1]=1
beta[1]=1

#likelihood function
f.eval = function(alpha, beta){
  p = a^alpha*t^beta
  d = dbinom(c, 1, p, log=TRUE)
  if (any(is.nan(d))) {
    print(list(alpha=alpha, beta=beta, p=p))
  }
  sum(d)
}

#chain
for (m in 2:M){
  #print(m)
  
  #sample alpha given all else
  cand=rnorm(1,alpha[m-1],a.sdTune)
  f.cur=f.eval(alpha[m-1], beta[m-1]) + dexp(alpha[m-1], a.lambda, log=TRUE)
  prob = 0
  if (cand >= 0) {
    f.cand=f.eval(cand, beta[m-1]) + dexp(cand, a.lambda, log=TRUE)
    prob=min(exp(f.cand-f.cur),1)
  }
  if(runif(1)<prob){
    alpha[m]=cand
    a.counter=a.counter+1
  }else{
    alpha[m]=alpha[m-1]
  }
  
  #sample beta given all else
  cand=rnorm(1,beta[m-1],b.sdTune)
  f.cur=f.eval(alpha[m], beta[m-1]) + dexp(beta[m-1], b.lambda, log=TRUE)
  prob = 0
  if (cand >= 0) {
    f.cand=f.eval(alpha[m], cand) + dexp(cand, b.lambda, log=TRUE)
    prob=min(exp(f.cand-f.cur),1)
  }
  if(runif(1)<prob){
    beta[m]=cand
    b.counter=b.counter+1
  }else{
    beta[m]=beta[m-1]
  }	
}

# print(paste("Acceptence Rate for alpha  ",a.counter/(M)))
# print(paste("Acceptence Rate for beta ",b.counter/(M)))
# ma=mean(alpha)
# mb=mean(beta)
# print(paste("alpha:", ma))
# print(paste("beta:", mb))
# thinned_ma=mean(alpha[seq(1000,10000,10)]) #checked param means for thinned samples
# thinned_mb=mean(beta[seq(1000,10000,10)])
# print(paste("alpha[seq(1000,10000,10)]:", thinned_ma))
# print(paste("beta[sea(1000,10000,10)]:", thinned_mb))
# hist(alpha)
# hist(beta)
# 
# par(mfrow=c(2,2))
# plot(alpha[seq(1000, 10000, 10)],typ='l')
# acf(alpha[seq(1000, 10000, 10)])
# plot(beta[seq(1000, 10000, 10)],typ='l')
# acf(beta[seq(1000, 10000, 10)])

out=cbind(alpha,beta)
write.table(file="T2.chain",out,quote=F,row.names=F)
