#Compare mult model predictions against logit regression model predictions at time 2

require(aod)
require(ggplot2)
setwd("~/Desktop/Analysis/Cardinality")
data=read.csv(file="Y1_ABC_clean.csv")
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


########################################################################################################
#mult model mean parameter estimates from MCMC at time 2
alpha=1.25
beta=0.36

#mult model predictions
my.mult=function(a, t){
  (a^alpha)*(t^beta)  
}
mult.p=my.mult(a,t)
X=data.frame(cbind(X,mult.p))

#need to clean OTS data for plots
#table(X$t) shows shows 0.7143 x 1 and 0.25 x 1
#simplify OTS intervals a bit by subsuming these four points into nearest OTS "level"
Xmod<-X #first, make copy of data
# Xmod$t[match(0.25,Xmod$t)]<-0.375 #actually, leave .25 alone; it's too far from .375
Xmod$t[match(0.7143,Xmod$t)]<-0.75

d=c+1 #for geom_point color
#prediction plots for mult model
#plot ANS as x-axis across OTS levels
ggplot(Xmod, aes(x=a, y=mult.p, colour=factor(t), group=factor(t))) + geom_point(color=d, position=position_jitter(w=0.0075, h=0.0075)) + geom_line() + ylim(0,1) + xlab("ANS") + ylab("Predicted CP-knower Probability") + scale_colour_brewer(name="OTS",palette="Blues") + theme_bw()
#plot OTS as x-axis across ANS levels
ggplot(Xmod, aes(x=t, y=mult.p, colour=a, group=a)) + geom_point(color=d, position=position_jitter(w=0.0075, h=0.0075)) + geom_line() + ylim(0,1) + xlab("OTS") + ylab("Predicted CP-knower Probability") + scale_colour_continuous(name="ANS") + theme_bw()


########################################################################################################
#run logit model
m1 <- glm(c ~ a + t , family = binomial(link = "logit"), data=X)
# summary(m1)

#logit model prediction
my.reg=function(a,t){
  plogis(m1$coef[1]+m1$coef[2]*a+m1$coef[3]*t)
}
reg.p=my.reg(a,t)
Xmod=cbind(Xmod,reg.p)

#prediction plots for logit model
#plot ANS as x-axis across OTS levels
ggplot(Xmod, aes(x = a, y = reg.p, colour = factor(t), group=factor(t))) + geom_point(color=d, position = position_jitter(w = 0.0075, h = 0.0075)) + geom_line() + ylim(0,1) + xlab("ANS") + ylab("Predticed CP-knower Probability") + scale_colour_brewer(name="OTS",palette="Blues") + theme_bw()
#plot OTS as x-axis across ANS levels
ggplot(Xmod, aes(x = t, y = reg.p, colour = a, group=a)) + geom_point(color=d, position = position_jitter(w = 0.0075, h = 0.0075)) + geom_line() + ylim(0,1) + xlab("OTS") + ylab("Predicted CP-knower Probability") + scale_colour_continuous(name="ANS") + theme_bw()


########################################################################################################
#CP-knower vs. non-CP-knower binned plots for mult and logit models
#save hist counts to get mids for x-lab on barplots (output for these plots can be ignored)
xrange=seq(0,1, by =.1)
c0_reg_hist=hist(reg.p[c==0], breaks=xrange, ylim=c(0,30)) 
c0_reg=hist(reg.p[c==0], breaks=xrange, ylim=c(0,30))$counts 
c0_mult_hist=hist(mult.p[c==0], breaks=xrange, ylim=c(0,30))
c0_mult=hist(mult.p[c==0], breaks=xrange, ylim=c(0,30))$counts
c1_reg_hist=hist(reg.p[c==1], breaks=xrange, ylim=c(0,20))
c1_reg=hist(reg.p[c==1], breaks=xrange, ylim=c(0,20))$counts
c1_mult_hist=hist(mult.p[c==1], breaks=xrange, ylim=c(0,20))
c1_mult=hist(mult.p[c==1], breaks=xrange, ylim=c(0,20))$counts
total_reg=c0_reg+c1_reg
total_mult=c0_mult+c1_mult

#finally draw the binned plots; first two are raw, next two are proportional
barplot(rbind(c1_reg,c0_reg),col=c("red3", "white"), xlab="Predicted CP-knower Probability (Logistic Model)", ylab="CP-knower vs. non-CP-knower Frequency", names.arg=c0_reg_hist$mids, las=1)
barplot(rbind(c1_mult,c0_mult),col=c("red3", "white"),  xlab="Predicted CP-knower Probability (Multiplicative Model)", ylab="CP-knower vs. non-CP-knower Frequency",names.arg=c0_mult_hist$mids, las=1)
barplot(rbind(c1_reg/total_reg,c0_reg/total_reg),col=c("red3", "white"), xlab="Predicted CP-knower Probability (Logistic Model)", ylab="CP-knower vs. non-CP-knower Proportion",names.arg=c1_reg_hist$mids, las=1)
lines(seq(0.5, 11.5,1), seq(0,1,0.0909), lty=2, lwd=2)
barplot(rbind(c1_mult/total_mult,c0_mult/total_mult),col=c("red3", "white"), xlab="Predicted CP-knower Probability (Multiplicative Model)", ylab="CP-knower vs. non-CP-knower Proportion",names.arg=c1_mult_hist$mids, las=1)
lines(seq(0.5, 11.5,1), seq(0,1,0.0909), lty=2, lwd=2)


########################################################################################################
#compute ePCP as per Herron (1999): results are nearly identical for both models 

my.ePCP=function(c,mod.p){
  one_sum=0
  zero_sum=0
  for(i in 1:N){
    if(c[i]==1){
      one_sum=one_sum+mod.p[i]
    }
    else{
      zero_sum=zero_sum+(1-mod.p[i])
    }
  }
  ePCP=1/N*(one_sum+zero_sum)
  ePCP
}

my.ePCP(c,mult.p)
my.ePCP(c,reg.p)
