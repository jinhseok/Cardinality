Sys.setenv(RGL_USE_NULL=TRUE)
library(shinyRGL)
library(shiny)
library(rgl)
library(grid)

data=read.csv(file="~/Desktop/Analysis/PS_Data/ABCY1Y2MI.csv")
a1=data$y1dqdaccT1/100
t1=data$y1mbT1/100
c1=ifelse(data$y1gnT1>3,1,0)
X1=cbind(a1,t1,c1)
X1=na.omit(X1)
m1=c(1.53,0.64,length(X1[,1]))
X1=cbind(X1,m1)

a2=data$y1dqdaccT2/100
t2=data$y1mbT2/100
c2=ifelse(data$y1gnT2>3,1,0)
X2=cbind(a2,t2,c2)
X2=na.omit(X2)
m2=c(1.25,.37,length(X2[,1]),0)
X2=cbind(X2,m2)

shinyServer(function(input, output, session) {
  output$threedplot<-renderWebGL({
    data=switch(input$radio, "T1"=X1, "T2"=X2)
    a=data[,1]
    t=data[,2]
    c=data[,3]
    ma=data[,4][1]
    mb=data[,4][2]
    N=data[,4][3]
    #3d plot
    cp.prob=1:N #store vector of all p(a^alpha*b^beta)
    color=1:N #assign color based on knower status
    for (i in 1:N) {  
      cp.prob[i]=(a[i]^ma)*(t[i]^mb)
      if (c[i]>0){
        color[i]="springgreen2" #CP-knower
      }else {
        color[i]="coral" #non-CP-knower
      }
    }
    #prediction plane
    surface=function(f, n=10, ...){ 
      ranges=rgl:::.getRanges()
      x=seq(ranges$xlim[1], ranges$xlim[2], length=n)
      y=seq(ranges$ylim[1], ranges$ylim[2], length=n)
      z=outer(x, y, f)
      surface3d(x, y, z, ...)
    }
    #mult model predicted values
    f=function(a,t){
      a^ma*t^mb
    }
    par3d(windowRect=c(0,0,500,500))
    plot3d(a, t, jitter(cp.prob), col=color, xlab="ANS Accuracy", ylab="OTS Accuracy", zlab="Predicted CP-Knower Probability", aspect=TRUE)
    surface(f, alpha=.35)
    #stop R thread upon closing of browser window
    session$onSessionEnded(function(){ 
      stopApp()  
    })
  }) #webGL
}) #server