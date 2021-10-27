# set working directory
setwd("/Users/ognyansimeonov/Desktop/Humans vs zombies/Data & Code/R Code")

# load libraries
library(ggplot2) #library for plotting
library(reshape2) # library for reshaping data (tall-narrow <-> short-wide)
library(deSolve) # library for solving differential equations
library(minpack.lm) # library for least squares fit using levenberg-marquart algorithm

#load concentration data
df=read.csv("/Users/ognyansimeonov/Desktop/Humans vs zombies/HvZ Data & Code/R Code/2012datacsv.csv")
names(df)=c("time","H","Z","C")
df = df[,c(1,2,3,4)]
df = df[c(75:165),]
# plot data
tmp=melt(df,id.vars=c("time"),variable.name="species",value.name="conc")
ggplot(data=tmp,aes(x=time,y=conc,color=species))+geom_point(size=1)



# rate function
rxnrate=function(t,c,parms){
  
  # rate constant passed through a list called parms
  k1=parms$k1
  b=parms$b
  
  # c is the concentration of species
  
  # derivatives dc/dt are computed below
  r=rep(0,length(c))
  r[1]=-k1*c["A"]*c["B"] #dcA/dt
  r[2]= k1*c["A"]*c["B"]-b*c["B"] #dcB/dt
  r[3]= b*c["B"] #dcC/dt
  
  # the computed derivatives are returned as a list
  # order of derivatives needs to be the same as the order of species in c
  return(list(r))
  
}

# predicted concentration for a given parameter set
cinit=c(A=69,B=1,C=2)
t=df$time
parms=list(k1=0.05,b=0.05)
out=ode(y=cinit,times=t,func=rxnrate,parms=parms)

ssq=function(parms){
  
  # inital concentration
  cinit=c(A=69,B=1,C=2)
  # time points for which conc is reported
  # include the points where data is available
  t=c(seq(75,165,1),df$time)
  t=sort(unique(t))
  # parameters from the parameter estimation routine
  k1=parms[1]
  b=parms[2]
  # solve ODE for a given set of parameters
  out=ode(y=cinit,times=t,func=rxnrate,parms=list(k1=k1,b=b))
  
  # Filter data that contains time points where data is available
  outdf=data.frame(out)
  outdf=outdf[outdf$time %in% df$time,]
  # Evaluate predicted vs experimental residual
  preddf=melt(outdf,id.var="time",variable.name="species",value.name="conc")
  expdf=melt(df,id.var="time",variable.name="species",value.name="conc")
  ssqres=preddf$conc-expdf$conc
  
  # return predicted vs experimental residual
  return(ssqres)
  
}
# parameter fitting using levenberg marquart algorithm
# initial guess for parameters
parms=c(k1=0.05,b=0.05)
# fitting
fitval=nls.lm(par=parms,fn=ssq)

# Summary of fit
summary(fitval)

# Estimated parameter
parest=as.list(coef(fitval))
# degrees of freedom: # data points - # parameters
dof=3*nrow(df)-2
# mean error
ms=sqrt(deviance(fitval)/dof)
# variance Covariance Matrix
S=vcov(fitval)

# plot of predicted vs experimental data

# simulated predicted profile at estimated parameter values
cinit=c(A=69,B=1,C=2)
t=seq(75,165,1)
parms=as.list(parest)
out=ode(y=cinit,times=t,func=rxnrate,parms=parms)
outdf=data.frame(out)
names(outdf)=c("time","H_pred","Z_pred","C_pred")

# Overlay predicted profile with experimental data
tmppred=melt(outdf,id.var=c("time"),variable.name="species",value.name="conc")
tmpexp=melt(df,id.var=c("time"),variable.name="species",value.name="conc")
p=ggplot(data=tmppred,aes(x=time,y=conc,color=species,linetype=species))+geom_line()
p=p+geom_line(data=tmpexp,aes(x=time,y=conc,color=species,linetype=species))
p=p+geom_point(data=tmpexp,aes(x=time,y=conc,color=species))
p=p+scale_linetype_manual(values=c(0,1,0,1,0,1))
p=p+scale_color_manual(values=rep(c("red","blue","green"),each=2))+theme_bw()
print(p)
