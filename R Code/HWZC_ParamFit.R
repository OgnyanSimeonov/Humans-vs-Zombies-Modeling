# set working directory
setwd("/Users/ognyansimeonov/Desktop/Humans vs zombies/Data & Code/R Code")

# load libraries
library(ggplot2) #library for plotting
library(reshape2) # library for reshaping data (tall-narrow <-> short-wide)
library(deSolve) # library for solving differential equations
library(minpack.lm) # library for least squares fit using levenberg-marquart algorithm

#load concentration data
df=read.csv("/Users/ognyansimeonov/Desktop/Humans vs zombies/Data & Code/R Code/2012datacsv.csv")
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
  k2=parms$k2
  b=parms$b
  
  # c is the concentration of species
  
  # derivatives dc/dt are computed below
  r=rep(0,length(c))
  r[1]=-k1*c["H"]*c["Z"]
  r[2]=-k2*c["W"]*c["Z"] #dcA/dt
  r[3]= k1*c["H"]*c["Z"]+k2*c["W"]*c["Z"]-b*c["Z"] #dcB/dt
  r[4]= b*c["Z"] #dcC/dt
  
  
  # the computed derivatives are returned as a list
  # order of derivatives needs to be the same as the order of species in c
  return(list(r))
  
}

# predicted concentration for a given parameter set
cinit=c(H=0.40644866 + 12, W= 9.97929516,Z=52.86935408 - 12,C=6.7449021)
t=df$time
parms=list(k1=5.6276e-04,k2=5.6292e-04,b=0.04687684)
out=ode(y=cinit,times=t,func=rxnrate,parms=parms)

ssq=function(parms){
  
  # inital concentration
  cinit=c(H=0.40644866 + 12, W= 9.97929516,Z=52.86935408 - 12,C=6.7449021)
  # time points for which conc is reported
  # include the points where data is available
  t=c(seq(75,165,1),df$time)
  t=sort(unique(t))
  # parameters from the parameter estimation routine
  k1=parms[1]
  k2=parms[2]
  b=parms[3]
  # solve ODE for a given set of parameters
  out=ode(y=cinit,times=t,func=rxnrate,parms=list(k1=k1,k2=k2,b=b))
  out1 = data.frame(out[,1], out[,2]+out[,3], out[,4], out[,5])
  names(out1)=c("time","H","Z","C")
  # Filter data that contains time points where data is available
  out1df=data.frame(out1)
  out1df=out1df[out1df$time %in% df$time,]
  # Evaluate predicted vs experimental residual
  preddf=melt(out1df,id.var="time",variable.name="species",value.name="conc")
  expdf=melt(df,id.var="time",variable.name="species",value.name="conc")
  ssqres=preddf$conc-expdf$conc
  
  # return predicted vs experimental residual
  return(ssqres)
  
}
# parameter fitting using levenberg marquart algorithm
# initial guess for parameters
parms=c(k1=5.6276e-04,k2=5.6292e-04,b=0.04687684)
# fitting
fitval=nls.lm(par=parms,fn=ssq)

# Summary of fit
summary(fitval)

# Estimated parameter
parest=as.list(coef(fitval))
# degrees of freedom: # data points - # parameters
dof=3*nrow(df)-3
# mean error
ms=sqrt(deviance(fitval)/dof)
# variance Covariance Matrix
S=vcov(fitval)

# plot of predicted vs experimental data

# simulated predicted profile at estimated parameter values
cinit=c(H=0.40644866 + 12, W= 9.97929516,Z=52.86935408 - 12,C=6.7449021)
t=seq(75,165,1)
parms=as.list(parest)
out=ode(y=cinit,times=t,func=rxnrate,parms=parms)
out1 = data.frame(out[,1], out[,2]+out[,3], out[,4], out[,5])
names(out1)=c("time","H","Z","C")
out1df=data.frame(out1)
names(out1df)=c("time","H_pred","Z_pred","C_pred")

# Overlay predicted profile with experimental data
tmppred=melt(out1df,id.var=c("time"),variable.name="species",value.name="conc")
tmpexp=melt(df,id.var=c("time"),variable.name="species",value.name="conc")
p=ggplot(data=tmppred,aes(x=time,y=conc,color=species,linetype=species))+geom_line()
p=p+geom_line(data=tmpexp,aes(x=time,y=conc,color=species,linetype=species))
p=p+geom_point(data=tmpexp,aes(x=time,y=conc,color=species))
p=p+scale_linetype_manual(values=c(0,1,0,1,0,1))
p=p+scale_color_manual(values=rep(c("red","blue","green"),each=2))+theme_bw()
print(p)
