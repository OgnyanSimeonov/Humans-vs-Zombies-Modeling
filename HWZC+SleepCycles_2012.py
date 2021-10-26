#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Humans vs Zombies Modeling
Code for the HWZC+Sleep Cycles Model
-------------------------------------------------------------------------------
Created by: Ognyan Simeonov (Bates College)
Date: 09/24/2021
"""
#Include Needed Packages
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
from scipy.integrate import odeint
from lmfit import minimize, Parameters, report_fit

#CHOOSE PARAMETER OPTIMIZATION METHOD
#Uncomment the chosen method and comment out the rest
#------------------------------------
method = "leastsq"
#method = "differential_evolution"
#method = "dual_annealing"
#------------------------------------

#Plot the initial data. Here we use the 2012 HvZ data. 
dataI = pd.read_csv("2012datacsv.csv")
plt.plot(dataI['Hour'], dataI['Humans'], 'b.', label="Humans Data")
plt.plot(dataI['Hour'], dataI['Zombies'], 'g.', label="Zombies Data")
plt.plot(dataI['Hour'], dataI['Corpses'], 'r.', label="Corpses Data")

#Label the x and y axis.
plt.xlabel('Time (h)')
plt.ylabel('Population')

#------------------------------------------------------------------------------
#DEFINE THE MODELS FOR THE FIRST TIME PERIOD (0-62h)
#First time period from t=0h to t=62h. We keep b=0, thus it is not included in the 
#parameter estimation. Here we use a combination of the warrior and the sleep cycle models.

def ode_model(z, t, k1, k2, b = 0): #define the ODE model
    H, W, Z, C = z
    #N = H + Z + C
    if H > 0 and Z > 0: #System of ODEs to model the change in populations
        dHdt = -k1*0.5*(math.sin(t*math.pi/12+2)+1)*H*Z
        dWdt = -k2*0.5*(math.sin(t*math.pi/12+2)+1)*W*Z
        dZdt = k1*0.5*(math.sin(t*math.pi/12+2)+1)*H*Z + k2*0.5*(math.sin(t*math.pi/12+2)+1)*W*Z - b*Z
        dCdt = b*Z
    else:
        dHdt = 0
        dWdt = 0
        dZdt = - b*Z
        dCdt = b*Z
    return [dHdt, dWdt, dZdt, dCdt]

#DEFINE THE ODE SOLVER
# The ODE solver takes as inputs the time span, initial conditions, and parameter
#and uses them to solve the model of ODEs we defined
def ode_solver(t, initial_conditions, params):
    initH, initW, initZ, initC = initial_conditions
    k1, k2 = params['k1'].value, params['k2'].value
    res = odeint(ode_model, [initH, initW, initZ, initC], t, args=(k1, k2))
    return res

#DEFNE THE ERROR TERM
#This is the difference between what  the ODE model predicts 
#and what we observe in the data.
def error(params, initial_conditions, tspan, data):
    sol0 = ode_solver(tspan, initial_conditions, params)
    return np.matrix([(sol0[:, 1] + sol0[:, 0])-data[:, 0], sol0[:, 2] - data[:, 1], sol0[:, 3] - data[:, 2]]).ravel()

#------------------------------------------------------------------------------
#Set up Initial Conditions:
initH = 55 #Initial Human Polpulation
initW = 14 #Initial Warrior Polpulation
initZ = 1 #Initial Zombie Polpulation
initC = 2 #Initial Corpse Polpulation
initial_conditions = [initH, initW, initZ, initC] #Unite all initial conditions

#Set parameter initial estimates
k1 = 0.0024
k2 = 0.0001
b = 0

#Create a variable for the paprmeters and save the initial estimates in it
params = Parameters()
params.add('k1', value=k1, min=0, max=5)
params.add('k2', value=k2, min=0, max=5)


#Set the time span for prediction
days = 62
tspan = np.arange(1, days, 1)

#Exctrace the necessary data for Human, Zombies, and Corpses
data = dataI.loc[1:(days-1), ['Humans', 'Zombies', 'Corpses']].values

#Fit model and find predicted values
result = minimize(error, params, args=(initial_conditions, tspan, data), method=method)

#Display fitted statistics
print("\nFit report for time period 0-62 h, where b = 0\n")
report_fit(result)

#Set the parameters equal to whatever the best-fit parameters the model found
b = 0 #b=o because we are in the first time period while the zombies don't di of starvation yet
k1 = result.params['k1'].value
k2 = result.params['k2'].value

#Solve the ODE system, we are going to use these results for plotting
sol0 = odeint(ode_model, initial_conditions, tspan,args=(k1,k2))

#Plotting the results for time t = 0-62h 
plt.plot(tspan,sol0[:,0]+sol0[:,1],'b-',label="Humans Fit") #Plot for Humans
plt.plot(tspan,sol0[:,2],'g-',label="Zombies Fit")#Plot for Zombies
plt.plot(tspan,sol0[:,3],'r-',label="Corpses Fit") #Plot for Corpses


#------------------------------------------------------------------------------
#DEFINE THE MODELS FOR THE REST OF THE TIME (t=62-165h)
#Second and third time periods from t=62h to t=75h and 7=75h to 165h.
#We try to find best-fit estimates for all parameters( k1, k2, and b)
#Here we use a combination of the warrior and the sleep cycle models.

def ode_model1(z, tspan1, k1, k2, b):
    H, W, Z, C = z
    #N = H + Z + C
    if H > 0 and Z > 0:
        dHdt = -k1*0.5*(math.sin(tspan1*math.pi/12 +2)+1)*H*Z
        dWdt = -k2*0.5*(math.sin(tspan1*math.pi/12 +2)+1)*W*Z
        dZdt = k1*0.5*(math.sin(tspan1*math.pi/12 +2)+1)*H*Z + k2*0.5*(math.sin(tspan1*math.pi/12 +2)+1)*W*Z - b*Z
        dCdt = b*Z
    else:
        dHdt = 0
        dWdt = 0
        dZdt = - b*Z
        dCdt = b*Z
    return [dHdt, dWdt, dZdt, dCdt]

#DEFINE THE ODE SOLVER
# The ODE solver takes as inputs the time span, initial conditions, and parameter
#and uses them to solve the model of ODEs we defined

def ode_solver1(tspan1, initial_conditions1, params1):
    initH1, initW1, initZ1, initC1 = initial_conditions1
    k1, k2, b = params1['k1'].value, params1['k2'].value, params1['b'].value
    res = odeint(ode_model1, [initH1, initW1, initZ1, initC1], tspan1, args=(k1, k2, b))
    return res

#DEFNE THE ERROR TERM
#This is the difference between what  the ODE model predicts 
#and what we observe in the data.
    
def error1(params1, initial_conditions1, tspan1, data1):
    sol1 = ode_solver1(tspan1, initial_conditions1, params1)
    return np.matrix([(sol1[:, 1] + sol1[:, 0])-data1[:, 0], sol1[:, 2] - data1[:, 1], sol1[:, 3] - data1[:, 2]]).ravel()

#The Initial conditions are set automatically. The code looks for the last 
#values of the ODE solution in the first time period and takes those as the new 
#initial conditions for the second time period.
initH1 =  sol0[days-2, 0]
initW1 = sol0[days-2, 1]
initZ1 =  sol0[days-2, 2]
initC1 = sol0[days-2, 3]
initial_conditions1 = [initH1, initW1 , initZ1, initC1]

#Set parameter initial estimates
k1 = 0.024
k2 = 0.001
b = 0.04

#Create a variable for the paprmeters and save the initial estimates in it
params1 = Parameters()
params1.add('k1', value=k1, min=0, max=5)
params1.add('k2', value=k2, min=0, max=5)
params1.add('b', value=b, min=0, max=5)


#Set the time span for the second period
days1 = 75
tspan1 = np.arange(days, days1, 1)

#Extract the necessary data for the second time priod for Humans, Zombies, and Corpses
data1 = dataI.loc[days:(days1-1), ['Humans', 'Zombies', 'Corpses']].values

#Fit model and find predicted values
result1 = minimize(error1, params1, args=(initial_conditions1, tspan1, data1), method=method)

#Display fitted statistics
print("\nFit report for time period 65-75 h\n")
report_fit(result1)

#Set the parameters equal to whatever the best-fit parameters the model found
b =  result1.params['b'].value
k1 = result1.params['k1'].value
k2 = result1.params['k2'].value
 
#Solve the ODE system, we are going to use these results for plotting   
sol1 = odeint(ode_model1, initial_conditions1, tspan1, args=(k1,k2,b))

#Plotting the results for time t = 62-75h 
plt.plot(tspan1,sol1[:,0]+sol1[:,1],'b-')
plt.plot(tspan1,sol1[:,2],'g-')
plt.plot(tspan1,sol1[:,3],'r-')


#----------------------------------------------------------------------------------
#TIME PERIOD t=75h-165h. THIS TIME PERIOD USES MODEL 1 AS WELL 
#BECAUSE WE NEED TO ESTIMATE ALL 3 PAREMETERS

#The Initial conditions are set automatically. The code looks for the last 
#values of the ODE solution in the second time period (tspan1) and takes those as the new 
#initial conditions for the second time period.
initH1 =  sol1[len(tspan1)-1, 0]
initW1 = sol1[len(tspan1)-1, 1]
initZ1 =  sol1[len(tspan1)-1, 2]
initC1 = sol1[len(tspan1)-1, 3]
initial_conditions1 = [initH1+12, initW1, initZ1-12, initC1]

#Parameter initial estimates from
k1 = 0.024
k2 = 0.001
b = 0.04

#Create a variable for the paprmeters and save the initial estimates in it
params1 = Parameters() #we redefine params1 because the model expects to receive a variable called params1
params1.add('k1', value=k1, min=0, max=5)
params1.add('k2', value=k2, min=0, max=5)
params1.add('b', value=b, min=0, max=5)


#Time span for prediction
days1 = 165
tspan1 = np.arange(75, days1, 1) #We redefine tspan1 because the model expects to receive a var called tspan1 

#Extract the necessary data for the second time priod for Humans, Zombies, and Corpses
data1 = dataI.loc[75:(days1-1), ['Humans', 'Zombies', 'Corpses']].values

#Fit model and find predicted values
result2 = minimize(error1, params1, args=(initial_conditions1, tspan1, data1), method=method)

#Display fitted statistics
print("\nFit report for time period 75-165 h\n")
report_fit(result2)

#Set the parameters equal to whatever the best-fit parameters the model found
b = result2.params['b']
k1 = result2.params['k1']
k2 = result2.params['k2']
 
#Solve the ODE system, we are going to use these results for plotting     
sol2 = odeint(ode_model1, initial_conditions1, tspan1, args=(k1, k2, b))

#Plotting the results for time t = 75-165h 
plt.plot(tspan1,sol2[:,0]+sol2[:,1],'b-')
plt.plot(tspan1,sol2[:,2],'g-')
plt.plot(tspan1,sol2[:,3],'r-')

#Add a legend to the final plot and save as'TIMEW2012_Plot.png' in the working directory
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig('TIMEW2012_Plot.png',bbox_inches='tight', dpi=300)