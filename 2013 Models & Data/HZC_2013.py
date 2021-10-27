#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Humans vs Zombies 2013 Modeling
Code for the HZC Model
-------------------------------------------------------------------------------
Created by: Ognyan Simeonov (Bates College)
Date: 09/24/2021
"""

#Include Needed Packages
import matplotlib.pyplot as plt
import numpy as np
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
dataI = pd.read_csv("HvZData2013.csv")
plt.plot(dataI['Hour'], dataI['Humans'], 'b.', label="Humans Data")
plt.plot(dataI['Hour'], dataI['Zombies'], 'g.', label="Zombies Data")
plt.plot(dataI['Hour'], dataI['Corpses'], 'r.', label="Corpses Data")

#Label the x and y axis.
plt.xlabel('Time (h)')
plt.ylabel('Population')

#------------------------------------------------------------------------------
#DEFINE THE MODELS FOR THE FIRST TIME PERIOD (0-62h)
#First time period from t=0h to t=65h. We keep b=0, thus it is not included in the 
#parameter estimation. Here we use a combination of the warrior and the sleep cycle models.

def ode_model(z, t, k1, b = 0): #define the ODE model
    H, Z, C = z
    #N = H + Z + C
    if H > 0 and Z > 0: #System of ODEs to model the change in populations
        dHdt = -k1*H*Z
        dZdt = k1*H*Z - b*Z
        dCdt = b*Z
    else:
        dHdt = 0
        dZdt = - b*Z
        dCdt = b*Z
    return [dHdt, dZdt, dCdt]

#DEFINE THE ODE SOLVER
# The ODE solver takes as inputs the time span, initial conditions, and parameter
#and uses them to solve the model of ODEs we defined
def ode_solver(t, initial_conditions, params):
    initH, initZ, initC = initial_conditions
    k1, b = params['k1'].value, 0
    res = odeint(ode_model, [initH, initZ, initC], t, args=(k1, b))
    return res

#DEFNE THE ERROR TERM
#This is the difference between what  the ODE model predicts 
#and what we observe in the data.
def error(params, initial_conditions, tspan, data):
    sol0 = ode_solver(tspan, initial_conditions, params)
    return np.matrix([sol0[:, 0]-data[:, 0], sol0[:, 1] - data[:, 1], sol0[:, 2] - data[:, 2]]).ravel()

#------------------------------------------------------------------------------
initH = 115
initZ = 1
initC = 0
initial_conditions = [initH, initZ, initC]

#parameter initial estimates from simulation
k1 = 0.024
k2 = 0.01
b = 0

params = Parameters()
params.add('k1', value=k1, min=0, max=5)

#time span for prediction
days = 62
tspan = np.arange(1, days, 1)

#data of IRD - Infected, Recovered, Dead
data = dataI.loc[1:(days-1), ['Humans', 'Zombies', 'Corpses']].values

# fit model and find predicted values
result = minimize(error, params, args=(initial_conditions, tspan, data), method=method)
# display fitted statistics
print("\nFit report for time period 0-62 h, where b = 0\n")
report_fit(result)

b = 0
k1 = result.params['k1'].value
    
t0 = np.arange(1, 62, 1)
int_cond0 = [initH, initZ, initC]
sol0 = odeint(ode_model, int_cond0, t0, args=(k1, b))

plt.plot(t0,sol0[:,0],'b-',label="Humans Fit")
plt.plot(t0,sol0[:,1],'g-',label="Zombies Fit")
plt.plot(t0,sol0[:,2],'r-',label="Corpses Fit")


#----------------------------------------------------------------------

def ode_model1(z, tspan1, k1, b):
    H, Z, C = z
    #N = H + Z + C
    if H > 0 and Z > 0:
        dHdt = -k1*H*Z
        dZdt = k1*H*Z - b*Z
        dCdt = b*Z
    else:
        dHdt = 0
        dZdt = - b*Z
        dCdt = b*Z
    return [dHdt, dZdt, dCdt]

def ode_solver1(tspan1, initial_conditions1, params1):
    initH1, initZ1, initC1 = initial_conditions1
    k1, b = params1['k1'].value, params1['b'].value
    res = odeint(ode_model1, [initH1, initZ1, initC1], tspan1, args=(k1, b))
    return res

def error1(params1, initial_conditions1, tspan1, data1):
    sol1 = ode_solver1(tspan1, initial_conditions1, params1)
    return np.matrix([(sol1[:, 0])-data1[:, 0], sol1[:, 1] - data1[:, 1], sol1[:, 2] - data1[:, 2]]).ravel()

#initial conditions
initH1 =  sol0[days-2, 0]
initZ1 =  sol0[days-2, 1]
initC1 = sol0[days-2, 2]
initial_conditions1 = [initH1 , initZ1, initC1]

#parameter initial estimates from simulation
k1 = 0.024
b = 0.04

params1 = Parameters()
params1.add('k1', value=k1, min=0, max=5)
params1.add('b', value=b, min=0, max=5)


#time span for prediction
days1 = 115
tspan1 = np.arange(days, days1, 1)

#data of IRD - Infected, Recovered, Dead
data1 = dataI.loc[days:(days1-1), ['Humans', 'Zombies', 'Corpses']].values

# fit model and find predicted values
result1 = minimize(error1, params1, args=(initial_conditions1, tspan1, data1), method=method)

# display fitted statistics
print("\nFit report for time period 62-114 h\n")
report_fit(result1)


t1 = np.arange(62, 115, 1)
initH1 =  sol0[days-2, 0]
initZ1 =  sol0[days-2, 1]
initC1 = sol0[days-2, 2]
int_cond1 = [initH1 , initZ1, initC1]

b =  result1.params['b'].value
k1 = result1.params['k1'].value
    
sol1 = odeint(ode_model1, int_cond1, t1, args=(k1,b))


plt.plot(t1,sol1[:,0],'b-')
plt.plot(t1,sol1[:,1],'g-')
plt.plot(t1,sol1[:,2],'r-')


#----------------------------------------------------------------------------------

plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig('HZC2013_Plot.png',bbox_inches='tight', dpi=300)