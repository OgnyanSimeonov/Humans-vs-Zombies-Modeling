# Humans-vs-Zombies-Modeling
Code for the article *Humans vs Zombies: Data-driven Modeling of Disease Spread on College Campuses*, submitted for publication in the Bulletin of Mathematical Biology.

We present the code for the five ODE models we developed in the article: HWZ (SIR), Predator-Prey II, Predator-Prey III, HWZC (multiple susceptibility classes), and HWZC + sinusoidal sleep cycles. The data and code are available for 2012 and 2013. In 2012, we split all models in three time periods from t=0h to t=62h, t=62h to t=75h, and t=75h to t=165h. In 2013, we split all models in two time periods from t=0h to t=62h and t=62h to t=115h. All the code for these models is in Python and can be modified to account for different characteristics in the susceptible and the infected classes. 

We also present the code for the simplest model in R as well. This code gives the same results an the Python model. This code optimizes the parameters for each time period invividually and the results should be combined by te user. The code can be modified to fit the other models as well.
