# Humans vs Zombies Modeling
Code accompanying the article *Humans vs Zombies: Data-driven Modeling of Disease Spread on College Campuses*, published in SPORA: A Journal of Biomathematics. (https://doi.org/10.1186/s12889-021-12275-6)

We present the code for the five ODE models we developed and discussed in the article: HWZ (SIR), Predator-Prey II, Predator-Prey III, HWZC (multiple susceptibility classes), and HWZC + sinusoidal sleep cycles. The data and code are available for 2012 and 2013. In 2012, we split all models in three time periods from t=0h to t=62h, t=62h to t=75h, and t=75h to t=164h. In 2013, we split all models in two time periods from t=0h to t=62h and t=62h to t=115h. All the code for these models is in Python and can be further modified to account for different characteristics in the susceptible and the infected classes. The code for the best-fit model (HWZC + sleep cycles 2012) is commented in detail and all the comments apply to the other models as well.

We also present the code for the simplest model (HZC) in R which gives the same results an the Python model. The program optimizes the parameters for each time period individually and the results should be combined by the user. The R code can be further developed to fit the other models as well.
