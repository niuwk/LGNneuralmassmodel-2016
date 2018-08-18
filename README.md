# LGNneuralmassmodel-2016
This is the matlab code for the neural mass model used in the article:
“Causal Role of Thalamic Interneurons in Brain State Transitions: A Study Using a Neural Mass Model Implementing Synaptic Kinetics”, 
Frontiers in Computational Neuroscience, vol. 10 (115), 2016.   doi: 10.3389/fncom.2016.00115. 
Authors: Basabdatta Sen Bhattacharya, Tom Bond, Louise O'Hare, Daniel Turner, Simon J. Durrant
Please cite the work when using any or part of the code.

ThalamicmoduleWithIN_2016.m is the main script that defines the model, and also contains code for visualising the time series and power spectral plots shown in Figure 3 of the paper. The model differential equations are solved using Runge-Kutta 4/5th order method, which is manually coded in the function rk45func_thalmodwithin.m, which is called from within the main script.

The time series and power spectra for the case 'Without IN', i.e. disconnecting the IN population, as shown in Figure 4 of the paper, is obtained by putting the parameter Ctii=0 in the main script.
