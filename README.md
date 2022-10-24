# Emmanoulopoulos-Light-Curve-Simulations
Hi, I'm Paolo Cristarella Orestano, a PhD student in Perugia, Italy.

Instructions for light Curve simulations from Emmanoulopoulos algorithm, as according to Emmanoulopoulos et al 2013, Monthly Notices of the Royal Astronomical Society, 433, 907.
The script simulate_lc.py comes from https://github.com/ttshimiz/simulate_lc, for the first step of the algorithm: Timmer Koenig simulation.
While the scripts I personally wrote are named with "pco" at the end.

######################################################

INSTRUCTIONS:
Put all the scripts in the same folder.
data_from_LCR_pco.py and simulate_lc.py are imported in the main script: LC_simulations_pco.py

The guided version of the main script is Help_LC_simulations_pco.py
Try this one to be guided in the simulations.
Launch the script, then follow the instructions that appear on the terminal.

When you will become familiar with the code I suggest you to use the script LC_simulations_pco.py, or import it and then use the sim() function wherever you need it with your data or with data from the Light Curve Repository.

######################################################

Packages used in the script:
numpy
math
matplotlib
scipy
cmath
random
statistics
csv
datetime
requests
json
re

######################################################

Main function:

sim(plot_condition = False, condition_guidance = False, example = False, final_plot = False, time = None, flux = None, flux_error = None, temporalbin = None, PSDmodel = None, PSDparams = None, source_name_title = None, spec = None)
    
    # plot_condition     (boolean)   to show plots at each step of the algorithm.
    # condition_guidance (boolean)   to be guided during the script with descriptions at each step of the algorithm.
    # example            (boolean)   to see an example with PG 1553+113.
    # final_plot         (boolean)   to show just the plot of Emmanoulopoulos simulated Light Curve.
    # time                (array)    time of observations, it is considered in Modified Julian Date (MJD).
    # flux                (array)    flux of the source.
    # flux_error          (array)    flux uncertainty of the source.
    # temporalbin          (int)     sampling time, required with your own data.
    # PSDmodel            (string)   Power Spectral Density (PSD) model for Timmer Konig simulation, if None it is considered a simple unbroken law. Possibilities: unbroken, sharp, slow.
    # PSDparams      (list of float) Parameters for PSD, the numbers and the type depends on the model (4 unbroken, 5 sharp, 5 slow). See below for more details.
    # source_name_title   (string)   Title name of the source for plots if you use your own data.
    # spec          (list of string) Parameters used when you want data from the Light Curve Repository
    ## spec = [ 4FGL name, sampling, flux type, index, TS]

return sim_time, x_Final

    # sim_time (array)  time array (in days, it starts from 0)
    # x_Final  (array)  simulated flux from Emmanoulopoulos algorithm
    
######################################################

Algorithm Steps

------------ First step of the algorithm ------------
Timmer Koenig LC simulation, with the same PSD model and parameters of the source.
From the TK simulation we calculate the Discrete Fourier Function (DFT) and from it we extract the amplitude: A_TK

------------ Second step of the algorithm ------------
From the PDF of the source we reproduce white noise data: x_WN.
Estimate the DFT and hence the phase: phi_loop (it will change at each loop)
(At the fifth step a loop will start from here, taking the last simulated curve, the DFT is estimated and the phase phi_loop taken)

------------ Third step of the algorithm ------------
We take A_TK and phi_loop, from these two we obtain an Adjust DFT, then we estimate x_Adj through the Inverse DFT

------------ Fourth step of the algorithm ------------
We create a new time series with the values of x_WN, at first, ordering them as the order of x_Adj, so obtaining x_loop.
(This means that the highest value of x_Adj is replaced with the highest value of x_WN, the second highest value of x_Adj is replaced with the second highest value of x_WN, and so on, so we keep the same time distribution of x_Adj and the PDF of X_WN)
    
------------ Fifth step of the algorithm ------------
We repeat all the processes from the second step until the two light curves: x_Adj and x_loop, converge.
(second step) we take the phase phi_loop from the DFT of x_loop, (third step) we generate a new x_Adj from phi_loop and A_TK, (fourth step) we sort and replace and at the end we check for the convergence.
when the two series converge it means we have obtained the Emmanoulopoulos simulated light curve x_Final, with the same PSD and PDF of the original source.
In the practice the convergence of the two series happens when the PSD slope converges, usually between 30 and 50 cycles are required for this convergence. (Following the example of Emmanoulopoulos et al 2013)



((((    Sixth step of the algorithm     ))))
This step is designed to work only with photon flux from the LCR.
The last thing to do is the addition of Poisson noise, in order to consider the light curve as a product of a counting detector-process.
Each point of x_Final is replaced by: Pois[mu = x_Final * dt * effective_area] / (dt * effective_area)
As effective_area we use a fixed value of 9000cm^2.
You can change the value or the way it is calculated, or if add or not Poisson noise through cond_poisson (True or False).


######################################################

Contact me if you have questions or to give me advice on the code to make it better, simpler or faster.

paolo.cristarellaorestano@studenti.unipg.it
