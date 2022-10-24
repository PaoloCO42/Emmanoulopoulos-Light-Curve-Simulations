# Create artificial light curves with the algorithm from Emmanoulopoulos et al., 2013, Monthly Notices of the Royal Astronomical Society, 433, 907.
import simulate_lc    # to create Timmer Koenig light curves
import numpy as np  
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import PercentFormatter
from matplotlib import colors
from matplotlib.colors import LogNorm
import scipy
import scipy.stats as st
from scipy.optimize import curve_fit
from scipy import signal
from scipy import interpolate
import math
from math import e
import cmath
import random
from statistics import mode
import data_from_LCR_pco
params  =  {'xtick.labelsize':18,'ytick.labelsize':18}
pylab.rcParams.update(params)

def answer_function(answer):
    if answer == 'no' or answer == 'not' or answer == 'n' or answer == 'no/n' or answer == 'nope':
        risposta = False
    elif answer == 'yes' or answer == 'y' or answer == 'yes/y':
        risposta = True
    return risposta

condition_guidance = False
example = False
#''' uncomment this line if you don't want this first helpful interaction
print('Do you need guidance? (yes/y or no/n)')
answer_input = input()
condition_guidance = answer_function(answer_input)

if not condition_guidance:
    print(' ')
    print('Use the function sim( )')
else:
    guidance_cycle_condition = True
    while True:
        print(' ')
        print('Do you want to see an example with PG 1553+113 from the Fermi Light Curve Repository (LCR)?')
        answer_input = input()
        if answer_function(answer_input):
            example = True
            condition_spec = False
            break
        example = False
        print(' ')
        print('You can choose another source from the LCR.')
        print('If you want use your data answer no.')
        print('Otherwise answer yes')
        answer_input = input()
        if not answer_function(answer_input):
            condition_spec = False
            print(' ')
            print('I suggest you to comment this help section and work importing the script and then using the function sim(...)')
            print(' ')
            print('You can continue writing: ')
            print(' ')
            print('sim( plot_condition , condition_guidance = False ,example = False , time = None , flux = None , flux_error  =  None, temporalbin = None , PSDmodel = None , PSDparams = None , source_name_title = None , spec = None )')
            print(' ')
            print('Do you need informations on the parameters of the function?')
            answer_input = input()
            if answer_function(answer_input):
                print(' ')
                print('(i)  plot_condition: True if you want to see the plots for every step of the algorithm')
                print(' ')
                print('(ii)  condition_guidance: True if you need help during the steps of the algorithm')
                print(' ')
                print('(iii)  example: True if you want the preset example with PG 1553+113')
                print(' ')
                print('(iv)  time: numpy array of the time of the original light curve (LC). None if you use the Light Curve Repository (see the last parameter)')
                print(' ')
                print('(v)  flux: numpy array of the flux of the original LC. None if you use the LCR (see the last parameter)')
                print(' ')
                print('(vi)  flux_error: numpy array of the flux uncertainty of the original LC. None if you use the LCR (see the last parameter)')
                print(' ')
                print('(vii)  temporalbin: float, LC sampling, the delta time between two observations (LC points), in days. None if you use the LCR (see the last parameter)')
                print(' ')
                print('(viii)  PSDmodel: string, choose the model you want to fit the PSD of the LC. (More info during the algorithm if condition_guidance = True) ')
                print(' ')
                print('(ix)  PSDparams: list [ ... , ...   ... ] (float inside) PSD parameters for the chosen model. (More info during the algorithm if condition_guidance = True) ')
                print(' ')
                print('(x)  source_name_title: string, the name that will appear on the original LC plot')
                print(' ')
                print('(xi)  spec: list [ ... , ...   ... ] (string inside), source specification needed if you want data from the LCR. In the list there must be, in this order: 4FGL name, sampling (monthly, weekly or daily), flux type (photon or energy), index type (fixed or free), minimum TS (1, 2, 3 or 4)')
            break
        else:
            condition_spec = True
            print(' ')
            print('4FGL name of the source? (for example: 4FGL J1555.7+1111)')
            nomesorgente = input()
            print(' ')
            print('Sampling? (answer one: monthly, weekly, daily)')
            Cadenza = input()
            print(' ')
            print('Type of flux? (answer one: photon, energy)')
            TipoFlusso = input()
            print(' ')
            print('Index type? (answer one: fixed, free)')
            TipoIndice = input()
            print(' ')
            print('Minimum TS? (answer one: 4, 3, 2, 1)')
            TSminimo = input()
            spec = [nomesorgente,Cadenza,TipoFlusso,TipoIndice,TSminimo]
            break
#'''

def sim(plot_condition = False, condition_guidance = False, example = False, final_plot = False, time = None, flux = None, flux_error = None, temporalbin = None, PSDmodel = None, PSDparams = None, source_name_title = None, spec = None):
    
    # plot_condition     (boolean)   to show plots at each step of the algorithm.
    # condition_guidance (boolean)   to be guided during the script with descriptions at each step of the algorithm.
    # example            (boolean)   to see an example with PG 1553+113.
    # final_plot         (boolean)   to show just the plot of Emmanoulopoulos simulated Light Curve.
    # time             (numpy array) time of observations, it is considered in Modified Julian Date (MJD).
    # flux             (numpy array) flux of the source.
    # flux_error       (numpy array) flux uncertainty of the source.
    # temporalbin          (int)     sampling time, required with your own data.
    # PSDmodel            (string)   Power Spectral Density (PSD) model for Timmer Konig simulation, if None it is considered a simple unbroken law. Possibilities: unbroken, sharp, slow.
    # PSDparams      (list of float) Parameters for PSD, the numbers and the type depends on the model (4 unbroken, 5 sharp, 5 slow). See below for more details.
    # source_name_title   (string)   Title name of the source for plots if you use your own data.
    # spec          (list of string) Parameters used when you want data from the Light Curve Repository
    ## spec = [ 4FGL name, sampling, flux type, index, TS]
    
    #######################
    
    #  return  sim_time, x_Final
    
    ULflux = None
    ULtime = None
    
    if flux != None:
        owndata = True
        cond_poisson = False
        label_flux = 'Flux'
        if source_name_title == None:
            nometitolo = ' '
    else:
        owndata = False
        
    print(' ')
    if example:
        flux,flux_error,time,ULflux,ULtime = data_from_LCR_pco.data_from_LCR()
        nometitolo = ' of PG 1553+113'
        label_flux = 'Photon Flux $(ph$ $cm^{-2}$$s^{-1})$'
        temporalbin = 30
        cond_poisson = True
        print('Light Curve of PG 1553+113, Photon Flux, fixed index, monthly sampling, minimum TS  =  4')
    
    if spec != None:
        condition_spec = True
        flux,flux_error,time,ULflux,ULtime = data_from_LCR_pco.data_from_LCR(spec[0],spec[1],spec[2],spec[3],spec[4])
        nometitolo = ' of '+spec[0]
        if spec[2] == 'photon':
            label_flux = 'Photon Flux $(ph$ $cm^{-2}$$s^{-1})$'
            cond_poisson = True
        else:
            label_flux =  'Energy Flux $(MeV$ $cm^{-2}$$s^{-1})$'
            cond_poisson = False
        print('Light Curve of '+spec[0]+', '+spec[2]+' Flux, '+spec[3]+' index, '+spec[1]+' sampling, minimum TS  =  '+spec[4])
        if spec[1] == 'monthly':
            temporalbin = 30
        elif spec[1] == 'weekly':
            temporalbin = 7
        elif spec[1] == 'daily':
            temporalbin = 3
    else:
        condition_spec = False
        
    if owndata and temporalbin == None:
        delta_time = []
        for i in range(0,(len(time)-1)):
            delta_time.append(time[i+1]-time[i])
        sampling_days = mode(delta_time)
        temporalbin = sampling_days
    elif owndata and time == None:
        time = np.arange(0,len(flux)*temporalbin,temporalbin)
        temporalbin = sampling_days
    
    if condition_guidance and PSDmodel == None:
        print('Do you need info on Power Spectral Density models and parameters?')
        answer_input = input()
        if answer_function(answer_input):
            print(' ')
            print('Models and parameters summary')
            print(' ')
            print('For unbroken, S  =  N*(nu/nu_0)**(-beta) + noise:')
            print('params[0]  =  N, normalization of power spectrum')
            print('params[1]  =  nu_0, the frequency at which the power  ==  the normalization')
            print('params[2]  =  beta, the slope of the power spectrum')
            print('params[3]  =  white noise level to add in')
            print(' ')
            print('For sharp, S  =  N*(nu/nu_c)^(-gamma) + noise      nu < nu_c')
            print('                   S  =  N*(nu/nu_c)^(-beta)  + noise      nu > nu_c')
            print('params[0]  =  N, normalization of power spectrum')
            print('params[1]  =  nu_c, the pivot frequency')
            print('params[2]  =  gamma, the low frequency slope')
            print('params[3]  =  beta, the high frequency slope')
            print('params[4]  =  white noise level')
            print(' ')
            print('For slow, S  =  N*nu^(alpha_lo)/(1+(nu/nu_k))^(alpha_hi - alpha_lo)')
            print('params[0]  =  N, normalization of power spectrum')
            print('params[1]  =  nu_k, the knee frequency where the spectrum rolls over')
            print('params[2]  =  alpha_lo, the low frequency slope')
            print('params[3]  =  alpha_high, the high frequency slope')
            print('params[4]  =  white noise level')
        print(' ')
        print('Choose a PSD model (unbroken power law, sharply broken powerlaw, slow bended knee), from model and parameters, a simulated LC is created based on paper Timmer, J., & Koenig, M. 1995, A&A, 300, 707 ')
        print(' ')
        print('PSDmodel: unbroken, sharp, slow.')
        print(' ')
        print('Which PSD model do you want?')
        PSDmodel=input()
            
    samplingfreq = 1/temporalbin
    simpleFreq, simplePSD = scipy.signal.periodogram(flux, samplingfreq, scaling = 'density')
    maskFreq = simpleFreq > 0
    simpleFreq = simpleFreq[maskFreq]
    simplePSD = simplePSD[maskFreq]
    
    z = np.polyfit(np.log10(simpleFreq),np.log10(simplePSD),1)
    p1 = z[0]   # angular coefficient
    p0 = z[1]   # intercept
    
    # Source plots: Light Curve, Probability Density Function, Power Spectral Density
    if plot_condition:
    
        a_fig = plt.figure(figsize = (16,9),tight_layout = True)  
        a_ax1 = plt.subplot(211)  # LC
        a_ax2 = plt.subplot(223)  # PDF
        a_ax3 = plt.subplot(224)  # PSD
        
        a_ax1.plot(time,flux,'b.')
        a_ax1.plot(ULtime,ULflux,'rv')
        a_ax1.errorbar(time,flux,yerr = flux_error,xerr = None,fmt = 'b.',ecolor = 'b')
        a_ax1.ticklabel_format(style = 'sci',axis = 'y',scilimits = (0,0))
        a_ax1.ticklabel_format(style = 'sci',axis = 'x',scilimits = (0,0))
        a_ax1.set_ylim(0,max(flux)*1.25) 
        a_ax1.set_xlabel('$MJD$',fontsize = 20)
        a_ax1.set_ylabel(label_flux,fontsize = 20) 
        if not source_name_title == None:
            nometitolo=' of '+source_name_title
        a_ax1.set_title('Light Curve'+nometitolo,fontsize = 25,pad = 12)
        
        binPDF = int(round(math.sqrt(len(flux)),0))
        a_ax2.hist(flux,bins = binPDF,density = True,histtype = u'step',alpha = 1,color = 'royalblue',linewidth = 3) # alpha = 0.7 histtype = u'step',
        a_ax2.set_title('Probability Density Function',fontsize = 20,pad = 12)
        a_ax2.set_xlabel(label_flux,fontsize = 20)
        a_ax2.set_ylabel('PDF',fontsize = 20)
        
        a_txt='$\delta$='+str(round(p1,2))
        a_coeffp=np.poly1d(z)
        a_ax3.plot(simpleFreq,simplePSD,'k.')
        a_ax3.plot(simpleFreq,10**(a_coeffp(np.log10(simpleFreq))),'r-', linewidth=2)
        a_ax3.set_xlabel('$\\nu$ $(days^{-1})$',fontsize = 20)
        a_ax3.set_ylabel('PSD',fontsize = 20)
        a_ax3.set_title('Power Spectral Density ('+a_txt+')',fontsize = 20,pad = 12)
        a_ax3.set_xscale('log')
        a_ax3.set_yscale('log')
        
    # Set n, the number of timesteps in the simulated light curve, and sampling_seconds, the length of each timestep in seconds.
    delta_time = []
    for i in range(0,(len(time)-1)):
        delta_time.append(time[i+1]-time[i])
        
    sampling_days = mode(delta_time)   # 30, 7, 3 days if it comes from the Fermi Light Curve Repository (LCR)
    sampling_seconds  = sampling_days*24*3600
    print('sampling_seconds: ',sampling_seconds)
    
    n = (time[-1]-time[0])/sampling_days 
    
    # Time array for simulations
    sim_time = np.arange(0,sampling_days*n,sampling_days)
    
    # Set the mean for the simulated light curve 
    mean_lc = np.mean(flux)
    
    # Due to the steps of the algorithm the slope is the only parameter that matters
    if PSDmodel == None:
        PSDmodel = 'unbroken'
        PSDparams = [1,1,-(p1),0.0]
        print('Unbroken PSD, parameters: ', PSDparams)
    elif PSDmodel == 'unbroken':
        if PSDparams == None:
            PSDparams = [1,1,-(p1),0.0]
            print('Unbroken PSD, parameters: ', PSDparams)
    elif PSDmodel == 'sharp':
        if PSDparams == None:
            PSDparams = [1,1,-(p1*0.75),-(p1),0.0]
            print('Sharp PSD, parameters: ', PSDparams)
    elif PSDmodel == 'slow':
        if PSDparams == None:
            PSDparams = [1,1,-(p1*0.75),-(p1),0.0]
            print('Slow PSD, parameters: ', PSDparams)
    
    ########################################################
    
    print(' ')
    print('------------ First step of the algorithm ------------')
    print('Timmer Koenig LC simulation, with the same PSD model and parameters of the source')
    # Simulate the light curve as in Timmer Koenig 1995
    x_TK  =  simulate_lc.lc_sim(int(n), sampling_days, mean_lc, PSDmodel, PSDparams)

    # TK simulation plots: LC, PDF, PSD
    if plot_condition:
        samplingfreq = 1/temporalbin
        simpleFreq, simplePSD = scipy.signal.periodogram(x_TK, samplingfreq, scaling = 'density')
        maskFreq = simpleFreq > 0
        simpleFreq = simpleFreq[maskFreq]
        simplePSD = simplePSD[maskFreq]
        z = np.polyfit(np.log10(simpleFreq),np.log10(simplePSD),1)
        p1 = z[0]   # angular coefficient
        p0 = z[1]   # intercept
            
        b_fig = plt.figure(figsize = (16,9),tight_layout = True)  
        b_ax1 = plt.subplot(211)  # LC
        b_ax2 = plt.subplot(223)  # PDF
        b_ax3 = plt.subplot(224)  # PSD
        
        b_ax1.plot(sim_time,x_TK,'b.')
        b_ax1.plot(sim_time,x_TK,'b-',linewidth = 1,alpha = 0.25)
        b_ax1.ticklabel_format(style = 'sci',axis = 'y',scilimits = (0,0))
        b_ax1.ticklabel_format(style = 'sci',axis = 'x',scilimits = (0,0))
        if min(x_TK)>0:
            min_norm = 0
        else:
            min_norm = 1.25
        b_ax1.set_ylim(min(x_TK+mean_lc)*min_norm,max(x_TK+mean_lc)*1.25) 
        b_ax1.set_xlabel('$t$ $(days)$',fontsize = 20)
        b_ax1.set_ylabel('Flux',fontsize = 20) #fluxname
        b_ax1.set_title('Timmer König Simulated Light Curve',fontsize = 25,pad = 12)
        
        binPDF = int(round(math.sqrt(len(x_TK)),0))
        b_ax2.hist(x_TK,bins = binPDF,density = True,histtype = u'step',alpha = 1,color = 'royalblue',linewidth = 3) # alpha = 0.7 histtype = u'step',
        b_ax2.ticklabel_format(style = 'sci',axis = 'y',scilimits = (0,0))
        b_ax2.set_title('Probability Density Function',fontsize = 20,pad = 12)
        b_ax2.set_xlabel('Flux',fontsize = 20)
        b_ax2.set_ylabel('PDF',fontsize = 20)
        
        
        b_txt='$\delta$='+str(round(p1,2))
        b_coeffp=np.poly1d(z)
        b_ax3.plot(simpleFreq,simplePSD,'k.')
        b_ax3.plot(simpleFreq,10**(b_coeffp(np.log10(simpleFreq))),'r-', linewidth=2)
        b_ax3.set_xlabel('$\\nu$ $(days^{-1})$',fontsize = 20)
        b_ax3.set_ylabel('PSD',fontsize = 20)
        b_ax3.set_title('Power Spectral Density ('+b_txt+')',fontsize = 20,pad = 12)
        b_ax3.set_xscale('log')
        b_ax3.set_yscale('log')
        
        
    print(' ')
    print('From the TK simulation we calculate the Discrete Fourier Function (DFT) and from it we extract the amplitude: A_TK')
    
    DFTxtk = np.fft.fft(x_TK)
    phi_TK = np.ones(len(DFTxtk.real))
    A_TK = np.ones(len(DFTxtk.real))
    for i in range(0,len(DFTxtk.real)):
        A_TK[i], phi_TK[i] = cmath.polar(DFTxtk[i])
    
    Counts, bins_counts = np.histogram(flux, bins = int(round(math.sqrt(len(flux)),0)))
    pdf_source = Counts / sum(Counts)
    cdf_source = np.cumsum(pdf_source)
    cdf_source = np.insert(cdf_source, 0, 0.0)
    inverseCDF = interpolate.interp1d(cdf_source,bins_counts)
        
    x_loop = 0
    j_index = 0
    index_converge = 0
    converge = False    
    while not converge:
        j_index += 1
        
        ########################################################
        
        if j_index == 1:
            print(' ')
            print('------------ Second step of the algorithm ------------')
            print('From the PDF of the source we reproduce white noise data: x_WN.')
            print('Estimate the DFT and hence the phase: phi_loop (it will change at each loop)')
        
            zero_one = np.random.uniform(0, 1, len(sim_time))
            x_WN = inverseCDF(zero_one)
            DFTxwn = np.fft.fft(x_WN)
            phi_loop = np.ones(len(DFTxwn.real))
            A_loop = np.ones(len(DFTxwn.real))
            for i in range(0,len(DFTxwn.real)):
                A_loop[i], phi_loop[i] = cmath.polar(DFTxwn[i])
            
        else: # after the first step we will take the final simulation of the algorithm (5th step)
            DFTloop = np.fft.fft(x_loop)
            phi_loop = np.ones(len(DFTloop.real))
            A_loop = np.ones(len(DFTloop.real))
            for i in range(0,len(DFTloop.real)):
                A_loop[i], phi_loop[i] = cmath.polar(DFTloop[i])
            
        
        if j_index == 1 and plot_condition:    
            if plot_condition:
                samplingfreq = 1/temporalbin
                simpleFreq, simplePSD = scipy.signal.periodogram(x_WN, samplingfreq, scaling = 'density')
                maskFreq = simpleFreq > 0
                simpleFreq = simpleFreq[maskFreq]
                simplePSD = simplePSD[maskFreq]
                z = np.polyfit(np.log10(simpleFreq),np.log10(simplePSD),1)
                p1 = z[0]   # angular coefficient
                p0 = z[1]   # intercept
                    
                c_fig = plt.figure(figsize = (16,9),tight_layout = True)  
                c_ax1 = plt.subplot(211)  # LC
                c_ax2 = plt.subplot(223)  # PDF
                c_ax3 = plt.subplot(224)  # PSD
                
                c_ax1.plot(sim_time,x_WN,'b.')
                c_ax1.plot(sim_time,x_WN,'b-',linewidth = 1,alpha = 0.25)
                c_ax1.ticklabel_format(style = 'sci',axis = 'y',scilimits = (0,0))
                c_ax1.ticklabel_format(style = 'sci',axis = 'x',scilimits = (0,0))
                if min(x_WN)>0:
                    min_norm = 0
                else:
                    min_norm = 1.25
                c_ax1.set_ylim(0,max(x_WN)*1.25) 
                c_ax1.set_xlabel('$t$ $(days)$',fontsize = 20)
                c_ax1.set_ylabel('Flux',fontsize = 20) #fluxname
                c_ax1.set_title('First White Noise Simulated Light Curve',fontsize = 25,pad = 12)
                
                binPDF = int(round(math.sqrt(len(x_WN)),0))
                c_ax2.hist(x_WN,bins = binPDF,density = True,histtype = u'step',alpha = 1,color = 'royalblue',linewidth = 3) # alpha = 0.7 histtype = u'step',
                c_ax2.ticklabel_format(style = 'sci',axis = 'y',scilimits = (0,0))
                c_ax2.set_title('Probability Density Function',fontsize = 20,pad = 12)
                c_ax2.set_xlabel('Flux',fontsize = 20)
                c_ax2.set_ylabel('PDF',fontsize = 20)
                
                
                c_txt='$\delta$='+str(round(p1,2))
                c_coeffp=np.poly1d(z)
                c_ax3.plot(simpleFreq,simplePSD,'k.')
                c_ax3.plot(simpleFreq,10**(c_coeffp(np.log10(simpleFreq))),'r-', linewidth=2)
                c_ax3.set_xlabel('$\\nu$ $(days^{-1})$',fontsize = 20)
                c_ax3.set_ylabel('PSD',fontsize = 20)
                c_ax3.set_title('Power Spectral Density ('+c_txt+')',fontsize = 20,pad = 12)
                c_ax3.set_xscale('log')
                c_ax3.set_yscale('log')
        
        
        ########################################################
        
        if j_index == 1:
            print(' ')
            print('------------ Third step of the algorithm ------------')
            print('We take A_TK and phi_loop, from these two we obtain an Adjust DFT, then we estimate x_Adj through the Inverse DFT')
        
        DFT_Adj = np.zeros(len(A_TK))*1j
        for i in range(0,len(A_TK)):
            DFT_Adj[i] = cmath.rect(A_TK[i],phi_loop[i])
        
        x_Adj = np.fft.ifft(DFT_Adj)
        x_Adj=x_Adj.real
        
        if j_index == 1 and plot_condition:
            samplingfreq = 1/temporalbin
            simpleFreq, simplePSD = scipy.signal.periodogram(x_Adj, samplingfreq, scaling = 'density')
            maskFreq = simpleFreq > 0
            simpleFreq = simpleFreq[maskFreq]
            simplePSD = simplePSD[maskFreq]
            z = np.polyfit(np.log10(simpleFreq),np.log10(simplePSD),1)
            p1 = z[0]   # angular coefficient
            p0 = z[1]   # intercept
                
            d_fig = plt.figure(figsize = (16,9),tight_layout = True)  
            d_ax1 = plt.subplot(211)  # LC
            d_ax2 = plt.subplot(223)  # PDF
            d_ax3 = plt.subplot(224)  # PSD
            
            d_ax1.plot(sim_time,x_Adj,'b.')
            d_ax1.plot(sim_time,x_Adj,'b-',linewidth = 1,alpha = 0.25)
            d_ax1.ticklabel_format(style = 'sci',axis = 'y',scilimits = (0,0))
            d_ax1.ticklabel_format(style = 'sci',axis = 'x',scilimits = (0,0))
            if min(x_Adj)>0:
                min_norm = 0
            else:
                min_norm = 1.25
            d_ax1.set_ylim(min(x_Adj)*min_norm,max(x_Adj)*1.25) 
            d_ax1.set_xlabel('$t$ $(days)$',fontsize = 20)
            d_ax1.set_ylabel('Flux',fontsize = 20) #fluxname
            d_ax1.set_title('First Adjusted ($A_{TK}$, $\phi_{loop}$) Simulated Light Curve',fontsize = 25,pad = 12)
            
            binPDF = int(round(math.sqrt(len(x_Adj)),0))
            d_ax2.hist(x_Adj,bins = binPDF,density = True,histtype = u'step',alpha = 1,color = 'royalblue',linewidth = 3) # alpha = 0.7 histtype = u'step',
            d_ax2.ticklabel_format(style = 'sci',axis = 'y',scilimits = (0,0))
            d_ax2.set_title('Probability Density Function',fontsize = 20,pad = 12)
            d_ax2.set_xlabel('Flux',fontsize = 20)
            d_ax2.set_ylabel('PDF',fontsize = 20)
            
            d_txt='$\delta$='+str(round(p1,2))
            d_coeffp=np.poly1d(z)
            d_ax3.plot(simpleFreq,simplePSD,'k.')
            d_ax3.plot(simpleFreq,10**(d_coeffp(np.log10(simpleFreq))),'r-', linewidth=2)
            d_ax3.set_xlabel('$\\nu$ $(days^{-1})$',fontsize = 20)
            d_ax3.set_ylabel('PSD',fontsize = 20)
            d_ax3.set_title('Power Spectral Density ('+d_txt+')',fontsize = 20,pad = 12)
            d_ax3.set_xscale('log')
            d_ax3.set_yscale('log')
    
        ########################################################

       
        if j_index == 1:
            print(' ')
            print('------------ Fourth step of the algorithm ------------')
            print('We create a new time series with the values of x_WN, at first, ordering them as the order of x_Adj. Obtaining x_loop.')
            print('(This means that the highest value of x_Adj is replaced with the highest value of x_WN, the second highest value of x_Adj is replaced with the second highest value of x_WN, and so on, so we keep the same time distribution)')
        
        
        perm_Adj = x_Adj.argsort()
        sim_time = sim_time[perm_Adj]
        perm_WN = x_WN.argsort()
        x_WN_sorted = x_WN[perm_WN]
        perm_time = sim_time.argsort()
        sim_time = sim_time[perm_time]
        x_loop = x_WN_sorted[perm_time]
        
        samplingfreq = 1/temporalbin
        simpleFreq, simplePSD = scipy.signal.periodogram(x_loop, samplingfreq, scaling = 'density')
        maskFreq = simpleFreq > 0
        simpleFreq = simpleFreq[maskFreq]
        simplePSD = simplePSD[maskFreq]
        z = np.polyfit(np.log10(simpleFreq),np.log10(simplePSD),1)
        p1 = z[0]   # angular coefficient
        p0 = z[1]   # intercept
        p1_loop = p1
        
        if j_index == 1 and plot_condition:
        
            e_fig = plt.figure(figsize = (16,9),tight_layout = True)  
            e_ax1 = plt.subplot(211)  # LC
            e_ax2 = plt.subplot(223)  # PDF
            e_ax3 = plt.subplot(224)  # PSD
            
            e_ax1.plot(sim_time,x_loop,'b.')
            e_ax1.plot(sim_time,x_loop,'b-',linewidth = 1,alpha = 0.25)
            e_ax1.ticklabel_format(style = 'sci',axis = 'y',scilimits = (0,0))
            e_ax1.ticklabel_format(style = 'sci',axis = 'x',scilimits = (0,0))
            if min(x_loop)>0:
                min_norm = 0
            else:
                min_norm = 1.25
            e_ax1.set_ylim(min(x_loop)*min_norm,max(x_loop)*1.25) 
            e_ax1.set_xlabel('$t$ $(days)$',fontsize = 20)
            e_ax1.set_ylabel('Flux',fontsize = 20) #fluxname
            e_ax1.set_title('First Loop Simulated Light Curve',fontsize = 25,pad = 12)
            
            binPDF = int(round(math.sqrt(len(x_loop)),0))
            e_ax2.hist(x_loop,bins = binPDF,density = True,histtype = u'step',alpha = 1,color = 'royalblue',linewidth = 3) # alpha = 0.7 histtype = u'step',
            e_ax2.ticklabel_format(style = 'sci',axis = 'y',scilimits = (0,0))
            e_ax2.set_title('Probability Density Function',fontsize = 20,pad = 12)
            e_ax2.set_xlabel('Flux',fontsize = 20)
            e_ax2.set_ylabel('PDF',fontsize = 20)
            
            e_txt='$\delta$='+str(round(p1,2))
            e_coeffp=np.poly1d(z)
            e_ax3.plot(simpleFreq,simplePSD,'k.')
            e_ax3.plot(simpleFreq,10**(e_coeffp(np.log10(simpleFreq))),'r-', linewidth=2)
            e_ax3.set_xlabel('$\\nu$ $(days^{-1})$',fontsize = 20)
            e_ax3.set_ylabel('PSD',fontsize = 20)
            e_ax3.set_title('Power Spectral Density ('+e_txt+')',fontsize = 20,pad = 12)
            e_ax3.set_xscale('log')
            e_ax3.set_yscale('log')

        ########################################################
    
        if j_index == 1:
            p1_pre = 0
            print(' ')
            print('------------ Fifth step of the algorithm ------------')
            print('We repeat all the processes from the second step until the two light curves: x_Adj and x_loop, converge.')
            print('(second step) we take the phase phi_loop from the DFT of x_loop, (third step) we generate a new x_Adj from phi_loop and A_TK, (fourth step) we sort and replace and at the end we check for the convergence.')
            print('when the two series converge it means we have obtained the Emmanoulopoulos simulated light curve, with the same PSD and PDF of the original source.')
            print(' ')
            print('Looking for convergence...')
            
        # converging condition based on the PSD
        if p1_pre == p1_loop: 
            index_converge += 1
            
        if index_converge >= 5:
            converge = True
            
        # "extreme" condition, the convergence should be already reached
        if j_index>100: 
            if round(p1_pre,8) == round(p1_loop,8):
                converge = True
             
        p1_pre = p1_loop
            
        if converge:
            x_Final = x_loop
    
    if condition_guidance and cond_poisson:
        print(' ')
        print('------------ Last (sixth) step of the algorithm ------------')
        print('The last thing to do is the addition of Poisson noise, in order to consider the light curve as a product of a counting detector-process.')
        print('Each point of x_Final is replaced by: Pois[mu = x_Final * dt * effective_area] / (dt * effective_area)')
        print('(effective_area =  9000 cm**2 was used)')
        print(' ')    
    # The last step should consist in the add of the poissonian noise. Commented below how to add it.
    if cond_poisson:
        effective_area=9000
        for i in range(0,len(x_Final)):
            x_Final[i] = np.random.poisson(x_Final[i]*sampling_seconds*effective_area)/(sampling_seconds*effective_area)
        
    if plot_condition:
        final_plot = True
    if final_plot:
        samplingfreq = 1/temporalbin
        simpleFreq, simplePSD = scipy.signal.periodogram(x_Final, samplingfreq, scaling = 'density')
        maskFreq = simpleFreq > 0
        simpleFreq = simpleFreq[maskFreq]
        simplePSD = simplePSD[maskFreq]
        z = np.polyfit(np.log10(simpleFreq),np.log10(simplePSD),1)
        p1 = z[0]   # angular coefficient
        p0 = z[1]   # intercept
            
        f_fig = plt.figure(figsize = (16,9),tight_layout = True)  
        f_ax1 = plt.subplot(211)  # LC
        f_ax2 = plt.subplot(223)  # PDF
        f_ax3 = plt.subplot(224)  # PSD
        
        f_ax1.plot(sim_time,x_Final,'b.')
        #f_ax1.plot(sim_time,x_Final,'b-',linewidth = 1,alpha = 0.25)
        f_ax1.ticklabel_format(style = 'sci',axis = 'y',scilimits = (0,0))
        f_ax1.ticklabel_format(style = 'sci',axis = 'x',scilimits = (0,0))
        if min(x_Final)>0:
            min_norm = 0
        else:
            min_norm = 1.25
        f_ax1.set_ylim(min(x_Final)*min_norm,max(x_Final)*1.25) 
        f_ax1.set_xlabel('$t$ $(days)$',fontsize = 20)
        f_ax1.set_ylabel('Flux',fontsize = 20) #fluxname
        f_ax1.set_title('Emmanoulopoulos Simulated Light Curve',fontsize = 25,pad = 12)
        
        binPDF = int(round(math.sqrt(len(x_Final)),0))
        f_ax2.hist(x_Final,bins = binPDF,density = True,histtype = u'step',alpha = 1,color = 'royalblue',linewidth = 3) # alpha = 0.7 histtype = u'step',
        f_ax2.ticklabel_format(style = 'sci',axis = 'y',scilimits = (0,0))
        f_ax2.set_title('Probability Density Function',fontsize = 20,pad = 12)
        f_ax2.set_xlabel('Flux',fontsize = 20)
        f_ax2.set_ylabel('PDF',fontsize = 20)
        
        f_txt='$\delta$='+str(round(p1,2))
        f_coeffp=np.poly1d(z)
        f_ax3.plot(simpleFreq,simplePSD,'k.')
        f_ax3.plot(simpleFreq,10**(f_coeffp(np.log10(simpleFreq))),'r-', linewidth=2)
        f_ax3.set_xlabel('$\\nu$ $(days^{-1})$',fontsize = 20)
        f_ax3.set_ylabel('PSD',fontsize = 20)
        f_ax3.set_title('Power Spectral Density ('+f_txt+')',fontsize = 20,pad = 12)
        f_ax3.set_xscale('log')
        f_ax3.set_yscale('log')
        plt.show()
    
    print('Number of iterations to Convergence: ', j_index)
    if plot_condition:
        plt.show()    
        
    return  sim_time, x_Final

if condition_guidance:
    example = example
    if not condition_spec:
        spec = None
    if example or condition_spec:
        sim(plot_condition = True, condition_guidance=True, example = example, spec = spec)   
    