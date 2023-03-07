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
#from statistics import mode

params  =  {'xtick.labelsize':18,'ytick.labelsize':18}
pylab.rcParams.update(params)

def graphs(Flux,Time,Title, sim_step = False, Flux_error = None, ULtime = None, ULflux = None, label_flux = None):
    SamplingFreq = 1/(Time[2] - Time[1])
    Freq, PSD = scipy.signal.periodogram(Flux, SamplingFreq, scaling = 'density')
    mask = Freq > 0
    Freq = Freq[mask]
    PSD = PSD[mask]
    Z = np.polyfit(np.log10(Freq),np.log10(PSD),1)
    p1 = Z[0]   # angular coefficient
    p0 = Z[1]   # intercept
    
    Fig = plt.figure(figsize = (16,9),tight_layout = True)
    ax1 = plt.subplot(211)  # LC
    ax2 = plt.subplot(223)  # PDF
    ax3 = plt.subplot(224)  # PSD

    ax1.plot(Time,Flux,'b.')
    if sim_step:
        ax1.plot(Time,Flux,'b-',linewidth = 1,alpha = 0.25)
    if isinstance(Flux_error,np.ndarray):
        ax1.errorbar(Time,Flux,yerr = Flux_error,xerr = None,fmt = 'b.',ecolor = 'b')
    if isinstance(ULflux,np.ndarray) and isinstance(ULtime,np.ndarray):
        ax1.plot(ULtime,ULflux,'rv')
    ax1.ticklabel_format(style = 'sci',axis = 'y',scilimits = (0,0))
    ax1.ticklabel_format(style = 'sci',axis = 'x',scilimits = (0,0))
    ax1.set_title(Title,fontsize = 25,pad = 12)
    
    if min(Flux) > 0:
        min_norm = 0
    else:
        min_norm = 1.25
    ax1.set_ylim(min(Flux)*min_norm,max(Flux)*1.25)
    
    if Time[0] > 1000:
        ax1.set_xlabel('$MJD$',fontsize = 20)
    else:
        ax1.set_xlabel('$t$ $(days)$',fontsize = 20)
        
    if label_flux == None:
        ax1.set_ylabel('Flux',fontsize = 20)
    else:
        ax1.set_ylabel(label_flux,fontsize = 20)
        

    binPDF = int(round(math.sqrt(len(Flux)),0))
    ax2.hist(Flux,bins = binPDF,density = True,histtype = u'step',alpha = 1,color = 'royalblue',linewidth = 3)
    ax2.ticklabel_format(style = 'sci',axis = 'y',scilimits = (0,0))
    ax2.set_title('Probability Density Function',fontsize = 20,pad = 12)
    if label_flux == None:
        ax2.set_xlabel('Flux',fontsize = 20)
    else:
        ax2.set_xlabel(label_flux,fontsize = 20)
    ax2.set_ylabel('PDF',fontsize = 20)

    Text='$\delta$='+str(round(p1,2))
    Coeff=np.poly1d(Z)
    ax3.plot(Freq,PSD,'k.')
    ax3.plot(Freq,10**(Coeff(np.log10(Freq))),'r-', linewidth=2)
    ax3.set_xlabel('$\\nu$ $(days^{-1})$',fontsize = 20)
    ax3.set_ylabel('PSD',fontsize = 20)
    ax3.set_title('Power Spectral Density ('+Text+')',fontsize = 20,pad = 12)
    ax3.set_xscale('log')
    ax3.set_yscale('log')

def sim(plot_condition = False, final_plot = False, time = None, flux = None, flux_error = None, temporalbin = None, PSDmodel = None, PSDparams = None, source_name_title = None, label_flux = None):
    
    # plot_condition     (boolean)   to show plots at each step of the algorithm.
    # final_plot         (boolean)   to show just the plot of Emmanoulopoulos simulated Light Curve.
    # time             (numpy array) time of observations, it is considered in Modified Julian Date (MJD).
    # flux             (numpy array) flux of the source.
    # flux_error       (numpy array) flux uncertainty of the source.
    # temporalbin          (int)     sampling time, required with your own data.
    # PSDmodel            (string)   Power Spectral Density (PSD) model for Timmer Konig simulation, if None it is considered a simple unbroken law. Possibilities: unbroken, sharp, slow.
    # PSDparams      (list of float) Parameters for PSD, the numbers and the type depends on the model (4 unbroken, 5 sharp, 5 slow). See below for more details.
    # source_name_title   (string)   Title name of the source for plots if you use your own data.
    # label_flux   (string)   label for the flux
        
    #  return  sim_time, x_Final
        
    #######################
        
    ULflux = None
    ULtime = None
    
    if label_flux == None:
        label_flux = 'Flux'
    if source_name_title == None:
        nometitolo = ' '
        
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
        if not source_name_title == None:
            title = 'Light Curve of '+source_name_title
        else:
            title = 'Light Curve of '+nometitolo
        graphs(Flux = flux,Time = time,Title = title, sim_step = False, Flux_error = flux_error, ULtime = ULtime, ULflux = ULflux, label_flux = label_flux)
    
        
    # Set n, the number of timesteps in the simulated light curve, the length of each timestep in seconds.
    delta_time = []
    for i in range(0,(len(time)-1)):
        delta_time.append(time[i+1]-time[i])
        
    sampling_days = temporalbin #mode(delta_time)   # 30, 7, 3 days if it comes from the Fermi Light Curve Repository (LCR)
    sampling_seconds  = sampling_days*24*3600
    
    n = (time[-1]-time[0])/sampling_days
    
    # Time array for simulations
    sim_time = np.arange(0,sampling_days*n,sampling_days)
    
    # Set the mean for the simulated light curve 
    mean_lc = np.mean(flux)
    
    # Due to the steps of the algorithm the slope is the only parameter that matters
    if PSDmodel == None:
        PSDmodel = 'unbroken'
        if PSDparams == None:
            PSDparams = [1,1,-(p1),0.0]
    elif PSDmodel == 'unbroken':
        if PSDparams == None:
            PSDparams = [1,1,-(p1),0.0]
    elif PSDmodel == 'sharp':
        if PSDparams == None:
            PSDparams = [1,1,-(p1*0.75),-(p1),0.0]
    elif PSDmodel == 'slow':
        if PSDparams == None:
            PSDparams = [1,1,-(p1*0.75),-(p1),0.0]
            
    
    ########################################################
    # '------------ First step of the algorithm ------------'
    # 'Timmer Koenig LC simulation, with the same PSD model and parameters of the source')
    # Simulate the light curve as in Timmer Koenig 1995
    x_TK  =  simulate_lc.lc_sim(int(n), sampling_days, mean_lc, PSDmodel, PSDparams)

    # TK simulation plots: LC, PDF, PSD
    if plot_condition:
        title='Timmer Koenig Simulated Light Curve'
        graphs(Flux = x_TK,Time = sim_time,Title = title, sim_step = True)
        
        
    # 'From the TK simulation we calculate the Discrete Fourier Function (DFT) and from it we extract the amplitude: A_TK'
    
    
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
    
    # we look for the convergence in the next steps
    converge = False    
    while not converge:
        j_index += 1
        
        ########################################################
        
        if j_index == 1:
            
            # '------------ Second step of the algorithm ------------'
            # 'From the PDF of the source we reproduce white noise data: x_WN.'
            # 'Estimate the DFT and hence the phase: phi_loop (it will change at each loop)'
        
            zero_one = np.random.uniform(0, 1, len(sim_time))
            x_WN = inverseCDF(zero_one)
            DFTxwn = np.fft.fft(x_WN)
            phi_loop = np.ones(len(DFTxwn.real))
            A_loop = np.ones(len(DFTxwn.real))
            for i in range(0,len(DFTxwn.real)):
                A_loop[i], phi_loop[i] = cmath.polar(DFTxwn[i])
            
        else: # after the first step we will take the final simulation of the algorithm (5th step)
            x_WN = x_loop
            DFTloop = np.fft.fft(x_loop)
            phi_loop = np.ones(len(DFTloop.real))
            A_loop = np.ones(len(DFTloop.real))
            for i in range(0,len(DFTloop.real)):
                A_loop[i], phi_loop[i] = cmath.polar(DFTloop[i])
            
        
        if j_index == 1 and plot_condition:
            title='First White Noise Simulated Light Curve'
            graphs(Flux = x_WN,Time = sim_time,Title = title, sim_step = True)
                
        
        
        ########################################################
        
        #if j_index == 1:
            # '------------ Third step of the algorithm ------------'
            # 'We take A_TK and phi_loop, from these two we obtain an Adjust DFT, then we estimate x_Adj through the Inverse DFT'
        
        DFT_Adj = np.zeros(len(A_TK))*1j
        for i in range(0,len(A_TK)):
            DFT_Adj[i] = cmath.rect(A_TK[i],phi_loop[i])
        
        x_Adj = np.fft.ifft(DFT_Adj)
        x_Adj=x_Adj.real
        
        if j_index == 1 and plot_condition:
            title='First Adjusted ($A_{TK}$, $\phi_{loop}$) Simulated Light Curve'
            graphs(Flux = x_Adj,Time = sim_time,Title = title, sim_step = True)
        
    
        ########################################################

       
        #if j_index == 1:
            # '------------ Fourth step of the algorithm ------------')
            # 'We create a new time series with the values of x_WN, at first, ordering them as the order of x_Adj. Obtaining x_loop.'
            
            # '(This means that the highest value of x_Adj is replaced with the highest value of x_WN, the second highest value of x_Adj is replaced with the second highest value of x_WN, and so on, so we keep the same time distribution)'
        
        
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
            title='First Loop Simulated Light Curve'
            graphs(Flux = x_loop,Time = sim_time,Title = title, sim_step = True)

        ########################################################
    
        if j_index == 1:
            p1_pre = 0
            
            
                # '------------ Fifth step of the algorithm ------------'
                
                # 'We repeat all the processes from the second step until the two light curves: x_Adj and x_loop, converge.'
                
                # '(second step) we take the phase phi_loop from the DFT of x_loop, (third step) we generate a new x_Adj from phi_loop and A_TK, (fourth step) we sort and replace and at the end we check for the convergence.'
                # 'when the two series converge it means we have obtained the Emmanoulopoulos simulated light curve x_Final, with the same PSD and PDF of the original source.'
                
                
            
        # converging condition based on the PSD
        if p1_pre == p1_loop: 
            index_converge += 1
            
        if index_converge >= 30:
            converge = True
        
        # "extreme" condition, the convergence should be already reached
        if j_index>100: 
            if round(p1_pre,8) == round(p1_loop,8):
                converge = True
             
        p1_pre = p1_loop
            
        if converge:
            x_Final = x_loop
            
    
    ############################################################
    
    
    # '------------ Last (sixth) step of the algorithm ------------'
    # 'The last thing to do is the addition of Poisson noise, in order to consider the light curve as a product of a counting detector-process.'
    # 'Each point of x_Final is replaced by: Pois[mu = x_Final * dt * effective_area] / (dt * effective_area)')
    # '(effective_area =  9000 cm**2 was used)' 
    
    
    
    # The last step should consist in the add of the poissonian noise. Commented below how to add it.
    # effective_area=9000
    # for i in range(0,len(x_Final)):
        # x_Final[i] = np.random.poisson(x_Final[i]*sampling_seconds*effective_area)/(sampling_seconds*effective_area)
    
    if plot_condition:
        final_plot = True
    if final_plot:
        title='Emmanoulopoulos Simulated Light Curve'
        graphs(Flux = x_Final,Time = sim_time,Title = title)
        
    if plot_condition or final_plot:
        plt.show()    
        
    return  sim_time, x_Final
    