def lc_sim(nn, delt, mean_lc, model, params):
	
	'''
	Python function to simulate a light curve given a model power spectrum
	Based on paper Timmer, J., & Koenig, M. 1995, A&A, 300, 707
	nn = number of points in simulated light curve
	delt = time sampling of simulated light curve in seconds
	mean_lc = mean value of simulated light curve
	model = "unbroken", "sharp", or "slow" to designate what type of power spectrum the simulated
	         light curve will be based on
	params = 4 (for 'unbroken' model) or 5 (for the other two) element array containing the
	         parameters of the power spectrum model
	         For 'unbroken', S = N*(nu/nu_0)**(-beta) + noise:
	         	params[0] = N, normalization of power spectrum
	         	params[1] = nu_0, the frequency at which the power = the normalization
	         	params[2] = beta, the slope of the power spectrum
	         	params[3] = white noise level to add in
	         For 'sharp', S = N*(nu/nu_c)^(-gamma) + noise      nu < nu_c
	                      S = N*(nu/nu_c)^(-beta)  + noise      nu > nu_c
	            params[0] = N, normalization of power spectrum
	            params[1] = nu_c, the pivot frequency
	            params[2] = gamma, the low frequency slope
	            params[3] = beta, the high frequency slope
	            params[4] = white noise level
	        For 'slow', S = N*nu^(alpha_lo)/(1+(nu/nu_k))^(alpha_hi - alpha_lo)
	        	params[0] = N, normalization of power spectrum
	        	params[1] = nu_k, the 'knee' frequency where the spectrum rolls over
	        	params[2] = alpha_lo, the low frequency slope
	        	params[3] = alpha_high, the high frequency slope
	        	params[4] = white noise level
	'''
	
	import numpy as np
	
	#Seed the random number generator
	#np.random.seed(6543289)
	
	#Fourier frequencies
	fi = np.arange(1, nn/2+1, dtype='float64')/nn/delt
	#print fi
	
	'''
	Power at each frequency depending on power spectrum model
	Unbroken model: S = N*(nu/nu0)^(-beta)
	Sharply broken model: S = N*(nu/nu_c)^(-gamma)       nu < nu_c
	                      S = N*(nu/nu_c)^(-beta)        nu > nu_c
	Slowly bending "knee" model: S = N*nu^(alpha_low)/(1+(nu/nu_k))^(alpha_hi - alpha_lo)
	'''
	if model == 'unbroken':
	
		a = params[0]
		nu_0 = params[1]
		beta = params[2]
		noise = params[3]
		
		s = a*(fi/nu_0)**(-beta) + noise
		
	elif model == 'sharp':
		
		a = params[0]
		nu_c = params[1]
		gamma = params[2]
		beta = params[3]
		noise = params[4]
		
		s = np.zeros(len(fi))
		s[fi <= nu_c] = a*(fi[fi <= nu_c]/nu_c)**(-gamma) + noise
		s[fi > nu_c] = a*(fi[fi > nu_c]/nu_c)**(-beta) + noise
		
	elif model == 'slow':
	
		a = params[0]
		nu_knee = params[1]
		alpha_lo = params[2]
		alpha_hi = params[3]
		noise = params[4]
		
		s = (a*fi**(alpha_lo))/(1+(fi/nu_knee))**(alpha_hi - alpha_lo) + noise
	
	#Generate two sets of normally distributed random numbers 	
	aa = np.random.randn(len(fi))
	bb = np.random.randn(len(fi))
	
	#Fourier transform of light curve  = SQRT(S/2) * (A + B*i)
	flc = np.sqrt(.5*s)*(aa + bb*1.j)
	#plt.plot(np.real(flc))
	
	if np.mod(nn, 2) == 0:
		flc[-1] = np.sqrt(.5*s[-1])*1
	
	del aa, bb, s
	
	#Put the mean of the light curve at frequency = 0
	flc = np.hstack([mean_lc, flc])
	
	#Take the inverse fourier transform to generate synthetic light curve
	lc = np.fft.irfft(flc, n=nn)
	
	return lc