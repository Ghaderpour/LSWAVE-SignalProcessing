"""
Antileakage Least-Squares Spectral Analysis (ALLSSA)
This module iteratively finds an optimal set of sinusoids that along with trend 
fit best to the time series. LSSA is called during each iteration.

The input time series does not have to be equally spaced, 
and the uncertainties in the time series values can also be considered.

Reference: 
Ghaderpour, E., Liao, W., Lamoureux, M.P., 2018. Antileakage least-squares 
       spectral analysis for seismic data regularization and random noise attenuation. 
       Geophysics, 83, V157--V170. doi: 10.1190/geo2017-0284.1

Author: Ebrahim Ghaderpour 
Email:  ebrahim.ghaderpour@ucalgary.ca
Copyright (c) 2021
"""
#==========================================================================================
import numpy as np
import time
from matplotlib import pyplot
from LSSA import LSSA
#==========================================================================================
def round_half_up(n, decimals = 0):
    multiplier = 10 ** decimals
    return np.floor(n*multiplier + 0.5) / multiplier
#------------------------------------------------------------------------------------------
def ALLSSA(t, f, P = 1, Omega = [], ind = [], level = 0.01, 
                                    trend = 'linear', slope = False, decimal = 1):
    """
    Summary: This function iteratively finds an optimal set of sinusoids that 
             along with trend fit best to the time series 
             Reference: 
             Ghaderpour, E., Liao, W., Lamoureux, M.P., 2018. 
                 Antileakage least-squares spectral analysis for seismic data 
                 regularization and random noise attenuation. 
                 Geophysics, 83, V157-V170. https://doi.org/10.1190/geo2017-0284.1
    
    Note: In the following descriptions, time and cyclic frequency can be
          distance and wavenumber, respectively
    Inputs:
    1) t           - Vector of order n. Times of observations or measurements 
    2) f           - Vector of order n. The time series values
    3) P           - Vector or Matrix of order n. The weights associated with
                     the observations. Default is 1 (a scale)
    4) Omega       - Vector. The initial cyclic frequencies 
    5) ind         - Vector. The indices of jumps in the trend component 
    6) level       - Numeric. The significance level, usually 0.01 or 0.05
                     Default is 0.01. 
    7) trend       - String. 'none','constant','linear','quadratic','cubic'
    8) slope       - Logical. Select true to estimate a unique slope for the
                     linear trend with multiple pieces
                     Default is false allowing the linear pieces to have
                     different slopes
    9) decimal     - Numeric. The number of decimals for frequency estimation
   
    Outputs:
    1) tr_coeff    - The estimated trend coefficients ordered as follows: 
                     For constant trend with multiple pieces: shifts
                     For linear trend with multiple pieces: intercepts,slopes 
                     For the quadratic trend (a+b*t+c*t^2): a, b, c
                     For cubic trend (a+b*t+c*t^2+d*t^3): a, b, c, d
    2) CS_coeff    - Estimated coefficients of cosine and sine functions of
                     the estimated cyclic frequencies (freq)
    3) res         - The residual series
    4) norm_res    - The squared weighted L2 norm of the residual series 
    5) cov         - Covariance matrix of the estimated coefficients
                     ordered as: intercepts, slopes, cosine and sine of
                     1st frequency, cosine and sine of 2nd frequency, etc.
    6) freq        - Estimated cyclic frequencies
    
    """
    #--------------------------------------------------------------------------------------
    Lt = len(t)    # Number of observations or measurements
    if not (Lt == len(f)): raise ValueError("t and f must have the same size")    
    #--------- Weights (Vector or Matrix)--------------------------------------------------
    SP = np.shape(P)
    CheckP = len(SP)
    if CheckP == 1 and (not len(P) == Lt):
        raise ValueError("The weights must be a vector of order n") 
    elif CheckP > 1 and not (SP[0] == Lt and SP[1] == Lt):
        raise ValueError("Weight matrix must be a square matrix of order n") 
    #--------- Set of Cyclic Frequencies --------------------------------------------------
    if len(Omega) == 0:
        M = (Lt-1)/(t[Lt-1]-t[0])        # A simple estimate for sampling rate
        if int(M/2) == M/2: Omega = [k for k in range(1, int(M/2))]
        else: Omega = [k for k in range(1, int(M/2)+1)]
    Omega = np.sort(Omega)
    #--------- Jump Indices ---------------------------------------------------------------
    if len(ind) > 0:
        if not (np.floor(ind)==ind).all():
            raise ValueError("The jump indices must be positive integers") 
        elif min(ind) < 2 or max(ind) > Lt-2:
            raise ValueError("The jump indices must be integers between 1 and n-1")
        ind = np.sort(ind)
    #--------- Significance Level ---------------------------------------------------------
    if level > 1 or level < 0:
        raise ValueError("The significance level must be between 0 and 1")
    #--------- decimal --------------------------------------------------------------------
    if (not int(decimal) == decimal) or decimal < 0 or decimal > 4:
        raise ValueError("The decimal must be a nonnegative integer less than five")    
    #======================================================================================
    check = True; freq1 = [];
    inc = pow(10,decimal)                  # For the cyclic frequency increment
    while check:
        Results = LSSA(t, f, P = P, Omega = Omega, ind = ind, level = level, 
                       trend = trend, slope = slope, freq = freq1)
        freq2 = freq1[:]
        Om_ind = np.argwhere(Results[0] == max(Results[0]))[0][0]
        Om = round_half_up(Omega[Om_ind])  # the selected integer cyclic frequency
        #----------------------------------------------------------------------------------
        # Find and remove the highest frequency from freq1 if it is close enough to Om  
        ZeroInd = np.nonzero([round_half_up(fr) for fr in freq1] == Om)[0]
        LZ = len(ZeroInd) 
        if LZ > 0: del freq1[ZeroInd[LZ-1]]            
        #------- Estimate LSS for a small neighbourhood around Om -------------------------
        if inc > 1: Om_range = [Om + k /inc - 0.5  for k in range(inc+1)]
        else: Om_range = [Om]; print(Om_range)
        Resultsnew = LSSA(t, f, P = P, Omega = Om_range, ind = ind, 
                          level = level, trend = trend, slope = slope, freq = freq1)
        Om_new_ind = np.argwhere(Resultsnew[0] == max(Resultsnew[0]))[0][0]
        #------- Check if the peak at Om_new is statistically significant -----------------
        if Resultsnew[0][Om_new_ind] > Resultsnew[1]:
            if inc > 1: Om_new = Om + Om_new_ind/inc - 0.5
            else: Om_new = Om 
            if Om_new in freq2: freq1 = freq2 # If peaks at the former peak (rare)
            else: freq1.append(round_half_up(Om_new, decimal)); freq1 = sorted(freq1)
        if freq1 == freq2:  check = False     # The termination condition     
    #-------------------------------------------------------------------------------------
    # coeff contains the estimated coefficients of trend and sinusoidal 
    # functions. The following command will separate them
    #-------------------------------------------------------------------------------------
    coeff = Results[5]; Lcoeff = len(coeff); Lfreq = len(freq1)
    if Lfreq > 0:
        CS_coeff = coeff[Lcoeff - 2*Lfreq : Lcoeff]
        tr_coeff = coeff[0 : Lcoeff - 2*Lfreq]
    else:
        CS_coeff = []
        tr_coeff = coeff
        
    return tr_coeff, CS_coeff, Results[3], Results[4], Results[6], freq1   
#==========================================================================================
if __name__ == '__main__':
    import argparse
    import sys
    # initialize the parser
    parser = argparse.ArgumentParser()
    # Add the parameters positional/optional 
    parser.add_argument('--datafile', help = "The path of data file. File should contain "    
                        "at least two columns of same order n: First column is the times, "
                        "second column is the time series values, third column should be "
                        "the weights OR third to n+2 columns should be the weight matrix "
                        "or the inverse of the covariance matrix associated with "
                        "the time series", type = str, default = 0)     
    parser.add_argument('--Lfreq', help = "Numeric. Lower cyclic frequency", 
                        type = float, default = 1) 
    parser.add_argument('--Ufreq', help = "Numeric. Upper cyclic frequency", 
                        type = float, default = 0)
    parser.add_argument('--Numfreq', help = "Numeric. The number of cyclic frequencies between "
                        "Lfreq and Ufreq", type = float, default = 0) 
    parser.add_argument('--ind', help = "Jump indices", nargs ='+', type = int, default = []) 
    parser.add_argument('--level', help = "Numeric. Significance level", 
                        type = float, default = 0.01) 
    parser.add_argument('--trend', help = "constant or linear or quadratic or cubic ",
                        type = str, default = "Linear") 
    parser.add_argument('--decimal', help = "Numeric. The number of decimals for frequency estimation: " 
                        "1, 2, 3, 4", type = int, default = 1)     
    parser.add_argument('--slope', help = "True estimates a unique slope, and "
                        "False allows estimation of multiple slopes",
                        type = bool, default = False)      
    parser.add_argument('--save', help = "Save the ALLSSA results as .xlsx", 
                        type = bool, default = False)    
    parser.add_argument('--display', help = "Display the amplitude spectrum ",
                        type = bool, default = False)     
    # parse the argument
    args = parser.parse_args()
    print(args)
    if args.datafile == 0:
        print("ERROR! missing input data file")
        print("Usage: {} datafile".format("JUSTmonitor.py"))
        sys.exit()
    
    filepath = args.datafile  
    # filepath = "D:\JUST\Sim_Geophysics.dat"                       
    # The same simulated example in https://doi.org/10.1190/geo2017-0284.1
    # True wavenumbers are 4.0744, 20.3718, 22.2818 by setting decimal = 4 in ALLSSA 
    
    
    if filepath[-3:] == 'dat' or filepath[-3:] == 'txt':
        with open(filepath) as TS:
            tfP = np.array([line.split() for line in TS], dtype = float)
    elif filepath[-3:] == 'csv': 
        tfP = np.loadtxt(open(filepath, "rb"), delimiter=",")
    else: 
        print("ERROR! only .csv and .dat and .txt are supported!")
        sys.exit()
    
    Dim = np.shape(tfP)   
    if Dim[1] == 1:
        print("The data set should at least have two columns!")
        sys.exit()  
    elif Dim[1] == 2 or Dim[1] == 3 or Dim[1] == Dim[0] + 2:
        t = tfP[:,0]; f = tfP[:,1]
    if Dim[1] == 3: P = tfP[:,2]
    elif Dim[1] == Dim[0] + 2: P = tfP[:,2:]
    else: P = 1
    
    if args.Ufreq > 0 and args.Numfreq > 0:
        Omega = np.linspace(args.Lfreq, args.Ufreq, args.Numfreq)
    else: Omega = []  
    
    t0  = t - t[0]   # For the sake of computational efficiency
    
    tic = time.clock()     
    Results = ALLSSA(t0, f, P = P, Omega = Omega, ind = args.ind, level = args.level, 
                   trend = args.trend, slope = args.slope, decimal = args.decimal)
    toc = time.clock()
    print("Computational Time: ", round(toc-tic,2), "s") 
    freq = Results[5]
    print("The ALLSSA estimated cyclic frequencies:", freq)
    
    if args.save:
        import os
        ExcelPath = filepath[:-4]
        if not os.path.exists(ExcelPath): os.makedirs(ExcelPath)         
        Path = ExcelPath + '\ALLSSAResults.xlsx' 
        LSSA(t0, f, P = P, ind = args.ind, trend = args.trend, 
                   slope = args.slope, freq = freq, saveto = Path)      
    
    if args.display:
        CS_coeff = Results[1]
        # ALLSSA Amplitude Spectrum 
        AmpSpectrum = [np.sqrt(CS_coeff[2*k]**2+CS_coeff[2*k+1]**2) for k in range(len(freq))]          
        fig, (ax1, ax2) = pyplot.subplots(2, 1, gridspec_kw={'hspace': 0.5})
        ax1.plot(t, f, 'blue', label = 'Time Series')
        ax1.plot(t, Results[2], 'red', label = 'Estimated Residual')
        ax1.set(ylabel='Value')
        ax1.set(xlabel='Time')
        ax1.grid(True) 
        ax1.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 3, ncol= 4, 
                   mode = "expand", borderaxespad = 0)  
        ax2.set_title('ALLSSA Amplitude Spectrum')
        ax2.stem(freq, AmpSpectrum,'blue')
        ax2.set(ylabel='Amplitude')
        ax2.set(xlabel='Cyclic frequency')
        ax2.grid(True)    
        pyplot.show()        
    
        
 
    '''
    # ----- Simulated Example -------------------------------------------------------------
    t = [k/23.0 for  k in range(70)]
    for kk in range(20):
        indt = int((70-kk)*np.random.rand(1)[0])
        t.remove(t[indt])
    Lt = len(t)
    t = np.array(t) 
    trend = 0.5 - 0.01*np.array(t)
    season = 0.1*np.cos(2*np.pi*1*np.array(t)) + 0.05*np.sin(2*np.pi*2.3*np.array(t))
    noise1 = 0.01*np.random.normal(0,1,Lt)
    noise2 = 0.03*np.random.normal(0,1,Lt)
    f = trend + season + noise1 - abs(noise2)
    P = 1/(0.01+(noise2)**2)
    Omega = [Om/5.0 for Om in range(4,20)]
    '''
    #--------------------------------------------------------------------------------------
 
