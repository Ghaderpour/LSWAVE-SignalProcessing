"""
Least-Squares Wavelet Analysis (LSWA)

This module computes a time-frequency spectrogram for a given time series based on 
the least-squares fit of sinusoids of different cyclic frequencies to the 
time series segments. The constitunets of known forms, such as trends and sinusoids 
of known frequencies are fitted first to each segment. 
Then, the spectra of the residual series segments will be estimated considering 
the effects of the removed known constituents that will ultimately obtain a spectrogram.

The input time series does not have to be equally spaced, 
and the uncertainties in the time series values can also be considered.

Reference:
Ghaderpour, E., Pagiatakis, S.D., 2017. Least-squares wavelet analysis of 
       unequally spaced and non-stationary time series and its applications. 
       Mathematical Geosciences, 49, 819-844.doi: 10.1007/s11004-017-9691-0

Author: Ebrahim Ghaderpour 
Email:  ebrahim.ghaderpour@ucalgary.ca
Copyright (c) 2021
"""
#==========================================================================================
import numpy as np
import time
from matplotlib import pyplot
from LSSA import LSSA
#------------------------------------------------------------------------------------------
def LSWA(t, f, P = 1, tt = [], rate = None, Omega = [], ind = [], level = 0.01, 
         trend = 'linear', slope = False, freq = [], L1 = None, L0 = None, morlet = 0.0125):
    """
    Summary: This function computes a frequency spectrum for a given time
             series. The constitunets of known forms, such as trends and 
             sinusoids of known frequencies are fitted first 
             Reference:
             Ghaderpour, E., Pagiatakis, S.D., 2017. Least-squares wavelet analysis 
                  of unequally spaced and non-stationary time series and its applications. 
                  Mathematical Geosciences 49, 819-844 doi: 10.1007/s11004-017-9691-0
             
    Inputs:
    1) t           - Vector of order n. Times of observations or measurements 
    2) f           - Vector of order n. The time series values
    3) P           - Vector or Matrix of order n. The weights associated with
                     the observations. Default is 1 (a scale)
    4) tt          - Vector. Equally spaced times for regularization purposes
    5) rate        - Numeric. The (average) sampling rate (M)
    6) Omega       - Vector. The initial cyclic frequencies 
    7) ind         - Vector. The indices of jumps in the trend component 
    8) level       - Numeric. The significance level (alpha), usually 0.01 or 0.05
                     Default is 0.01. 
    9) trend       - String. 'none','constant','linear','quadratic','cubic'
    10)slope       - Logical. Select true to estimate a unique slope for the
                     linear trend with multiple pieces
                     Default is false allowing the linear pieces to have
                     different slopes
    11)freq        - Vector. The known cyclic frequencies of the sinusoids
                     to be estimated/removed from the time series while their
                     removal effect is considered in the spectrum of residual
    12) L1         - Numeric. The number of cycles of sinusoids being fitted 
                     to the segments of the time series. 
                     Recommeded: L1 = 2 when morlet = 0
    13) L0         - Numeric. The number of additional observations to increase 
                     the translating window size (e.g., L0 = 5, L0 = 20)
    14) morlet     - Numeric. The Morlet coefficient. 
                     Selecting 0 will fit the sinusoids to the segments 
                     without any Gaussian weights
                     Default is 0.0125 with L1 = 6 cycles and L0 = 0. Thus,
                     the sinusoids will be adapted to the Morlet wavelet
                     in the least-squares sense via Gaussian values. 
                     Note: The Gaussian values may not be attenuated to zero 
                           at both ends of the windows for certain L1 and L0
   
    Outputs:
    1) spectrogram        - The normalized least-squares wavelet spectrogram
    2) stoch_surf         - The stochastic surface at (1-alpha) confidence 
                            level obtained for the normalized spectrogram
    3) amp_spectrogram    - The amplitude spectrogram 
    """
    #--------------------------------------------------------------------------------------
    Lt = len(t)    # Number of observations or measurements
    if not (Lt == len(f)):
        raise ValueError("t and f must have the same size")    
    #--------- Weights (Vector or Matrix)--------------------------------------------------
    SP = np.shape(P)
    CheckP = len(SP)
    if CheckP == 1 and (not len(P) == Lt):
        raise ValueError("The weights must be a vector of order n") 
    elif CheckP > 1 and not (SP[0] == Lt and SP[1] == Lt):
        raise ValueError("Weight matrix must be a square matrix of order n") 
    #--------- Equally Spaced Times -------------------------------------------------------
    Ltt = len(tt) 
    if Ltt > 0:
        if (t[Lt-1] < tt[Ltt-1]) or (tt[0] < t[0]):
            raise ValueError("The range of time must be within the range of t") 
        check = True
    else:
        tt = t[:]
        check = False
    Ltt = len(tt) 
    #--------- Sampling Rate ---------------------------------------------------------------
    # M is a simple estimate for sampling rate. Note: this estimate may not be 
    # valid for certain segments in an unequally spaced time series that may 
    # cause aliasing when selecting the range of frequencies based on it
    if rate == None: M = (Ltt-1)/(tt[Ltt-1]-tt[0])    # A simple estimate for sampling rate
    else: M = rate
    #--------- Set of Cyclic Frequencies --------------------------------------------------
    if len(Omega) == 0:
        if int(M/2) == M/2: Omega = [k for k in range(1, int(M/2))]
        else: Omega = [k for k in range(1, int(M/2)+1)]
    else: Omega = np.sort(Omega)    
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
    #--------- Window Dilation (L1) -------------------------------------------------------
    if L1 == None and morlet: L1 = 6
    elif L1 == None and not morlet: L1 = 2
    #--------- Additional Data Points (L0) ------------------------------------------------
    if L0 == None and morlet: L0 = 0
    elif L0 == None and not morlet: L0 = 20   
    #======================================================================================
    LOm = len(Omega)
    spectrogram = []; stoch_surf = []; amp_spectrogram = []
    start_time = tt[0]
    tt = tt - start_time;  # For the sake of computational efficiency 
    t  = t  - start_time;  # For the sake of computational efficiency   
    #--------------------------------------------------------------------------------------
    # Remove the trend with multiple pieces first using LSSA if ind is given 
    if len(ind) > 0:
        Results = LSSA(t, f, P = P, ind = ind, level = level, trend = trend, slope = slope)
        f = Results[3]
    #--------------------------------------------------------------------------------------
    for k in range(LOm):                         #loop over the frequencies
        L_k = L0 + int((L1*M)/Omega[k])    
        if not L_k % 2: L_k += 1
        spectrogram_kj = []; stoch_surf_kj  = []; amp_spectrogram_kj = []      
        for j in range(Ltt):                     #loop over the times
            if j >= (L_k-1)/2 and j < Ltt-(L_k-1)/2:         #non-marginal windows
                lb = int(j - (L_k-1)/2); ub = int(j + (L_k+1)/2)
            elif j < (L_k-1)/2 and j < Ltt-(L_k-1)/2:        #left marginal windows
                lb = 0; ub = int(j + (L_k+1)/2)
            elif j >= Ltt-(L_k-1)/2 and j >= (L_k-1)/2:      #right marginal windows
                lb = int(j - (L_k-1)/2); ub = Ltt
            else:                                            #for the entire series
                lb = 0; ub = Ltt
            if check:
                lb = np.where(t >= tt[lb])[0][0]        
                ub = np.where(t >= tt[ub-1])[0][0] + 1   
            ty = t[lb:ub]                             # time segment
            y  = f[lb:ub]                             # time series segment
            # The Gaussian weights based on the morlet value to define Py
            if morlet != 0: MorlW = np.exp(-morlet*((2*np.pi*Omega[k])**2)*(ty-tt[j])**2)
            else: MorlW = []  
            LM = len(MorlW)
            if CheckP > 1 and LM: Py = np.diag(P[lb:ub, lb:ub])*MorlW  # extract submatrix  
            elif CheckP > 1 and not LM: Py = P[lb:ub, lb:ub]           # extract submatrix 
            elif CheckP == 1 and LM: Py = P[lb:ub]*MorlW               # extract submatrix
            elif CheckP == 1 and not LM: Py = P[lb:ub]                 # extract submatrix
            elif CheckP == 0 and LM: Py = MorlW
            else: Py = 1          # If the time series is equally weighted and morlet == 0
            #----- Implement LSSA to segments ty, y, Py ------------------------------------
            Results = LSSA(ty, y, P = Py, Omega = [Omega[k]], level = level, 
                           trend = trend, slope = slope, freq = freq)   
            spectrogram_kj.append(Results[0][0])    # spectral peak  for (Omega(k), t(j))
            stoch_surf_kj.append(Results[1])        # critical value for (Omega(k), t(j))
            CS_coeff = Results[2][0]
            amp_spectrogram_kj.append(np.sqrt(CS_coeff[0]**2 + CS_coeff[1]**2))
        spectrogram.append(spectrogram_kj)
        stoch_surf.append(stoch_surf_kj)
        amp_spectrogram.append(amp_spectrogram_kj)
    return np.array(spectrogram), np.array(stoch_surf), np.array(amp_spectrogram)
#==========================================================================================
def PlotSpectrogram (t, f, tt, Omega, spectrogram, stoch_surf, amp_spectrogram, display):
    '''
    This function displays and saves the least-squares spectrograms.
    Inputs:
    1) t                - Vector of order n. Times of observations or measurements 
    2) f                - Vector of order n. The time series values
    3) tt               - Vector. Equally spaced times for regularization purposes
    4) Omega            - Vector. The initial cyclic frequencies 
    5) spectrogram      - The normalized least-squares wavelet spectrogram
    6) stoch_surf       - The stochastic surface at (1-alpha) confidence 
                          level obtained for the normalized spectrogram
    7) amp_spectrogram  - The amplitude spectrogram 
    8) display          - Numeric. Diplays the spectrogram: 
                          1 displays the normalized spectrogram
                          2 displays the statistically significant spectrogram
                          3 displays the amplitude spectrogram
    '''

    font = {'family' : 'Arial','size' : 12}
    pyplot.rc('font', **font)     
    fig, (ax1, ax2) = pyplot.subplots(2, 1, sharex=True, figsize=(7,4),
                                     gridspec_kw={'height_ratios': [3.5, 5.5]}, 
                                     constrained_layout=True)  
    ax1.set_title('Time Series Decomposition')
    ax1.plot(t,f, '--gD', linewidth=1, markersize=3)
    ax1.set(ylabel='Value')
    ax1.grid(True)  
    if args.display == 1:
        im = ax2.pcolor(tt, Omega, spectrogram*100, cmap='jet') 
        fig.colorbar(im, label = 'Percentage variance', orientation='vertical') 
    elif args.display == 2:
        spectrogram = np.where(spectrogram < stoch_surf, np.nan, spectrogram)
        im = ax2.pcolor(tt, Omega, spectrogram*100, cmap='jet') 
        fig.colorbar(im, label = 'Percentage variance', orientation='vertical') 
    elif args.display == 3:
        im = ax2.pcolor(tt, Omega, amp_spectrogram, cmap='jet') 
        fig.colorbar(im, label = 'Amplitude', orientation='vertical')     
    ax2.set(ylabel = 'Cyclic frequency')
    ax2.set(xlabel = 'Time')   
    ax2.axis([t.min(), t.max(), np.min(Omega), np.max(Omega)])      
    
    fig.savefig("LSWAdecomp.png", dpi = 300, bbox_inches = 'tight')
    pyplot.show()
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
    parser.add_argument('--Ltime', help = "Numeric. The start time for regularization", 
                        type = float, default = 0) 
    parser.add_argument('--Utime', help = "Numeric. The end time for regularization",
                        type = float, default = 0) 
    parser.add_argument('--Numtime', help = "Numeric. The number of equally spaced times "
                        "between Ltime and Utime", type = float, default = 0)     
    parser.add_argument('--rate', help = "Numeric. M: An estimate for the sampling rate", 
                        type = int, default = None)
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
    parser.add_argument('--slope', help = "True estimates a unique slope, and "
                        "False allows estimation of multiple slopes",
                        type = bool, default = False)      
    parser.add_argument('--freq', help = "Known cyclic frequencies", 
                        nargs ='+', type = float, default = []) 
    parser.add_argument('--L1', help = "Numeric. Window dilation parameter", 
                        type = int, default = 6) 
    parser.add_argument('--L0', help = "Numeric. Additional number of data points to enlarge " 
                        "the window size", type = int, default = 0)                          
    parser.add_argument('--morlet', help = "Numeric. The morlet coefficient. "
                        "Zero mean no morlet wavelet applies", type = float, default = 0.0125)
    parser.add_argument('--save', help = "Save the spectrograms", 
                        type = bool, default = False)    
    parser.add_argument('--display', help = "Display the spectrogram: "
                        "1 displays the normalized spectrogram "
                        "2 displays the statistically significant normalized spectrogram "
                        "3 displays the amplitude spectrogram", 
                        type = int, default = 1)     
    # parse the argument
    args = parser.parse_args()
    print(args)
    if args.datafile == 0:
        print("ERROR! missing input data file")
        print("Usage: {} datafile".format("JUSTmonitor.py"))
        sys.exit()
    
    filepath = args.datafile  
    #filepath = "D:\JUST\Sim_JumpInd_48_104.dat"           # Jump indices: [47,103]
    #filepath = "D:\JUST\EVIA_CA.dat"                      # Jump indices: [13, 111, 144]
    
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
     
    if args.Numtime > 0:
        if args.Ltime == 0: Ltim = t[0]
        else: Ltim = args.Ltime
        if args.Utime == 0: Utim = t[len(t)-1]
        else: Utim = args.Utime
        tt = np.linspace(Ltim, Utim,int(args.Numtime))
    else: tt = []
    
    if args.Ufreq > 0 and args.Numfreq > 0:
        Omega = np.linspace(args.Lfreq, args.Ufreq,int(args.Numfreq))
    else: Omega = []  
    
    t0  = t - t[0]   # For the sake of computational efficiency
    if len(tt) > 0: tt0 = tt - t[0]   # For the sake of computational efficiency 
    else: tt0 = []
    
    tic = time.time()     
    Results = LSWA(t0, f, P = P, tt = tt0, rate = args.rate, Omega = Omega, ind = args.ind, 
                   level = args.level, trend = args.trend, slope = args.slope, freq = args.freq, 
                   L1 = args.L1, L0 = args.L0, morlet = args.morlet)
    toc = time.time()
    print("Computational Time: ", round(toc-tic,2), "s") 
    # Select larger window size parameters (e.g., L0 = 5) if you got the following error: 
    # 'The size of the time series (segment) is too small!'
    #  Or reduce the upper bound for Omega      

    spectrogram     = Results[0]
    stoch_surf      = Results[1]
    amp_spectrogram = Results[2]
    if args.save:
        import os
        CSVPath = filepath[:-4]
        if not os.path.exists(CSVPath): os.makedirs(CSVPath)          
        np.savetxt(CSVPath + '\Spectrogram.csv', spectrogram, delimiter=',')
        np.savetxt(CSVPath + '\StochasticSurface.csv', stoch_surf, delimiter=',') 
        np.savetxt(CSVPath + '\AmplitudeSpectrogram.csv', amp_spectrogram, delimiter=',') 
    
    if len(tt) == 0: tt = t[:]
    if len(Omega) == 0: Omega = [k for k in range(1,np.shape(spectrogram)[0]+1)]        
    if args.display:
        PlotSpectrogram (t, f, tt, Omega, spectrogram, stoch_surf, 
                         amp_spectrogram, args.display) 
   

    
