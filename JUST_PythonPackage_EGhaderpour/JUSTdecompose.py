"""
JUST Decomposition into Trend, Seasonal, and Remainder Components

This module uses the jump indices estimated by JUST to decompose the time series into 
trend, seasonal, and remainder components via ALLSSA or OLS    

The input time series does not have to be equally spaced, 
and the uncertainties in the time series values can also be considered.

Reference: 
Ghaderpour, E., Vujadinovic, T. Change Detection within Remotely Sensed 
      Satellite Image Time Series via Spectral Analysis, 
      Remote Sensing, 2020, 12, 4001; doi:10.3390/rs12234001

Author: Ebrahim Ghaderpour 
Email:  ebrahim.ghaderpour@ucalgary.ca
Copyright (c) 2021
"""
#==========================================================================================
import numpy as np
import time
from matplotlib import pyplot
from LSSA import LSSA
from ALLSSA import ALLSSA
#------------------------------------------------------------------------------------------
def TrendComponent(t, tt, ind, tr_coeff):
    """
    This function gets the estimated intercepts and slopes of linear trend
    with multiple pieces and jump indices to create the trend component.

    Inputs:
    1) t            - Times of the observations or measurements as a vector
    2) tt           - Another time vector for regularization purposes
    3) ind          - The indices of jumps in the trend component 
    4) tr_coeff     - Estimated intercepts and then slopes of linear pieces

    Output:
    1) trend        - The trend series
    """
    #--------------------------------------------------------------------------------------
    Ltt = len(tt); Lind = len(ind); JumpInd = []
    for k in range(Lind): JumpInd.append(np.where(tt >= t[ind[k]])[0][0] + 1)     
    # Create a design matrix D for the linear trend with multiple pieces
    D = np.zeros([Ltt,2*Lind+2]); q = 0
    if Lind == 0:  
        D[:,0] = np.ones([1,Ltt]) # constituent of known form for the intercept
        D[:,1] = tt               # constituent of known form for slope
    else:
        # constituents of known forms for intercepts
        D[0:JumpInd[0], q] = np.ones([1,JumpInd[0]])
        D[0:JumpInd[0], q+Lind+1] = tt[0:JumpInd[0]]
        q += 1        
        for j in range(Lind):
            if j < Lind-1: 
                D[JumpInd[j]:JumpInd[j+1], q] = np.ones([1,JumpInd[j+1]-JumpInd[j]])                
                D[JumpInd[j]:JumpInd[j+1], q+Lind+1] = tt[JumpInd[j]:JumpInd[j+1]]
            else: 
                D[JumpInd[j]:Ltt, q] = np.ones([1,Ltt-JumpInd[j]])
                D[JumpInd[j]:Ltt, q+Lind+1] = tt[JumpInd[j]:Ltt]
            q += 1          
    return D.dot(tr_coeff) 
#==========================================================================================
def SeasonalComponent(t, freq, CS_coeff):
    """
    This function gets the estimated coefficients for sinusoids and their 
    frequencies to create the seasonal time series.
    
    Inputs:
    1) t            - Times of the observations or measurements as a vector
    2) freq         - Cyclic frequencies (enter as a vector) 
    3) CSCoeff      - A vector that contains the coefficients of cosine and
                      sine functions of each cyclic frequency consecutively
    Output:
    1) seasonal     - The seasonal series

    Note: This function can also be used for regularization, i.e., one may 
          input any equally spaced vector t to generate equally spaced series
    """
    #-------------------------------------------------------------------------------------
    k = 0; freq_ind = 0; seasonal = np.zeros(len(t))
    while freq_ind < len(freq):
        seasonal += CS_coeff[k]   * np.cos(2*np.pi*freq[freq_ind]*np.array(t))
        seasonal += CS_coeff[k+1] * np.sin(2*np.pi*freq[freq_ind]*np.array(t))
        freq_ind += 1; k += 2    
    return seasonal    
#==========================================================================================
def JUSTdecompose(t, f, P = 1, tt = [], size = None, step = None, season = 'ALLSSA', 
                  Omega = [], ind = [], level = 0.01, display = False):
    """
    Summary: This function uses the jump indices estimated by JUST to
             decompose the time series into trend, seasonal, and remainder 
             components via ALLSSA or OLS
    
             Reference: 
             Ghaderpour, E., Vujadinovic, T. Change Detection within Remotely 
                    Sensed Satellite Image Time Series via Spectral Analysis, 
                    Remote Sensing, 2020, 12, 4001; doi:10.3390/rs12234001
                    
    Inputs:
    1) t           - Vector of order n. Times of observations or measurements 
    2) f           - Vector of order n. The time series values
    3) P           - Vector or Matrix of order n. The weights associated with
                     the observations. Default is 1 (a scale)
    4) tt          - Vector. Equally spaced times for regularization purposes
    5) size        - Numeric. The window or segment size (R)
                     Default is the average sampling rate tripled (3M)
    6) step        - Numeric. The translation step that is less than size
                     Default is the average sampling rate (M)
    7) season      - String. Three options: 
                             ALLSSA season-trend model ('ALLSSA')
                             OLS season-trend model    ('OLS')
                             Only the trend model      ('none')
    8) Omega       - Vector. The initial cyclic frequencies 
                     Default is Omega = [1,2,3,4]   for OLS and 
                     Omega = [Om/5.0 for Om in range(4,20)] for ALLSSA as the initial set   
    9) ind         - Vector. The indices of jumps in the trend component                 
    10)level       - Numeric. The significance level, usually 0.01 or 0.05
                     Default is 0.01. Applies only for ALLSSA not OLS
    11)display     - Logical. To display the result
   
    Outputs:
    1) trend        - The estimated trend component of the time series
    2) seasonal     - The estimated seasonal component of the time series 
    3) remainder    - The estimated remainder component of the time series
    The following outputs will be generated only when tt is given for
    regularization: 
    4) trend_reg    - The regularized trend component of the time series
    5) seasonal_reg - The regularized seasonal component of the time series 

    NOTE: The estimated trend component may have curvature due to the 
          smoothing operation, i.e., when the Gaussian moving average is used
          for all the estimated linear trends within the translating windows.
          Also, at the end of the process, the components will be diplayed 
          in such a way that the seasonal component has zero mean.
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
    #--------- Window Size ----------------------------------------------------------------
    M = int(Ltt/(tt[Ltt-1]-tt[0]))        # A simple estimate for sampling rate
    if size != None:
        if (size > Ltt or size < M):
            raise ValueError("The window size must be an integer between M and n")
    else: size = 3*M
    #--------- Translation Step -----------------------------------------------------------
    if step != None:
        if (step >= size or step < 1):
            raise ValueError("Translation step must be a positive integer less than size")
    else: step = M
    #--------- Set of Cyclic Frequencies --------------------------------------------------
    if len(Omega) == 0:
        if season.upper() == 'ALLSSA': Omega = [Om/5.0 for Om in range(4,20)]
        else: Omega = [1,2,3,4]  
    LOm = len(Omega)
    #--------- Significance Level ---------------------------------------------------------
    if level > 1 or level < 0:
        raise ValueError("The significance level must be between 0 and 1")
    #--------- Jump Indices ---------------------------------------------------------------
    if len(ind) > 0:
        if not (np.floor(ind)==ind).all():
            raise ValueError("The jump indices must be positive integers") 
        elif min(ind) < 2 or max(ind) > Lt-2:
            raise ValueError("The jump indices must be integers between 1 and n-1")
        ind = np.sort(ind)    
    #======================================================================================
    lb = 0; ub = lb + size
    DilFac = np.log(0.1)/((size/(6*M))**2)
    trends0 = np.zeros(Lt); seasonals0 = np.zeros(Lt); remainders0 = np.zeros(Lt)
    trends = np.zeros(Ltt); seasonals = np.zeros(Ltt)   
    sumW0 = np.zeros(Lt); sumW = np.zeros(Ltt)
    stat = True
    start_time = tt[0]
    tt = tt - start_time;  # For the sake of computational efficiency 
    t  = t  - start_time;  # For the sake of computational efficiency     
    while ub <= Ltt + size and stat:
        # The condition for the last segment to also have the same size
        if ub > Ltt: lb = Ltt - size; ub = Ltt; stat = False
        tty = tt[lb:ub]
        if check:           
            lb0 = np.where(t >= tt[lb])[0][0]        
            ub0 = np.where(t >= tt[ub-1])[0][0] + 1                  
            size0  = ub0 - lb0 
        else: lb0 = lb; ub0 = ub; size0  = size
        ty = t[lb0:ub0]
        y  = f[lb0:ub0]
        if CheckP == 0: Py = 1               # If P is a scale 
        elif CheckP == 1: Py = P[lb0:ub0]    # If P is a column vector  
        else: Py = P[lb0:ub0, lb0:ub0]       # If P is a square matrix  
        # Locate the jump index within each segment
        JumpInd = [J-lb0 for J in ind if lb0+2 <= J and J < ub0-1]
        if season.upper() == 'ALLSSA':
            Results = ALLSSA(ty, y, P = Py, Omega = Omega, ind = JumpInd, level = level)
            tr_coeff = Results[0]; CS_coeff = Results[1]; remainder = Results[2]
            freq = Results[5]
        elif season.upper() == 'OLS':
            Results = LSSA(ty, y, P = Py, ind = JumpInd, freq = Omega) 
            remainder = Results[3]; coeff = Results[5]; Lcoeff = len(coeff)
            tr_coeff = coeff[0:len(coeff)-2*LOm]
            CS_coeff = coeff[len(coeff)-2*LOm:Lcoeff]
            freq = Omega
        else:
            Results = LSSA(ty, y, P = Py, ind = JumpInd) 
            remainder = Results[3]; tr_coeff = Results[5]; CS_coeff = []; freq = []
        season_t = SeasonalComponent(ty, freq, CS_coeff)
        trend_t  = y - season_t - remainder
        #------------- Enforce the seasonal component to have zero mean -------------------
        # NOTE: The trend and seasonal components are simultaneously estimated 
        ms = np.mean(season_t); trend_t = trend_t + ms; season_t = season_t - ms
        #----------------------------------------------------------------------------------        
        if size0 % 2: tau0 = ty[int(size0/2)] 
        else: tau0 = (ty[int(size0/2)] + ty[int(size0/2)-1])/2
        GaussW = np.exp(DilFac * (ty - tau0)**2)   # np.ones(len(t0))
        trends0 [lb0:ub0]     = trends0 [lb0:ub0]     + trend_t * GaussW
        seasonals0 [lb0:ub0]  = seasonals0 [lb0:ub0]  + season_t * GaussW
        remainders0 [lb0:ub0] = remainders0 [lb0:ub0] + remainder * GaussW
        sumW0 [lb0:ub0] = sumW0 [lb0:ub0] + GaussW
        if check:
            trend_tt  =  TrendComponent(ty, tty, JumpInd, tr_coeff)
            season_tt =  SeasonalComponent(tty, freq, CS_coeff)
            ms1 = np.mean(season_tt); trend_tt = trend_tt+ms1; season_tt = season_tt-ms1
            if size % 2: tau = tty[int(size/2)] 
            else: tau = (tty[int(size/2)] + tty[int(size/2) - 1])/2
            GW = np.exp(DilFac * (tty - tau)**2) # The Gaussian weight function
            trends[lb:ub]     = trends[lb:ub]     + trend_tt * GW
            seasonals[lb:ub]  = seasonals[lb:ub]  + season_tt * GW
            sumW[lb:ub] = sumW[lb:ub] + GW                
        lb = lb + step; ub = lb + size
    np.seterr(divide='ignore', invalid='ignore')
    trend     = trends0/sumW0                # The weighted average for trend
    seasonal  = seasonals0/sumW0             # The weighted average for seasonal
    remainder = remainders0/sumW0            # The weighted average for remainder 
    if check:
        trend_reg    = trends/sumW           # The weighted average for trend
        seasonal_reg = seasonals/sumW        # The weighted average for seasonal
    else: trend_reg = []; seasonal_reg = [] 

    return trend, seasonal, remainder, trend_reg, seasonal_reg

#==========================================================================================
def PlotDecompose (t, f, P, trend, seasonal, remainder, 
                   tt = [], trend_reg = [], seasonal_reg = []):
    """
    This function displays the decomposition result and saves the result as .png
    Inputs:
    1) t            - Vector of order n. Times of observations/measurements 
    2) f            - Vector of order n. The time series values
    3) P            - Vector or Matrix of order n. The weights associated with
                      the observations. Default is 1 (a scale)
    4) trend        - Vector of order n. The estimated trend component 
    5) seasonal     - Vector of order n. The estimated seasonal component 
    6) remainder    - Vector of order n. The estimated remainder component
    7) tt           - Vector of order m. Equally spaced times 
    8) trend_reg    - Vector of order m. The regularized trend component 
    9) seasonal_reg - Vector of order m. The regularized seasonal component  
    """
    #-------------------------------------------------------------------------------------
    vs = 0.1 * (max(f)-min(f))
    font = {'family' : 'Arial','size'   : 12}
    pyplot.rc('font', **font) 
    SP = np.shape(P); CheckP = len(SP)        
    if CheckP:
        if CheckP > 1: Weight = np.diag(P)
        else: Weight = P
        fig, (ax1, ax, ax2, ax3, ax4) = pyplot.subplots(5, 1, figsize=(7,9), 
                                        sharex='all', gridspec_kw={'hspace': 0.01})
        markerline, stemlines, baseline = ax.stem(t, Weight)
        pyplot.setp(stemlines, color = 'k', linewidth = 0.2)
        pyplot.setp(markerline, color = 'k', markersize = 3)
        pyplot.setp(baseline, color = 'k', linewidth = 0.2, linestyle  = 'solid')
        ax.set(ylabel='Weight')
        ax.grid(True)
    else:
        fig, (ax1, ax2, ax3, ax4) = pyplot.subplots(4, 1, figsize=(7,7), 
                                    sharex='all', gridspec_kw={'hspace': 0.01})
    #-------------------------------------------------------------------------------------- 
    ax1.set_title('Time Series Decomposition')
    ax1.plot(t, f, '--gD', linewidth=1, markersize=3)
    ax1.set(ylabel='Value')
    ax1.grid(True)
    ax1.axis(ymin=min(f)-vs,ymax=max(f)+vs) 
    #--------------------------------------------------------------------------------------
    if len(trend_reg): ax2.plot(tt, trend_reg,'b', linewidth=1)
    else: ax2.plot(t, trend,'b', linewidth=1)
    ax2.set(ylabel='Trend')
    ax2.grid(True)
    ax2.axis(ymin=min(f)-vs,ymax=max(f)+vs)  
    #--------------------------------------------------------------------------------------
    if len(seasonal_reg): ax3.plot(tt, seasonal_reg, '-b', linewidth=1)
    else: ax3.plot(t, seasonal, '-bo', linewidth=1, markersize=3)
    ax3.set(ylabel='Seasonal')
    ax3.grid(True)
    #--------------------------------------------------------------------------------------
    ax4.plot(t, remainder, '--bo', linewidth=1, markersize=3)
    ax4.set(ylabel='Remainder')
    ax4.set(xlabel='Time')
    ax4.grid(True)
    #--------------------------------------------------------------------------------------
    # You may change the directory here to save the results
    fig.savefig("TSdecomp.png",dpi=300, bbox_inches='tight')
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
    parser.add_argument('--Numtime', help = "Numeric. The number of equally spaced times"
                        "for regularization purposes", type = float, default = 0)     
    parser.add_argument('--size', help = "Numeric. Segment size", type = int, default = None)
    parser.add_argument('--step', help = "Numeric. Translation step", type = int, default = None)
    parser.add_argument('--season', help = "Season-trend model: ALLSSA or OLS or None", 
                        type = str, default = "ALLSSA")
    parser.add_argument('--level', help = "Numeric. Significance level", 
                        type = float, default = 0.01) 
    parser.add_argument('--ind', help = "Jump indices", nargs ='+', type = int, default = [])  
    parser.add_argument('--Lfreq', help = "Numeric. Lower cyclic frequency", 
                        type = float, default = 1) 
    parser.add_argument('--Ufreq', help = "Numeric. Upper cyclic frequency", 
                        type = float, default = 0)
    parser.add_argument('--Numfreq', help = "Numeric. The number of cyclic frequencies between "
                        "Lfreq and Ufreq", type = float, default = 0) 
    parser.add_argument('--save', help = "Save the JUST decomposition results as .csv", 
                        type = bool, default = False)
    parser.add_argument('--display', help = "Display the JUST decomposition results", 
                        type = bool, default = False)    
    # parse the argument
    args = parser.parse_args()
    print(args)
    if args.datafile == 0:
        print("ERROR! missing input data file")
        print("Usage: {} datafile".format("JUSTdecompose.py"))
        sys.exit()
        
    filepath = args.datafile
    # Source of the following examples:   https://doi.org/10.3390/rs12234001
    #filepath = "D:\JUST\Sim_JumpInd_48_104.dat"                # [47,103]
    #filepath = "D:\JUST\EVIA_CA.dat"                           # [13, 111, 144]
    
    # May also try the examples above for Omega = [1,2,3,4] when season = ALLSSA 
    # (can be faster but less accurate) 
    
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
    
    if args.Numtime > 0: tt = np.linspace(t[0], t[len(t)-1], args.Numtime)
    else: tt = [] 
    
    if args.Ufreq > 0 and args.Numfreq > 0:
        Omega = np.linspace(args.Lfreq, args.Ufreq, args.Numfreq)
    else: Omega = []
    
    tic = time.clock()    
    Results = JUSTdecompose(t, f, P = P, tt = tt, size = args.size, step = args.step, 
                        season = args.season, Omega = Omega, ind = args.ind,  
                        level = args.level)        
    toc = time.clock()
    print("Computational Time: ", round(toc-tic,2), "s")     
    
    # Save the results as .csv 
    if args.save: 
        import os
        CSVPath = filepath[:-4]
        if not os.path.exists(CSVPath): os.makedirs(CSVPath)                 
        np.savetxt(CSVPath + '\JUSTdecomposition.csv', [r for r in zip(t, Results[0], Results[1], Results[2])], 
                   delimiter=',', header = 'Time, Trend, Seasonal, Remainder', comments="", fmt = '%.4f')        
        if len(tt) > 0:
            np.savetxt(CSVPath + '\JUSTdecompositionRegularized.csv', [r for r in zip(tt, Results[3], Results[4])], 
                       delimiter=',', header = 'Time, Trend, Seasonal', comments="", fmt = '%.4f')  
    # Display the JUST decomposition results       
    if args.display:  
        PlotDecompose (t, f, P, Results[0], Results[1], Results[2], 
                       tt = tt, trend_reg = Results[3], seasonal_reg = Results[4])      
    
  


    
    