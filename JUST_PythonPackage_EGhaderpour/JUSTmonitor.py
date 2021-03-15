"""
JUST Decomposition into Trend, Seasonal, and Remainder Components

This module detects and classifies a possible disturbance within a time series 
in Near-Real-Time (NRT). It requires a stable historical period. 
The time index of start of the stable historical period can be obtained by using JUST 
OR entered by user via other information and/or methods. 
The time index of the start time of the monitoring period needs to be entered as well.

The input time series does not have to be equally spaced, 
and the uncertainties in the time series values can also be considered.

References: 
Ghaderpour, E., and Vujadinovic, T. The potential of the least-squares spectral 
      and cross-wavelet analyses for ear-real-time disturbance detection within 
      unequally spaced satellite image time series. 
      Remote Sensing, 2020, 12, 2446; doi: 10.3390/rs12152446

Ghaderpour, E., Vujadinovic, T. Change Detection within Remotely Sensed 
      Satellite Image Time Series via Spectral Analysis, 
      Remote Sensing, 2020, 12, 4001; doi:10.3390/rs12234001

Author: Ebrahim Ghaderpour 
Email:  ebrahim.ghaderpour@ucalgary.ca
Copyright (c) 2021
"""
#==========================================================================================
import numpy as np
from LSSA import LSSA
from ALLSSA import ALLSSA
import time
from matplotlib import pyplot
from matplotlib.patches import Rectangle
#==========================================================================================
def PlotNRT(t, f, signal, P = 1, start = None, history = 0, Type = 0, Ind = None):
    """
    This function plots the Near-Real-Time (NRT) monitoring results                                            
    Inputs:
    1) t           - Vector of order n. Times of observations or measurements 
    2) f           - Vector of order n. The time series values
    3) signal      - Vector. The least-squares fit to the stable historical time series 
                     extending to the monitoring period
    4) P           - Vector or Matrix of order n. The weights associated with
                     the observations.                    
    5) start       - Numeric. The time index of the start time of monitoring 
    6) history     - Numeric. The time index of start of the stable 
                     historical period.  
    7) Type        - Numeric. The disturbance Type: 0 means no disturbance, 
                     1 means abrupt change, and 2 means phenological change                 
    8) Ind         - Numeric. The time index of detected disturbance  
    """  
    #--------------------------------------------------------------------------------------
    Lt = len(t)
    if start == None: start = Lt - 2       
    Dif = max(f)-min(f)
    xpos = t[start-1] + (t[start]-t[start-1])/2
    SP = np.shape(P); CheckP = len(SP)
    if CheckP > 1:    Std = np.sqrt(np.diag(np.linalg.inv(P)))   # If P is a square matrix    
    elif CheckP == 1: Std = 1/np.sqrt(P)                         # If P is a column vector  
    font = {'family' : 'Arial','size' : 12}
    pyplot.rc('font', **font)             
    fig, ax = pyplot.subplots(figsize = (8,4))         
    if CheckP >= 1:
        MStd = 0.5*Dif * (Std - min(Std))/(max(Std)-min(Std)) + 0.1*Dif
        ax.errorbar(t,f, yerr = MStd, fmt = '--gD', ecolor = 'g', capsize = 2.0,
                    linewidth = 0.5, markersize = 1.5, 
                    markeredgewidth = 0.5, label = 'Time Series & Error bar')   
    else: ax.plot(t,f, '--gD',linewidth = 0.2,markersize = 1.5, label = 'Time Series')
    ax.plot(t[history:Lt], signal, '--ob', linewidth = 0.2, 
            markersize = 1.5, label = 'Signal (Forecast)')
    if Type == 1:
        ax.plot(t[Ind], f[Ind], '*r', markersize = 10, label  = 'Abrupt Change')
    elif Type == 2:
        ax.plot(t[Ind], f[Ind], '*c', markersize = 10, label = 'Intra-annual Change')
    ax.set_title('Near-Real-Time Monitoring') 
    ax.legend(bbox_to_anchor = (0.,.75,1.,.102), loc = 3, ncol= 1, borderaxespad = 0)         
    ax.add_patch(Rectangle((xpos, -2*Dif), 3, 6, color = 'lightgray'))
    pyplot.xlabel('Time')
    pyplot.ylabel('Value')
    pyplot.grid(True)             
    fig.savefig("NRTmonitor.png", dpi = 300, bbox_inches = 'tight')
    pyplot.show()
#==========================================================================================
def JUSTmonitor(t, f, P = 1, start = None, history = 0, h = 0.25, season = 'ALLSSA', 
                Omega = [], trend = 'linear', level = 0.01, display = False):
    """
    Summary: This function detects and classifies a disturbance within a 
             time series in Near-Real-Time (NRT)

             References: 
             Ghaderpour, E., and Vujadinovic, T. The potential of the 
                    least-squares spectral and cross-wavelet analyses for 
                    near-real-time disturbance detection within unequally 
                    spaced satellite image time series. 
                    Remote Sensing, 2020, 12, 2446; doi: 10.3390/rs12152446
             Ghaderpour, E., Vujadinovic, T. Change Detection within Remotely 
                    Sensed Satellite Image Time Series via Spectral Analysis, 
                    Remote Sensing, 2020, 12, 4001; doi:10.3390/rs12234001
                                               
    Inputs:
    1) t           - Vector of order n. Times of observations or measurements 
    2) f           - Vector of order n. The time series values
    3) P           - Vector or Matrix of order n. The weights associated with
                     the observations. Default is none
    4) start       - Numeric. The time index of the start time of monitoring 
    5) history     - Numeric. The time index of start of the stable 
                     historical period. This index can be obtained by JUST OR
                     entered by user via other information and/or methods.
                     Default is 0 which corresponds to t[0]  
    6) h           - Numeric. The fraction of the historical time series 
                     used during monitoring. Default. 0.25                
    7) season      - String. Three options: 
                                 ALLSSA season-trend model ('ALLSSA')
                                 OLS season-trend model    ('OLS')
                                 Only the trend model      ('none')
    8) Omega       - Vector. The initial cyclic frequencies 
                     Default is Omega = [1,2,3,4]   for OLS and 
                     Omega = [Om/5.0 for Om in range(4,20)] for ALLSSA as the initial set                           
    9) trend       - String. 'constant', 'linear'. Default. 'linear'
                     Choose 'constant' to only estimate the jump magnitude
                     Choose 'linear' to estimate the jump magnitude and direction                             
    10)level       - Numeric. The significance level, usually 0.01 or 0.05
                     Default is 0.01. 
    11)display     - Logical. To display the result
     
    Outputs:
    1) Loc         - Detected disturbance location
    2) Mag         - Estimated magnitude of the detected jump
    3) Dir         - Estimated direction of the detected jump only provided
                     when 'trend' is selected to be 'linear'
    4) Type        - The disturbance Type: 0 means no disturbance, 
                     1 means abrupt change, and 2 means phenological change   
    """
    #--------------------------------------------------------------------------------------
    Lt = len(t)    # Number of observations or measurements
    if not (Lt == len(f)): raise ValueError("t and f must have the same size")  
    start_time = t[0]
    t  = t - start_time      # For the sake of computational efficiency
    #--------- Weights (Vector or Matrix)--------------------------------------------------
    SP = np.shape(P); CheckP = len(SP)
    if CheckP == 1 and (not len(P) == Lt):
        raise ValueError("The weights must be a vector of order n") 
    elif CheckP > 1 and not (SP[0] == Lt and SP[1] == Lt):
        raise ValueError("Weight matrix must be a square matrix of order n") 
    #--------- Set of Cyclic Frequencies --------------------------------------------------
    if len(Omega) == 0:
        if season.upper() == 'ALLSSA': Omega = [Om/5.0 for Om in range(4,20)]
        else: Omega = [1,2,3,4]        
    #--------- Time Index of the Start of Monitoring --------------------------------------
    if start == None: start = Lt - 2 
    #--------- Check the Start Time of the Monitoring Period ------------------------------
    Lhis = start - history
    if Lhis < Lt/(t[Lt-1]-t[0]):
        raise ValueError("Insufficient size for the stable historical time series") 
    elif start >= Lt: 
        raise ValueError("The start time index of the monitoring period must be < n") 
    #--------------------------------------------------------------------------------------
    # Size of Series in the Monitoring Process Divided By Size of the Historical Series 
    if h < 0 or h > 1: raise ValueError("h must be in interval [0,1]")  
    #--------- Significance Level ---------------------------------------------------------
    if level > 1 or level < 0:
        raise ValueError("The significance level must be between 0 and 1")
    #======================================================================================
    t1 = t[history:Lt]            
    f1 = f[history:Lt]
    if CheckP == 0:                # If the series is equally weighted 
        P0 = 1
        P1 = 1
    elif CheckP == 1:              # If P is a column vector 
        P0 = P[history:Lt]
        P1 = P0[0:Lhis]
    else:                          # If P is a square matrix     
        P0 = P[history:Lt, history:Lt]
        P1 = P0[0:Lhis, 0:Lhis]
    #---- Estimate the Trend and Seasonal Coefficients for the Stable History Period ------
    if season.upper() == 'ALLSSA':
        Results = ALLSSA(t1[0:Lhis], f1[0:Lhis], P = P1, Omega = Omega, level = level)
        tr_coeff = Results[0]
        CS_coeff = Results[1]
        freq = Results[5]        
        coeff = np.hstack((tr_coeff, CS_coeff))                                  
    elif season.upper() == 'OLS': 
        Results = LSSA(t1[0:Lhis], f1[0:Lhis], P = P1, freq = Omega)
        coeff = Results[5]
        freq = Omega[:]  
    else: 
        Results = LSSA(t1[0:Lhis], f1[0:Lhis], P = P1)
        coeff = Results[5]
        freq = []        
    Lfreq = len(freq)
    #---- Signal (Forecast) ---------------------------------------------------------------
    Loc = 0; Mag = 0; Dir = 0
    nr = Lt - history
    D = np.zeros([nr,2*Lfreq+2])
    D[:,0] = np.ones([1,nr])               # constituent of known form for the intercept
    D[:,1] = t1                            # constituent of known form for slope    
    for q in range(Lfreq):
        D[:,2*q+2] = np.cos(2*np.pi*freq[q]*t1)
        D[:,2*q+3] = np.sin(2*np.pi*freq[q]*t1)        
    signal = np.matmul(D, coeff)           # The history and forecast signal
    #--------------------------------------------------------------------------------------
    #  Near-Real-Time Monitoring
    #--------------------------------------------------------------------------------------
    LSSAOmega = [Om/10.0 for Om in range(1,36)]  # Cyclic frequencies for monitoring 
    # The size of the time series during the monitoring process
    Lh = Lhis - int(h*Lhis) - 1
    if Lh > Lhis - 3: Lh -= 3
    elif Lh == -1: Lh += 1
    for d in range(1, Lt-start+1):
        t2 = t1[Lh:Lhis+d]                     # t2 = t1[Lh+d:Lhis+d] is another option!
        f2 = f1[Lh:Lhis+d] - signal[Lh:Lhis+d]
        if CheckP == 0:   P2 = 1               # If the series is equally weighted 
        elif CheckP == 1: P2 = P0[Lh:Lhis+d]   # If P is a column vector 
        else: P2 = P0[Lh:Lhis+d, Lh:Lhis+d]    # If P is a square matrix 
        Results = LSSA(t2, f2, P = P2, Omega = LSSAOmega, level = level, trend ='constant')        
        spectrum = Results[0]
        CritVal  = Results[1]  
        Ind = start + d - 1
        Loc = t[Ind] + start_time 
        MaxSpec = max(spectrum[0:9])
        if MaxSpec > CritVal: Type = 1
        elif max(spectrum) > CritVal and MaxSpec <= CritVal: Type = 2; break
        else: Type = 0  
        if Type == 1:
            JInd = Lhis + d - 1
            if JInd < len(t1)-2:
                if season.upper() == 'ALLSSA':
                    Result = ALLSSA(t1, f1, P = P0, Omega = Omega, ind = [JInd], 
                                    level = level, trend = trend)
                    tr_coeff = Result[0]
                elif season.upper() == 'OLS': 
                    Result = LSSA(t1,f1,P = P0, ind = [JInd], trend = trend, freq = Omega) 
                    tr_coeff = Result[5]
                else:
                    Result = LSSA(t1, f1, P = P0, ind = [JInd], trend = trend) 
                    tr_coeff = Result[5] 
                if trend.lower() == 'linear':
                    Mag = (tr_coeff[1]-tr_coeff[0]) + (tr_coeff[3]-tr_coeff[2])*t1[JInd]
                    Dir = tr_coeff[3]-tr_coeff[2]
                else:
                    Mag = tr_coeff[1]-tr_coeff[0] 
            break
    #--------------------------------------------------------------------------------------   
    if display: PlotNRT(t+start_time, f, signal, P = P, start = start, 
                        history = history, Type = Type, Ind = Ind)
    #--------------------------------------------------------------------------------------
    return Loc, Mag, Dir, Type
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
    parser.add_argument('--start', help = "Numeric. The start index of monitoring", 
                        type = int, default = None) 
    parser.add_argument('--history', help = "Numeric. The start index of the stable history", 
                        type = int, default = 0) 
    parser.add_argument('--h', help = "Numeric. The fraction of the historical time series " 
                        "used during monitoring", type = float, default = 0.25)                          
    parser.add_argument('--season', help = "Season-trend model: ALLSSA or OLS or None", 
                        type = str, default = "ALLSSA")
    parser.add_argument('--Lfreq', help = "Numeric. Lower cyclic frequency", 
                        type = float, default = 1) 
    parser.add_argument('--Ufreq', help = "Numeric. Upper cyclic frequency", 
                        type = float, default = 0)
    parser.add_argument('--Numfreq', help = "Numeric. The number of cyclic frequencies between "
                        "Lfreq and Ufreq", type = float, default = 0) 
    parser.add_argument('--trend', help = "constant or linear to estimate the magnitude "
                        "and direction of a jump", type = str, default = "Linear")     
    parser.add_argument('--level', help = "Numeric. Significance level", 
                        type = float, default = 0.01) 
    parser.add_argument('--display', help = "Display the JUST decomposition results", 
                        type = bool, default = False)     
    # parse the argument
    args = parser.parse_args()
    print(args)
    if args.datafile == 0:
        print("ERROR! missing input data file")
        print("Usage: {} datafile".format("JUSTmonitor.py"))
        sys.exit()
    
    filepath = args.datafile  
    # Source of the following examples: https://doi.org/10.3390/rs12152446
    #filepath = "D:\JUST\EVIA_AB.dat"                        # start = 72 history = 17
    #filepath = "D:\JUST\EVIB_AB.dat"                        # start = 71 history = 30
    #filepath = "D:\JUST\EVIC_AB.dat"                        # start = 71 history = 15
    
    # Source of the following examples:   https://doi.org/10.3390/rs12234001
    #filepath = "D:\JUST\Sim_JumpInd_48_104.dat"             # Jump Indices: [47,103]
    #filepath = "D:\JUST\EVIA_CA.dat"                        # start = 140, history = 13
    #The detected jump indices for EVIA_CA via JUSTjumps are [13, 111, 144]     

        
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
    
    if args.Ufreq > 0 and args.Numfreq > 0:
        Omega = np.linspace(args.Lfreq, args.Ufreq, args.Numfreq)
    else: Omega = []
    
    tic = time.clock()  
    Res = JUSTmonitor(t, f, P = P, start = args.start, history = args.history, h = args.h,
                      season = args.season, Omega = Omega, trend = args.trend,  
                      level = args.level, display = args.display)        
    toc = time.clock()
    print('Loc = {:0.3f}, Mag = {:0.3f}, Dir = {:0.3f}, Type = {:0.0f}'.format(Res[0], Res[1], Res[2], Res[3]))
    print("Computational Time: ", round(toc-tic,2), "s") 

    
    '''
    # ----- Simulated Example -------------------------------------------------------------    
    t = [k/23.0 for  k in range(95)]
    for kk in range(20):
        indt = int((95-kk)*np.random.rand(1)[0])
        t.remove(t[indt])
    Lt = len(t)
    t = np.array(t) 
    #t.shape = (Lt, 1)
    trend = 0.5 - 0.01*t
    season = 0.1*np.cos(2*np.pi*1*t) + 0.05*np.sin(2*np.pi*2.3*t)
    noise1 = 0.01*np.random.normal(0,1,Lt)
    noise2 = 0.03*np.random.normal(0,1,Lt)
    f = trend + season + noise1 - abs(noise2)
    f[30:Lt] = f[30:Lt] + 0.2
    f[Lt-5:Lt] = f[Lt-5:Lt] - 0.15
    P = 1/(0.01+(noise2)**2)
    '''

  