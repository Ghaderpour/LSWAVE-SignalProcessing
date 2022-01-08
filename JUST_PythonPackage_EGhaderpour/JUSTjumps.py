"""
Jumps Upon Spectrum and Trend (JUST)

This module searches a time series and finds potential jumps within the time series 
and estimates the magnitudes and directions of jumps.

There are three search methods in this module: 
1) Anti-Leakage Least-Squares Spectral Analysis (ALLSSA)
2) OLS season-trend model (OLS)
3) Only the trend model      

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
from collections import Counter
from LSSA import LSSA
from ALLSSA import ALLSSA
import time
#------------------------------------------------------------------------------------------
def JumpDetect(t, f, P = 1, season = 'ALLSSA', Omega = [], level = 0.01):
    """
    The sequential approach in JUST
    Summary: A linear trend with two pieces will be fitted sequentially 
             along with the seasonal component based on the user input 
             (i.e., 'ALLSSA', 'OLS', 'none') to estimate the potential jumps 
             in the trend component of time series f. 
                                               
    Inputs:
    1) t           - Vector of order n. Times of observations or measurements 
    2) f           - Vector of order n. The time series values
    3) P           - Vector or Matrix of order n. The weights associated with
                     the observations. Default is none
    4) season      - String. Three options: 
                                 ALLSSA season-trend model ('ALLSSA')
                                 OLS season-trend model    ('OLS')
                                 Only the trend model      ('none')
    5) Omega       - Vector. The initial cyclic frequencies 
                     Default is Omega = [1,2,3,4]   for OLS and 
                     Omega = [Om/5.0 for Om in range(4,20)] for ALLSSA as the initial set                              
    6) level       - Numeric. The significance level, usually 0.01 or 0.05
                     Default is 0.01. Applies only for ALLSSA not OLS
     
    Outputs:
    1) Ind         - Time index of the potential jump (jump index)
    2) Mag         - Estimated magnitude of the potential jump
    3) Dir         - Estimated direction of the potential jump   
    """
    #--------------------------------------------------------------------------------------
    Lt = len(t)    # Number of observations or measurements
    if not (Lt == len(f)): raise ValueError("Arrays must have the same size")    
    #--------- Weights (Vector or Matrix)--------------------------------------------------
    SP = np.shape(P)
    CheckP = len(SP)
    if CheckP == 1 and (not len(P) == Lt):
        raise ValueError("The weights must be a vector of order n") 
    elif CheckP > 1 and not (SP[0] == Lt and SP[1] == Lt):
        raise ValueError("Weight matrix must be a square matrix of order n") 
    #--------- Set of Cyclic Frequencies --------------------------------------------------
    if len(Omega) == 0:
        if season.upper() == 'ALLSSA': Omega = [Om/5.0 for Om in range(4,20)]
        else: Omega = [1,2,3,4] 
    LOm = len(Omega)
    #--------- Significance Level ---------------------------------------------------------
    if level > 1 or level < 0:
        raise ValueError("The significance level must be between 0 and 1")
    #======================================================================================
    Attributes = []      
    for Ind in range(3,Lt-2):
        if season.upper() == 'ALLSSA':
            Results = ALLSSA(t, f, P = P, Omega = Omega, ind = [Ind], level = level)
            norm_res = Results[3]
            tr_coeff = Results[0]
        elif season.upper() == 'OLS':
            Results = LSSA(t, f, P = P, ind = [Ind], freq = Omega) 
            norm_res = Results[4]
            coeff = Results[5]
            tr_coeff = coeff[0:len(coeff)-2*LOm]
        else:
            Results = LSSA(t, f, P = P, ind = [Ind]) 
            norm_res = Results[4]
            tr_coeff = Results[5]
        Mag = (tr_coeff[1]-tr_coeff[0])+(tr_coeff[3]-tr_coeff[2])*t[Ind]
        Dir = tr_coeff[3]-tr_coeff[2]
        Attributes.append([Ind, Mag, Dir, norm_res])
    # Ascending sort according to norm_res
    Attributes = sorted(Attributes, key = lambda d: d[3])
    return [Attributes[0][0], Attributes[0][1], Attributes[0][2]]
#==========================================================================================
def AllJumps(t, f, P = 1, size = None, step = None, season = 'ALLSSA', Omega = [], 
             level = 0.01):
    """
    The segmentation or windowing approach in JUST
    Summary: A window of fixed size translates by step over time.
             For each segment within a translating window, 
             the JumpDetect function (the sequential approach) will be called 
             to find a jump candidate along with its magnitude and direction.
                    
    Inputs:
    1) t           - Vector of order n. Times of observations or measurements 
    2) f           - Vector of order n. The time series values
    3) P           - Vector or Matrix of order n. The weights associated with
                     the observations. Default is 1 (a scale)
    4) size        - Numeric. The window or segment size (R)
                     Default is the average sampling rate tripled (3M)
    5) step        - Numeric. The translation step that is less than R 
                     Default is the average sampling rate (M)
    6) season      - String. Three options: 
                             ALLSSA season-trend model ('ALLSSA')
                             OLS season-trend model    ('OLS')
                             Only the trend model      ('none')
    7) Omega       - Vector. The initial cyclic frequencies 
                     Default is Omega = [1,2,3,4]   for OLS and 
                     Omega = [Om/5.0 for Om in range(4,20)] for ALLSSA as the initial set                         
    8) level       - Numeric. The significance level, usually 0.01 or 0.05
                     Default is 0.01. Applies only for ALLSSA not OLS

    Output:
    IndMagDirDis   - A matrix with four columns
                     Each row contains the attributes of a potential jump:
                          1st column: Time index of the jump
                          2nd column: Estimated magnitude of the jump
                          3rd column: Estimated direction of the jump
                          4th column: Distance between the time index of the  
                                      jump and the time index of its 
                                      corresponding window location 
    """
    #--------------------------------------------------------------------------------------
    Lt = len(t)    # Number of observations or measurements
    if not (Lt == len(f)): raise ValueError("t and f must have the same size")    
    #--------- Weights (Vector or Matrix)--------------------------------------------------
    SP = np.shape(P); CheckP = len(SP)
    if CheckP == 1 and (not len(P) == Lt):
        raise ValueError("The weights must be a vector of order n") 
    elif CheckP > 1 and not (SP[0] == Lt and SP[1] == Lt):
        raise ValueError("Weight matrix must be a square matrix of order n") 
    #--------- Window Size ----------------------------------------------------------------
    M = int(Lt/(t[Lt-1]-t[0]))        # A simple estimate for sampling rate
    if size != None:
        if (size > Lt or size < M):
            raise ValueError("The window size must be an integer between M and n")
    elif size == None: size = 3*M
    #--------- Translation Step -----------------------------------------------------------
    if step != None:
        if (step >= size or step < 1):
            raise ValueError("Translation step must be a positive integer less than size")
    else: step = M
    #--------- Set of Cyclic Frequencies --------------------------------------------------
    if len(Omega) == 0:
        if season.upper() == 'ALLSSA': Omega = [Om/5.0 for Om in range(4,20)]
        else: Omega = [1,2,3,4] 
    #--------- Significance Level ----------------------------------------------------------
    if level > 1 or level < 0:
        raise ValueError("The significance level must be between 0 and 1")
    #======================================================================================
    inda = 0; indb = inda + size
    IndMagDirDis = []; check = True
    while indb <= Lt + size and check:
        if indb > Lt:   # The condition for the last segment to also have the same size
            inda = Lt - size; indb = Lt; check = False
        t0 = t[inda:indb]
        f0 = f[inda:indb] 
        if CheckP == 0:   P0 = 1             # If P is a scale
        elif CheckP == 1: P0 = P[inda:indb]  # If P is a column vector  
        else: P0 = P[inda:indb, inda:indb]   # If P is a square matrix   
        IndMagDir = JumpDetect(t0, f0, P = P0, Omega = Omega, 
                                   season = season, level = level)        
        IndMagDirDis.append([IndMagDir[0]+inda, IndMagDir[1], 
                              IndMagDir[2], abs(IndMagDir[0]-(size-1)/2.0)])
        inda = inda + step; indb = inda + size      
    IndMagDirDis = sorted(IndMagDirDis, key = lambda d: d[0])
    return IndMagDirDis     
#==========================================================================================
def JUSTjumps(t, f, P = 1, size = None, step = None, season = 'ALLSSA', Omega = [], 
              level = 0.01, mag_th = 0.05, dir_th = 0.01, jump_th = 0.75, saveto = 0):
    """
    Jumps Upon Spectrum and Trend (JUST)
    Summary: This function calls the AllJumps function (the windowing or
             segmentation approach). For each segment within a translating 
             window, the AllJumps function calls the JumpDetect function 
             (the sequential approach) to find a jump candidate along with 
             its magnitude and direction. 
             Among the detected jumps within the translating windows,
             certain jumps will be selected based on the following criteria:
                1) Maximum occurrence 
                2) Minimum distance to the window locations
                3) Magnitudes and Directions
             Thresholds for the minimum time difference between
             the potential jumps and magnitude and direction will be applied 
             at this stage.
    
             Reference: 
             Ghaderpour, E., Vujadinovic, T. Change Detection within Remotely 
                    Sensed Satellite Image Time Series via Spectral Analysis, 
                    Remote Sensing, 2020, 12, 4001; doi:10.3390/rs12234001
                    
    Inputs:
    1) t           - Vector of order n. Times of observations or measurements 
    2) f           - Vector of order n. The time series values
    3) P           - Vector or Matrix of order n. The weights associated with
                     the observations. Default is 1 (a scale)
    4) size        - Numeric. The window or segment size (R)
                     Default is the average sampling rate tripled (3M)
    5) step        - Numeric. The translation step that is less than R 
                     Default is the average sampling rate (M)
    6) season      - String. Three options: 
                             ALLSSA season-trend model ('ALLSSA')
                             OLS season-trend model    ('OLS')
                             Only the trend model      ('none')
    7) Omega       - Vector. The initial cyclic frequencies 
                     Default is Omega = [1,2,3,4] for OLS and 
                     Omega = [Om/5.0 for Om in range(4,20)] for ALLSSA as the initial set                        
    8) level       - Numeric. The significance level, usually 0.01 or 0.05
                     Default is 0.01. Applies only for ALLSSA not OLS
    9) mag_th      - Numeric. The threshold for the absolute value of jump magnitude
                     Default is 0.05 (|Mag| >= 0.05)
    10)dir_th      - Numeric. The threshold for the absolute value of jump direction
                     Default is 0.01 (|Dir| >= 0.01)
    11)jump_th     - Numeric. The minimum time difference between jumps
                     Default is 0.75 (e.g., jumps are at least 9 months apart)
    12)saveto      - String. Save the detected jumps and their attributes; Default: 0

    Output:
    LocIndMagDir   - A matrix with four columns
                     Each row contains the attributes of a potential jump:
                          1st column: Jump location
                          2nd column: Jump index
                          3rd column: Estimated magnitude of the jump
                          4th column: Estimated direction of the jump
    """
    #--------------------------------------------------------------------------------------
    Lt = len(t)    # Number of observations or measurements
    if not (Lt == len(f)): raise ValueError('t and f must have the same size')   
    start_time = t[0]
    t = t - start_time              # For the sake of computational efficiency     
    #--------- Weights (Vector or Matrix)--------------------------------------------------
    SP = np.shape(P); CheckP = len(SP)
    if CheckP == 1 and (not len(P) == Lt):
        raise ValueError("The weights must be a vector of order n") 
    elif CheckP > 1 and not (SP[0] == Lt and SP[1] == Lt):
        raise ValueError("Weight matrix must be a square matrix of order n") 
    #--------- Window Size ----------------------------------------------------------------
    M = int(Lt/(t[Lt-1]-t[0]))        # A simple estimate for sampling rate
    if size != None:
        if (size > Lt or size < M):
            raise ValueError("The window size must be an integer between M and n")
    elif size == None: size = 3*M
    #--------- Translation Step -----------------------------------------------------------
    if step != None:
        if (step >= size or step < 1):
            raise ValueError("Translation step must be a positive integer less than size")
    else: step = M
    #--------- Set of Cyclic Frequencies --------------------------------------------------
    if len(Omega) == 0:
        if season.upper() == 'ALLSSA': Omega = [Om/5.0 for Om in range(4,20)]
        else: Omega = [1,2,3,4] 
    #--------- Significance Level ---------------------------------------------------------
    if level > 1 or level < 0:
        raise ValueError("The significance level must be between 0 and 1")
    #======================================================================================  
    IndMagDirDis = AllJumps(t, f, P = P, size = size, step = step, Omega = Omega, 
                             season = season, level = level) 
    t = t + start_time
    LJ = len(IndMagDirDis)
    JumpInd = [IndMagDirDis[ind][0] for ind in range(LJ)]
    JumpDict = dict(Counter(JumpInd))
    IndOccMagDirDis= []; results = []; val = 0
    # The following loop will select a window whose location is closest to the 
    # jump location with occurence greater than one in order to assign
    # the attributes of the jump detected within that window.    
    for key, value in sorted(JumpDict.items()):
        MinDist = sorted(IndMagDirDis[val:val+value][:], key = lambda d: d[3])
        Mag = round(MinDist[0][1],2)
        Dir = round(MinDist[0][2],2)
        Dis = int(MinDist[0][3])
        val += value        
        IndOccMagDirDis.append([key, value, Mag, Dir, Dis])   
        if saveto: results.append([round(t[key],3), key, value, Mag, Dir, Dis])
    LJ = len(IndOccMagDirDis)   
    for k in range(LJ): IndOccMagDirDis[k][4] = int(IndOccMagDirDis[k][4]/(size/6.0))
    # The following while loop will select the optimal breakpoints whose
    # distances are greater than or equal to jump_th
    r = 1
    while r < LJ:
        if t[IndOccMagDirDis[r][0]] - t[IndOccMagDirDis[r-1][0]] < jump_th:
            group = IndOccMagDirDis[r-1:r+1][:]
            # Sort the jump occurrences descendingly, then sort ascendingly the 
            # scaled jump distances to the window locations where they were 
            # estimated, then sort abs magnitudes descendingly
            group = sorted(group, key = lambda d: d[1], reverse = True)
            if group[0][1] == group[1][1]: 
                group = sorted(group, key = lambda d: d[4]) 
                if group[0][4] == group[1][4]: 
                    group = sorted(group, key = lambda d: abs(d[2]), reverse=True)  
            IndOccMagDirDis [r-1][:] = group [0][:]
            del IndOccMagDirDis[r]
            LJ -= 1 
        else: r += 1
    # The following while loop applies the magnitude and direction thresholds
    r = 0
    while r < LJ:
        if (abs(IndOccMagDirDis[r][2]) < mag_th and abs(IndOccMagDirDis[r][3]) < dir_th): 
            del IndOccMagDirDis[r]
            LJ -= 1 
        else: r += 1 
    LocIndMagDir = ([[round(t[IndOccMagDirDis[k][0]],3), IndOccMagDirDis[k][0], 
                      IndOccMagDirDis[k][2], IndOccMagDirDis[k][3]] for k in range(LJ)]) 
    if saveto:
        log = open(saveto, "w")  
        log.write("Jump Location    Jump Index    Occurrence    ")
        log.write("Magnitude    Direction    Distance")
        [log.write("\n  {a1:<10}{a2:>11}{a3:>13}{a4:>16}{a5:>13}{a6:>11}".format(
            a1 = results[k][0], a2 = results[k][1], a3 = results[k][2], a4 = results[k][3], 
            a5 = results[k][4], a6 = results[k][5])) for k in range(len(results))]
        log.write("\n\nSelected jumps after applying the thresholds:\n")
        log.write("Jump Location    Jump Index    Magnitude     Direction")
        [log.write("\n  {a1:<10}{a2:>11}{a3:>15}{a4:>14}".format(
            a1 = LocIndMagDir[k][0], a2 = LocIndMagDir[k][1], 
            a3 = LocIndMagDir[k][2], a4 = LocIndMagDir[k][3])) for k in range(LJ)]    
        log.close() 
    return LocIndMagDir  
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
    parser.add_argument('--size', help = "Numeric. Segment size", type = int, default = None)
    parser.add_argument('--step', help = "Numeric. Translation step", type = int, default = None)
    parser.add_argument('--season', help = "Season-trend model: ALLSSA or OLS or None", 
                        type = str, default = "ALLSSA")
    parser.add_argument('--level', help = "Numeric. Significance level", 
                        type = float, default = 0.01) 
    parser.add_argument('--mag_th', help = "Numeric. Jump magnitude threshold", 
                        type = float, default = 0.05) 
    parser.add_argument('--dir_th', help = "Numeric. Jump direction threshold", 
                        type = float, default = 0.01)  
    parser.add_argument('--jump_th', help = "Numeric. Minimum time difference between jumps as "
                        "a multiple of unit time", type = float, default = 0.75) 
    parser.add_argument('--save', help = "Save the jump results", 
                        type = bool, default = False) 
    parser.add_argument('--Lfreq', help = "Numeric. Lower cyclic frequency", 
                        type = float, default = 1) 
    parser.add_argument('--Ufreq', help = "Numeric. Upper cyclic frequency", 
                        type = float, default = 0)
    parser.add_argument('--Numfreq', help = "Numeric. The number of cyclic frequencies between "
                        "Lfreq and Ufreq", type = float, default = 0) 
    # parse the argument
    args = parser.parse_args()
    print(args)
    if args.datafile == 0:
        print("ERROR! missing input data file")
        print("Usage: {} datafile".format("JUSTjumps.py"))
        sys.exit()
    
    filepath = args.datafile
    #Source of the following examples:   https://doi.org/10.3390/rs12234001
    #filepath = "D:\JUST\EVIA_CA.dat"
    #filepath = "D:\JUST\EVIB_CA.dat"
    #filepath = "D:\JUST\EVIC_CA.dat"
    #filepath = "D:\JUST\ExA_NDVI_Australia.dat"   
    #filepath = "D:\JUST\ExB_NDVI_Australia.dat"
    #filepath = "D:\JUST\Sim_JumpInd_48_104.dat"
    
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
    
    if args.save:
        import os
        TxtPath = filepath[:-4]
        if not os.path.exists(TxtPath): os.makedirs(TxtPath) 
        logname = TxtPath + '\JumpResultsPython.txt'
    else: logname = 0
    
    if args.Ufreq > 0 and args.Numfreq > 0:
        Omega = np.linspace(args.Lfreq, args.Ufreq, args.Numfreq)
    else: Omega = []
    
    tic = time.clock()
    LocIndMagDir  = JUSTjumps(t, f, P = P, size = args.size, step = args.step, 
                              season = args.season, Omega = Omega, level = args.level, 
                              mag_th = args.mag_th, dir_th = args.dir_th, 
                              jump_th = args.jump_th, saveto = logname)
    print("Jump Location, Jump Index, Magnitude, Direction:\n", LocIndMagDir)
    toc = time.clock()
    print("Computational Time: ", round(toc-tic,2), "s")     
   