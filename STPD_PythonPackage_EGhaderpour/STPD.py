"""
Sequential Turning Point Detection (STPD)

This module finds the potential linear trend turning points (TPs) for a given time series.
It uses a windowing strategy where a TP and its direction are estimated within each window.
The TPs are selected based on the following criteria:
1) Signal-to-noise ratio (SNR) 
2) Normalized difference residual index (NDRI)
3) Confidence level
4) TP direction 
5) Minimum distance to the window locations

The input time series does not have to be equally spaced, but for this algorithm,
it is suggested to resample the input time series on regular time intervals, e.g., monthly

Reference: 
Ghaderpour, E.; Benedetta, A.; Bozzano, F.; Scarascia Mugnozza, G.; Mazzanti, P.
A Fast and Robust Method for Detecting Trend Turning Points in InSAR Displacement
Time Series, Computers & Geosciences, 2023.

Author: Dr. Ebrahim Ghaderpour 
Email:  ebrahim.ghaderpour@uniroma1.it
Copyright (c) 2023
"""
#============================================================================================
import numpy as np
from numpy.linalg import inv 
import scipy.stats as ss
#--------------------------------------------------------------------------------------------
def pvalue(t, f, ind, intercept1, slope1, slope2):
    """
    p-value corresponding to TP (t-test)
    Summary: This function gets a segment t and f(t) and its estimated 
             connected linear trend with two pieces and returns its p-value.
    
    Inputs:
    1) t           - Vector of order n. Times of observations or measurements 
    2) f           - Vector of order n. The time series values
    3) ind         - The TP
    4) intercept1  - The estimated intercept of the first linear piece
    5) slope1      - The estimated slope of the first linear piece
    6) slope2      - The estimated slope of the second linear piece
    
    Outputs:
    pval           - The p-value
    """
    #----------------------------------------------------------------------------------------
    Lt = len(t)    # Number of observations or measurements
    if not (Lt == len(f)): raise ValueError("Arrays must have the same size")
    t = np.array(t); f = np.array(f)
    # ---- For the first linear piece: 
    D1 = np.zeros([ind+1,2])
    D1[:,0] = np.ones(ind+1)
    D1[:,1] = t[0:ind+1]
    y1 = f[0:ind+1]
    g1 = y1-D1.dot([intercept1, slope1])
    # --- Two degrees of freedom are lost due to estimation of intercept and slope:
    SE1 = np.sqrt((g1.dot(g1))/(ind-1))/np.linalg.norm(D1[:,1]-np.mean(D1[:,1]))
    # --- For the second linear piece:
    D2 = t[ind:Lt]-t[ind];
    y2 = f[ind:Lt]-intercept1-slope1*t[ind]
    g2 = y2-D2*slope2
    # --- One degree of freedom is lost:
    SE2 = np.sqrt((g2.dot(g2))/(Lt-ind-1))/np.linalg.norm(D2-np.mean(D2)) 
    # --- Student t-test:
    ttest = abs(slope1-slope2)/np.sqrt(SE1**2+SE2**2)
    pval = 2*(1-ss.t.cdf(ttest,Lt-2))
    return pval
#============================================================================================
def TPDetect(t, f, margin):
    """
    Turning Point Detection (TPD)
    Summary: This function gets a segment t and f(t) and finds a potential TP by
             the ordinary least-squares method (OLS).

    Inputs:
    1) t        - Vector of order n. Times of observations or measurements 
    2) f        - Vector of order n. The time series values
    3) margin   - The size of marginal segment where no TP is estimated.
             
    Outputs:
    An array consisting of a potential TP, intercept and slope of the first linear piece,
    the slope of the second linear piece, and the L2 norm of the residual segment.
    """
    #----------------------------------------------------------------------------------------
    Lt = len(t)    # Number of observations or measurements
    if not (Lt == len(f)): raise ValueError("Arrays must have the same size")
    t = np.array(t); f = np.array(f)
    Attributes = []
    r = np.zeros(len(f))
    for ind in range(margin,Lt-margin):
        # Use OLS to estimate the first linear piece from the beginning of segment to ind:
        D1 = np.zeros([ind+1,2])
        D1[:,0]= np.ones(ind+1)
        D1[:,1] = t[0:ind+1]
        y1 = f[0:ind+1]
        DT = np.transpose(D1)        
        Ninv  = inv(DT.dot(D1))
        chat = Ninv.dot(DT.dot(y1))        
        r[0:ind+1] = y1-D1.dot(chat)
        # Use OLS to estimate the second linear piece from the end of the first estimated 
        # linear piece to the end of the segment:
        D2 = t[ind:Lt]-t[ind]
        y2 = f[ind:Lt]-chat[0]-chat[1]*t[ind]
        slope = (D2.dot(y2))/(D2.dot(D2))
        r[ind:Lt] = y2-slope*D2
        # ----------------
        Attributes.append([ind, chat[0], chat[1], slope, np.linalg.norm(r)])

    Attributes = sorted(Attributes, key = lambda d: d[4])
    return Attributes[0] # Choose the TP whose residual norm is minimum
#============================================================================================
def TPD(t, f, margin):
    """
    Turning Point Detection (TPD)
    Summary: This function gets a segment t and f(t) and finds a potential TP.
             It calls the TPDetect function for each of forward and backward
             TP estimations. Then it selects the TP whose residual norm is
             smaller and returns the TP along with its statistics.
    
    Inputs:
    1) t           - Vector of order n. Times of observations or measurements 
    2) f           - Vector of order n. The time series values
    3) margin      - The size of marginal segment where no TP is estimated.
    
    Outputs:
    1) TP          - The potential TP for segment t and f(t)
    2) SNR         - The estimated SNR corresponding to the TP
    3) pval        - The p-value corresponding to the TP
    4) Diff        - The difference between the TPs of forward and backward
    5) DIR         - The TP direction 
    6) SNR         - The absolute value of NDRI corresponding to the TP
    7) y           - The residual segment 
    """
    #----------------------------------------------------------------------------------------
    Lt = len(t)    # Number of observations or measurements
    if not (Lt == len(f)): raise ValueError("Arrays must have the same size")
    t = np.array(t); f = np.array(f)
    # The forward estimation of TP:
    t = t-t[0]
    TPStats = TPDetect(t, f, margin)
    TInd = TPStats[0]
    pval = pvalue(t,f,TPStats[0],TPStats[1],TPStats[2],TPStats[3])
    # The backward estimation of TP:
    t0 = np.flip(t); f0 = np.flip(f)
    t0 = t0[0]-t0
    TPStats0 = TPDetect(t0, f0, margin)
    TInd0 = Lt-TPStats0[0]-1
    pval0 = pvalue(t0,f0,TPStats0[0],TPStats0[1],TPStats0[2],TPStats0[3])
    # The difference between the TPs of forward and backward:
    Diff = abs(TInd0-TInd)
    y = np.zeros(Lt)
    # From the estimated TPs for forward and backward estimations, choose the TP and 
    # its statistics whose residual norm is smaller:
    if TPStats[4] <= TPStats0[4]:
        TP = TPStats[0]; intercept1 = TPStats[1]; slope1 = TPStats[2]; slope2 = TPStats[3]
        g1 = intercept1+slope1*t[0:TP+1]
        g2 = slope2*(t[TP:Lt]-t[TP])+intercept1+slope1*t[TP]
        y[0:TP+1] = g1
        y[TP:Lt]  = g2
        SNR = np.linalg.norm(y)/np.linalg.norm(f-y) 
        nr1 = np.linalg.norm(f[0:TP+1]-g1)
        nr2 = np.linalg.norm(f[TP:Lt]-g2)
        NDRI = abs(nr1-nr2)/(nr1+nr2)
    else:
        TP = TPStats0[0]; intercept1 = TPStats0[1]; slope1=TPStats0[2]; slope2=TPStats0[3]
        g1 = intercept1+slope1*t0[0:TP+1]
        g2 = slope2*(t0[TP:Lt]-t0[TP])+intercept1+slope1*t0[TP]
        y[0:TP+1] = g1
        y[TP:Lt]  = g2
        SNR = np.linalg.norm(y)/np.linalg.norm(f0-y)
        nr1 = np.linalg.norm(f0[0:TP+1]-g1)
        nr2 = np.linalg.norm(f0[TP:Lt]-g2) 
        NDRI = abs(nr1-nr2)/(nr1+nr2)
        pval = pval0
        y = np.flip(y)
        TP = Lt-TP-1
    DIR = slope2-slope1
    return [TP, SNR, pval, Diff, DIR, NDRI, y]
#============================================================================================
def STPD(t,f, size=None,step=None, SNR=1, NDRI=0.3, alpha=0.01, dir_th=0, tp_th=1, margin=1):
    """
    Sequential Turning Point Detection (STPD)
    Summary: This function takes the time series t and f(t) and finds the potential linear 
             trend turning points (TPs). It calls the TPD function for each segment within a 
             translating window to find a TP along with its direction, i.e., the slope of
             the linear piece after the TP minus the one before the TP. 
    The TPs are selected based on the following criteria:
    1) Signal-to-noise ratio (SNR) 
    2) Normalized difference residual index (NDRI)
    3) Confidence level
    4) TP Direction 
    5) Minimum distance to the window locations

    Inputs:
    1) t           - Vector of order n. Times of observations or measurements 
    2) f           - Vector of order n. The time series values
    3) size        - Numeric. The window or segment size (R)
                     Default is five times the average sampling rate (5M)
    4) step        - Numeric. The translation step should be less than size
                     Default is the average sampling rate (M)
    5) SNR         - Numeric. Threshold for the estimated SNR. Default is 1
    6) NDRI        - Numeric. Threshold for |NDRI| to seperate TPs from jumps
                     Default is 0.3
    7) alpha       - Numeric. The significance level, usually 0.01 or 0.05
                     Default is 0.01. 
    8) dir_th      - Numeric. Threshold for the abs of TP direction. 
                     Default is 0 (|Dir| > 0)
    9) tp_th       - Numeric. The minimum time difference between TPs
                     Default is 1 (TPs are at least one year apart)
    10)margin      - The size of marginal segment where no TP is estimated.
                     Default is 1
    
    Outputs:
    TPs            - The indices of the potential TPs

   """
    #----------------------------------------------------------------------------------------
    Lt = len(t)    # Number of observations or measurements
    if not (Lt == len(f)): raise ValueError("Arrays must have the same size") 
    t = np.array(t); f = np.array(f)
    #--------- Window Size ------------------------------------------------------------------
    M = int(Lt/(t[Lt-1]-t[0]))        # A simple estimate for sampling rate
    if size != None:
        if (size > Lt or size < 3):
            raise ValueError("The window size must be an integer between 3 and n")
    elif size == None: size = 5*M
    #--------- Translation Step -------------------------------------------------------------
    if step != None:
        if (step >= size or step < 1):
            raise ValueError("Translation step must be a positive integer less than size")
    else: step = M    
    #--------- Signal-to-Noise Ratio  -------------------------------------------------------
    if SNR <= 0: raise ValueError("SNR must be positive") 
    #--------- Normlized Difference Residual Index (NDRI)  ----------------------------------
    if NDRI > 1 or NDRI < 0: raise ValueError("Threshold for |NDRI| must be in [0,1]")
    #--------- Significance Level -----------------------------------------------------------
    if alpha>1 or alpha<0: raise ValueError("Significance level must be between 0 and 1")
    #--------- TP Direction Threshold -------------------------------------------------------
    if dir_th < 0: raise ValueError("Threshold for |DIR| must be positive")
    #--------- Minimum Time Difference Between Potential TPs --------------------------------
    if tp_th <= 0: 
        raise ValueError("Minimum time difference between potential TPs must be positive")    
    #-----The size of marginal segment where no TP is estimated------------------------------
    if margin >= size/2 or margin < 1:
        raise ValueError("margin should be an integer between 0 and half of window size")  
    #----------------------------------------------------------------------------------------    
    Lt = len(t)    # Number of observations or measurements
    if not (Lt == len(f)): raise ValueError("Arrays must have the same size") 
    inda = 0; indb = inda + size
    TIndDis = []; check = True
    while indb <= Lt + size and check:
        if indb >= Lt:   # The condition for the last segment to also have the same size
            inda = Lt - size; indb = Lt; check = False
        t0 = t[inda:indb]
        f0 = f[inda:indb] 
        # Call function TPD to find a potential TP for each segment:        
        [TInd, SNR0, pval, Diff, Dir0, NDRI0, y] = TPD(t0, f0, margin) 
        # Apply the input thresholds to the estimated TPs:
        if Diff<4 and SNR0>SNR and pval<alpha and abs(Dir0)>dir_th and NDRI0<NDRI:        
            TIndDis.append([TInd+inda, abs(TInd-(size-1)/2.0)])
        inda = inda + step; indb = inda + size      
    TIndDis = sorted(TIndDis, key = lambda d: d[0]) 
    nr = np.size(TIndDis, axis=0)
    if len(TIndDis)==0: TPs=[]
    else: 
        r = 1
        while r < nr:
            # Time difference between TPs should be >= tp_th
            if t[TIndDis[r][0]]- t[TIndDis[r-1][0]] < tp_th:
                group = TIndDis[r-1:r+1][:]
                group = sorted(group, key = lambda d: d[1])
                TIndDis [r-1][:] = group [0][:]
                del TIndDis[r]
                nr -= 1 
            else: r += 1    
        TPs = [TIndDis[k][0] for k in range(nr)]
    return np.array(TPs)
#=============================================================================================