"""
Turning Point Trend Estimation (TPTR)

This module gets a time series t and f(t) and its estimated turning points (TPs)
to estimate the final linear trend that contains the connected linear pieces.

Reference: 
Ghaderpour, E.; Benedetta, A.; Bozzano, F.; Scarascia Mugnozza, G.; Mazzanti, P.
A Fast and Robust Method for Detecting Trend Turning Points in InSAR Displacement
Time Series, Computers & Geosciences, 2023.

Author: Dr. Ebrahim Ghaderpour 
Email:  ebrahim.ghaderpour@uniroma1.it
Copyright (c) 2023
"""
#==============================================================================================
import numpy as np
from numpy.linalg import inv 
#----------------------------------------------------------------------------------------------
def TPTRF(t,f,TPs):
    """
    Turning Point Trend Forward Estimation (TPTRF)
    
    Summary: This function gets a time series t and f(t) and its estimated turning points (TPs)
             to estimate the final linear trend pieces using the ordinary least-squares (OLS).
             
    Inputs:
    1) t           - Vector of order n. Times of observations or measurements 
    2) f           - Vector of order n. The time series values
    3) TPs         - The potential TPs estimated by STPD
    
    Outputs:
    1) slopes      - The estimated slopes of linear pieces after TPs
    2) y           - The fitted linear trend whose turning points are the TPs

    """
    #--------------------------------------------------------------------------------------------
    Lt = len(t)    # Number of observations or measurements
    if not (Lt == len(f)): raise ValueError("Arrays must have the same size")
    t = np.array(t); f = np.array(f); TPs = np.array(TPs)
    LTPs = len(TPs)
    TPs = np.append(TPs, Lt-1)
    TPs = np.sort(TPs)
    y = np.zeros(Lt)
    if LTPs == 0:
        D1 = np.zeros([Lt,2])
        D1[:,0] = np.ones(Lt)
        D1[:,1] = t
        DT = np.transpose(D1)        
        Ninv  = inv(DT.dot(D1))
        chat = Ninv.dot(DT.dot(f))
        y = D1.dot(chat)
        slopes = chat[1]
    else:
        slopes = []
        TInd = TPs[0] + 1
        D1 = np.zeros([TInd,2])
        D1[:,0] = np.ones(TInd)
        D1[:,1] = t[0:TInd]
        g = f[0:TInd]
        DT = np.transpose(D1)        
        Ninv  = inv(DT.dot(D1))
        chat = Ninv.dot(DT.dot(g))        
        slopes.append(chat[1])
        y[0:TInd] = D1.dot(chat)
        for k in range(LTPs):
            D2 = t[TPs[k]:TPs[k+1]+1]-t[TPs[k]]
            y2 = f[TPs[k]:TPs[k+1]+1]-y[TPs[k]]
            slope = (D2.dot(y2))/(D2.dot(D2))
            y[TPs[k]:TPs[k+1]+1] = D2*slope + y[TPs[k]]
            slopes.append(slope)
    return slopes, y
#================================================================================================
def TPTR(t, f, TPs):
    """
    Turning Point Trend Estimation (TPTR)
    Summary: This function gets a time series t, f(t) and its estimated turning points (TPs)
    to estimate the final linear trend pieces.
    The forward and backward linear trends are estimated and the one
    whose residual series is mimimum will be selected.  
    Inputs:
    1) t           - Vector of order n. Times of observations or measurements 
    2) f           - Vector of order n. The time series values
    3) TPs         - The potential TPs estimated by STPD

    Outputs:
    1) stats       - A matrix whose rows contain TP, the new estimated slope of the second
                     linear piece after the TP, the new estimated TP direction, the new
                     estimated signal-to-noise ratio (SNR), and the new estimated 
                     normalized difference residual index (NDRI)for the TP
    2) y           - The fitted linear trend whose turning points are the TPs
    """
    #--------------------------------------------------------------------------------------------
    Lt = len(t)    # Number of observations or measurements
    if not (Lt == len(f)): raise ValueError("Arrays must have the same size")
    t = np.array(t); f = np.array(f); TPs = np.array(TPs)
    stats = []
    # The forward estimation of linear trend given TPs:
    [slopes1, y1] = TPTRF(t, f, TPs)
    # The backward estimation of linear trend given TPs:
    t0 = np.flip(t); f0 = np.flip(f)
    t0 = t0[0] - t0
    TPs0 = Lt - TPs - 1
    [slopes0, y0] = TPTRF(t0, f0, TPs0)
    slopes2 = -np.flip(slopes0)
    y2 = np.flip(y0)
    # Choose the linear trend (the forward or backward) whose norm is minimum: 
    if np.linalg.norm(f-y1) < np.linalg.norm(f0-y0):
        slopes = slopes1; y = y1;
    else:
        slopes = slopes2; y = y2;
    LTPs = len(TPs)
    TPs = np.append(TPs, [0, Lt - 1])
    TPs = np.sort(TPs)
    nr = np.linalg.norm(f-y); nsig = np.linalg.norm(y); SNR = nsig/nr
    for k in range(1, LTPs+1):
        r1 = f[TPs[k-1]:TPs[k]+1]-y[TPs[k-1]:TPs[k]+1]
        r2 = f[TPs[k]:TPs[k+1]+1]-y[TPs[k]:TPs[k+1]+1]
        r = f[TPs[k-1]:TPs[k+1]+1]-y[TPs[k-1]:TPs[k+1]+1]
        nr1 = np.linalg.norm(r1)
        nr2 = np.linalg.norm(r2)
        nr = np.linalg.norm(r)
        nsig = np.linalg.norm(y[TPs[k-1]:TPs[k+1]+1])
        SNR = nsig/nr
        NDRI = abs(nr2-nr1)/(nr1+nr2)
        stats.append(np.round([TPs[k], slopes[k], slopes[k]-slopes[k-1], SNR, NDRI], 2))
    return stats, y 
#================================================================================================