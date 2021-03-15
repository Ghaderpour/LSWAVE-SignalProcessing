"""
Least-Squares Spectral Analysis (LSSA)

This module computes a frequency spectrum for a given time series based on 
the least-squares fit of sinusoids of different cyclic frequencies to the time series. 
The constitunets of known forms, such as trends and sinusoids of known frequencies
are fitted first which is the Ordinary Least-Squares (OLS) estimation. 
Then, the spectrum of the residual series will be estimated considering 
the effect of removed known constituents.

The input time series does not have to be equally spaced, 
and the uncertainties in the time series values can also be considered.

Author: Ebrahim Ghaderpour 
Email:  ebrahim.ghaderpour@ucalgary.ca
Copyright (c) 2021
"""
#==========================================================================================
import numpy as np
from numpy.linalg import inv as inv
# The user may import pinv as inv for certain applications (special cases), 
# however, pinv is generally slower than inv
import time
from matplotlib import pyplot
#------------------------------------------------------------------------------------------
def LSSA(t, f, P = 1, Omega = [], ind = [], level = 0.01, 
        trend = 'linear', slope = False, freq = [], saveto = 0):
    """
    Summary: This function computes a frequency spectrum for a given time
             series. The constitunets of known forms, such as trends and 
             sinusoids of known frequencies are fitted first  
   
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
    9) freq        - Vector. The known cyclic frequencies of the sinusoids
                     to be estimated/removed from the time series along with 
                     the trend while their removal effect is considered 
                     in the spectrum of residual
    10) saveto     - String. Save the OLS results as .xlsx. Default: 0
    
    Outputs:
    1) Spectrum    - The normalized least-squares spectrum (LSS)
    2) CritVal     - The Critical Value at (1-level) confidence level
    3) CS_coeff    - Estimated coefficients of cosine and sine functions of
                     input cyclic frequencies (Omega) to determine the
                     amplitude and phase for each frequency (out-of-context)
    4) res         - The residual series
    5) norm_res    - The squared weighted L2 norm of the residual series 
    6) coeff       - The ordinary least-squares estimated coefficients of 
                     the trend and sinusoids of the known cyclic frequencies 
                     ordered as follows: 
                     For constant trend with multiple pieces: shifts
                     For linear trend with multiple pieces: intercepts,slopes 
                     For the quadratic trend (a+b*t+c*t^2): a, b, c
                     For cubic trend (a+b*t+c*t^2+d*t^3): a, b, c, d
                     Followed by the estimated coefficients of the cosine and 
                     sine functions of the known cyclic frequencies if given
    7) cov         - Covariance matrix of the estimated coefficients (coeff)
                     with the same order as coeff 
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
    #======================================================================================
    Omega = np.sort(Omega); LOm = len(Omega); Lfreq = len(freq); Lind  = len(ind);
    CScoeff = []; spectrum = []; coeff = []
    #--------------------------------------------------------------------------------------    
    # Determine the number of constituents of known forms
    nc = 2*Lfreq;
    if   trend.lower() == 'constant':                 nc = nc + Lind + 1
    elif trend.lower() == 'linear' and slope:         nc = nc + Lind + 2
    elif trend.lower() == 'linear' and not slope:     nc = nc + 2*Lind + 2
    elif trend.lower() == 'quadratic' and slope:      nc = nc + Lind + 3
    elif trend.lower() == 'quadratic' and not slope:  nc = nc + 2*Lind + 3
    elif trend.lower() == 'cubic' and slope:          nc = nc + Lind + 4
    elif trend.lower() == 'cubic' and not slope:      nc = nc + 2*Lind + 4
    if Lt < nc+2: raise ValueError("The size of the time series (segment) is too small!")
    # -------------------------------------------------------------------------------------
    #     Create a design matrix (D) based on the constituents of known forms 
    #--------------------------------------------------------------------------------------
    # ----Constituent of intercept (Datum shift)
    D = np.zeros([Lt,nc]); q = 0;
    if (trend.lower() == 'constant' or trend.lower() == 'linear' or 
        trend.lower() == 'quadratic' or trend.lower() == 'cubic') and Lind == 0:                           
        D[:,q] = np.ones([1,Lt]) # constituent of known form for the intercept
        q += 1
    elif (trend.lower() == 'constant' or trend.lower() == 'linear' or 
          trend.lower() == 'quadratic' or trend.lower() == 'cubic') and Lind > 0: 
        # constituents of known forms for intercepts
        D[0:ind[0],q] = np.ones([1,ind[0]])
        q += 1
        for j in range(Lind):
            if j < Lind-1: 
                D[ind[j]:ind[j+1],q] = np.ones([1,ind[j+1]-ind[j]])
                q += 1
            else: 
                D[ind[j]:Lt,q] = np.ones([1, Lt-ind[j]])
                q += 1
    # ----Constituent of slope
    if (trend.lower() == 'linear' or trend.lower() == 'quadratic' or 
        trend.lower() == 'cubic') and (Lind == 0 or slope): 
        D[:,q] = t
        q += 1          
    elif (trend.lower() == 'linear' or trend.lower() == 'quadratic' or 
        trend.lower() == 'cubic') and (Lind > 0 and not slope): 
        D[0:ind[0],q] = t[0:ind[0]]
        q += 1        
        for j in range(Lind):
            if j < Lind-1: 
                D[ind[j]:ind[j+1],q] = t[ind[j]:ind[j+1]]
                q += 1
            else: 
                D[ind[j]:Lt,q] = t[ind[j]:Lt]
                q += 1           
    #-----Constituent of quadratic term
    if trend.lower() == 'quadratic' or trend.lower() == 'cubic':
        D[:,q] = t**2
        q += 1  
    #------Constituent of cubic term
    if trend.lower() == 'cubic':
        D[:,q] = t**3
        q += 1     
    #------Constituent of sinusoids of known frequencies
    for jj in range(Lfreq):
        D[:,q]   = np.cos(2*np.pi*freq[jj]*t)
        D[:,q+1] = np.sin(2*np.pi*freq[jj]*t)
        q += 2
    D  = np.array(D)
    DT = np.transpose(D)
    #--------------------------------------------------------------------------------------
    #  Estimate the coefficients (coeff) for the constituents of known forms
    #--------------------------------------------------------------------------------------
    if nc > 0:
        if CheckP == 0:    # If the time series is equally weighted
            Ninv  = inv(DT.dot(D))
            coeff = Ninv.dot(DT.dot(f))
        elif CheckP == 1:  # If P is a column vector
            KronP = np.repeat(P,nc).reshape(SP[0], nc)
            Ninv  = inv(DT.dot(D*KronP))
            coeff = Ninv.dot(DT.dot(f*P))
        else:              # If P is a square matrix   
            Ninv  = inv(DT.dot(P.dot(D)))
            coeff = Ninv.dot(DT.dot(P.dot(f)))
    #--------------------------------------------------------------------------------------
    # Find the residual series: after removing the constituents of known forms
    if nc == 0: res = f
    else: res = f - D.dot(coeff)
    #--------------------------------------------------------------------------------------
    if CheckP == 0: resP = res             # If the time series is equally weighted 
    elif CheckP == 1: resP = res * P       # If P is a column vector  
    else: resP = P.dot(res)                # If P is a square matrix
    norm_res = np.dot(res,resP)
    #---------------------------------------------------------------------------------------
    # Estimate the covariance matrix of coeff
    if nc > 0: cov = (norm_res/(Lt-nc)) * Ninv
    else: cov = []
    # --------------------------------------------------------------------------------------
    # Save the OLS results as .xlsx or .mat
    # Source: https://doi.org/10.1007/s10291-019-0841-3
    if saveto: 
        import xlsxwriter
        LTr = len(coeff)-2*len(freq); Tr_coeff = coeff[:LTr]; CS_coeff = coeff[LTr:];
        R1 = Tr_coeff[:Lind+1]; R2 = []; R3 = None; R4 = None
        R5 = []; R6 = []; R7 = []; R8 = []; R9 = []
        diagCov = np.diag(cov); Tr_std = np.sqrt(diagCov[0:LTr]); 
        CS_cov = cov[LTr:,LTr:]; CS_var = np.diag(CS_cov)
        S1 = Tr_std[:Lind+1]; S2 = []; S3 = None; S4 = None
        S6 = []; S7 = []; S8 = []; S9 = []
        
        if   trend.lower() == 'linear':  
            R2 = Tr_coeff[Lind+1:]; S2 = Tr_std[Lind+1:]
        elif trend.lower() == 'quadratic' and slope:     
            R2 = Tr_coeff[Lind+1:Lind+2]; R3 = Tr_coeff[Lind+2]
            S2 = Tr_std[Lind+1:Lind+2]; S3 = Tr_std[Lind+2]
        elif trend.lower() == 'quadratic' and not slope: 
            R2 = Tr_coeff[Lind+1:2*Lind+2]; R3 = Tr_coeff[2*Lind+2]
            S2 = Tr_std[Lind+1:2*Lind+2]; S3 = Tr_std[2*Lind+2]
        elif trend.lower() == 'cubic' and slope:    
            R2 = Tr_coeff[Lind+1:Lind+2]; R3 = Tr_coeff[Lind+2]; R4 = Tr_coeff[Lind+3]
            S2 = Tr_std[Lind+1:Lind+2]; S3 = Tr_std[Lind+2]
            S4 = Tr_std[Lind+3] 
        elif trend.lower() == 'cubic' and not slope:     
            R2 = Tr_coeff[Lind+1:2*Lind+2]; R3 = Tr_coeff[2*Lind+2]; R4 = Tr_coeff[2*Lind+3]
            S2 = Tr_std[Lind+1:2*Lind+2]; S3 = Tr_std[2*Lind+2]
            S4 = Tr_std[2*Lind+3]            
        R5 = freq; R6 = CS_coeff[0::2]; R7 = CS_coeff[1::2]; 
        R8 = [np.sqrt(CS_coeff[2*k]**2+CS_coeff[2*k+1]**2) for k in range(len(freq))] 
        R9 = 2*np.arctan((R8-R7)/R6)
        S6 = np.sqrt(CS_var[0::2]); S7 = np.sqrt(CS_var[1::2]);
        for k in range(len(freq)):
            c = R6[k]; s = R7[k];  amp = R8[k]
            ce = S6[k]; se = S7[k]; cse = CS_cov[2*k,2*k+1]
            AmpError = np.sqrt((c*ce)**2+(s*se)**2 + 2*c*s*cse)/amp
            if c < 1e-10: c = 1e-10
            PhaseError = 2*(np.sqrt(((s/c)*(amp-s))**2*ce**2+(s-amp)**2*se**2-2*(s/c)*(amp-s)**2*cse)
                            /(abs(amp*c)*(1+(amp-s)**2/c**2)))
            S8.append(AmpError)
            S9.append(PhaseError)
        
        workbook = xlsxwriter.Workbook(saveto)
        ws = workbook.add_worksheet()
        ResList = ["Intercept", "Slope", "Quad_Coeff", "Cub_Coeff", 
                   "Cos_Coeff", "Sin_Coeff", "Amplitude", "Phase"]
        for k in range(4): ws.write(0,2*k,ResList[k]); ws.write(0,2*k+1,"Error")
        ws.write(0,8,"Frequency"); ws.write(0,17,"Residual Norm"); ws.write(1,17,round(norm_res,4));
        for k in range(4,8): ws.write(0,2*k+1,ResList[k]); ws.write(0,2*k+2,"Error")
        for k in range(len(R1)): ws.write(k+1,0,round(R1[k],4));  ws.write(k+1,1,round(S1[k],6))
        for k in range(len(R2)): ws.write(k+1,2,round(R2[k],4));  ws.write(k+1,3,round(S2[k],6))
        if not R3 == None: ws.write(1,4,round(R3,4));  ws.write(1,5,round(S3,6))
        if not R4 == None: ws.write(1,6,round(R4,4));  ws.write(1,7,round(S4,6))
        for k in range(len(R5)): 
            ws.write(k+1,8,round(R5[k],4))  
            ws.write(k+1,9,round(R6[k],4)); ws.write(k+1,10,round(S6[k],6))
            ws.write(k+1,11,round(R7[k],4)); ws.write(k+1,12,round(S7[k],6))
            ws.write(k+1,13,round(R8[k],4)); ws.write(k+1,14,round(S8[k],6)) 
            ws.write(k+1,15,round(R9[k],4)); ws.write(k+1,16,round(S9[k],6))     
        workbook.close()     
    #--------------------------------------------------------------------------------------
    # Critical value at (1-level) confidence level for the normalized LSS
    CritVal = 1 - pow(level, (2/(Lt-nc-2)))
    #--------------------------------------------------------------------------------------
    #             Calculate the least-squares spectrum (LSS)  
    #--------------------------------------------------------------------------------------
    for k in range(LOm):
        if Omega[k] not in freq:
            Phi1 = np.cos(2*np.pi*Omega[k]*t)
            Phi2 = np.sin(2*np.pi*Omega[k]*t)
            Phi  = np.array([Phi1,Phi2]) 
            PhiT = np.transpose(Phi)
            if CheckP == 0:    # If the time series is equally weighted
                N22 = Phi.dot(PhiT)
                if nc > 0: N12 = DT.dot(PhiT)
            elif CheckP == 1:  # If P is a column vector
                KronP = np.repeat(P,2).reshape(SP[0], 2)
                N22 = Phi.dot(PhiT*KronP)
                if nc > 0: N12 = DT.dot(PhiT*KronP)
            else:              # If P is a square matrix
                N22 = Phi.dot(P.dot(PhiT))
                if nc > 0: N12 = DT.dot(P.dot(PhiT))
            N12T = np.transpose(N12)
            if nc > 0: chat = inv(N22-(N12T.dot(Ninv)).dot(N12)).dot(Phi.dot(resP))
            else:      chat = inv(N22).dot(Phi.dot(resP))
            spec = ((Phi.dot(resP)).dot(chat))/norm_res    # The normalized LSS
            if   spec < 0.0:  spectrum.append(0.0)   # It may happen due to singularity   
            elif spec > 1.0:  spectrum.append(1.0)   # It may happen due to singularity
            else:             spectrum.append(spec)
            CScoeff.append(chat)  # The estimated cosine and sine coeff
        else:
            spectrum.append(0.0) 
            CScoeff.append([0.0,0.0])

    return np.array(spectrum), CritVal, CScoeff, res, norm_res, coeff, cov
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
    parser.add_argument('--slope', help = "True estimates a unique slope, and "
                        "False allows estimation of multiple slopes",
                        type = bool, default = False)      
    parser.add_argument('--freq', help = "Known cyclic frequencies", 
                        nargs ='+', type = float, default = []) 
    parser.add_argument('--save', help = "Save the OLS results as .xlsx", 
                        type = bool, default = False)    
    parser.add_argument('--display', help = "Display the spectrum ",
                        type = bool, default = False)     
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
    
    if args.save:
        import os
        ExcelPath = filepath[:-4]
        if not os.path.exists(ExcelPath): os.makedirs(ExcelPath)         
        Path = ExcelPath + '\OLSResults.xlsx' 
    else: Path = 0
    
    if args.Ufreq > 0 and args.Numfreq > 0:
        Omega = np.linspace(args.Lfreq, args.Ufreq, args.Numfreq)
    else: Omega = []  
    
    t0  = t - t[0]   # For the sake of computational efficiency
    
    tic = time.clock()     
    Results = LSSA(t0, f, P = P, Omega = Omega, ind = args.ind, level = args.level, 
                   trend = args.trend, slope = args.slope, freq = args.freq, saveto = Path)
    toc = time.clock()
    print("Computational Time: ", round(toc-tic,2), "s") 
    
    if args.display and len(Omega) == 0:
        print("Enter a set of cyclic frequencies, e.g., --Lfreq 0.8 --Ufreq 3.8 --Numfreq 16")
        sys.exit()        
    
    if args.display:
        fig, (ax1, ax2) = pyplot.subplots(2, 1, gridspec_kw={'hspace': 0.5})
        ax1.plot(t, f, 'blue', label = 'Time Series')
        ax1.plot(t, Results[3], 'red', label = 'Estimated Residual')    
        ax1.set(ylabel='Value')
        ax1.set(xlabel='Time')
        ax1.grid(True) 
        ax1.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 3, ncol= 4, 
                   mode = "expand", borderaxespad = 0)     
        ax2.plot(Omega,Results[0]*100,'blue', label = 'Least-Square Spectrum of Residual')
        ax2.plot(Omega,Results[1]*100*np.ones([len(Omega),1]),'red', label ='Critical Line')
        ax2.set(ylabel='Percentage variance')
        ax2.set(xlabel='Cyclic frequency')
        ax2.grid(True)   
        ax2.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 3, ncol= 4, 
                   mode = "expand", borderaxespad = 0) 
        fig.savefig("LSSA.png", dpi = 300, bbox_inches = 'tight')
        pyplot.show()

           

    '''
    # ----- Simulated Example ------------------------------------------------------------- 
    t = [k/23.0 for  k in range(65)]
    for kk in range(20):
        indt = int((65-kk)*np.random.rand(1)[0])
        t.remove(t[indt])
    Lt = len(t)
    t = np.array(t) 
    #t.shape = (Lt, 1)
    trend = 0.5 - 0.1*t
    season = 0.1*np.cos(2*np.pi*1*t) + 0.05*np.sin(2*np.pi*2.3*t)
    noise1 = 0.01*np.random.normal(0,1,Lt)
    noise2 = 0.03*np.random.normal(0,1,Lt)
    f = trend + season + noise1 - abs(noise2)
    P = 1/(0.01+(noise2)**2)
    Omega = [Om/5.0 for Om in range(1,20)]
    '''