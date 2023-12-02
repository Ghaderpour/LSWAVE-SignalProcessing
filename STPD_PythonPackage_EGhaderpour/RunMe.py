"""
This script calls the Sequential Turning Point Detection (STPD) function
and nicely plots a given time series with its linear trend results which may have
several connected linear pieces.

Note: The input time series does not have to be equally spaced, but for STPD,
it is suggested to resample the input time series on regular time intervals, e.g., monthly

Reference: 
Ghaderpour, E.; Benedetta, A.; Bozzano, F.; Scarascia Mugnozza, G.; Mazzanti, P.
A Fast and Robust Method for Detecting Trend Turning Points in InSAR Displacement
Time Series, Computers & Geosciences, 2023.

Author: Dr. Ebrahim Ghaderpour 
Email:  ebrahim.ghaderpour@uniroma1.it
Copyright (c) 2023
"""
#==================================================================================================
import csv
import numpy as np
from matplotlib import pyplot as plt
from STPD import STPD
from TPTR import TPTR
#--------------------------------------------------------------------------------------------------
if __name__== "__main__":
    # Enter the csv path containing two columns time and series values:
    #csvpath = "D:\STPD_PythonPackage_EGhaderpour\SimulatedExample.csv"    
    csvpath = "D:\STPD_PythonPackage_EGhaderpour\DESC-PSInSAR-MonthlyResampled.csv"
    with open(csvpath, 'r') as csvfile:
        reader = csv.reader(csvfile, skipinitialspace=True)
        next(reader)
        t=[]; f=[]
        for row in reader: t.append(float(row[0])); f.append(float(row[1]))
    t = np.array(t); f = np.array(f)    
    TPs = STPD(t,f, size=60, step=12, SNR=1, NDRI=0.3, dir_th=0, tp_th=1, margin=12, alpha=0.01)
    stats, y = TPTR(t, f, TPs)
    #==============================================================================================
    # The following command lines save the TP results in a .txt
    # Note that slope, direction, SNR, |NDRI| may differ from the ones estimated within each
    # window because these values correspond to the entire time series    
    log = open('Result.txt', "w")  
    log.write("TP         slope        direction         SNR           |NDRI|")
    [log.write("\n{a1:<6}{a2:>10}{a3:>15}{a4:>15}{a5:>15}".format(
        a1 = int(stats[k][0]), a2 = stats[k][1], a3 = stats[k][2], a4 = stats[k][3], 
        a5 = stats[k][4])) for k in range(len(TPs))]  
    log.close()    
    #==============================================================================================
    # The following command lines just nicely display and save the results 
    plt.plot(t, f, '-ok', label = 'Time Series', linewidth=1, markersize=3)
    plt.plot(t, y, 'b', label = 'Linear trend')
    if len(TPs) > 0:
        plt.plot(t[np.array(TPs)], y[np.array(TPs)], 'db', label = 'Turning point')
    for k in range(len(TPs)):
        C = "SNR = "+ str(np.round(stats[k][3],2)) + "\nNDRI = "+ str(np.round(stats[k][4],2))
        plt.text(t[TPs[k]], f[TPs[k]], C, bbox=dict(facecolor='b', alpha=0.1))   
        
    plt.xlabel('Time (year)')    
    plt.ylabel('Displacement (mm)')
    plt.grid(True)
    # This is to re-write the x-axis value to values you like nicely:
    plt.xticks([0, 1, 2, 3, 4, 5, 6, 7], ['2015','2016','2017','2018','2019','2020','2021','2022'])
    plt.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 3, ncol= 4, 
               mode = "expand", borderaxespad = 0)
    fig = plt.gcf()
    fig.set_size_inches(8, 4)    
    #----------------------------------------------------------------------------------------------
    # You may change the directory here to save the results
    plt.savefig("AnExample.png",dpi=300, bbox_inches='tight')     
    plt.show()
   

