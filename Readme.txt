Name of Software: Sequential Turning Point Detection (STPD)

Developed by: Dr. Ebrahim Ghaderpour 
Email:  ebrahim.ghaderpour@uniroma1.it

This package contains:
1) RunMe
2) Sequential Turning Point Detection (STPD)
3) Turning Point Trend Estimation (TPTR)

The main script for detecing trend turning points (TPs) in a time series is STPD. 

TPTR takes a time series t and f(t) and its estimated TPs to estimate the final 
linear trend that contains the connected linear pieces.

To begin, users may run "RunMe" which gets a time series and plots it along with 
the fitted linear trend with multiple connected pieces. It also saves the results.

Note: The input time series does not have to be equally spaced, but for STPD, it is 
suggested to resample the input time series on regular time intervals, e.g., monthly. 

Please pay attention to the format of input time series:
Example 1: A simulated example: SimulatedExample.csv that is Figure 3 of the following article.
Example 2: A PS-InSAR example: DESC-PSInSAR-MonthlyResampled.csv that is Figure 10-D of the 
following article.


Please cite the following article if you use this software for research purposes:

Reference: 
Ghaderpour, E.; Benedetta, A.; Bozzano, F.; Scarascia Mugnozza, G.; Mazzanti, P.  
A Fast and Robust Method for Detecting Trend Turning Points in InSAR Displacement 
Time Series, Computers & Geosciences, 2023.