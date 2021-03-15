JUST software package written by Ebrahim Ghaderpour
emails: ebrahim.ghaderpour@ucalgary.ca   
        ebrahimg2@gmail.com

There are several time series included in this package in .dat and one in .csv and one in .txt format
whose first, second, and third columns are the times, time series values, and weights.

----------------------------------------------------------------------------------------------------------------------

In Python: 

The Python code can be run individually and the examples with their descriptions 
are given in the scripts (after if __name__ == '__main__':)

The supported format for the data sets are .csv, .dat, and .txt
Each module imports argparse and can be easily run from the Command Prompt (Fully commented)

Use two consecutive hyphens -- for each input parameter

Example: Try 

D:\JUST>python JUSTdecompose.py --data "D:\JUST\EVIA_CA.csv" --ind 13 111 144 --Numtime 300 --display True --save True

This should display the JUST decomposition results and save the results in a new folder 
called EVIA_CA in the same directory as the input data 

Example: Try

D:\JUST>python ALLSSA.py --datafile "D:\JUST\Sim_Geophysics.dat" --decimal 2 --save True

This will print:

The ALLSSA estimated cyclic frequencies: [4.07, 20.37, 22.29]

And will save the ALLSSA results as .xlsx in a new folder called Sim_Geophysics in the same directory as the input data 

Example:

The input frequencies can be modified by --Lfreq, --Ufreq, and --Numfreq

The input times for regularization in LSWA.py can also be modified by --Ltime, --Utime, and --Numtime

Try:

D:\JUST>python LSWA.py --data "D:\JUST\EVIA_CA.txt" --ind 13 111 144 --Lfreq 0.2 --Ufreq 4.1 --Numfreq 40 --Numtime 300 --display 3

This should display the amplitude spectrogram. The cyclic frequencies are in cycles per year (c/y).
Users may also suppress the annual peaks in the spectrogram (1 c/y) by adding --freq 1 in the command above.
In cases where the matrix of normal equations becomes singular (singularity error), users may change/alter the window size parameters
i.e., by adding --L0 5 to the command as an example 

The full description of parameters is available from the Command Prompt using pydoc:

Example:

D:\JUST>python -m pydoc JUSTmonitor

-------------------------------------------------------------------------------------------

Sources of the examples included in this package:

Sim_JumpInd_48_104.dat, ExA_NDVI_Australia.dat, ExB_NDVI_Australia.dat, 
EVIA_CA.dat, EVIB_CA.dat, and EVIC_CA.dat are in https://doi.org/10.3390/rs12234001


For near-real-time monitoring, the following three time series are shown in https://doi.org/10.3390/rs12152446

 
You may run JUSTmonitor by inputting the values for 'history' and 'start' as listed below:
  Time series            EVIA_AB.dat    EVIB_AB.dat    EVIC_AB.dat 
    history                  17             30             15
     start                   72             71             71



sim_Geophysics.dat is the same simulated example in 
https://doi.org/10.1190/geo2017-0284.1
True wavenumbers are 4.0744, 20.3718, 22.2818 achieved by setting 'decimal' to 4 in ALLSSA 


Sim_JumpInd_101_321.dat is the same simulated example in 
https://doi.org/10.1007/s10291-019-0841-3

-------------------------------------------------------------------------------------------

To compute the cross-spectrogam for two time series: (t1, f1, P1) and (t2, f2, P2)
Users may select a common time vector tt for both time series and run LSWA for 
each time series by inputting 'tt' = tt to compute a normalized spectrogram. 
Then multiply the normalized spectrograms entrywise to obtain the cross-spectrogram!
(Reference: https://doi.org/10.1007/s00190-018-1156-9)

Example:

A common time vector tt (equally spaced) for both time series may be obtained by:
a = max (t1(1), t2(1));
b = min (t1(length(t1)), t2(length(t2)));
tt = linspace(a, b, num);  

-------------------------------------------------------------------------------------------

While the software has been tested and debugged using thousands of time series, 
users may still see some bugs for certain time series and inputs. 
I appreciate your feedback that can be posted either on Github or emailed to me. 

Thank you!
