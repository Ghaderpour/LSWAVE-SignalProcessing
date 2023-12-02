%{ 
This script runs the Sequential Turning Point Detection (STPD) on a 
displacement InSAR time series and plot it along with its fitted linear 
trend with multiple connected pieces. Note that STPD also works on 
unequally spaced time series, but it is recommended to resample the 
time series, e.g., to monthly intervals 

Reference: 
Ghaderpour, E.; Benedetta, A.; Bozzano, F.; Scarascia Mugnozza, G.; 
Mazzanti, P.  A Fast and Robust Method for Detecting Trend Turning Points 
in InSAR Displacement Time Series, Computers & Geosciences, 2023.

Author: Dr. Ebrahim Ghaderpour 
Email:  ebrahim.ghaderpour@uniroma1.it
Copyright (c) 2023
%}
%TS = readmatrix("D:\STPD_MATLABPackage_EGhaderpour\SimulatedExample.csv");
TS = readmatrix("D:\STPD_MATLABPackage_EGhaderpour\DESC-PSInSAR-MonthlyResampled.csv");
t = TS(:,1);
f = TS(:,2);
TPs = STPD(t,f,'size',60,'step',12,'SNR',1,'NDRI',0.3,'margin',12,'alpha',0.01);
[stats, y] = TPTR(t,f,TPs);
%==========================================================================
%-- The following command lines save the TP results in a .csv -------------
% Note that slope, direction, SNR,|NDRI| may differ from the ones estimated 
% for each window because these values correspond to the entire time series 
Tb = array2table(round(stats,2));
Tb.Properties.VariableNames(1:5)={'TP','slope','direction','SNR','|NDRI|'};
writetable(Tb,'Result.csv')  
%==========================================================================
%--The following command lines just nicely display and save the results ---
figure
hold on
plot(t,f,'-ok','LineWidth',0.5,'MarkerFaceColor','k','MarkerSize',5)
plot(t,y,'-b','LineWidth',1)
plot(t(TPs),y(TPs),'db','MarkerFaceColor','b','MarkerSize',5)
legend('Time series','Linear trend','Turning point','Location','southwest') 
for k=1:length(TPs)
    A = "SNR = "+num2str(round(stats(k,4),2)); 
    B = "NDRI = "+num2str(round(stats(k,5),2));
    text(t(TPs(k)), f(TPs(k)), [A;B], 'EdgeColor', 'b')
end
hold off
grid on
ylabel("Displacement (mm)")
xlabel("Time (year)")
% This is to re-write the x-axis value to values you like nicely:
set(gca, 'XTick', 0:7, ...
   'XTickLabel', {'2015','2016','2017','2018','2019','2020','2021','2022'}) 
set(gcf, 'Color', 'w');     % To set background white
set(gca, 'fontsize' ,14) ;  % To adjust the font size
width = 8;                  % To set the width of figure
height = 4;                 % To set the height of figure
set(gcf, 'PaperUnits', 'inches','PaperPosition',[0 0 width height]);       
print( '-dtiff',  'AnExample', '-r300');                                                      