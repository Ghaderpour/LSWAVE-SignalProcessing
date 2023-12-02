function [stats, y] = TPTR(t, f, TPs)
%% Turning Point Trend Estimation (TPTR)
%{
function [stats, y] = TPTR(t, f, TPs)

Summary: This function gets a time series t and f(t) and its estimated 
         turning points (TPs) to estimate the final linear trend pieces.
         The forward and backward linear trends are estimated and the one
         whose residual series is mimimum will be selected.  
Inputs:
1) t           - Vector of order n. Times of observations or measurements 
2) f           - Vector of order n. The time series values
3) TPs         - The potential TPs estimated by STPD

Outputs:
1) stats       - A matrix whose rows contain TP, the new estimated slope 
                 of the second linear piece after the TP, the new 
                 estimated TP direction, the new estimated 
                 signal-to-noise ratio (SNR), and the new estimated 
                 normalized difference residual index (NDRI)for the TP
2) y           - The fitted linear trend whose turning points are the TPs
--------------------------------------------------------------------------
Author: Dr. Ebrahim Ghaderpour 
Email:  ebrahim.ghaderpour@uniroma1.it
Copyright (c) 2023
%}
%==========================================================================
Lt = length(t);
% The forward estimation of linear trend given TPs:
[slopes1, y1] = TPTRF(t,f,TPs);
% The backward estimation of linear trend given TPs:
t0 = flip(t); f0 = flip(f);
t0 = t0(1)-t0;
TPs0 = 1+Lt-TPs;
[slopes0, y0] = TPTRF(t0,f0,TPs0);
slopes2 = -flip(slopes0);
y2 = flip(y0);
% Choose the linear trend (the forward or backward) whose norm is minimum: 
if norm(f-y1)<norm(f0-y0)
    slopes = slopes1; y = y1;
else
    slopes = slopes2; y = y2;
end
LTPs = length(TPs);
TPs = [1,TPs,Lt];
nr = norm(f-y); nsig = norm(y); SNR = nsig/nr;
stats = [0,slopes(1),0,SNR,0];
for k=2:LTPs+1
    r1 = f(TPs(k-1):TPs(k))-y(TPs(k-1):TPs(k)); 
    r2 = f(TPs(k):TPs(k+1))-y(TPs(k):TPs(k+1));
    r = f(TPs(k-1):TPs(k+1))-y(TPs(k-1):TPs(k+1));
    nr1 = norm(r1); nr2 = norm(r2); 
    nr = norm(r); nsig = norm(y(TPs(k-1):TPs(k+1))); SNR = nsig/nr;
    NDRI = abs(nr2-nr1)/(nr1+nr2);
    stats(k-1,:) = [TPs(k),slopes(k),slopes(k)-slopes(k-1),SNR,NDRI];
end
end    
%==========================================================================
function [slopes, y] = TPTRF(t,f,TPs)
%% Turning Point Trend Forward Estimation (TPTRF)
%{
function [slopes, y] = TPTRF(t,f,TPs)

Summary: This function gets a time series t and f(t) and its estimated 
         turning points (TPs) to estimate the final linear trend pieces
         using the ordinary least-squares method (OLS).
Inputs:
1) t           - Vector of order n. Times of observations or measurements 
2) f           - Vector of order n. The time series values
3) TPs         - The potential TPs estimated by STPD

Outputs:
1) slopes      - The estimated slopes of linear pieces after TPs
2) y           - The fitted linear trend whose turning points are the TPs
%}
%--------------------------------------------------------------------------
LTPs = length(TPs);
Lt = length(t);
TPs = sort([TPs, Lt]);
slopes = zeros(LTPs+1,1);
y = zeros(Lt,1);
if LTPs == 0 
    D1 = zeros(Lt,2);
    D1(:,1) = ones(Lt,1);
    D1(:,2) = t(1:Lt);
    chat = (D1'*D1)\D1'*f;
    y = D1*chat;
    slopes = chat(2);
else
    TInd = TPs(1);
    D1 = zeros(TInd,2);
    D1(:,1) = ones(TInd,1);
    D1(:,2) = t(1:TInd);
    g = f(1:TInd);
    chat = (D1'*D1)\D1'*g;
    y(1:TInd) = D1*chat;
    slopes(1) = chat(2);
    for k = 1:LTPs
        D2 = t(TPs(k):TPs(k+1))-t(TPs(k));
        y2 = f(TPs(k):TPs(k+1))-y(TPs(k));
        slope = (D2'*D2)\D2'*y2;
        slopes(k+1) = slope;
        y(TPs(k):TPs(k+1)) = D2*slope+y(TPs(k));
    end
end
end