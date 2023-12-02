function TPs = STPD(t,f,varargin)
%% Sequential Turning Point Detection (STPD)
%{
function TPs = STPD(t,f,varargin)

Summary: This function takes the time series t and f(t) and finds the
         potential linear trend turning points (TPs).
         It calls the TPD function for each segment within a translating 
         window to find a TP along with its direction, i.e.,
         the slope of the linear piece after the TP minus the one before
         the TP. The TPs are selected based on the following criteria:
            1) Signal-to-noise ratio (SNR) 
            2) Normalized difference residual index (NDRI)
            3) Confidence level
            4) TP Direction 
            5) Minimum distance to the window locations

Inputs:
1) t           - Vector of order n. Times of observations or measurements 
2) f           - Vector of order n. The time series values
3) 'size'      - Numeric. The window or segment size (R)
                 Default is five times the average sampling rate (5M)
4) 'step'      - Numeric. The translation step should be less than size 
                 Default is the average sampling rate (M)
5) 'SNR'       - Numeric. Threshold for the estimated SNR
                 Default is 1
6) 'NDRI'      - Numeric. Threshold for |NDRI| to seperate TPs from jumps
                 Default is 0.3
7) 'alpha'     - Numeric. The significance level, usually 0.01 or 0.05
                 Default is 0.01. 
8) 'dir_th'    - Numeric. Threshold for the abs of TP direction. 
                 Default is 0 (|Dir| > 0)
9) 'tp_th'     - Numeric. The minimum time difference between TPs
                 Default is 1 (TPs are at least one year apart)
10)'margin'    - The size of marginal segment where no TP is estimated.
                 Default is 1

Outputs:
TPs            - The indices of the potential TPs

Reference: 
Ghaderpour, E.; Benedetta, A.; Bozzano, F.; Scarascia Mugnozza, G.; 
Mazzanti, P.  A Fast and Robust Method for Detecting Trend Turning Points 
in InSAR Displacement Time Series, Computers & Geosciences, 2023.
--------------------------------------------------------------------------
Author: Dr. Ebrahim Ghaderpour 
Email:  ebrahim.ghaderpour@uniroma1.it
Copyright (c) 2023
%}
%==========================================================================
%% Check the input arguments
%-------------------------------------------------------------------------
Lt = length(t);
%--------- Window Size ----------------------------------------------------
size_ind = find(strcmpi(varargin,'size'));
M = floor(Lt/(t(Lt)-t(1)));   % A simple estimate for sampling rate
if ~isempty(size_ind)
    R = varargin{size_ind+1};
    if R > Lt || R < 3
        error('The window size must be an integer between 3 and n')
    end
else
    R = 5*M;
end
%--------- Translation Step -----------------------------------------------
margin = find(strcmpi(varargin,'step'));
if ~isempty(margin)
    delta = varargin{margin+1};
    if delta >= R || delta < 1
        error('Translation step must be a positive integer less than size')
    end
else
    delta = M;
end
%--------- Signal-to-Noise Ratio  -----------------------------------------
SNR = find(strcmpi(varargin,'SNR'));
if ~isempty(SNR)
    SNR = varargin{SNR+1};
    if SNR <= 0; error('SNR must be positive'); end
else
    SNR = 1;
end
%--------- Normlized Difference Residual Index (NDRI)  --------------------
NDRI = find(strcmpi(varargin,'NDRI'));
if ~isempty(NDRI)
    NDRI = varargin{NDRI+1};
    if NDRI > 1 || NDRI < 0
        error('Threshold for |NDRI| must be in [0,1]')
    end
else
    NDRI = 0.3;
end
%--------- Significance Level ---------------------------------------------
alpha = find(strcmpi(varargin,'alpha'));
if ~isempty(alpha)
    alpha = varargin{alpha+1};
    if alpha > 1 || alpha < 0
        error('The significance level must be between 0 and 1')
    end
else
    alpha = 0.01;
end
%--------- TP Direction Threshold -----------------------------------------
dir_ind = find(strcmpi(varargin,'dir_th'));
if ~isempty(dir_ind); DIR = varargin{dir_ind+1};
else; DIR = 0; end
%--------- Minimum Time Difference Between Potential TPs ------------------
tp_th = find(strcmpi(varargin,'tp_th'));
if ~isempty(tp_th); tp_th = varargin{tp_th+1};
else; tp_th = 1; end
%--------- The size of marginal segment where no TP is estimated  ---------
bnd = find(strcmpi(varargin,'margin'));
if ~isempty(bnd) 
    bnd = varargin{bnd+1};
    if bnd >= R/2 || bnd < 1
        error('margin should be an integer between 0 and half of the size')
    end
else
    bnd = 1; 
end
%--------------------------------------------------------------------------                                     
R1 = R - 1;
inda = 1; indb = 1 + R1;
TIndDis = zeros(1,2);
check = true;  nr = 1;
while indb <= Lt + R1 && check 
    if indb >= Lt  % The condition for the last segment to also have size R
        inda = Lt - R1; indb = Lt; check = false;
    end
    t0=t(inda:indb);
    f0=f(inda:indb);
    % Call function TPD to find a potential TP for each segment:
    [TInd, SNR0, pval, Diff, Dir0, NDRI0, ~] = TPD(t0, f0, bnd);
    % Apply the input thresholds to the estimated TPs:
    if Diff<4 && SNR0>SNR && pval<alpha && abs(Dir0)>DIR && NDRI0<NDRI
        TIndDis(nr,:) = [TInd+inda-1, abs(TInd-1-R1/2)];
        nr = nr + 1; 
    end
    inda = inda + delta; indb = inda + R1;
end
TIndDis = sortrows(TIndDis, 1);
if TIndDis(1,1)==0
    TPs=[];
else
    r = 2;
    while r < nr
        if t(TIndDis(r,1))- t(TIndDis(r-1,1)) < tp_th
            group = TIndDis(r-1:r,:);
            group = sortrows(group,2);
            TIndDis (r-1,:) = group(1,:);
            TIndDis (r,:) = [];
            nr = nr - 1;
        else
            r = r + 1;
        end
    end
    TPs = TIndDis(:,1)';
end
end
%==========================================================================
function [TP, SNR, pval, Diff, DIR, NDRI, y] = TPD(t, f, bnd)
%% Turning Point Detection (TPD)
%{
function [TP, SNR, pval, Diff, DIR, NDRI, y] = TPD(t, f, bnd)

Summary: This function gets a segment t and f(t)and finds a potential TP.
         It calls the TPDetect function for each of forward and backward
         TP estimations. Then it selects the TP whose residual norm is
         smaller and returns the TP along with its statistics.

Inputs:
1) t           - Vector of order n. Times of observations or measurements 
2) f           - Vector of order n. The time series values
3) bnd         - The size of marginal segment where no TP is estimated.

Outputs:
1) TP          - The potential TP for segment t and f(t)
2) SNR         - The estimated SNR corresponding to the TP
3) pval        - The p-value corresponding to the TP
4) Diff        - The difference between the TPs of forward and backward
5) DIR         - The TP direction 
6) SNR         - The absolute value of NDRI corresponding to the TP
7) y           - The residual segment 
%}
%==========================================================================
Lt = length(t);
% The forward estimation of TP 
t = t - t(1);
TPStats = TPDetect(t, f, bnd);
TInd = TPStats(1);
pval = pvalue(t,f,TPStats(1),TPStats(2),TPStats(3),TPStats(4));
% The backward estimation of TP
t0 = flip(t); f0 = flip(f);
t0 = t0(1)-t0;
TPStats0 = TPDetect(t0, f0, bnd);
TInd0 = Lt-TPStats0(1)+1;
pval0 = pvalue(t0,f0,TPStats0(1),TPStats0(2),TPStats0(3),TPStats0(4));
% The difference between the TPs of forward and backward
Diff = abs(TInd0-TInd);
y = zeros(Lt,1); 
% From the estimated TPs for forward and backward estimations, choose 
% the TP and its statistics whose residual norm is smaller:
if TPStats(5) <= TPStats0(5)
    TP = TPStats(1); intercept1 = TPStats(2); 
    slope1 = TPStats(3); slope2 = TPStats(4);
    g1 = intercept1+slope1*t(1:TP);
    g2 = slope2*(t(TP:Lt)-t(TP))+intercept1+slope1*t(TP);
    y(1:TP) = g1; y(TP:Lt) = g2;
    SNR = norm(y)/norm(f-y); 
    nr1 = norm(f(1:TP)-g1); nr2 = norm(f(TP:Lt)-g2); 
    NDRI = abs(nr1-nr2)/(nr1+nr2);
else
    TP = TPStats0(1); intercept1 = TPStats0(2); 
    slope1 = TPStats0(3); slope2 = TPStats0(4);
    g1 = intercept1+slope1*t0(1:TP);
    g2 = slope2*(t0(TP:Lt)-t0(TP))+intercept1+slope1*t0(TP);
    y(1:TP) = g1; y(TP:Lt) = g2;
    SNR = norm(y)/norm(f0-y); 
    nr1 = norm(f0(1:TP)-g1); nr2 = norm(f0(TP:Lt)-g2); 
    NDRI = abs(nr1-nr2)/(nr1+nr2);
    pval = pval0;
    y = flip(y);
    TP = Lt-TP+1;
end
DIR = slope2-slope1;
end
%==========================================================================
function TurnStats = TPDetect(t, f, bnd)
%% Turning Point Detection (TPD)
%{
function TurnStats = TPDetect(t, f, bnd)

Summary: This function gets a segment t and f(t) and finds a potential
         TP by the ordinary least-squares method (OLS).

Inputs:
1) t           - Vector of order n. Times of observations or measurements 
2) f           - Vector of order n. The time series values
3) bnd         - The size of marginal segment where no TP is estimated.

Outputs:
TurnStats      - An array consisting of TP, intercept and slope of 
                 the first linear piece and the slope of the second
                 linear piece and L2 norm of the residual segment.
%}
%==========================================================================
Lt = length(t);
Attributes = zeros(Lt-2*bnd,5); 
k = 1; r = f;
for ind = bnd+1:Lt-bnd  % From bnd+1
    % Use OLS to estimate the first linear piece from the beginning of the
    % segment to ind:
    D1=zeros(ind,2); D1(:,1)=ones(ind,1); D1(:,2)=t(1:ind); y1=f(1:ind);
    chat = (D1'*D1)\D1'*y1;
    r(1:ind) = y1-D1*chat;
    % Use OLS to estimate the second linear piece from the end of the
    % first estimated linear piece to the end of the segment:
    D2=t(ind:Lt)-t(ind); y2=f(ind:Lt)-chat(1)-chat(2)*t(ind);
    slope = (D2'*D2)\D2'*y2;
    r(ind:Lt) = y2-D2*slope;
    %----------------
    Attributes(k,:) = [ind, chat(1), chat(2), slope, norm(r)];
    k = k + 1;
end
Attributes = sortrows(Attributes, 5);
TurnStats = Attributes(1,:); % Choose the TP whose residual norm is minimum
end

function pval = pvalue(t, f, ind, intercept1, slope1, slope2)
%% p-value corresponding to TP (t-test)
%{
function pval = pvalue(t, f, ind, intercept1, slope1, slope2)

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
%}
%==========================================================================
Lt = length(t);
% For the first linear piece:
D1(:,1) = ones(ind,1);
D1(:,2) = t(1:ind);
y1 = f(1:ind);
g1 = y1-D1*[intercept1, slope1]';
% Two degrees of freedom are lost due to estimation of intercept and slope:
SE1 = sqrt((g1'*g1)/(ind-2))/norm(D1(:,2)-mean(D1(:,2))); 
% For the second linear piece:
D2 = t(ind:Lt)-t(ind);
y2 = f(ind:Lt)-intercept1-slope1*t(ind);
g2 = y2-D2*slope2;
SE2 = sqrt((g2'*g2)/(Lt-ind))/norm(D2-mean(D2)); % One degree of freedom is lost
ttest = abs(slope1-slope2)/sqrt(SE1^2+SE2^2);
pval = 2*(1-tcdf(ttest,Lt-2));
end
