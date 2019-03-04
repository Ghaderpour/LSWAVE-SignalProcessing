function [TrendsCoeff, SpecCosSin, g, Q, COVcoeff]=...
                               ALLSSA(t, f, P, Omega, Q, IndexDatum, alpha)
%% The Antileakage Least-Squares Spectral Analysis (ALLSSA)
% Author: Ebrahim Ghaderpour www.ghader.org
% Copyright (c) 2019 
%
% Please acknowledge the use of this software in any publications:
% "Ghaderpour, E., Liao, W., Lamoureux, M.P., 2018b. 
% Antileakage least-squares spectral analysis for seismic data 
% regularization and random noise attenuation. Geophysics 83, V157–V170."
% 
% Note: In the following descriptions time and cyclic frequency 
%       can be replaced by distance and wavenumber, respectively
% Inputs:
% 1) t           - Time column vector of the series 
% 2) f           - Column vector of the time series values 
% 3) P           - The weight matrix (inverse of the covariance matrix 
%                  associated with the series if exists) 
%                  Note: if only standard deviations are provided, then 
%                        P is a column vector whose elements are the 
%                        inverse of variances of series values, and set
%                        P=1 if the series is equally weighted
% 4) Omega       - The vector of cyclic frequencies  
% 5) Q           - The vector of form Q = [q1,q2,q3,q4,q5,...,qn] 
%                  to remove constituents of known forms. To find Q, let
%                  p(t)=m0+m1*t+m2*t^2+m3*t^3, where "m0" is intercept, 
%                  "m1" is slope, "m2" and "m3" are the quadratic and
%                  cubic coefficients, respectively
%                  Valid options for q1, q2, q3, and q4 are 0 or 1
%                  Set q1=1 to remove shift (m0) 
%                  Set q2=1 to remove linear trend (m1*t) 
%                  Set q3=1 to remove quadratic trends (m2*t^2) 
%                  Set q4=1 to remove cubic trends (m3*t^3) 
%                  from the entire f
%                  Set q5 to qn, the frequenncies of the sinusoids to be 
%                  removed from the entire f, o.w., Q must have size 4
%                  Example1: Q=[0,0,0,0], no constituents of known forms
%                  Example2: Q=[1,1,0,0,5,7.5] for linear trend and 
%                            frequencies 5, 7.5
% 6) IndexDatum  - The indices for the starts of datum shifts 
%                  Note: To account for datum shifts, set q1=1 in Q
%                  Example1: IndexDatum=[] for no datum shifts
%                  Example2: IndexDatum=[101, 321] for two datum breaks 
%                            at indices 101 and 321
% 7) alpha       - The significance level, usually 0.01 or 0.05
%
% Outputs:
%
% 1) TrendsCoeff - The estimated coefficients for datum shifts and trends
% 2) SpecCosSin  - The estimated cosine and sine coefficients of the 
%                  antileakage spectrum
% 3) g           - The residual series after simultaneously removing the 
%                  constituents of known forms
% 4) Q           - The updated Q that contains the estimated frequencies
% 5) COVcoeff    - The covariance matrix of TrendsCoeff and SpecCosSin
%--------------------------------------------------------------------------
%% Check the input arguments
if nargin < 2
    error('Too few input arguments');
elseif nargin==2
    P=1; Omega=1:0.5*(length(t)-1)/(t(length(t))-t(1));
    Q=[1,0,0,0]; IndexDatum=[]; alpha= 0.01; 
elseif nargin==3
    Omega=1:0.5*(length(t)-1)/(t(length(t))-t(1));
    Q=[1,0,0,0]; IndexDatum=[]; alpha= 0.01; 
elseif nargin==4
    Q=[1,0,0,0]; IndexDatum=[]; alpha= 0.01; 
elseif nargin==5
    IndexDatum=[]; alpha= 0.01; 
elseif nargin==6
    alpha= 0.01; 
elseif nargin > 7
    error('Too many input arguments');
end
%--------------------------------------------------------------------------
TrendsCoeff=zeros;
SpecCosSin=zeros;
Q=Q(1:4);
jj=1;
check=true;
g=f;
LOm=length(Omega);
% The following while loop estimates the cyclic frequencies (corresponding 
% to statistically significant peaks) after one round of iterations
while check && norm(g)>1e-6
    [spectrum, ~, ~, g, xhat, COVcoeff, ~]=...
        LSSA(t, f, P, Omega, Q, IndexDatum, alpha, true);
    Q1=Q;
    Om_inds=find(spectrum==max(spectrum(:)));
    Om_ind=Om_inds(1);   % the index of peak with maximum energy/power 
    if Om_ind>1 && Om_ind<LOm
        LB=Omega(Om_ind-1);
        UB=Omega(Om_ind+1);
    elseif Om_ind==1 && Om_ind<LOm
        LB=Omega(Om_ind);
        UB=Omega(Om_ind+1);
    elseif Om_ind>1 && Om_ind==LOm
        LB=Omega(Om_ind-1);
        UB=Omega(Om_ind);  
    else
        return
    end
    Om0=Omega(Om_ind); % selected frequency from the preselected set
    d1=Om0-LB;         % for the small neighbourhood around Om0
    d2=UB-Om0;         % for the small neighbourhood around Om0
    %----------------------------------------------------------------------
    % Find and remove a frequency from Q if it is close enough to Om0 
    loc=[];
    for bb=5:length(Q)
        if (0<=(Q(bb)-Om0) && (Q(bb)-Om0)<=d2/2) ||...
                (0<=(Om0-Q(bb)) && (Om0-Q(bb))<=d1/2)
            loc=bb;
        end
    end
    if ~isempty(loc) 
        if loc>4
            Q(loc)=[];
            if jj>1
                jj=jj-1;
            end
        end
    end
    %----------------------------------------------------------------------
    inc=(d1+d2)/200;                 % cyclic frequency increment
    Om_range=Om0-d1/2:inc:Om0+d2/2;  % small neighbourhood around Om0 
    [spectrum, ~, CritVal, ~, ~, ~, ~]=...
        LSSA(t, f, P, Om_range, Q, IndexDatum, alpha, true);
    Om_new_inds=find(spectrum==max(spectrum(:)));
    Om_new_ind=Om_new_inds(1);
    %-------check if the peak at Om_new is statistically significant ------
    if spectrum(Om_new_ind) > CritVal 
        Om_new=Om0-d1/2 -inc + inc*Om_new_ind;
        Q(4+jj)=Om_new;
        jj=jj+1;
        Q=[Q(1:4),sort(Q(5:length(Q)))]; % sort frequencies in Q
    end
    if isequal(Q,Q1)  % The termination condition
        check=false;
    end

end
%--------------------------------------------------------------------------
% separate the sinusoidal coefficients from trend coefficients (optional)
Lid=length(IndexDatum);
if Q(1)==1
    SpecCosSin=xhat(Lid+nnz(Q(1:4))+1:length(xhat));
    TrendsCoeff=xhat(1:Lid+nnz(Q(1:4)));
elseif nnz(Q(2:4))>1
    SpecCosSin=xhat(nnz(Q(2:4))+1:length(xhat));
    TrendsCoeff=xhat(1:nnz(Q(2:4)));
else
    SpecCosSin=xhat;
end
