function [spectrogram, stoch_surf, window_size, chat_cos, chat_sin, num_cons]=...
    LSWA(t, f, P_f, Omega, M, L1, L0, Q, IndexDatum, alpha, MorlCoeff, normalized)
%% The Least-Squares Wavelet Analysis (LSWA)
% Author: Ebrahim Ghaderpour www.ghader.org
% Copyright (c) 2019
%
% Please acknowledge the use of this software in any publications:
% "Ghaderpour, E., Pagiatakis, S.D., 2017. 
% Least-squares wavelet analysis of unequally spaced and non-stationary 
% time series and its applications. Mathematical Geosciences 49, 819–844."
%
% Note: In the following descriptions time and cyclic frequency 
%       can be replaced by distance and wavenumber, respectively
% Inputs:
% 1) t            - Time column vector of the series 
% 2) f            - Column vector of the time series values 
% 3) P_f          - The weight matrix (inverse of the covariance matrix 
%                   associated with the series if exists). 
%                   Note: if only standard deviations are provided, then 
%                         P_f is a column vector whose elements are the 
%                         inverse of variances of series values, and set
%                         P_f=1 if the series is equally weighted
% 4) Omega        - The vector of cyclic frequencies   
% 5) M            - The average sampling rate  
%                   Recommended:  M = length(t)/(t(length(t))-t(1))
% 6) L1           - The number of cycles of sinusoids being fitted to the 
%                   segments of the time series. Recommeded: L1=2
% 7) L0           - The number of additional samples to enlarge the 
%                   translating window size. Recommended: L0=20
% 8) Q            - The vector of form Q = [q1,q2,q3,q4,q5,...,qn] 
%                   to remove constituents of known forms. To find Q, let
%                   p(t)=m0+m1*t+m2*t^2+m3*t^3, where "m0" is intercept, 
%                   "m1" is slope, "m2" and "m3" are the quadratic and
%                   cubic coefficients, respectively
%                   Valid options for q1, q2, q3, and q4 are 0 or 1
%                   Set q1=1 to remove shift (m0) 
%                   Set q2=1 to remove linear trend (m1*t) 
%                   Set q3=1 to remove quadratic trends (m2*t^2) 
%                   Set q4=1 to remove cubic trends (m3*t^3) 
%                   from the segments of f
%                   Set q5 to qn, the frequenncies of the sinusoids to be 
%                   removed from segments of f, o.w., Q must have size 4
%                   Example1: Q=[0,0,0,0], no constituents of known forms
%                   Example2: Q=[1,1,0,0,5,7.5] for linear trend and 
%                             frequencies 5, 7.5
% 9) IndexDatum   - The indices for the starts of datum shifts 
%                   Note: To account for datum shifts, set q1=1 in Q
%                   Example1: IndexDatum=[] for no datum shifts
%                   Example2: IndexDatum=[101, 321] for two datum breaks 
%                             at indices 101 and 321
% 10) alpha       - The significance level, usually 0.01 or 0.05
% 11) MorlCoeff   - The Morlet coefficient. Recommended: MorlCoeff=0
%                   o.w., we recommend MorlCoeff=0.0125 and L1=6 cycles, 
%                   Morlet wavelet applied in the least-squares sense. 
%                   Note: the Gaussian values may not be attenuated to zero 
%                         at both ends of the windows for certain L1 and L0
% 12) normalized  - Set ture for the normalized spectrogram = Q_s/(Q_s+Q_n)
%                   Set false for spectrogram = Q_s (signals)
%                   Note: We only define the stochastic surface for 
%                         the normalized spectrogram
% 
% Outputs:
% 1) spectrogram  - The least-squares wavelet spectrogram (LSWS)
% 2) stoch_surf   - The stochastis surface at (1-alpha) confidence level
%                   NOTE: The stochastic is valid when normalized = true
% 3) window_size  - The sizes of translating windows
% 4) chat_cos     - The estimated cosine coefficients for segments 
% 5) chat_sin     - The estimated sine coefficients for segments
% 6) num_const    - The numbers of constituents of known forms within 
%                   translating windows (vary due to datum shifts)
%--------------------------------------------------------------------------
%% Check the input arguments
if nargin < 2
    error('Too few input arguments');
elseif nargin==2
    P_f=1; M=(length(t)-1)/(t(length(t))-t(1));
    Omega=1:0.5*M; L1=2; L0=20;
    Q=[1,0,0,0]; IndexDatum=[]; alpha= 0.01; 
    MorlCoeff=0; normalized=true;
elseif nargin==3
    M=(length(t)-1)/(t(length(t))-t(1));
    Omega=1:0.5*M; L1=2; L0=20;
    Q=[1,0,0,0]; IndexDatum=[]; alpha= 0.01; 
    MorlCoeff=0; normalized=true;
elseif nargin==4
    M=(length(t)-1)/(t(length(t))-t(1)); L1=2; L0=20;
    Q=[1,0,0,0]; IndexDatum=[]; alpha= 0.01; 
    MorlCoeff=0; normalized=true;
elseif nargin==5
    L1=2; L0=20; Q=[1,0,0,0]; IndexDatum=[]; alpha= 0.01; 
    MorlCoeff=0; normalized=true;
elseif nargin==6
    L0=20; Q=[1,0,0,0]; IndexDatum=[]; alpha= 0.01; 
    MorlCoeff=0; normalized=true;
elseif nargin==7
    Q=[1,0,0,0]; IndexDatum=[]; alpha= 0.01; 
    MorlCoeff=0; normalized=true;
elseif nargin==8
    IndexDatum=[]; alpha= 0.01; 
    MorlCoeff=0; normalized=true;
elseif nargin==9
    alpha= 0.01; MorlCoeff=0; normalized=true;
elseif nargin==10
    MorlCoeff=0; normalized=true;
elseif nargin==11
    normalized=true;
elseif nargin > 12
    error('Too many input arguments');
end
%--------------------------------------------------------------------------
Lt=length(t);
LOm=length(Omega);
LQ=length(Q);
% Sort out frequencies in Q
if LQ>4
   Q1=sort(Q(5:LQ));
   Q= [Q(1:4), Q1];
end
spectrogram=zeros(LOm,Lt); window_size=zeros(LOm,Lt);
num_cons=zeros(LOm,Lt);
stoch_surf=zeros(LOm,Lt); chat_cos=zeros(LOm,Lt); chat_sin=zeros(LOm,Lt);
[nrow,ncol]=size(P_f);
%--------------------------------------------------------------------------
if MorlCoeff~=0
    if ~isempty(IndexDatum)
        % For Morlet wavelet remove datum shifts from the entire series (if any)
        [~, ~, ~, f, ~, ~]=LSSA(t, f, P_f, [], Q(1:4), IndexDatum, alpha, normalized); 
    end
end
%--------------------------------------------------------------------------
for k=1:LOm      %loop over the cyclic frequencies
     L_k=L0+floor((L1*M)/Omega(k)); % L(omega_k) 
     if (mod(L_k,2)==0) % check if L(omega_k) is odd, o.w., add one to it
         L_k=L_k+1;
     end
     for j=1:Lt     %loop over the times
        if j >=(L_k+1)/2 && j<= Lt-(L_k-1)/2   % non-marginal windows
            lb=j-(L_k-1)/2; ub=j+(L_k-1)/2;
            tt = t(lb:ub);
            y = f(lb:ub);                      % segment
            Rj=L_k;                            % window size (segment size)
            if MorlCoeff~=0
                sub_ind_datum = [];
                MorlW=exp(-MorlCoeff*((2*pi*Omega(k))^2)*...
                    (tt-t(j)*ones(1,length(tt))').^2);
            else
                % Extract those datum shift indices within the window
                sub_ind_datum = intersect (lb:ub,IndexDatum)- lb+1;
                MorlW=1;
            end
            if ncol>1 && nrow>1 
                P_y = P_f(lb:ub,lb:ub).*MorlW;  % extract submatrix 
            elseif ncol==1 && nrow>1
                P_y = P_f(lb:ub).*MorlW;        % extract subvector
            else               
                P_y = MorlW;
            end
        elseif j < (L_k+1)/2 && j <= Lt-(L_k-1)/2  % left marginal windows
            ub=j+(L_k-1)/2;
            tt = t(1:ub);
            y = f(1:ub);
            Rj=ub;
            if MorlCoeff~=0
                sub_ind_datum = [];
                MorlW=exp(-MorlCoeff*((2*pi*Omega(k))^2)*...
                    (tt-t(j)*ones(1,length(tt))').^2);
            else
                sub_ind_datum = intersect (1:ub,IndexDatum);
                MorlW=1;
            end
            if ncol>1 && nrow>1
                P_y = P_f(1:ub,1:ub).*MorlW;
            elseif ncol==1 && nrow>1
                P_y = P_f(1:ub).*MorlW;
            else
                P_y=MorlW;
            end
        elseif j > Lt-(L_k-1)/2 && j>=(L_k+1)/2   % right marginal windows
            lb=j-(L_k-1)/2;
            tt = t(lb:Lt);
            y = f(lb:Lt);
            Rj=Lt-lb+1;
            if MorlCoeff~=0
                sub_ind_datum =[];
                MorlW=exp(-MorlCoeff*((2*pi*Omega(k))^2)*...
                    (tt-t(j)*ones(1,length(tt))').^2);
            else
                sub_ind_datum = intersect (lb:Lt,IndexDatum)- lb+1;
                MorlW=1;
            end
            if ncol>1 && nrow>1
                P_y = P_f(lb:Lt,lb:Lt).*MorlW;
            elseif ncol==1 && nrow>1
                P_y = P_f(lb:Lt).*MorlW;
            else
                P_y=MorlW;
            end
        else         % for entire series
            tt = t;
            y = f;
            Rj=Lt;
            if MorlCoeff~=0
                sub_ind_datum =[];
                MorlW=exp(-MorlCoeff*((2*pi*Omega(k))^2)*...
                    (tt-t(j)*ones(1,length(tt))').^2);
            else
                sub_ind_datum = IndexDatum;
                MorlW=1;
            end
            P_y = P_f.*MorlW;
        end  
        %----- Implement the LSSA to segments tt, y, P_y------------------
        [spectrum, chats, CritVal, ~, ~, ~, nc]=...
            LSSA(tt, y, P_y, Omega(k), Q, sub_ind_datum, alpha, normalized);  
        num_cons(k,j)=nc;  % number of constituents of known forms within the window
        chat_cos(k,j)=chats(1);    % estimated coefficient for cosine
        chat_sin(k,j)=chats(2);    % estimated coefficient for sine       
        spectrogram(k,j)=spectrum; % spectral peak corresponding to Omega(k)
        window_size(k,j)=Rj;       %window size for each pair (k,j) 
        stoch_surf(k,j)=CritVal;   % critical value 
     end
 end