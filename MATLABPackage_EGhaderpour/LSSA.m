function [spectrum, chats, CritVal, g, xhat, COVcoeff, nc]=...
                     LSSA(t, f, P, Omega, Q, IndexDatum, alpha, normalized)
%% The Least-Squares Spectral Analysis (LSSA)
% Author: Ebrahim Ghaderpour www.ghader.org
% Copyright (c) 2019
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
% 8) normalized  - Set ture for the normalized spectrum = Q_s/(Q_s+Q_n)
%                  Set false for spectrum = Q_s (signals)
%                  Note: We only define the critical value for 
%                        the normalized spectrum
%
% Outputs:
%
% 1) spectrum    - The least-squares spectrum (LSS)
% 2) chats       - The simultaneously estimated coefficients of sinusoids
% 3) CritVal     - The critical value at (1-alpha) confidence level 
%                  NOTE: The critical value is valid if normalized = true
% 4) g           - The residual series after simultaneously removing the 
%                  constituents of known forms
% 5) xhat        - The estimated coefficients of constituents of known forms
% 6) COVcoeff    - The estimated covariance matrix of xhat
% 7) nc          - The number of constituents of known forms
%--------------------------------------------------------------------------
%% Check the input arguments
if nargin < 2
    error('Too few input arguments');
elseif nargin==2
    P=1; Omega=1:0.5*(length(t)-1)/(t(length(t))-t(1));
    Q=[1,0,0,0]; IndexDatum=[]; alpha= 0.01; normalized= true;
elseif nargin==3
    Omega=1:0.5*(length(t)-1)/(t(length(t))-t(1));
    Q=[1,0,0,0]; IndexDatum=[]; alpha= 0.01; normalized= true;
elseif nargin==4
    Q=[1,0,0,0]; IndexDatum=[]; alpha= 0.01; normalized= true;
elseif nargin==5
    IndexDatum=[]; alpha= 0.01; normalized= true;
elseif nargin==6
    alpha= 0.01; normalized= true;
elseif nargin==7
    normalized= true;
elseif nargin > 8
    error('Too many input arguments');
end
%--------------------------------------------------------------------------
spectrum=zeros; xhat=zeros; chat=zeros;                 
Lt=length(t);
LOm=length(Omega);
chats=zeros(1,2*LOm);
IndexDatum=sort(IndexDatum);
Lid=length(IndexDatum);
if Lid>0
    if IndexDatum(1)==1
        % the first datum break index cannot be 1 
        IndexDatum(1)=[];
        Lid = Lid -1;
    end
end
% sort Q
LQ=length(Q);
if LQ>4
   Q1=sort(Q(5:LQ));
   Q= [Q(1:4), Q1];
end
% -------------------------------------------------------------------------
% Calculate the number of constituents of known forms (nc)
nc=0; % counter for the number of constituents of known forms
for i=1:LQ
    if Q(i)~=0 && i<5
        nc=nc+1;
    elseif Q(i)~=0 && i>=5
        nc=nc+2;
    end
end
if Q(1)~=0 && Lid>0
    nc=nc+Lid;
end
% -------------------------------------------------------------------------
% Make an array "QQ" for removing the sinusoids entered by users:
if LQ>4
    QQ=Q(5:LQ);
else
    QQ=[];
end
% -------------------------------------------------------------------------
% Create the design matrix B based on the constituents of known forms 
Phi_=zeros(Lt,nc);
q=1;
if Q(1)==1 && Lid==0                           
    Phi_(:,q)=ones(1,Lt);  % constituent of known form for the intercept
    q=q+1;
elseif Q(1)==1 && Lid>0
    % constituents of known forms for datum shifts
    ind=(1:IndexDatum(1)-1)';
    Phi_(ind,q)=ones(length(ind),1);
    q=q+1;
    for j=1:Lid
        if j<Lid
            ind=(IndexDatum(j):IndexDatum(j+1)-1)';
        else
            ind=(IndexDatum(j):Lt)';
        end
        Phi_(ind,q)=ones(length(ind),1);
        q=q+1;
    end
end
if Q(2)==1   %constituent of known form for linear trend (slope)
 Phi_(:,q) = t;
 q=q+1; 
end
if Q(3)==1   %constituent of known form for quadratic trend
 Phi_(:,q) = t.^2;
 q=q+1;
end
if Q(4)==1    %constituent of known form for cubic trend
 Phi_(:,q) = t.^3;
 q=q+1;
end
for qq=5:LQ   %constituent of known forms for sinusoids
  Phi_(:,q) = cos(2*pi*Q(qq)*t);
  q=q+1;
  Phi_(:,q) = sin(2*pi*Q(qq)*t); 
  q=q+1;
end
% -------------------------------------------------------------------------
% Estimate the coefficients (xhat) for the constituents of known forms
[nrow,ncol]=size(P);
if nc>0
    if ncol==1 && nrow==1 % If the time series is equally weighted
        Ninv=inv(Phi_'*Phi_);
        xhat=Ninv*(Phi_'*f);
    elseif nrow>1 && ncol==1 % If P is a column vector
        Ninv=inv(Phi_'*(Phi_.* kron(P,ones(1,nc))));
        xhat=Ninv*(Phi_'*(f.*P));
    elseif nrow>1 && ncol>1 % If P is a square matrix
        Ninv=inv(Phi_'*P*Phi_);
        xhat=Ninv*(Phi_'*P*f); 
    end
end
% -------------------------------------------------------------------------
% Find the residual series: after removing the constituents of known forms
if (nc==0)
   g=f;
else
   g= f- Phi_*xhat; % Residual time series
end
% -------------------------------------------------------------------------
if nrow==1 && ncol==1    % If the time series is equally weighted
    gP=g;
elseif nrow>1 && ncol==1 % If P is a column vector
    gP=g.*P;
elseif nrow>1 && ncol>1  % If P is a square matrix     
    gP=P*g;    
end
norm_g=g'*gP;
% -------------------------------------------------------------------------
% Estimate the covariance matrix of xhat
if nc > 0
    COVcoeff=(norm_g/(Lt-nc))*Ninv;
else
    COVcoeff = [];
end
% -------------------------------------------------------------------------
% Calculate the least-squares spectrum (LSS)
for k=1:LOm
      if ~any(Omega(k)==QQ)
          Phi(:,1) = cos(2*pi*Omega(k)*t);
          Phi(:,2) = sin(2*pi*Omega(k)*t);
          if nrow==1 && ncol==1 % If the series is equally weighted
              if nc>0
                    N12=Phi_'*Phi;
                    N22=Phi'*Phi;
                    chat=(N22-N12'*Ninv*N12)\(Phi'*gP);
              else
                    N22=Phi'*Phi;
                    chat=N22\(Phi'*gP);
              end
          elseif nrow>1 && ncol==1 % If P is a column vector
              if nc>0
                    N12=Phi_'*(Phi.*[P P]);
                    N22=Phi'*(Phi.*[P P]);
                    chat=(N22-N12'*Ninv*N12)\(Phi'*gP);
              else
                    N22=Phi'*(Phi.*[P P]);
                    chat=N22\(Phi'*gP);
              end
          elseif nrow>1 && ncol>1 % If P is a square matrix
              if nc>0
                    N12=Phi_'*P*Phi;
                    N22=Phi'*P*Phi;
                    chat=(N22-N12'*Ninv*N12)\(Phi'*gP);
              else
                    N22=Phi'*P*Phi;
                    chat=N22\(Phi'*gP);
              end
          end
          spec=gP'*Phi*chat;    % This is Q_s (signal): Not normalized
          if normalized
              spec_k=spec/norm_g;
              if spec_k<0     % It may happen due to singularity
                  spectrum(k)=0;
              elseif spec_k>1 % It may happen due to singularity
                  spectrum(k)=1;
              else
                  spectrum(k)=spec_k;
              end
          else
              if spec<0        % It may happen due to singularity
                  spectrum(k)=0;
              else
                  spectrum(k)=spec;
              end
          end
          chats(2*k-1:2*k)=chat; % The estimated cosine and sine coeff
      else
          spectrum(k)=0;
      end 
end
% -------------------------------------------------------------------------
% Critical value at (1-alpha) confidence level for the LSS
% NOTE: WE ONLY DEFINE THE CRITICAL VALUE FOR NORMALIZED SPECTRUM!
CritVal=ones(1,LOm)*(1-alpha^(2/(Lt-nc-2))); 


