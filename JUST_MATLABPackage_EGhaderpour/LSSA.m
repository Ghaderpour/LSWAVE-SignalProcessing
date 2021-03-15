function [spectrum, CritVal, CScoeff, res, norm_res, coeff, cov] = ...
                                                       LSSA(t, f, varargin)
%% Least-Squares Spectral Analysis (LSSA)
%
% function [spectrum, CritVal, CScoeff, res, norm_res, coeff, cov] = ...
%                                                      LSSA(t, f, varargin)
%
% Summary: This function computes a frequency spectrum for a given 
%          time series based on the least-squares fit of sinusoids of 
%          different cyclic frequencies to the time series. 
%          The constitunets of known forms, such as trends and sinusoids 
%          of known frequencies are fitted first which is the 
%          Ordinary Least-Squares (OLS) estimation. Then, the spectrum of 
%          the residual series will be estimated considering the effect of 
%          removed known constituents.  
%
% Inputs:
% 1) t           - Vector of order n. Times of observations or measurements 
% 2) f           - Vector of order n. The time series values
% 3) 'P'         - Vector or Matrix of order n. The weights associated with
%                  the observations. Default is none
% 4) 'Omega'     - Vector. The initial cyclic frequencies 
% 5) 'ind'       - Vector. The indices of jumps in the trend component 
% 6) 'level'     - Numeric. The significance level, usually 0.01 or 0.05
%                  Default is 0.01. 
% 7) 'trend'     - String. 'none','constant','linear','quadratic','cubic'
% 8) 'slope'     - Logical. Select true to estimate a unique slope for the
%                  linear trend with multiple pieces
%                  Default is false allowing the linear pieces to have
%                  different slopes
% 9)'freq'       - Vector. The known cyclic frequencies of the sinusoids
%                  to be estimated/removed from the time series along with 
%                  the trend while their removal effect is considered 
%                  in the spectrum of residual 
% 10) 'save'     - Logical. Save the OLS results as .xlsx 
%
% Outputs:
% 1) Spectrum    - The normalized least-squares spectrum (LSS)
% 2) CritVal     - The Critical Value at (1-alpha) confidence level
% 3) CS_coeff    - Estimated coefficients of cosine and sine functions of
%                  input cyclic frequencies (Omega) to determine the
%                  amplitude and phase for each frequency (out-of-context)
% 4) res         - The residual series
% 5) norm_res    - The squared weighted L2 norm of the residual series 
% 6) coeff       - The ordinary least-squares estimated coefficients of 
%                  the trend and sinusoids of the known cyclic frequencies 
%                  ordered as follows: 
%                  For constant trend with multiple pieces: shifts
%                  For linear trend with multiple pieces: intercepts,slopes 
%                  For the quadratic trend (a+b*t+c*t^2): a, b, c
%                  For cubic trend (a+b*t+c*t^2+d*t^3): a, b, c, d
%                  Followed by the estimated coefficients of the cosine and 
%                  sine functions of the known cyclic frequencies if given
% 7) cov         - Covariance matrix of the estimated coefficients (coeff)
%                  with the same order as coeff  

%--------------------------------------------------------------------------
% Author: Ebrahim Ghaderpour
% Email:  ebrahim.ghaderpour@ucalgary.ca
% Copyright (c) 2021
%==========================================================================
%% Check the input arguments
Lt = length(t);  % Number of observations or measurements
St = size(t); if St(1) == 1; t = t'; end
Sf = size(f); if Sf(1) == 1; f = f'; end
if ~(Lt == length(f)); error('t and f must have the same size'); end
%--------- Weights (Vector or Matrix)--------------------------------------
P_ind = find(strcmpi(varargin,'P'));
if ~isempty(P_ind)
    P = varargin{P_ind+1};
    SP = size(P);
    if ~((SP(1) == Lt && SP(2) == Lt) || (SP(1) == 1 && SP(2) == Lt) || ...
            (SP(1) == Lt && SP(2) == 1) || (SP(1) == 1 && SP(2) == 1))
        error('P must be a vector or a matrix of order n')
    end
    if SP(1) == 1; P = P'; end
else
    P = 1;
end
%--------- Set of Cyclic Frequencies --------------------------------------
Omega_ind = find(strcmpi(varargin,'Omega'));
if ~isempty(Omega_ind)
    Omega = sort(varargin{Omega_ind+1});
    if ~isempty(Omega) && ((~isvector(Omega)) || ischar(Omega))
        error('Cyclic frequencies must be entered as a vector')
    end
else
    M = (Lt-1)/(t(Lt)-t(1));        % A simple estimate for sampling rate
    if isequal(floor(M/2),M/2); Omega = 1:floor(M/2)-1;
    else; Omega = 1:floor(M/2); end
end
%--------- Jump Indices ---------------------------------------------------
ind_ind = find(strcmpi(varargin,'ind'));
if ~isempty(ind_ind)
    ind = sort(varargin{ind_ind+1});
    if ~isempty(ind)
        if ~isequal(floor(ind),ind)
             error('The jump indices must be positive integers')
        elseif max(ind) > Lt-1 || min(ind) < 3
            error('The jump indices must be integers between 2 and n')
        end
    end
else
    ind = [];
end
%--------- Significance Level ---------------------------------------------
alpha_ind = find(strcmpi(varargin,'level'));
if ~isempty(alpha_ind)
    alpha = varargin{alpha_ind+1};
    if alpha > 1 || alpha < 0
        error('The significance level must be between 0 and 1')
    end
else
    alpha = 0.01;
end
%--------- Trend ----------------------------------------------------------
trend_ind = find(strcmpi(varargin,'trend'));
if ~isempty(trend_ind); trend = varargin{trend_ind+1};
else; trend = 'linear'; end
%--------- Slope ----------------------------------------------------------
slope_ind = find(strcmpi(varargin,'slope'));
if ~isempty(slope_ind); slope = varargin{slope_ind+1};
else; slope = false; end
%--------- Known Cyclic Frequencies ---------------------------------------
freq_ind = find(strcmpi(varargin,'freq'));
if ~isempty(freq_ind)
    freq = varargin{freq_ind+1};
    if (~isempty(freq)) && ((~isvector(freq)) || ischar(freq))
        error('The known cyclic frequencies must be entered as a vector')
    end
else
    freq = [];
end
%--------- Save -----------------------------------------------------------
save_ind = find(strcmpi(varargin,'save'));
if ~isempty(save_ind); save = varargin{save_ind+1};
else; save = false; end
%==========================================================================
Lt = length(t); LOm = length(Omega); Lfreq = length(freq);
ind = sort(ind); Lind = length(ind); 
CScoeff = zeros(1,2*LOm); spectrum = zeros; coeff = zeros;
% -------------------------------------------------------------------------
% Determine the number of constituents of known forms
nc = 2*Lfreq;
if     strcmpi(trend,'constant');             nc = nc + Lind + 1;  
elseif strcmpi(trend,'linear') && slope;      nc = nc + Lind + 2;
elseif strcmpi(trend,'linear') && ~slope;     nc = nc + 2*Lind + 2;
elseif strcmpi(trend,'quadratic') && slope;   nc = nc + Lind + 3;
elseif strcmpi(trend,'quadratic') && ~slope;  nc = nc + 2*Lind + 3;
elseif strcmpi(trend,'cubic') && slope;       nc = nc + Lind + 4;
elseif strcmpi(trend,'cubic') && ~slope;      nc = nc + 2*Lind + 4; 
end
if Lt < nc + 2
    error('The size of the time series (segment) is too small!')
end
%--------------------------------------------------------------------------
% Create a design matrix (D) based on the constituents of known forms 
% ---------- Constituent of intercept (Datum shift)
D = zeros(Lt,nc); q = 1;
if (strcmpi(trend,'constant') || strcmpi(trend,'linear') || ...
        strcmpi(trend,'quadratic') || strcmpi(trend,'cubic')) && Lind == 0                           
    D(:,q) = ones(1,Lt);  % constituent of known form for the intercept
    q = q + 1;
elseif (strcmpi(trend,'constant') || strcmpi(trend,'linear') || ...
        strcmpi(trend,'quadratic') || strcmpi(trend,'cubic')) && Lind > 0
    % constituents of known forms for intercepts
    range = (1:ind(1)-1)';
    D(range,q) = ones(length(range),1);
    q = q + 1;
    for j = 1:Lind
        if j < Lind; range = (ind(j):ind(j+1)-1)';
        else; range = (ind(j):Lt)'; end
        D(range,q) = ones(length(range),1);
        q = q + 1;
    end
end
%---------- Constituent of slope 
if (strcmpi(trend,'linear') || strcmpi(trend,'quadratic') || ...
        strcmpi(trend,'cubic')) && (Lind == 0 || slope)                        
    D(:,q) = t;  
    q = q + 1;
elseif (strcmpi(trend,'linear') || strcmpi(trend,'quadratic') || ...
        strcmpi(trend,'cubic')) && (Lind > 0 && ~slope)
    range = (1:ind(1)-1)';
    D(range,q )= t(range);
    q = q + 1;
    for j = 1:Lind
        if j < Lind; range = (ind(j):ind(j+1)-1)';
        else;        range = (ind(j):Lt)'; end
        D(range,q) = t(range);
        q = q + 1;
    end
end
%--------- Constituent of quadratic term
if strcmpi(trend,'quadratic') || strcmpi(trend,'cubic')
    D(:,q) = t.^2; q = q + 1; 
end
%--------- Constituent of cubic term
if strcmpi(trend,'cubic'); D(:,q) = t.^3; q = q + 1; end
% -------- Constituent of sinusoids of known frequencies
for jj = 1:Lfreq  
  D(:,q) = cos(2*pi*freq(jj)*t); q = q + 1;
  D(:,q) = sin(2*pi*freq(jj)*t); q = q + 1;
end
%--------------------------------------------------------------------------
%  Estimate the coefficients (coeff) for the constituents of known forms
%--------------------------------------------------------------------------
[nrow, ncol] = size(P);
if nc > 0
    if ncol == 1 && nrow == 1    % If the time series is equally weighted
        Ninv  = inv(D'*D);
        coeff = Ninv * (D'*f);
    elseif nrow > 1 && ncol == 1 % If P is a column vector
        Ninv  = inv(D'*(D.* kron(P,ones(1,nc))));
        coeff = Ninv * (D'*(f.*P));  
    else                         % If P is a square matrix
        Ninv  = inv(D'*P*D);
        coeff = Ninv * (D'*P*f); 
    end
end
% -------------------------------------------------------------------------
% Find the residual series: after removing the constituents of known forms
if (nc == 0); res = f;
else; res = f - D*coeff;  end    % Residual time series
% -------------------------------------------------------------------------
if nrow == 1 && ncol == 1    % If the time series is equally weighted
    resP = res;
elseif nrow > 1 && ncol == 1 % If P is a column vector
    resP = res.*P;
elseif nrow > 1 && ncol > 1  % If P is a square matrix     
    resP = P*res;    
end
norm_res = res'*resP;
% -------------------------------------------------------------------------
% Estimate the covariance matrix of coeff
if nc > 0; cov = (norm_res/(Lt-nc))*Ninv;
else; cov = []; end
% -------------------------------------------------------------------------
% Save the OLS results as .xlsx
% Source: https://doi.org/10.1007/s10291-019-0841-3
if save
    OLS_Est{1,1}  = 'Intercept';     OLS_Est{1,2}  = 'Error';
    OLS_Est{1,3}  = 'Slope';         OLS_Est{1,4}  = 'Error';
    OLS_Est{1,5}  = 'Quad_Coeff';    OLS_Est{1,6}  = 'Error';
    OLS_Est{1,7}  = 'Cub_Coeff';     OLS_Est{1,8}  = 'Error';
    OLS_Est{1,9}  = 'Frequency'; 
    OLS_Est{1,10} = 'Cos_Coeff';     OLS_Est{1,11} = 'Error'; 
    OLS_Est{1,12} = 'Sin_Coeff';     OLS_Est{1,13} = 'Error'; 
    OLS_Est{1,14} = 'Amplitude';     OLS_Est{1,15} = 'Error'; 
    OLS_Est{1,16} = 'Phase';         OLS_Est{1,17} = 'Error';
    OLS_Est{1,18} = 'Residual Norm'; OLS_Est{2,18} = round(norm_res,4);
    LTr = nc - 2*Lfreq;
    Tr_coeff = coeff(1:LTr); CS_coeff = coeff(LTr+1:nc);
    diagCov = diag(cov); Tr_std = sqrt(diagCov(1:LTr)); 
    CS_cov = cov(LTr+1:nc,LTr+1:nc); CS_var = diag(CS_cov);
    R1 = Tr_coeff(1:Lind+1); S1 = Tr_std(1:Lind+1);
    R2 = []; R3 = []; R4 = []; S2 = []; S3 = []; S4 = [];
    if strcmpi(trend,'linear')      
        R2 = Tr_coeff(Lind+2:LTr); S2 = Tr_std(Lind+2:LTr);
    elseif strcmpi(trend,'quadratic') && slope   
        R2 = Tr_coeff(Lind+2); S2 = Tr_std(Lind+2);
        R3 = Tr_coeff(Lind+3); S3 = Tr_std(Lind+3);
    elseif strcmpi(trend,'quadratic') && ~slope
        R2 = Tr_coeff(Lind+2:2*Lind+2); S2 = Tr_std(Lind+2:2*Lind+2);
        R3 = Tr_coeff(2*Lind+3); S3 = Tr_std(2*Lind+3);
    elseif strcmpi(trend,'cubic') && slope
        R2 = Tr_coeff(Lind+2); S2 = Tr_std(Lind+2);
        R3 = Tr_coeff(Lind+3); S3 = Tr_std(Lind+3);
        R4 = Tr_coeff(Lind+4); S4 = Tr_std(Lind+4);
    elseif strcmpi(trend,'cubic') && ~slope   
        R2 = Tr_coeff(Lind+2:2*Lind+2); S2 = Tr_std(Lind+2:2*Lind+2);
        R3 = Tr_coeff(2*Lind+3); S3 = Tr_std(2*Lind+3);
        R4 = Tr_coeff(2*Lind+4); S4 = Tr_std(2*Lind+4);
    end
    for j=1:length(R1)
        OLS_Est{j+1,1} = round(R1(j),4); OLS_Est{j+1,2} = round(S1(j),6); 
    end
    for j=1:length(R2)
        OLS_Est{j+1,3} = round(R2(j),4); OLS_Est{j+1,4} = round(S2(j),6); 
    end
    if ~isempty(R3)    
        OLS_Est{2,5} = round(R3,4); OLS_Est{2,6} = round(S3,6); 
    end
    if ~isempty(R4)
        OLS_Est{2,7} = round(R4,4); OLS_Est{2,8} = round(S4,6); 
    end
    for j = 1:Lfreq
        OLS_Est{j+1,9} = round(freq(j),4);
        c = CS_coeff(2*j-1); s = CS_coeff(2*j); 
        ce = sqrt(CS_var(2*j-1)); se = sqrt(CS_var(2*j));
        cse = CS_cov(2*j-1,2*j);
        OLS_Est{j+1,10} = round(c,4); OLS_Est{j+1,11} = round(ce,6);
        OLS_Est{j+1,12} = round(s,4); OLS_Est{j+1,13} = round(se,6);
        amp = sqrt(c^2+s^2); ampE = sqrt((c*ce)^2+(s*se)^2+2*c*s*cse)/amp;
        OLS_Est{j+1,14} = round(amp,4); OLS_Est{j+1,15} = round(ampE,6);
        OLS_Est{j+1,16} = round(2*atan((amp-s)/c),4);
        if c < 1e-10; c = 1e-10; end
        phaseE = 2*sqrt(((s/c)*(amp-s))^2*ce^2+(s-amp)^2*se^2-2*(s/c)*(amp-s)^2*cse)...
                        /(abs(amp*c)*(1+(amp-s)^2/c^2));
        OLS_Est{j+1,17} = round(phaseE,6);
    end
    filter = {'*.xlsx'};
    [file, path] = uiputfile(filter,'File Selection','LS_Estimates');
    FullFileName = fullfile(path,file);
    xlswrite(FullFileName, OLS_Est)
end
%--------------------------------------------------------------------------
% Critical value at (1-alpha) confidence level for the normalized LSS
CritVal = 1 - alpha^(2/(Lt - nc - 2)); 
%--------------------------------------------------------------------------
%             Calculate the least-squares spectrum (LSS)  
%--------------------------------------------------------------------------
for k = 1:LOm
  if ~any(round(Omega(k),5) == round(freq,5))
      Phi(:,1) = cos(2*pi*Omega(k)*t);
      Phi(:,2) = sin(2*pi*Omega(k)*t);
      if nrow == 1 && ncol == 1   % If the series is equally weighted
          N22 = Phi'*Phi;
          if nc > 0; N12 = D'*Phi; end
      elseif nrow > 1 && ncol == 1 % If P is a column vector
          N22 = Phi'*(Phi.*[P P]);
          if nc > 0; N12 = D'*(Phi.*[P P]); end
      else                         % If P is a square matrix
          N22 = Phi'*P*Phi;
          if nc > 0; N12 = D'*P*Phi; end
      end
      if nc > 0; chat = (N22 - N12'*Ninv*N12) \ (Phi'*resP);
      else; chat = (N22) \ (Phi'*resP); end
      spec = (resP'*Phi*chat)/norm_res;    % The normalized LSS
      if spec < 0; spectrum(k) = 0;   %It may happen due to singularity     
      elseif spec>1; spectrum(k) = 1; %It may happen due to singularity   
      else; spectrum(k) = spec;
      end
      % The estimated cosine and sine coefficients for Omega
      CScoeff(2*k-1 : 2*k) = chat; 
  else
      spectrum(k) = 0;
  end 
end