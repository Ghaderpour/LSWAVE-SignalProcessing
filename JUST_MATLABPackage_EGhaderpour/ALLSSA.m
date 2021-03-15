function [tr_coeff,CS_coeff,res,norm_res,cov,freq1] = ALLSSA(t,f,varargin)
%% Antileakage Least-Squares Spectral Analysis (ALLSSA)
%
% function [tr_coeff,CS_coeff,res,norm_res,cov,freq1] = ALLSSA(t,f,varargin)
%
% Summary: This function iteratively finds an optimal set of sinusoids 
%           that along with trend fit best to the time series 
%
%          Reference: 
%          Ghaderpour, E., Liao, W., Lamoureux, M.P., 2018. 
%                 Antileakage least-squares spectral analysis for seismic 
%                 data regularization and random noise attenuation. 
%                 Geophysics, 83, V157–V170.  doi: 10.1190/geo2017-0284.1
% 
% Note: In the following descriptions, time and cyclic frequency can be
%       distance and wavenumber as well, respectively
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
% 9) 'decimal'   - Numeric. The number of decimals for frequency estimation
% 10) 'save'     - Logical. Save the ALLSSA results as .xlsx 
%
% Outputs:
% 1) tr_coeff    - The estimated trend coefficients ordered as follows: 
%                  For constant trend with multiple pieces: shifts
%                  For linear trend with multiple pieces: intercepts,slopes 
%                  For the quadratic trend (a+b*t+c*t^2): a, b, c
%                  For cubic trend (a+b*t+c*t^2+d*t^3): a, b, c, d
% 2) CS_coeff    - Estimated coefficients of cosine and sine functions of
%                  the estimated cyclic frequencies (freq)
% 3) res         - The residual series
% 4) norm_res    - The squared weighted L2 norm of the residual series 
% 5) cov         - Covariance matrix of the estimated coefficients
%                  ordered as: intercepts, slopes, cosine and sine of
%                  1st frequency, cosine and sine of 2nd frequency, etc.
% 6) freq        - Estimated cyclic frequencies

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
    if (~isvector(Omega)) || ischar(Omega)
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
%--------- Decimal --------------------------------------------------------
dec_ind = find(strcmpi(varargin,'decimal'));
if ~isempty(dec_ind)
    decimal = varargin{dec_ind+1};
    if (~isequal(floor(decimal),decimal)) || decimal < 0 || decimal > 4
        error('The decimal must be a nonnegative integer less than five')
    end
else
    decimal = 1; 
end
%--------- Save -----------------------------------------------------------
save_ind = find(strcmpi(varargin,'save'));
if ~isempty(save_ind); save = varargin{save_ind+1};
else; save = false; end
%==========================================================================
k = 1; check = true; freq1 = [];
inc = 1/(10^decimal);                  % cyclic frequency increment
while check 
    [spectrum, ~, ~, res, norm_res, coeff, cov] = LSSA(t,f, 'P',P, ...
                       'Omega',Omega, 'ind', ind, 'level',alpha,  ...
                        'trend',trend, 'slope',slope, 'freq',freq1);
    freq2 = freq1;
    Om_inds = find(spectrum == max(spectrum(:)));
    Om_ind = Om_inds(1);   % the index of peak with maximum energy/power 
    Om = round(Omega(Om_ind)); % the selected integer cyclic frequency
    %----------------------------------------------------------------------
    % Find and remove the highest cyclic frequency from freq1 
    % if it is close enough to Om 
    ZeroInd = find(round(freq1) == Om); LZ = length(ZeroInd); 
    if LZ > 0; freq1(ZeroInd(LZ)) = []; k = k - 1; end
    %------ Estimate LSS for a small neighbourhood around Om --------------
    if inc < 1; Om_range = Om - 0.5 : inc : Om + 0.5;   
    else; Om_range = Om; end
    [spectrum, CritVal, ~, ~, ~, ~, ~] = LSSA (t,f, 'P',P, ...
                         'Omega',Om_range, 'ind', ind, 'level',alpha,  ...
                         'trend',trend, 'slope',slope, 'freq',freq1);
    Om_new_inds = find(spectrum == max(spectrum(:)));
    Om_new_ind = Om_new_inds(1);
    %----- Check if the peak at Om_new is statistically significant -------
    if spectrum(Om_new_ind) > CritVal 
        if inc < 1; Om_new = Om - 0.5 + inc * (Om_new_ind - 1); 
        else; Om_new = Om; end
        % If peaks at the previous peak (rarely happens)
        if any(round(freq2,5) == round(Om_new,5)); freq1 = freq2;
        else; freq1(k) = Om_new; freq1 = sort(freq1); k = k + 1; end
    end
    %------ The termination condition -------------------------------------
    if  isequal(freq1,freq2); check = false; end
end
%==========================================================================
% coeff contains the estimated coefficients of trend and sinusoidal 
% functions. The following command will separate them
%--------------------------------------------------------------------------
Lcoeff = length(coeff); Lfreq = length(freq1);
if Lfreq > 0
    CS_coeff = coeff(Lcoeff - 2*Lfreq + 1 : Lcoeff);
    tr_coeff = coeff(1 : Lcoeff - 2*Lfreq);
else
    CS_coeff = [];
    tr_coeff = coeff;
end
%----Save the ALLSSA results as .xlsx  ------------------------------------
if save
    LSSA (t,f, 'P',P, 'Omega',[], 'ind',ind, 'trend',trend, ...
        'slope',slope, 'freq',freq1, 'save',true);
end