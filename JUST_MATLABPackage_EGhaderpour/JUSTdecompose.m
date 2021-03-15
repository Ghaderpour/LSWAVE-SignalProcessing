function [trend, seasonal, remainder, trend_reg, seasonal_reg] = ...
                                              JUSTdecompose(t, f, varargin)
%% JUST Decomposition into Trend, Seasonal, and Remainder Components
%
% function [trend, seasonal, remainder, trend_reg, seasonal_reg] = ...
%                                             JUSTdecompose(t, f, varargin)
%
% Summary: This function uses the jump indices estimated by JUST to
%          decompose the time series into trend, seasonal, and remainder 
%          components via ALLSSA or OLS
%
%          Reference: 
%          Ghaderpour, E., Vujadinovic, T. Change Detection within Remotely 
%                 Sensed Satellite Image Time Series via Spectral Analysis, 
%                 Remote Sensing, 2020, 12, 4001; doi:10.3390/rs12234001
% 
% Inputs:
% 1) t           - Vector of order n. Times of observations or measurements 
% 2) f           - Vector of order n. The time series values
% 3) 'P'         - Vector or Matrix of order n. The weights associated with
%                  the observations. Default is none
% 4) 'tt'        - Vector. Equally spaced times for regularization purposes
% 5) 'size'      - Numeric. The window or segment size (R)
%                  Default is the average sampling rate tripled (3M)
% 6) 'step'      - Numeric. The translation step delta that is less than R 
%                  Default is the average sampling rate (M)
% 7) 'season'    - String. Three options: 
%                          ALLSSA season-trend model ('ALLSSA')
%                          OLS season-trend model    ('OLS')
%                          Only the trend model      ('none')
% 8) 'Omega'     - Vector. The initial cyclic frequencies 
%                  Default is Omega = [1,2,3,4]   for OLS and 
%                  Omega = 0.8:0.2:3.8 for ALLSSA as the initial set
% 9) 'ind'       - Vector. The indices of jumps in the trend component
% 10)'level'     - Numeric. The significance level, usually 0.01 or 0.05
%                  Default is 0.01. Applies only for ALLSSA not OLS 
% 11)'display'   - Logical. To display the decomposition result
% 
% Outputs:
% 1) trend        - The estimated trend component of the time series
% 2) seasonal     - The estimated seasonal component of the time series 
% 3) remainder    - The estimated remainder component of the time series
% The following outputs will be generated only when 'tt' is given for
% regularization: 
% 4) trend_reg    - The regularized trend component of the time series
% 5) seasonal_reg - The regularized seasonal component of the time series 
%
% NOTE: The estimated trend component may have curvature due to the 
%       smoothing operation, i.e., when the Gaussian moving average is used
%       for all the estimated linear trends within the translating windows.
%       Also, at the end of the process, the components will be diplayed 
%       in such a way that the seasonal component has zero mean.

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
%--------- Equally Spaced Times -------------------------------------------
tt_ind = find(strcmpi(varargin,'tt'));
if ~isempty(tt_ind)
    tt = varargin{tt_ind+1};
    check = true;
else
    tt = t;
    check = false;
end
Stt = size(tt); if Stt(1) == 1; tt = tt'; end
Ltt = length(tt);
if (t(Lt) < tt(Ltt)) || (tt(1) < t(1))
    error('The range of time must be within the range of t')
end
%--------- Window Size ----------------------------------------------------
size_ind = find(strcmpi(varargin,'size'));
M = floor(Ltt/(tt(Ltt)-tt(1))); % A simple estimate for sampling rate
if ~isempty(size_ind)
    R = varargin{size_ind+1};
    if R > Ltt || R < M
        error('The window size must be an integer between M and n')
    end
else
    R = 3*M;
end
%--------- Translation Step -----------------------------------------------
step_ind = find(strcmpi(varargin,'step'));
if ~isempty(step_ind)
    delta = varargin{step_ind+1};
    if delta >= R || delta < 1
        error('Translation step must be a positive integer less than R')
    end
else
    delta = M;
end
%--------- Seasonal Component ('ALLSSA', 'OLS', 'none')--------------------
season_ind = find(strcmpi(varargin,'season'));
if ~isempty(season_ind); season = varargin{season_ind+1};
else; season = 'ALLSSA'; end
%--------- Initial Set of Cyclic Frequencies ------------------------------
Omega_ind = find(strcmpi(varargin,'Omega'));
if ~isempty(Omega_ind)
    Omega = sort(varargin{Omega_ind+1});
    if (~isvector(Omega)) || ischar(Omega)
        error('Cyclic frequencies must be entered as a vector')
    end
else
    if strcmpi(season,'ALLSSA'); Omega = 0.8:0.2:3.8; 
    else; Omega = 1:4; end % Cyclic frequencies for OLS
end
LOm = length(Omega);
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
%--------- Display --------------------------------------------------------
display_ind = find(strcmpi(varargin,'display'));
if ~isempty(display_ind); display = varargin{display_ind+1};
else; display = true; end
%==========================================================================
% The dilation factor for the Gaussian moving average:
DilFac = log(0.1)/((R/(6*M))^2);
trends0 = zeros(Lt,1); seasonals0 = zeros(Lt,1); remainders0 = zeros(Lt,1); 
trends = zeros(Ltt,1); seasonals = zeros(Ltt, 1); 
sumW0 = zeros(Lt,1); sumW  = zeros(Ltt, 1); 
R1 = R - 1; lb = 1; ub = lb + R1;
[nrow, ncol] = size(P);
stat = true; 
start_time = tt(1);
tt = tt - start_time;  % For the sake of computational efficiency 
t  = t  - start_time;  % For the sake of computational efficiency 
while ub <= Ltt + R1 && stat
    % The condition for the last segment to also have size R
    if ub > Ltt; lb = Ltt - R1; ub = Ltt; stat = false; end
    tty = tt(lb:ub);
    if check
        lb0 = find (t >= tt(lb),1 ,'first');
        ub0 = find (t <= tt(ub),1 ,'last');
        R0  = ub0 - lb0 + 1;
    else
        lb0 = lb; ub0 = ub; R0  = R;
    end
    ty = t(lb0:ub0);
    y  = f(lb0:ub0);
    if ncol == 1 && nrow == 1;  Py = 1;           % If P is a scale    
    elseif nrow>1 && ncol == 1; Py = P(lb0:ub0);  % If P is a column vector        
    else; Py = P(lb0:ub0, lb0:ub0);               % If P is a square matrix   
    end
    % Locate the jump index within each segment
    JumpInd = intersect (lb0 + 2:ub0 - 1, ind) - lb0 + 1;
    if strcmpi(season,'ALLSSA')
        [tr_coeff, CS_coeff, remainder, ~, ~, freq] = ALLSSA(ty, y, ...
               'P',Py, 'Omega', Omega, 'ind', JumpInd, 'level', alpha);
    elseif strcmpi(season,'OLS')                  
        [~, ~, ~, remainder, ~, coeff, ~] = LSSA(ty, y, 'P',Py, ...
                           'Omega',[], 'ind', JumpInd, 'freq',Omega); 
        Lcoeff = length(coeff);
        tr_coeff = coeff(1 : Lcoeff-2*LOm);
        CS_coeff  = coeff(Lcoeff-2*LOm+1 : Lcoeff);
        freq = Omega;
    else
        [~, ~, ~, remainder, ~, coeff, ~] = LSSA(ty, y, 'P',Py, ...
                                               'Omega',[], 'ind', JumpInd);
        tr_coeff = coeff; CS_coeff  = []; freq = [];
    end
    season_t =  SeasonalComponent(ty, freq, CS_coeff);
    trend_t  =  y - season_t - remainder;
    %----- Enforce the seasonal component to have zero mean ---------------
    % NOTE: The trend and seasonal components are simultaneously estimated 
    ms = mean(season_t); trend_t = trend_t + ms; season_t = season_t - ms;
    %----------------------------------------------------------------------
    if mod(R0, 2); tau0 = ty(floor(R0/2)+1); 
    else; tau0 = (ty(R0/2) + ty(R0/2 + 1))/2; end  
    GaussW = exp(DilFac * (ty - tau0).^2);   % The Gaussian weight function
    trends0(lb0:ub0)     = trends0(lb0:ub0)     + trend_t.*GaussW;
    seasonals0(lb0:ub0)  = seasonals0(lb0:ub0)  + season_t.*GaussW;
    remainders0(lb0:ub0) = remainders0(lb0:ub0) + remainder.*GaussW;
    sumW0(lb0:ub0) = sumW0(lb0:ub0) + GaussW;
    if check
        trend_tt  =  TrendComponent(ty, tty, JumpInd, tr_coeff);
        season_tt =  SeasonalComponent(tty, freq, CS_coeff);
        ms1 = mean(season_tt);
        trend_tt  = trend_tt  + ms1;
        season_tt = season_tt - ms1;
        if mod(R, 2)
            tau = tty(floor(R/2)+1); 
        else 
            tau = (tty(R/2) + tty(R/2 + 1))/2;
        end
        GW = exp(DilFac * (tty - tau).^2); % The Gaussian weight function
        trends(lb:ub)    = trends(lb:ub)    + trend_tt.*GW;
        seasonals(lb:ub) = seasonals(lb:ub) + season_tt.*GW;
        sumW(lb:ub) = sumW(lb:ub) + GW;
    end
    lb = lb + delta; ub = lb + R1;
end
trend     = trends0./sumW0;          % The weighted average for trend
seasonal  = seasonals0./sumW0;       % The weighted average for seasonal
remainder = remainders0./sumW0;      % The weighted average for remainder
if check
    trend_reg    = trends./sumW;     % The weighted average for trend
    seasonal_reg = seasonals./sumW;  % The weighted average for seasonal
else
    trend_reg = []; seasonal_reg = [];
end
if display
    PlotDecompResults(t+start_time, f, P, trend, seasonal, remainder, ...
                      tt+start_time, trend_reg, seasonal_reg);
end
%==========================================================================
function trend = TrendComponent(t, tt, ind, tr_coeff)
%% Trend Component Creation
%
% function trend = TrendComponent(t, tt, ind, tr_coeff)
%
% Summary: This function gets the estimated intercepts and slopes of linear
%          trend with multiple pieces and jump indices to create the 
%          trend component.
%
% Inputs:
% 1) t            - Times of the observations or measurements as a vector
% 2) tt           - Another time vector for regularization purposes
% 3) ind          - The indices of jumps in the trend component 
% 4) tr_coeff     - Estimated intercepts and then slopes of linear pieces
%
% Outputs:
% 1) trend        - The trend series

%==========================================================================
%% Check the input arguments
Ltt = length(tt); Lind = length(ind); JumpInd = zeros(Lind, 1);
for k = 1:Lind; JumpInd(k) = find (tt <= t(ind(k)),1 ,'last'); end
% Create a design matrix D for the linear trend with multiple pieces
D = zeros(Ltt, 2*Lind+2); q = 1;
if Lind == 0                           
    D(:,1) = ones(1,Ltt);  % constituent of known form for the intercept
    D(:,2) = tt;           % constituent of known form for slope
else
    % constituents of known forms for intercepts
    ind = (1:JumpInd(1)-1)';
    D(ind, q)  = ones(length(ind),1);
    D(ind, q+Lind+1) = tt(ind);
    q = q + 1;
    for j = 1:Lind
        if j < Lind; ind = (JumpInd(j):JumpInd(j+1)-1)';
        else; ind = (JumpInd(j):Ltt)'; end
        D(ind, q) = ones(length(ind),1);
        D(ind, q+Lind+1) = tt(ind);
        q = q + 1;
    end
end
trend = D * tr_coeff; 
%==========================================================================
function seasonal = SeasonalComponent(t, freq, CS_Coeff)
%% Seasonal Component Creation
% 
% function seasonal = SeasonalComponent(t, freq, CS_Coeff)
%
% Summary: This function gets the estimated coefficients for sinusoids and 
%          their frequencies to create the seasonal component.
%
% Inputs:
% 1) t            - Times of the observations or measurements as a vector
% 2) freq         - Cyclic frequencies (enter as a vector) 
% 3) CSCoeff      - A vector that contains the coefficients of cosine and
%                   sine functions of each cyclic frequency consecutively
%
% Outputs:
% 1) seasonal     - The seasonal series
% Note: This function can also be used for regularization, i.e., one may 
%       input any equally spaced vector t to generate equally spaced series

%==========================================================================
%% Check the input arguments
seasonal = zeros(length(t),1); freq_ind = 1; k = 1;
while freq_ind <= length(freq) 
    seasonal = seasonal + CS_Coeff(k) * cos(2*pi*freq(freq_ind)*t) + ...
                          CS_Coeff(k + 1) * sin(2*pi*freq(freq_ind)*t); 
    freq_ind = freq_ind + 1;
    k = k + 2;
end
%==========================================================================
function PlotDecompResults(t, f, P, trend, seasonal, remainder, ...
                                               tt, trend_reg, seasonal_reg)
%% Display JUST Decomposition Results  
%
% function PlotDecompResults(t, f, P, trend, seasonal, remainder, ...
%                                              tt, trend_reg, seasonal_reg)
%
% Summary: This function displays the decomposition result and saves it
%
% Inputs:
% 1) t            - Vector of order n. Times of observations/measurements 
% 2) f            - Vector of order n. The time series values
% 3) P            - Vector or Matrix of order n. The weights associated 
%                   with the observations. Default is 1 (a scale)
% 4) trend        - Vector of order n. The estimated trend component 
% 5) seasonal     - Vector of order n. The estimated seasonal component  
% 6) remainder    - Vector of order n. The estimated remainder component 
% 7) tt           - Vector of order m. Equally spaced times 
% 8) trend_reg    - Vector of order m.The regularized trend component
% 9) seasonal_reg - Vector of order m .The regularized seasonal component  
%==========================================================================
%% Check the input arguments
Ltt = length(tt);
[nrow, ncol] = size(P);
if ncol == 1 && nrow == 1;    Wt = 1; pn = 1;  % If P is a scale   
elseif nrow > 1 && ncol == 1; Wt = P; pn = 0;  % If P is a column vector  
else; Wt = diag(P); pn = 0;                    % If P is a square matrix  
end
%--------------Plot the time series ---------------------------------------
figure;
subplot(5-pn,1,1);
plot(t, f,':d','color', [0,0.5,0],'MarkerFaceColor',[0,0.5,0],...
    'LineWidth', 0.1,'MarkerSize',4)
title('Time Series Decomposition')
ylabel('Value')
vs = 0.1 * (max(f)-min(f));
if vs == 0; vs = 0.01; end
axis([tt(1)  tt(Ltt)  min(f) - vs   max(f) + vs])
grid on
set(gca,'Position',[0.1 0.79-pn*0.045 0.8 0.18+pn*0.045]); %[X Y W H]
set(gca,'XTickLabel',[]);   %remove x-axis tick labels 
set(gca, 'fontsize' ,12) ;
%--------------Plot the weights -------------------------------------------
if ~pn
    subplot(5,1,2);
    stem(t,Wt,'sk','MarkerFaceColor','k', 'LineWidth',0.1, 'MarkerSize',3)
    set(gca,'XTickLabel',[]);      %remove x-axis tick labels 
    set(gca, 'YAxisLocation', 'right')
    ylabel('Weight')
    vw = 0.1 * (max(Wt)-min(Wt));
    if vw == 0; vw = 0.1; end
    axis([tt(1) tt(Ltt) min(Wt) - vw   max(Wt) + vw])
    grid on
    set(gca,'Position',[0.1 0.61 0.8 0.18] );  %[X Y W H]
    set(gca, 'fontsize' ,12) ;
end
%--------------Plot the trend component -----------------------------------
subplot(5-pn,1,3-pn);
if ~isempty(trend_reg)
%     plot(t, trend,'ob','MarkerFaceColor','b','MarkerSize',3)
%     hold on
    plot(tt, trend_reg,'-b','LineWidth', 0.5)
%     hold off
else
    plot(t, trend,'-b','LineWidth', 1)
end
set(gca,'XTickLabel',[]);       %remove x-axis tick labels 
ylabel('Trend')
axis([tt(1)  tt(Ltt)  min(f) - vs   max(f) + vs])
grid on
set(gca,'Position',[0.1 0.43+pn*0.09 0.8 0.18+pn*0.045]);
set(gca, 'fontsize' ,12) ;
%--------------Plot the seasonal component --------------------------------
subplot(5-pn,1,4-pn);
if ~isempty(seasonal_reg)
%     plot(t, seasonal,'ob','MarkerFaceColor','b','MarkerSize',3)
%     hold on
    plot(tt, seasonal_reg,'-b','LineWidth', 0.5)
%     hold off
else
    plot(t, seasonal,':ob','MarkerFaceColor','b',...
        'MarkerSize',3,'LineWidth', 0.1)
end
set(gca,'XTickLabel',[]);       %remove x-axis tick labels
set(gca, 'YAxisLocation', 'right')
ylabel('Seasonal')
vs = 0.1 * (max(seasonal)-min(seasonal));
if vs == 0; vs = 0.01; end
axis([tt(1)  tt(Ltt)  min(seasonal) - vs    max(seasonal) + vs])
grid on
set(gca,'Position',[0.1 0.25+pn*0.045 0.8 0.18+pn*0.045] );
set(gca, 'fontsize' ,12) ;
%--------------Plot the remainder component -------------------------------
subplot(5-pn,1,5-pn);
plot(t, remainder,':ob','MarkerFaceColor','b',...
    'MarkerSize',3,'LineWidth', 0.1)
xlabel('Time')
ylabel('Remainder')
vr = 0.1 * (max(remainder)-min(remainder));
if vr == 0; vr = 0.01; end
axis([tt(1)  tt(Ltt)  min(remainder) - vr    max(remainder) + vr])
grid on
set(gca,'Position',[0.1 0.07 0.8 0.18+pn*0.045] );
set(gca, 'fontsize' ,12) ;
%--------------Save the figure --------------------------------------------
set(gcf, 'Color', 'w'); width = 10; height = 8;
set(gcf, 'PaperUnits', 'inches','PaperPosition',[0 0 width height]); 
% print( '-dtiff',  'TS_Decomposition', '-r300');