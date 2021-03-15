function [Loc, Mag, Dir, Type] = JUSTmonitor(t,f,varargin)
%% Near-Real-Time (NRT) Monitoring via Jumps Upon Spectrum and Trend (JUST) 
%
% function [Loc, Mag, Dir, Type] = JUSTmonitor(t,f,varargin)
%
% Summary: This function detects and classifies a disturbance within a 
%          time series in Near-Real-Time (NRT)
%
%          References: 
%          Ghaderpour, E., and Vujadinovic, T. The potential of the 
%                 least-squares spectral and cross-wavelet analyses for 
%                 near-real-time disturbance detection within unequally 
%                 spaced satellite image time series. 
%                 Remote Sensing, 2020, 12, 2446; doi: 10.3390/rs12152446
%          Ghaderpour, E., Vujadinovic, T. Change Detection within Remotely 
%                 Sensed Satellite Image Time Series via Spectral Analysis, 
%                 Remote Sensing, 2020, 12, 4001; doi:10.3390/rs12234001
% 
% Inputs:
% 1) t           - Vector of order n. Times of observations or measurements 
% 2) f           - Vector of order n. The time series values
% 3) 'P'         - Vector or Matrix of order n. The weights associated with
%                  the observations. Default is none
% 4) 'start'     - Numeric. The time index of the start time of monitoring 
% 5) 'history'   - Numeric. The time index of start of the stable 
%                  historical period. This index can be obtained by JUST OR
%                  entered by user via other information and/or methods.
%                  Default is 1 which corresponds to t(1)   
% 6) 'h'         - Numeric. The fraction of the historical time series 
%                  used during monitoring. Default. 0.25
% 7) 'season'    - String. Two options:   
%                          ALLSSA season-trend model ('ALLSSA')
%                          OLS season-trend model    ('OLS')
%                          Only the trend model      ('none')
% 8) 'Omega'     - Vector. The initial cyclic frequencies 
%                  Default is Omega = [1,2,3,4]   for OLS and 
%                  Omega = 0.8:0.2:3.8 for ALLSSA as the initial set
% 9) 'trend'     - String. 'constant', 'linear'. Default. 'linear'
%                  Choose 'constant' to only estimate the jump magnitude
%                  Choose 'linear' to estimate the jump magnitude and
%                  direction                           
% 10) 'level'    - Numeric. The significance level; usually 0.01 or 0.05. 
%                  Default is 0.01
% 11)'display'   - Logical. To display the result
% 
% Outputs:
% 1) Loc         - Detected disturbance location
% 2) Mag         - Estimated magnitude of the detected jump
% 3) Dir         - Estimated direction of the detected jump only provided
%                  when 'trend' is selected to be 'linear'
% 4) Type        - The disturbance type: 0 means no disturbance, 
%                  1 means abrupt change, and 2 means phenological change

%--------------------------------------------------------------------------
% Author: Ebrahim Ghaderpour 
% Email:  ebrahim.ghaderpour@ucalgary.ca
% Copyright (c) 2021
%==========================================================================
%% Check the input arguments
Lt = length(t);        % Number of observations or measurements
start_time = t(1);
t  = t - start_time;   % For the sake of computational efficiency 
St = size(t); if St(1) == 1; t = t'; end
Sf = size(f); if Sf(1) == 1; f = f'; end
if ~(Lt == length(f)); error('t and f must have the same size'); end
%--------- Weights --------------------------------------------------------
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
%--------- Time Index of the Start of Monitoring --------------------------
start_ind = find(strcmpi(varargin,'start'));
if ~isempty(start_ind); start = varargin{start_ind+1};
else; start = Lt - 1; end
%--------- Time Index of the Stable Historical Period ---------------------
history_ind = find(strcmpi(varargin,'history'));
if ~isempty(history_ind); history = varargin{history_ind+1};
else; history = 1; end
%--------- Check the Start Time of the Monitoring Period ------------------
Lhis = start - history;
if Lhis < Lt/(t(Lt)-t(1))
    error('Insufficient size for the stable historical time series')
elseif start > Lt 
    error('The start time index of the monitoring period must be <= n')
end
%--------------------------------------------------------------------------
% Size of Time Series in the Monitoring Process Divided By Size of the 
% Historical Time Series --------------------------------------------------
h_ind = find(strcmpi(varargin,'h'));
if ~isempty(h_ind)
    h = varargin{h_ind+1};
    if h < 0 || h > 1; error('h must be in interval [0,1]'); end
else
    h = 0.25;
end
%--------- Seasonal Component ('ALLSSA', 'OLS')----------------------------
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
    else; Omega = 1:4; end
end
%--------- Trend Type to Estimate Jump Magnitude and Direction-------------
trend_ind = find(strcmpi(varargin,'trend'));
if ~isempty(trend_ind); trend = varargin{trend_ind+1};
else; trend = 'linear'; end
%--------- Significance Level for ALLSSA ----------------------------------
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
t1 = t(history:Lt);            
f1 = f(history:Lt);
[nrow,ncol] = size(P);
if ncol > 1 && nrow > 1 
    P0 = P(history:Lt, history:Lt); 
    P1 = P0(1:Lhis, 1:Lhis);
elseif ncol == 1 && nrow > 1
    P0 = P(history:Lt);
    P1 = P0(1:Lhis);
else               
    P0 = 1; P1 = 1;
end
%--------------------------------------------------------------------------
%Estimate the Trend and Seasonal Coefficients for the Stable History Period
if strcmpi(season,'ALLSSA')
    [tr_coeff, CS_coeff, ~, ~, ~, freq] = ALLSSA(t1(1:Lhis), ...
                      f1(1:Lhis), 'P', P1, 'Omega', Omega, 'level', alpha);
    coeff = cat(1, tr_coeff, CS_coeff);                                    
elseif strcmpi(season,'OLS')                  
    [~, ~, ~, ~, ~, coeff, ~] = LSSA(t1(1:Lhis), f1(1:Lhis), ...
                                   'P',P1, 'Omega',[], 'freq',Omega); 
    freq = Omega; 
else
    [~,~,~,~,~,coeff,~] = LSSA(t1(1:Lhis), f1(1:Lhis), 'P',P1, 'Omega',[]); 
    freq = [];  
end
Lfreq = length(freq);
%---- Signal (Forecast) ---------------------------------------------------
Loc = 0; Mag = 0; Dir = 0;
nr = Lt - history + 1; 
D = zeros(nr,2*Lfreq+2);
D(:,1) = ones(nr,1);
D(:,2) = t1;
for q =1:Lfreq
    D(:,2*q+1) = cos(2*pi*freq(q)*t1);
    D(:,2*q+2) = sin(2*pi*freq(q)*t1);
end
signal = D*coeff;        % The history and forecast signal
%--------------------------------------------------------------------------
%  Near-Real-Time Monitoring
%--------------------------------------------------------------------------
% The size of the time series during the monitoring process
Lh = Lhis - floor(h*Lhis);
if Lh > Lhis - 2; Lh = Lh - 2;
elseif Lh == 0; Lh = Lh + 1; end
for d = 1:Lt-start+1
    t2 = t1(Lh:Lhis+d);          % t2 = t1(Lh+d:Lhis+d) is another option!
    f2 = f1(Lh:Lhis+d) - signal(Lh:Lhis+d);
    if ncol > 1 && nrow > 1;  P2 = P0(Lh:Lhis+d, Lh:Lhis+d);
    elseif ncol == 1 && nrow > 1; P2 = P0(Lh:Lhis+d);
    else; P2 = 1; 
    end
    [spectrum, CritVal, ~, ~, ~, ~, ~] = LSSA (t2,f2, 'P',P2, ...
                     'Omega',0.1:0.1:3.5,'trend','constant','level',alpha);                               
    JumpInd = start + d - 1;
    Loc = t(JumpInd) + start_time; 
    MaxSpec = max(spectrum(1:9));
    if MaxSpec > CritVal; Type = 1; 
    elseif max(spectrum) > CritVal && MaxSpec <= CritVal; Type = 2; break;
    else; Type = 0;
    end      
    if Type == 1
        JInd = Lhis + d;
        if JInd < length(t1)-1
            if strcmpi(season,'ALLSSA')
                [tr_coeff, ~, ~, ~, ~, ~] = ALLSSA(t1, f1, 'P',P0, ...
                  'Omega',Omega, 'ind',JInd, 'trend',trend, 'level',alpha);
            elseif strcmpi(season,'OLS') 
                [~, ~, ~, ~, ~, coeff, ~] = LSSA(t1, f1, 'P',P0, ...
                      'Omega',[], 'ind',JInd, 'trend',trend, 'freq',Omega);
                tr_coeff = coeff;
            else
                [~, ~, ~, ~, ~, coeff, ~] = LSSA(t1, f1, 'P',P0, ...
                                    'Omega',[], 'ind',JInd, 'trend',trend);
                tr_coeff = coeff;
            end
            if strcmpi(trend,'linear')
                Mag = (tr_coeff(2)-tr_coeff(1))+...
                                     (tr_coeff(4)-tr_coeff(3))*t1(JInd);
                Dir = tr_coeff(4)-tr_coeff(3);
            else
                Mag = tr_coeff(2)-tr_coeff(1); 
            end
        end
        break
    end
end
%--------------------------------------------------------------------------
if display
    PlotNRT(t+start_time, f, P, start, history, Type, JumpInd, signal)
end
%==========================================================================
function PlotNRT(t, f, P, start, history, Type, JumpInd,  signal)
%% Display the Near-Real-Time Monitoring Results
%
% function PlotNRT(t, f, P, start, history, Type, JumpInd,  signal)
%
% Summary: This function displays the near-real-time monitoring results 
%
% Inputs:
% 1) t         - Vector of order n. Times of observations or measurements 
% 2) f         - Vector of order n. The time series values
% 3) P         - Vector or Matrix of order n. The weights associated with
%                the observations.                   
% 5) start     - Numeric. The time index of the start time of monitoring 
% 6) history   - Numeric. The time index of start of the stable 
%                historical period.  
% 7) Type      - Numeric. The disturbance Type: 0 means no disturbance, 
%                1 means abrupt change, and 2 means phenological change                 
% 8) JumpInd   - Numeric. The time index of detected disturbance 
% 9) signal    - Vector. The least-squares fit to the stable historical 
%                time series extending to the monitoring period

%==========================================================================
%% Check the input arguments
figure;
Lt = length(t);
Dif = max(f)-min(f);
xpos = t(start-1)+(t(start)-t(start-1))/2;
[nrow,ncol] = size(P);
if ncol > 1 && nrow > 1; Std = sqrt(diag(inv(P)));
elseif ncol == 1 && nrow > 1; Std = sqrt(1./P); end
hold on  
if length(P)>1
    MStd = 0.5*Dif * (Std - min(Std))/(max(Std)-min(Std)) + 0.1*Dif;                                                                                       
    errorbar(t, f, MStd, 'LineStyle','none','CapSize',3)
    x = xticks; Lx = length(xticks); y = yticks; Ly = length(yticks); 
    rectangle('Position', [xpos y(1) x(Lx)-xpos y(Ly)-y(1)],...
                                                 'FaceColor',[0.7,0.7,0.7])                                       
    axis([x(1) x(Lx)  y(1)  y(Ly)])
    [X1,Y1] = meshgrid(x,y) ;
    line(X1,Y1,'color','k','LineWidth', 0.1,'LineStyle',':')
    line(X1',Y1','color','k','LineWidth', 0.1,'LineStyle',':') 
    h1=errorbar(t, f, MStd,'Marker','d','MarkerSize',5,'LineStyle','--',...
                'MarkerFaceColor',[0,0.5,0],'Color',[0,0.5,0],'CapSize',3);
    lgd{1} = 'Value & Error Bar';
else                                            
    plot(t, f,'d','LineWidth',0.1,'MarkerSize',4)
    x = xticks; Lx = length(xticks); y = yticks; Ly = length(yticks); 
    rectangle('Position', [xpos y(1) x(Lx)-xpos y(Ly)-y(1)],...
                                                 'FaceColor',[0.7,0.7,0.7])                                       
    axis([x(1) x(Lx)  y(1)  y(Ly)])
    [X1,Y1] = meshgrid(x,y) ;
    line(X1,Y1,'color','k','LineWidth', 0.1,'LineStyle',':')
    line(X1',Y1','color','k','LineWidth', 0.1,'LineStyle',':') 
    h1=plot(t, f,'--d','color',[0,0.5,0],'LineWidth',0.1,'MarkerSize',4,...
            'MarkerFaceColor',[0,0.5,0]);
    lgd{1} = 'Time Series';
end
h2=plot(t(history:Lt), signal,'--ob','LineWidth', 0.5,...
                                     'MarkerFaceColor','b','MarkerSize',4);
lgd{2} = 'Signal (Forecast)';
if Type == 1
    h3=plot(t(JumpInd),f(JumpInd),'rp','MarkerFaceColor','r','MarkerSize',15);
    title('Near-Real-Time Monitoring')
    lgd{3} = 'Abrupt Change';
elseif Type == 2
    h3=plot(t(JumpInd),f(JumpInd),'gp','MarkerFaceColor','g','MarkerSize',15);
    title('Near-Real-Time Monitoring')
    lgd{3} = 'Intra-annual Change';
end
if length(lgd)==3; legend([h1,h2,h3],lgd,'Location','Best')
else; legend([h1,h2],lgd,'Location','Best'); end
hold off                              
ylabel('Value')
xlabel('Time')
set(gca, 'TickDir', 'out')
%--------------Save the figure --------------------------------------------
set(gcf, 'Color', 'w'); width = 10; height = 4;
set(gcf, 'PaperUnits', 'inches','PaperPosition',[0 0 width height]); 
% print( '-dtiff',  'Near-Real-Time', '-r300');