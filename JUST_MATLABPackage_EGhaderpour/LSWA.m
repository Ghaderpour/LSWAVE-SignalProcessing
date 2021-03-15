function [spectrogram, stoch_surf, amp_spectrogram] = LSWA(t, f, varargin)
%% Least-Squares Wavelet Analysis (LSWA)
% 
% function [spectrogram, stoch_surf, amp_spectrogram] = LSWA(t, f, varargin)
%
% Summary: This function computes a time-frequency spectrogram for a 
%          given time series based on the least-squares fit of sinusoids 
%          of different cyclic frequencies to the time series segments. 
%          The constitunets of known forms, such as trends and sinusoids 
%          of known frequencies are fitted first to each segment. 
%          Then, the spectra of the residual series segments will be 
%          estimated considering the effects of the removed known 
%          constituents that will ultimately obtain a spectrogram.  
% 
%          Reference:
%          Ghaderpour, E., Pagiatakis, S.D., 2017. Least-squares wavelet 
%                 analysis of unequally spaced and non-stationary time  
%                 series and its applications. Mathematical Geosciences,
%                 49, 819-844.doi: 10.1007/s11004-017-9691-0
%
% Inputs:
% 1) t           - Vector of order n. Times of observations or measurements 
% 2) f           - Vector of order n. The time series values
% 3) 'P'         - Vector or Matrix of order n. The weights associated with
%                  the observations. Default is none
% 4) 'tt'        - Vector. Equally spaced times for regularization purposes
% 5) 'rate'      - Numeric. The (average) sampling rate (M)
% 6) 'Omega'     - Vector. The cyclic frequencies 
% 7) 'ind'       - Vector. The indices of jumps in the trend component 
% 8) 'level'     - Numeric. The significance level (alpha), 
%                  usually 0.01 or 0.05. Default is 0.01. 
% 9) 'trend'     - String. 'constant', 'linear', 'quadratic', 'cubic'
% 10)'slope'     - Logical. Select true to estimate a unique slope for the
%                  linear trend with multiple pieces
%                  Default is false allowing the linear pieces to have
%                  different slopes
% 11)'freq'      - Vector. The known cyclic frequencies of the sinusoids
%                  to be estimated/removed from the time series while their
%                  removal effect is considered in the spectrum of residual
% 12) 'L1'       - Numeric. The number of cycles of sinusoids being fitted 
%                  to the segments of the time series. 
%                  Recommeded: L1=2 when morlet = 0
% 13) 'L0'       - Numeric. The number of additional observations to 
%                  increase the translating window size (e.g., L0=5, L0=20)
% 14) 'morlet'   - Numeric. The Morlet coefficient. 
%                  Selecting 0 will fit the sinusoids to the segments 
%                  without any Gaussian weights
%                  Default is 0.0125 with L1 = 6 cycles and L0 = 0. Thus,
%                  the sinusoids will be adapted to the Morlet wavelet
%                  in the least-squares sense via Gaussian values. 
%                  Note: The Gaussian values may not be attenuated to zero 
%                        at both ends of the windows for certain L1 and L0
% 15) 'display'  - Numeric. Diplays the spectrogram: 
%                  0 does not display 
%                  1 displays the normalized spectrogram
%                  2 displays the statistically significant spectrogram
%                  3 displays the amplitude spectrogram
%
% Outputs:
% 1) spectrogram        - The normalized least-squares wavelet spectrogram
% 2) stoch_surf         - The stochastic surface at (1-alpha) confidence 
%                         level obtained for the normalized spectrogram
% 3) amp_spectrogram    - The amplitude spectrogram

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
    tt = varargin{tt_ind+1}; check = true;
    Stt = size(tt); if Stt(1) == 1; tt = tt'; end
else
    tt = t; check = false;
end
Ltt = length(tt);
if (t(Lt) < tt(Ltt)) || (tt(1) < t(1))
    error('The range of time must be within the range of t')
end
%--------- Sampling Rate --------------------------------------------------
% M is a simple estimate for sampling rate. Note: this estimate may not be 
% valid for certain segments in an unequally spaced time series that may 
% cause aliasing when selecting the range of frequencies based on it
rate_ind = find(strcmpi(varargin,'rate'));
if ~isempty(rate_ind); M = varargin{rate_ind+1};
else; M = (Ltt-1)/(tt(Ltt)-tt(1)); end    %A simple estimate for sampling rate
%--------- Set of Cyclic Frequencies --------------------------------------
Omega_ind = find(strcmpi(varargin,'Omega'));
if ~isempty(Omega_ind)
    Omega = sort(varargin{Omega_ind+1});
    if (~isvector(Omega)) || ischar(Omega)
        error('Cyclic frequencies must be entered as a vector')
    end
else       
    if ~isequal(floor(M/2),M/2); Omega = 1:floor(M/2)-1;
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
if ~isempty(freq_ind); freq = sort(varargin{freq_ind+1});
else; freq = []; end
%--------- Morlet Coefficient ---------------------------------------------
morlet_ind = find(strcmpi(varargin,'morlet'));
if ~isempty(morlet_ind); morlet = varargin{morlet_ind+1};
else; morlet = 0.0125; end
MorlW = 1;
%--------- Window Dilation (L1) -------------------------------------------
L1_ind = find(strcmpi(varargin,'L1'));
if ~isempty(L1_ind); L1 = varargin{L1_ind+1};
elseif morlet; L1 = 6; 
else; L1 = 2; 
end
%--------- Additional Data Points (L0) ------------------------------------
L0_ind = find(strcmpi(varargin,'L0'));
if ~isempty(L0_ind); L0 = varargin{L0_ind+1};
elseif morlet; L0 = 0; 
else; L0 = 20; 
end
%--------- Display --------------------------------------------------------
display_ind = find(strcmpi(varargin,'display'));
if ~isempty(display_ind)
    display = varargin{display_ind+1};
    if ~any(display == [0,1,2,3])
        error('Valid options for the display are 0, 1, 2, 3')
    end
else
    display = 0;
end
%==========================================================================
LOm = length(Omega); 
spectrogram = zeros(LOm,Ltt); amp_spectrogram = zeros(LOm,Ltt);
stoch_surf = zeros(LOm,Ltt); 
[nrow,ncol] = size(P);
start_time = tt(1);
tt = tt - start_time;  % For the sake of computational efficiency 
t  = t  - start_time;  % For the sake of computational efficiency
%--------------------------------------------------------------------------
% Remove the trend with multiple pieces first using LSSA if 'ind' is given 
if ~isempty(ind)
    [~, ~, ~, res, ~, ~, ~] = LSSA(t,f, 'P',P, 'ind', ind, ...
                              'level',alpha, 'trend',trend, 'slope',slope);
else
    res = f;
end
%--------------------------------------------------------------------------
for k = 1:LOm                %loop over the frequencies
     L_k = L0 + floor((L1*M)/Omega(k));      
     if (mod(L_k,2)==0); L_k = L_k + 1; end
     for j = 1:Ltt           %loop over the times
        if j >= (L_k+1)/2 && j <= Ltt-(L_k-1)/2    % non-marginal windows
            lb = j - (L_k-1)/2;  ub = j + (L_k-1)/2;
        elseif j < (L_k+1)/2 && j <= Ltt-(L_k-1)/2 % left marginal windows
            lb = 1; ub = j + (L_k-1)/2;
        elseif j > Ltt-(L_k-1)/2 && j >= (L_k+1)/2 % right marginal windows
            lb = j - (L_k-1)/2; ub = Ltt;
        else                                       % for the entire series
            lb = 1; ub = Ltt;
        end 
        if check
            lb = find (t >= tt(lb),1 ,'first');
            ub = find (t <= tt(ub),1 ,'last');
        end
        ty = t(lb:ub);                             % time segment
        y  = res(lb:ub);                           % time series segment
        % The Gaussian weights based on the morlet value to define Py
        if morlet
            MorlW = exp(-morlet*((2*pi*Omega(k))^2)*...
                                    (ty-tt(j)*ones(1,length(ty))').^2); 
        end
        if ncol > 1 && nrow > 1 && morlet 
            Py = diag(P(lb:ub,lb:ub)).*MorlW;  % extract submatrix 
        elseif ncol > 1 && nrow > 1 && ~morlet 
            Py = P(lb:ub,lb:ub);               % extract submatrix
        elseif ncol == 1 && nrow > 1
            Py = P(lb:ub).*MorlW;              % extract subvector
        else               
            Py = MorlW;
        end 
        %----- Implement LSSA to segments ty, y, Py -----------------------
        [spectrum, CritVal, CScoeff, ~, ~, ~, ~] = LSSA(ty,y,'P',Py,...
                        'Omega',Omega(k), 'level',alpha, 'trend',trend, ...
                        'slope',slope, 'freq',freq);      
        spectrogram(k,j) = spectrum;  % spectral peak for (Omega(k), t(j))
        stoch_surf(k,j) = CritVal;    % critical value for (Omega(k), t(j))
        amp_spectrogram(k,j) = sqrt(CScoeff(1)^2 + CScoeff(2)^2);
     end
end
%--------------------------------------------------------------------------
if display
    PlotLSWS(t+start_time, tt+start_time, f, Omega, spectrogram, ...
                                      stoch_surf, amp_spectrogram, display)
end
%==========================================================================
function PlotLSWS(t, tt, f, Omega, spect, stoch_surf, amp_spect, display)
%% Display Spectrograms
%
% function PlotLSWS(t, tt, f, Omega, spect, stoch_surf, amp_spect, display)
%
% Summary: This function displays spectrograms
%
% Inputs:
% 1) t          - Vector of order n. Times of observations or measurements 
% 2) tt         - Vector. Equally spaced times for regularization purposes
% 3) f          - Vector of order n. The time series values
% 4) Omega      - Vector. The cyclic frequencies 
% 5) spect      - The normalized least-squares wavelet spectrogram
% 6) stoch_surf - The stochastic surface at (1-alpha) confidence level
%                 obtained for the normalized spectrogram
% 7) amp_spect  - The amplitude spectrogram
% 8) 'display'  - Numeric. Diplays the spectrogram: 
%                  1 displays the normalized spectrogram
%                  2 displays the statistically significant spectrogram
%                  3 displays the amplitude spectrogram

%==========================================================================
%% Check the input arguments
figure;
subplot(2,1,1);
plot(t, f,':d','color', [0,0.5,0],'MarkerFaceColor',[0,0.5,0],...
    'LineWidth', 0.1,'MarkerSize',4)
title('Time Series Decomposition via LSWA')
ylabel('Value')
vs = 0.1 * (max(f)-min(f));
if vs == 0; vs = 0.01; end
axis([tt(1)  tt(length(tt))  min(f) - vs   max(f) + vs])
grid on
set(gca,'Position',[0.1 0.7 0.72 0.25]); %[X Y W H]
set(gca,'XTickLabel',[]);   %remove x-axis tick labels 
set(gca, 'fontsize',10) ;
subplot(2,1,2);
if display == 1
    pcolor(tt, Omega, spect*100);
    % surf(tt, freq, spect*100);
    cb = colorbar; ylabel(cb, 'Percentage Variance');
elseif display == 2
    spect(spect < stoch_surf)= nan;
    pcolor(tt, Omega, spect*100);
    % surf(tt, freq, spect*100);
    cb = colorbar; ylabel(cb, 'Percentage Variance');
elseif display == 3
    pcolor(tt, Omega, amp_spect);
    % surf(tt, freq, amp_spect);
    cb = colorbar; ylabel(cb, 'Amplitude');
end
view(2); colormap jet; shading flat
xlabel('Time'); ylabel('Cyclic frequency');
set(gca, 'TickDir', 'out');
set(gca,'Position',[0.1 0.12 0.72 0.58]); %[X Y W H]
%--------------Save the figure --------------------------------------------
set(gcf, 'Color', 'w'); width = 7; height = 4;
set(gcf, 'PaperUnits', 'inches','PaperPosition',[0 0 width height]); 
%print( '-dtiff',  'LSWS', '-r300');

