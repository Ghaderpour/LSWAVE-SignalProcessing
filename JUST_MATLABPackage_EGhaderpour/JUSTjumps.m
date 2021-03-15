function LocIndMagDir = JUSTjumps(t,f,varargin)
%% Jumps Upon Spectrum and Trend (JUST)
%
% function LocIndMagDir = JUSTjumps(t,f,varargin)
%
% Summary: This function calls the AllJumps function (the windowing or
%          segmentation approach). For each segment within a translating 
%          window, the AllJumps function calls the JumpDetect function 
%          (the sequential approach) to find a jump candidate along with 
%          its magnitude and direction. 
%          Among the detected jumps within the translating windows,
%          certain jumps will be selected based on the following criteria:
%             1) Maximum occurrence 
%             2) Minimum distance to the window locations
%             3) Magnitudes and Directions
%          Thresholds for the minimum time difference between
%          the potential jumps and magnitude and direction will be applied 
%          at this stage.
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
% 4) 'size'      - Numeric. The window or segment size (R)
%                  Default is the average sampling rate tripled (3M)
% 5) 'step'      - Numeric. The translation step delta that is less than R 
%                  Default is the average sampling rate (M)
% 6) 'season'    - String. Three options: 
%                          ALLSSA season-trend model ('ALLSSA')
%                          OLS season-trend model    ('OLS')
%                          Only the trend model      ('none')
% 7) 'Omega'     - Vector. The initial cyclic frequencies 
%                  Default is Omega = [1,2,3,4] for OLS and 
%                  Omega = 0.8:0.2:3.8 for ALLSSA as the initial set
% 8) 'level'     - Numeric. The significance level, usually 0.01 or 0.05
%                  Default is 0.01. Applies only for ALLSSA not OLS
% 9) 'mag_th'    - Numeric. The threshold for the absolute value 
%                  of jump magnitude. Default is 0.05 (|Mag| >= 0.05)
% 10)'dir_th'    - Numeric. The threshold for the absolute value
%                  of jump direction. Default is 0.01 (|Dir| >= 0.01)
% 11)'jump_th'   - Numeric. The minimum time difference between jumps
%                  Default is 0.75 (e.g.,jumps are at least 9 months apart)
% 12)'save'      - Logical. Save the detected jumps and their attributes
% 
% Output:
% LocIndMagDir   - A matrix with four columns
%                  Each row contains the attributes of a potential jump:
%                       1st column: Jump location
%                       2nd column: Jump index
%                       3rd column: Estimated magnitude of the jump
%                       4th column: Estimated direction of the jump

%--------------------------------------------------------------------------
% Author: Ebrahim Ghaderpour 
% Email:  ebrahim.ghaderpour@ucalgary.ca
% Copyright (c) 2021
%==========================================================================
%% Check the input arguments
Lt = length(t);        % Number of observations or measurements
start_time = t(1);
t = t - start_time;  % For the sake of computational efficiency 
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
else
    P = 1;
end
%--------- Window Size ----------------------------------------------------
size_ind = find(strcmpi(varargin,'size'));
M = floor(Lt/(t(Lt)-t(1)));      % A simple estimate for sampling rate
if ~isempty(size_ind)
    R = varargin{size_ind+1};
    if R > Lt || R < M
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
    else; Omega = 1:4; end
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
%--------- Magnitude Threshold --------------------------------------------
mag_ind = find(strcmpi(varargin,'mag_th'));
if ~isempty(mag_ind); mag_th = varargin{mag_ind+1};
else; mag_th = 0.05; end
%--------- Direction Threshold --------------------------------------------
dir_ind = find(strcmpi(varargin,'dir_th'));
if ~isempty(dir_ind); dir_th = varargin{dir_ind+1};
else; dir_th = 0.01; end
%--------- Minimum Time Difference Between Potential Jumps ----------------
jump_ind = find(strcmpi(varargin,'jump_th'));
if ~isempty(jump_ind); jump_th = varargin{jump_ind+1};
else; jump_th = 0.75; end
%--------- Save -----------------------------------------------------------
save_ind = find(strcmpi(varargin,'save'));
if ~isempty(save_ind); save = varargin{save_ind+1};
else; save = false; end
%==========================================================================
IndMagDirDis = AllJumps(t, f, 'P',P, 'size',R, 'step',delta, ...
                         'Omega', Omega, 'season', season, 'level', alpha);
t = t + start_time;                     
[Occ, JumpInd] = hist(IndMagDirDis(:,1), IndMagDirDis(:,1));
NumJ = length(JumpInd);
IndOccMagDirDis = zeros(1, 5);
val = 0; kk = 1;
% The following loop will select a window whose location is closest to the 
% jump location with occurence greater than one in order to assign
% the attributes of the jump detected within that window.
for k = 1:NumJ
    if Occ(k)>0
        MinDist = sortrows(IndMagDirDis(val+1 : val+Occ(k),:),4);
        Mag =  round(MinDist(1,2),2);
        Dir =  round(MinDist(1,3),2);
        Dis =  floor(MinDist(1,4));
        val =  val + Occ(k);
        IndOccMagDirDis(kk,:) = [JumpInd(k),Occ(k),Mag,Dir,Dis];
        kk = kk + 1;
    end
end
if save
    results(:,1) = t(IndOccMagDirDis(:,1)); 
    results(:,2:6) = IndOccMagDirDis;
end
LJ = length(IndOccMagDirDis(:,1));
IndOccMagDirDis(:,5) = floor(IndOccMagDirDis(:,5)/(R/6));
% The following while loop will select the optimal jump locations whose
% time differences are greater than or equal to jumps_th
r = 2;
while r <= LJ
    if t(IndOccMagDirDis(r,1))-t(IndOccMagDirDis(r-1,1)) < jump_th
        group = IndOccMagDirDis(r-1:r,:);
        % Sort the jump occurrences descendingly, then sort ascendingly the 
        % scaled jump distances to the window locations where they were 
        % estimated, then sort abs magnitudes descendingly 
        group = sortrows(group,[2 5 3],{'descend' 'ascend' 'descend'},...
                                                 'ComparisonMethod','abs');
        IndOccMagDirDis (r-1,:) = group(1, 1:5);
        IndOccMagDirDis (r,:) = [];
        LJ = LJ - 1;
    else
        r = r + 1;
    end
end
% The following while loop applies the magnitude and direction thresholds
r = 1;
while r <= LJ
    if abs(IndOccMagDirDis(r,3)) < abs(mag_th) && ...
            abs(IndOccMagDirDis(r,4)) < abs(dir_th)
        IndOccMagDirDis (r,:) = [];
        LJ = LJ - 1;
    else
        r = r + 1;
    end
end
LocIndMagDir(:,1)   = t(IndOccMagDirDis(:,1));
LocIndMagDir(:,2)   = IndOccMagDirDis(:,1); 
LocIndMagDir(:,3:4) = IndOccMagDirDis(:,3:4);
if save
    filter = {'*.txt';'*.dat'};
    [file, path] = uiputfile(filter,'File Selection','JumpResults');
    FullFileName = fullfile(path,file);
    fileID = fopen(FullFileName,'w');
    fprintf(fileID,'Jump Location    Jump Index    Occurrence    ');     
    fprintf(fileID,'Magnitude    Direction    Distance\n');
    for k = 1:kk-1
        fprintf(fileID, '%10.3f %12g %12g %15.2f %12.2f %10g\n',...
                                                             results(k,:));
    end
    fprintf(fileID,'\nSelected jumps after applying the thresholds:\n');
    fprintf(fileID,'Jump Location    Jump Index    Magnitude    Direction\n');
    for k = 1:LJ
        fprintf(fileID, '%10.3f %12g %14.2f %12.2f\n', LocIndMagDir(k,:));
    end
    fclose(fileID);
end
%==========================================================================                         
function IndMagDirDis = AllJumps(t, f, varargin)
%% The Segmentation or Windowing Approach in JUST
%
% IndMagDirDis = AllJumps(t, f, varargin)
%
% Summary: A window of fixed size R translates by delta steps (observations)
%          over time. For each segment within a translating window, 
%          the JumpDetect function (the sequential approach) will be called 
%          to find a jump candidate along with its magnitude and direction. 
%                                            
% Inputs:
% 1) t           - Vector of order n. Times of observations or measurements 
% 2) f           - Vector of order n. The time series values
% 3) 'P'         - Vector or Matrix of order n. The weights associated with
%                  the observations. Default is none
% 4) 'size'      - Numeric. The window or segment size (R)
%                  Default is the average sampling rate tripled (3M)
% 5) 'step'      - Numeric. The translation step delta that is less than R 
%                  Default is the average sampling rate (M)
% 6) 'season'    - String. Three options: 
%                          ALLSSA season-trend model ('ALLSSA')
%                          OLS season-trend model    ('OLS')
%                          Only the trend model      ('none')
% 7) 'Omega'     - Vector. The initial cyclic frequencies 
%                  Default is Omega = [1,2,3,4]   for OLS and 
%                  Omega = 0.8:0.2:3.8 for ALLSSA as the initial set
% 8) 'level'     - Numeric. The significance level, usually 0.01 or 0.05
%                  Default is 0.01. Applies only for ALLSSA not OLS
% 
% Output:
% IndMagDirDis   - A matrix with four columns
%                  Each row contains the attributes of a potential jump:
%                       1st column: Time index of the jump
%                       2nd column: Estimated magnitude of the jump
%                       3rd column: Estimated direction of the jump
%                       4th column: Distance between the time index of the  
%                                   jump and the time index of its 
%                                   corresponding window location

%==========================================================================
%% Check the input arguments
Lt = length(t);
%--------- Weights (Vector or Matrix)--------------------------------------
P_ind = find(strcmpi(varargin,'P'));
if ~isempty(P_ind)
    P = varargin{P_ind+1};
    SP = size(P);
    if ~((SP(1) == Lt && SP(2) == Lt) || (SP(1) == 1 && SP(2) == Lt) || ...
            (SP(1) == Lt && SP(2) == 1) || (SP(1) == 1 && SP(2) == 1))
        error('P must be a vector or a matrix of order n')
    end
else
    P = 1;
end
%--------- Window Size ----------------------------------------------------
size_ind = find(strcmpi(varargin,'size'));
M = floor(Lt/(t(Lt)-t(1)));   % A simple estimate for sampling rate
if ~isempty(size_ind)
    R = varargin{size_ind+1};
    if R > Lt || R < M
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
    else; Omega = 1:4; end
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
%--------------------------------------------------------------------------                                     
R1 = R - 1;
inda = 1; indb = 1 + R1;
IndMagDirDis = zeros(1,4);
[nrow, ncol] = size(P);
check = true;  k = 1;
while indb <= Lt + R1 && check
    if indb > Lt  % The condition for the last segment to also have size R
        inda = Lt - R1; indb = Lt; check = false;
    end
    if ncol==1 && nrow==1; P0 = 1;              % If P is a scale  
    elseif nrow>1 && ncol==1; P0 = P(inda:indb);% If P is a column vector        
    else; P0 = P(inda:indb, inda:indb);         % If P is a square matrix     
    end
    t0=t(inda:indb);
    f0=f(inda:indb);
    [Ind, Mag, Dir] = JumpDetect(t0, f0, 'P', P0, ...
                         'Omega', Omega, 'season', season, 'level', alpha);
    IndMagDirDis(k,:) = [Ind+inda-1, Mag, Dir, abs(Ind-1-R1/2)];
    k = k + 1; inda = inda + delta; indb = inda + R1;
end
IndMagDirDis = sortrows(IndMagDirDis, 1);
%==========================================================================
function [Ind, Mag, Dir] = JumpDetect(t, f, varargin) 
%% The Sequential Approach in JUST
%
% function [Ind, Mag, Dir] = JumpDetect(t, f, varargin) 
%
% Summary: A linear trend with two pieces will be fitted sequentially 
%          along with the seasonal component based on the user input 
%          (i.e., 'ALLSSA', 'OLS', 'none') to estimate the potential jumps 
%          in the trend component of segment f. 
%                                            
% Inputs:
% 1) t           - Vector of order n. Times of observations or measurements 
% 2) f           - Vector of order n. The time series values
% 3) 'P'         - Vector or Matrix of order n. The weights associated with
%                  the observations. Default is none
% 4) 'season'    - String. Three options: 
%                             ALLSSA season-trend model ('ALLSSA')
%                             OLS season-trend model    ('OLS')
%                             Only the trend model      ('none')
% 5) 'Omega'     - Vector. The initial cyclic frequencies 
%                  Default is Omega = [1,2,3,4]   for OLS and 
%                  Omega = 0.8:0.2:3.8 for ALLSSA as the initial set
% 6) 'level'     - Numeric. The significance level, usually 0.01 or 0.05
%                  Default is 0.01. Applies only for ALLSSA not OLS
% 
% Outputs:
% 1) Ind         - Time index of the potential jump (jump index)
% 2) Mag         - Estimated magnitude of the potential jump
% 3) Dir         - Estimated direction of the potential jump

%==========================================================================
%% Check the input arguments
Lt = length(t);
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
else
    P = 1;
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
    else; Omega = 1:4; end
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
%--------------------------------------------------------------------------
Attributes = zeros(Lt-5,4);
k = 1;
for Ind = 4:Lt-2
    if strcmpi(season,'ALLSSA')
        [tr_coeff, ~, ~, norm_res, ~, ~] = ALLSSA(t, f, 'P', P, ...
                            'Omega', Omega, 'ind', Ind, 'level', alpha);                     
    elseif strcmpi(season,'OLS')                  
        [~, ~, ~, ~, norm_res, coeff, ~] = LSSA(t, f, 'P',P, ...
                            'Omega',[], 'ind', Ind, 'freq',Omega); 
        tr_coeff = coeff(1 : length(coeff)-2*length(Omega));
    else
        [~, ~, ~, ~, norm_res, coeff, ~] = LSSA(t, f, 'P',P, ...
                                                'Omega',[], 'ind', Ind);
        tr_coeff = coeff;
    end                  
    Mag = (tr_coeff(2)-tr_coeff(1))+(tr_coeff(4)-tr_coeff(3))*t(Ind);
    Dir = tr_coeff(4)-tr_coeff(3); 
    Attributes(k,:) = [Ind, Mag, Dir, norm_res];
    k = k + 1;
end
Attributes = sortrows(Attributes, 4);
Ind = Attributes(1,1); Mag = Attributes(1,2); Dir = Attributes(1,3);