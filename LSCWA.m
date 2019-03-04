function [t_union, cross_spectrogram, cross_stoch_surf, phase_diff]=...
    LSCWA(t1, t2, f1, f2, P_f1, P_f2, Omega, M1, L11, L01, Q1, IndexDatum1, ...
               M2, L12, L02, Q2, IndexDatum2, alpha, MorlCoeff, normalized)
%% The Least-Squares Cross-Wavelet Analysis (LSCWA)
% Author: Ebrahim Ghaderpour www.ghader.org
% Copyright (c) 2019
%
% Please acknowledge the use of this software in any publications:
% "Ghaderpour, E., Ince, E.S., Pagiatakis, S.D., 2018a. 
% Least-squares cross-wavelet analysis and its applications in 
% geophysical time series. Journal of Geodesy 92, 1223–1236."
%
% Inputs:
% 1) All the LSWA inputs for the first series
% 2) All the LSWA inputs for the second series
% Note: The set of cyclic frequencies or wavenumbers (Omega), 
%       significance level (alpha), the Morlet coefficient (MorlCoeff), 
%       and the normalized mode (normalized) are the same for both series
%
% Outputs:
% 1) t_union              - The union of both time vectors 
% 2) cross_spectrogram    - The least-squares cross-wavelet spectrogram
% 3) cross_stoch_surf     - The LSCWA stochastic surface
%                           NOTE: the cross_stoch_surface is valid only 
%                                 when normalized = true 
% 4) phase_diff           - The phase differences in the LSCWA
% -------------------------------------------------------------------------
%% Check the input arguments
if nargin < 4
    error('Too few input arguments');
elseif nargin==4
    P_f1=1; P_f2=1; 
    M1=(length(t1)-1)/(t1(length(t1))-t1(1));
    M2=(length(t2)-1)/(t2(length(t2))-t2(1));
    Omega=1:0.5*min(M1, M2); 
    L11=2; L01=20;  % For the first series
    L12=2; L02=20;  % For the second series
    Q1=[1,0,0,0];  Q2=[1,0,0,0];
    IndexDatum1=[]; IndexDatum2=[];
    alpha= 0.01; MorlCoeff=0; normalized=true;
elseif nargin==6
    M1=(length(t1)-1)/(t1(length(t1))-t1(1));
    M2=(length(t2)-1)/(t2(length(t2))-t2(1));
    Omega=1:0.5*min(M1, M2); 
    L11=2; L01=20;  % For the first series
    L12=2; L02=20;  % For the second series
    Q1=[1,0,0,0];  Q2=[1,0,0,0];
    IndexDatum1=[]; IndexDatum2=[];
    alpha= 0.01; MorlCoeff=0; normalized=true;
elseif 6 < nargin && nargin < 20
    error('Too few input arguments');
elseif nargin > 20
    error('Too many input arguments');
end
%--------------------------------------------------------------------------
% Run the LSWA code for the first series
[spectrogram1, ~, w_size1, chat1_cos, chat1_sin, num_cons1]=LSWA(t1, f1,...
    P_f1, Omega, M1, L11, L01, Q1, IndexDatum1, alpha, MorlCoeff, normalized);

% Run the LSWA code for the second series
[spectrogram2, ~, w_size2, chat2_cos, chat2_sin, num_cons2]=LSWA(t2, f2,...
    P_f2, Omega, M2, L12, L02, Q2, IndexDatum2, alpha, MorlCoeff, normalized);


t1=round(t1,9); t2=round(t2,9);
n1=length(t1); % number of samples in the first series
n2=length(t2); % number of samples in the second series
j=1;  % time index for the union of the two time vectors 
j1=1; % time index for the first series 
j2=1; % time index for the second series
LOm=length(Omega); % Number of frequencies, the same in both spectrograms

% Estimated amplitudes for segments of the first series
amp1=sqrt(chat1_cos.^2+chat1_sin.^2);
% Estimated amplitudes for segments of the second series
amp2=sqrt(chat2_cos.^2+chat2_sin.^2);

% Take the union of the times as the time vector of the cross-spectrogram
t_union=union(t1,t2);
LU=length(t_union);

cross_spectrogram=zeros(LOm,LU); phase_diff=zeros(LOm,LU); 
dof1=zeros(LOm,LU); dof2=zeros(LOm,LU);% degrees of freedom (dof)
% Note that for a given frequency, dof for each time segment may vary, 
% depending on the datum shift indices entered by users
% -------------------------------------------------------------------------
% The following while loop will adjust the locations of the windows
% in both series to compute the cross-spectrogram based on the union of t1
% and t2. NOTE: The time vector for cross-spectrogram may be selected 
% to be equally spaced
while j1<=n1 || j2<=n2
    if j1>n1 && j2<=n2
        cross_spectrogram(:,j)=spectrogram1(:,n1).*spectrogram2(:,j2);
        dof1(:,j)=w_size1(:,n1)-num_cons1(:,n1)-2; 
        dof2(:,j)=w_size2(:,j2)-num_cons2(:,j2)-2;
        phase_diff(:,j)=2*(atan((amp1(:,n1)-chat1_sin(:,n1))./(chat1_cos(:,n1))))...
            -2*(atan((amp2(:,j2)-chat2_sin(:,j2))./(chat2_cos(:,j2))));
        j2=j2+1;
    elseif j1<=n1 && j2>n2 
        cross_spectrogram(:,j)=spectrogram1(:,j1).*spectrogram2(:,n2);
        dof1(:,j)=w_size1(:,j1)-num_cons1(:,j1)-2;
        dof2(:,j)=w_size2(:,n2)-num_cons2(:,n2)-2;
        phase_diff(:,j)=2*(atan((amp1(:,j1)-chat1_sin(:,j1))./(chat1_cos(:,j1))))...
            -2*(atan((amp2(:,n2)-chat2_sin(:,n2))./(chat2_cos(:,n2))));
        j1=j1+1;
    elseif j1<=n1 && j2<=n2
        if t1(j1)<t2(j2)
            if j2>1 && abs(t1(j1)-t2(j2))>abs(t1(j1)-t2(j2-1))
                 cross_spectrogram(:,j)=spectrogram1(:,j1).*spectrogram2(:,j2-1);
                 dof1(:,j)=w_size1(:,j1)-num_cons1(:,j1)-2;
                 dof2(:,j)=w_size2(:,j2-1)-num_cons2(:,j2-1)-2;
                 phase_diff(:,j)=2*(atan((amp1(:,j1)-chat1_sin(:,j1))./(chat1_cos(:,j1))))...
                     -2*(atan((amp2(:,j2-1)-chat2_sin(:,j2-1))./(chat2_cos(:,j2-1))));
            else
                 cross_spectrogram(:,j)=spectrogram1(:,j1).*spectrogram2(:,j2);
                 dof1(:,j)=w_size1(:,j1)-num_cons1(:,j1)-2;
                 dof2(:,j)=w_size2(:,j2)-num_cons2(:,j2)-2;
                 phase_diff(:,j)=2*(atan((amp1(:,j1)-chat1_sin(:,j1))./(chat1_cos(:,j1))))...
                     -2*(atan((amp2(:,j2)-chat2_sin(:,j2))./(chat2_cos(:,j2))));
            end
            j1=j1+1;
        elseif t1(j1)>t2(j2)
            if j1>1 && abs(t1(j1)-t2(j2))>abs(t1(j1-1)-t2(j2))
                 cross_spectrogram(:,j)=spectrogram1(:,j1-1).*spectrogram2(:,j2);
                 dof1(:,j)=w_size1(:,j1-1)-num_cons1(:,j1-1)-2;
                 dof2(:,j)=w_size2(:,j2)-num_cons2(:,j2)-2;
                 phase_diff(:,j)=2*(atan((amp1(:,j1-1)-chat1_sin(:,j1-1))./(chat1_cos(:,j1-1))))...
                     -2*(atan((amp2(:,j2)-chat2_sin(:,j2))./(chat2_cos(:,j2))));
            else
                 cross_spectrogram(:,j)=spectrogram1(:,j1).*spectrogram2(:,j2);
                 dof1(:,j)=w_size1(:,j1)-num_cons1(:,j1)-2;
                 dof2(:,j)=w_size2(:,j2)-num_cons2(:,j2)-2;
                 phase_diff(:,j)=2*(atan((amp1(:,j1)-chat1_sin(:,j1))./(chat1_cos(:,j1))))...
                     -2*(atan((amp2(:,j2)-chat2_sin(:,j2))./(chat2_cos(:,j2))));
            end
            j2=j2+1;
        elseif t1(j1)==t2(j2)
            cross_spectrogram(:,j)=spectrogram1(:,j1).*spectrogram2(:,j2);
            dof1(:,j)=w_size1(:,j1)-num_cons1(:,j1)-2;
            dof2(:,j)=w_size2(:,j2)-num_cons2(:,j2)-2;
            phase_diff(:,j)=2*(atan((amp1(:,j1)-chat1_sin(:,j1))./(chat1_cos(:,j1))))...
            -2*(atan((amp2(:,j2)-chat2_sin(:,j2))./(chat2_cos(:,j2))));
            j1=j1+1;
            j2=j2+1;
        end
    end
    j=j+1;
end
% -------------------------------------------------------------------------
% Calculate the phase differences between the two series in radians
for r=1:LOm
    for c=1:LU
        if -2*pi<= phase_diff(r,c) && phase_diff(r,c)<-pi
            phase_diff(r,c)= phase_diff(r,c)+2*pi;
        elseif pi<= phase_diff(r,c) && phase_diff(r,c)<2*pi
            phase_diff(r,c)= phase_diff(r,c)-2*pi;
        end
    end
end
% -------------------------------------------------------------------------
% Calculate the stochastic surface at (1-alpha) confidence level
cross_stoch_surf=zeros(LOm,LU);

if alpha==0.1
    CV90=load('CriticalPV90.mat', 'CV90cl'); CV=CV90.CV90cl;
elseif alpha==0.05
    CV95=load('CriticalPV95.mat', 'CV95cl'); CV=CV95.CV95cl;
elseif alpha==0.01
    CV99=load('CriticalPV99.mat', 'CV99cl'); CV=CV99.CV99cl;
end 
for i=1:LOm
    for j=1:LU
        d1=dof1(i,j);
        d2=dof2(i,j);
        if d1<1 || d2<1
            cross_stoch_surf(i,j)=100;
        elseif (d1>300 || d2>300) 
            % CV beyond 300 samples is not calculated (approximately zero)
            cross_stoch_surf(i,j)=0;
        else
            cross_stoch_surf(i,j)=CV(d1,d2);
        end
        % Set  -pi<= phase_diff < pi
        if -2*pi<= phase_diff(i,j) && phase_diff(i,j)<-pi
            phase_diff(i,j)=phase_diff(i,j)+2*pi;
        elseif pi<= phase_diff(i,j) && phase_diff(i,j)<2*pi
            phase_diff(i,j)=phase_diff(i,j)-2*pi;
        end     
    end
end
