function [cross_spectrum, CrossCritVal, phase_diff]=LSCSA(t1, t2, f1, f2, ...
    P_f1, P_f2, Omega, Q1, Q2, IndexDatum1, IndexDatum2, alpha, normalized)
%% The Least-Squares Cross Spectral Analysis (LSCSA)
% Author: Ebrahim Ghaderpour www.ghader.org
% Copyright (c) 2019
%
% Please acknowledge the use of this software in any publications:
% "Ghaderpour, E., Ince, E.S., Pagiatakis, S.D., 2018a. 
% Least-squares cross-wavelet analysis and its applications in geophysical 
% time series. Journal of Geodesy 92, 1223–1236."
%
% Inputs:
% 1) All the LSSA inputs for the first series
% 2) All the LSSA inputs for the second series
%    Note: The set of cyclic frequencies or wavenumbers (Omega), 
%          significance level (alpha), and the normalized mode (normalized) 
%          are the same for both series
%
% Outputs:
% 1) cross_spectrum    - The least-squares cross-spectrum (LSCS)
% 2) CrossCritVal      - The critical value at (1-alpha) confidence level
%                        NOTE: This value is valid if normalized = true
% 3) phase_diff        - The phase differences in radian in the LSCSA
%--------------------------------------------------------------------------
%% Check the input arguments
if nargin < 4
    error('Too few input arguments');
elseif nargin==4
    P_f1=1; P_f2=1; 
    M1=(length(t1)-1)/(t1(length(t1))-t1(1));
    M2=(length(t2)-1)/(t2(length(t2))-t2(1));
    Omega=1:0.5*min(M1, M2); 
    Q1=[1,0,0,0];  Q2=[1,0,0,0];
    IndexDatum1=[]; IndexDatum2=[];
    alpha= 0.01; normalized=true;
elseif nargin==6
    M1=(length(t1)-1)/(t1(length(t1))-t1(1));
    M2=(length(t2)-1)/(t2(length(t2))-t2(1));
    Omega=1:0.5*min(M1, M2); 
    Q1=[1,0,0,0];  Q2=[1,0,0,0];
    IndexDatum1=[]; IndexDatum2=[];
    alpha= 0.01; normalized=true;
elseif 6 < nargin && nargin < 13
    error('Too few input arguments');
elseif nargin > 13
    error('Too many input arguments');
end
%--------------------------------------------------------------------------
% Run the LSSA code for the first series
[spectrum1, chats1, ~, ~, ~, ~, nc1]=...
    LSSA(t1, f1, P_f1, Omega, Q1, IndexDatum1, alpha, normalized);
% Run the LSSA code for the second series
[spectrum2, chats2, ~, ~, ~, ~, nc2]=...
    LSSA(t2, f2, P_f2, Omega, Q2, IndexDatum2, alpha, normalized);

% Calculate the cross-spectrum
cross_spectrum=spectrum1.*spectrum2;
%------------ Critical Variance for the LSCSA -----------------------------
LOm=length(Omega);
Lt1=length(t1);
Lt2=length(t2);
if alpha==0.1
    CV90=load('CriticalPV90.mat', 'CV90cl'); CV=CV90.CV90cl;
elseif alpha==0.05
    CV95=load('CriticalPV95.mat', 'CV95cl'); CV=CV95.CV95cl;
elseif alpha==0.01
    CV99=load('CriticalPV99.mat', 'CV99cl'); CV=CV99.CV99cl;
end
dof1=Lt1-nc1-2;    % the degrees of freedom of the first series
dof2=Lt2-nc2-2;    % the degrees of freedom of the second series
if dof1<1 || dof2<1
    Cross_CritVal=100;
elseif (dof1>300 || dof2>300) % CV beyond 300 data points is not calculated
    Cross_CritVal=0;
else
    Cross_CritVal=CV(dof1,dof2);
end
CrossCritVal=ones(1,LOm)*Cross_CritVal;
%--------------- Calculate phase differences in the LSCSA  ----------------
Lchats1=length(chats1);
ch1_cos=chats1(1:2:Lchats1);
ch1_sin=chats1(2:2:Lchats1);
amp1=sqrt(ch1_cos.^2+ch1_sin.^2);

ch2_cos=chats2(1:2:Lchats1);
ch2_sin=chats2(2:2:Lchats1);
amp2=sqrt(ch2_cos.^2+ch2_sin.^2);
     
phase_diff=zeros; 
for i=1:LOm
   phase_diff(i)= 2*(atan((amp1(i)-ch1_sin(i))./(ch1_cos(i))))...
       -2*(atan((amp2(i)-ch2_sin(i))./(ch2_cos(i)))); 
   if -2*pi<= phase_diff(i) && phase_diff(i)<-pi
       phase_diff(i)= phase_diff(i)+2*pi;
   elseif pi<= phase_diff(i) && phase_diff(i)<2*pi
       phase_diff(i)= phase_diff(i)-2*pi;
   end
end