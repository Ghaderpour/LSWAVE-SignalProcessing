% A GUI for time/data series analysis package (LSWAVE) 
% Author: Ebrahim Ghaderpour www.ghader.org
function varargout = LSWAVE_GUI(varargin)
%% LSWA MATLAB code for LSWA.fig
%      LSWA, by itself, creates a new LSWA or raises the existing
%      singleton*.
%
%      H = LSWA returns the handle to a new LSWA or the handle to
%      the existing singleton*.
%
%      LSWA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LSWA.SR with the given input arguments.
%
%      LSWA('Property','Value',...) creates a new LSWA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LSWA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LSWA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LSWA

% Last Modified by GUIDE v2.5 19-Jan-2019 09:28:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LSWA_OpeningFcn, ...
                   'gui_OutputFcn',  @LSWA_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before LSWA is made visible.
function LSWA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LSWA (see VARARGIN)

% Choose default command line output for LSWA
handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes LSWA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LSWA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SR_Callback(hObject, eventdata, handles)
% hObject    handle to SR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SR as text
%        str2double(get(hObject,'String')) returns contents of SR as a double


% --- Executes during object creation, after setting all properties.
function SR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function L1_Callback(hObject, eventdata, handles)
% hObject    handle to L1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L1 as text
%        str2double(get(hObject,'String')) returns contents of L1 as a double


% --- Executes during object creation, after setting all properties.
function L1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function L0_Callback(hObject, eventdata, handles)
% hObject    handle to L0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L0 as text
%        str2double(get(hObject,'String')) returns contents of L0 as a double


% --- Executes during object creation, after setting all properties.
function L0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function force_Callback(hObject, eventdata, handles)
% hObject    handle to force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of force as text
%        str2double(get(hObject,'String')) returns contents of force as a double

function DatumShifts1_Callback(hObject, eventdata, handles)
% hObject    handle to DatumShifts1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DatumShifts1 as text
%        str2double(get(hObject,'String')) returns contents of DatumShifts1 as a double

% --- Executes during object creation, after setting all properties.
function DatumShifts1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DatumShifts1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DatumShifts2_Callback(hObject, eventdata, handles)
% hObject    handle to DatumShifts2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DatumShifts2 as text
%        str2double(get(hObject,'String')) returns contents of DatumShifts2 as a double


% --- Executes during object creation, after setting all properties.
function DatumShifts2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DatumShifts2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function force_CreateFcn(hObject, eventdata, handles)
% hObject    handle to force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in datum.
function datum_Callback(hObject, eventdata, handles)
% hObject    handle to datum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of datum


% --- Executes on button press in lineartrend.
function lineartrend_Callback(hObject, eventdata, handles)
% hObject    handle to lineartrend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lineartrend


% --- Executes on button press in quadratictrend.
function quadratictrend_Callback(hObject, eventdata, handles)
% hObject    handle to quadratictrend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of quadratictrend


% --- Executes on button press in cubictrend.
function cubictrend_Callback(hObject, eventdata, handles)
% hObject    handle to cubictrend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cubictrend



% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Inputs_Callback(hObject, eventdata, handles)
% hObject    handle to Inputs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenFiles_Callback(hObject, eventdata, handles)
% hObject    handle to OpenFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%cla(handles.axes1,'reset');
%set(handles.axes1,'visible','off');
%cla(handles.axes2,'reset');
%set(handles.axes2,'visible','off');
set(handles.uipanelforce,'visible','off');
set(handles.uipanelfreq,'visible','off');
set(handles.uipanellevel,'visible','off');
set(handles.uipanelwindow,'visible','off');
set(handles.uipanelaxes,'visible','off');
set(handles.uipanelarrow,'visible','off');
guidata(hObject, handles);

% --------------------------------------------------------------------
function Exit_Callback(hObject, eventdata, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);

% --------------------------------------------------------------------
function Window_Callback(hObject, eventdata, handles)
% hObject    handle to Window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%cla(handles.axes1,'reset');
%set(handles.axes1,'visible','off');
set(handles.uipanelforce,'visible','off');
set(handles.uipanellevel,'visible','off');
set(handles.uipanelwindow,'visible','on');
set(handles.uipanelaxes,'visible','off');
set(handles.uipanelarrow,'visible','off');
set(handles.uipanelfreq,'visible','off');
guidata(hObject, handles);

% --------------------------------------------------------------------
function Slevel_Callback(hObject, eventdata, handles)
% hObject    handle to Slevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%cla(handles.axes1,'reset');
%set(handles.axes1,'visible','off');
set(handles.uipanelforce,'visible','off');
set(handles.uipanellevel,'visible','on');
set(handles.uipanelwindow,'visible','off');
set(handles.uipanelaxes,'visible','off');
set(handles.uipanelarrow,'visible','off');
set(handles.uipanelfreq,'visible','off');
guidata(hObject, handles);

% --- Executes on selection change in SignificanceLevel.
function SignificanceLevel_Callback(hObject, eventdata, handles)
% hObject    handle to SignificanceLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SignificanceLevel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SignificanceLevel


% --- Executes during object creation, after setting all properties.
function SignificanceLevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SignificanceLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Force_Callback(hObject, eventdata, handles)
% hObject    handle to Force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%cla(handles.axes1,'reset');
%set(handles.axes1,'visible','off');
set(handles.uipanelforce,'visible','on');
set(handles.uipanellevel,'visible','off');
set(handles.uipanelwindow,'visible','off');
set(handles.uipanelaxes,'visible','off');
set(handles.uipanelarrow,'visible','off');
set(handles.uipanelfreq,'visible','off');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function text15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function LSWA_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to LSWA_Spectrogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%cla(handles.axes1,'reset');
%axes(handles.axes1);
handles = guidata(hObject); 
if ~isfield(handles,'FirstTS')
    set(handles.Wait,'String','Please open a file containing time/data series');
    guidata(hObject, handles);
      pause(.5);
      return;
end
set(handles.uipanelforce,'visible','off');
set(handles.uipanellevel,'visible','off');
set(handles.uipanelwindow,'visible','off');
set(handles.uipanelaxes,'visible','off');
set(handles.uipanelarrow,'visible','off');
set(handles.uipanelfreq,'visible','off');
T_S=handles.FirstTS;
set(handles.Wait,'String','Please wait (it may take a while)...');
pause(0.05);
guidata(hObject, handles);
Omega=handles.Omega;
if contains(T_S,'.xlsx') || contains(T_S,'.xls')
    Data=xlsread(T_S);
else
    Data=load(T_S);
end
tt = Data(:,1);
t1=tt-Data(1,1);
f = Data(:,2);
C_f=zeros;
[~,n_c]=size(Data);
if  (n_c>3) 
    % formats the weight matrix (the inverse of the covariance matrix 
    % associated with the series)
     for i=1:length(t1)
        for j=1:length(t1)
            C_f(i,j)=Data(i,j+2);
        end
     end
     P_f=inv(C_f);
elseif n_c==3
     P_f =1./(Data(:,3).^2);
elseif n_c==2
     P_f =1;
end
%--------------------------------------------
M=str2double(get(handles.SR,'String'));
if isnan(M) || M<1 % Checks sampling rate M
    set(handles.SR,'String','');
    set(handles.Wait,'String','Error! Sampling rate must be a positive number');
    set(handles.uipanelwindow,'visible','on');
    return;
end
%--------------------------------------------
content=get(handles.SignificanceLevel,'Value');
if content==1
    alpha=0.01;
elseif content==2
    alpha=0.05;
else
    alpha=0.1;
end
%--------------------------------------------
L1=str2double(get(handles.L1,'String'));
if isnan(L1) || L1<0
    set(handles.L1,'String','');
    set(handles.Wait,'String',...
        'Error! L1 must be a small nonnegative integer');
    set(handles.uipanelwindow,'visible','on');
    return;
end
%--------------------------------------------
L0=str2double(get(handles.L0,'String'));
if isnan(L0) || L0<0 || rem(L0,1)~=0
    set(handles.L0,'String','');
    set(handles.Wait,'String',...
        'Error! L0 must be a nonnegative integer');
    set(handles.uipanelwindow,'visible','on');
    return;
end
%--------------------------------------------
% Adjust QQ, the vector of constituents of known forms
QQ(1)=get(handles.datum,'Value');
QQ(2)=get(handles.lineartrend,'Value');
QQ(3)=get(handles.quadratictrend,'Value');
QQ(4)=get(handles.cubictrend,'Value');
normalized=get(handles.Normalizedcheckbox,'Value');
force_sinusoids=get(handles.force,'String');
FS=str2num(force_sinusoids);
if isnan(FS)
    Q=QQ;
else
    Q=[QQ,FS];
end

IndexDatum1=get(handles.DatumShifts1,'String');
IndexDatum1=str2num(IndexDatum1);
if isnan(IndexDatum1)
   IndexDatum1=[];
end
Lt1=length(t1);
if any(floor(IndexDatum1)-IndexDatum1) || any(IndexDatum1>=Lt1) || any(IndexDatum1<1)
    set(handles.Wait,'String','Error! Datum shift indices must be natural numbers less than the size of series');
    set(handles.uipanelforce,'visible','on');
    return;
end
MorlCoeff=get(handles.MorletCheck,'Value');
MorlCoeff=MorlCoeff*0.0125;
%-----------Run the LSWA---------------------
[spectrogram, stoch_surf, window_size, chat1_cos, chat1_sin, ~]=...
    LSWA(t1, f, P_f, Omega, M, L1, L0, Q, IndexDatum1, alpha, MorlCoeff, normalized);
%--------------------------------------------
[~,v]=size(spectrogram); 
%dimension of LSWS (spec) for plotting, u the number of rows and v the number of columns
if v==1
    return;
end
 
set(handles.colorbarvariance,'state','off');
set(handles.colorbarstochastic,'state','off');
cla(handles.axes2,'reset');
set(handles.axes2,'visible','off');
cla(handles.axes3,'reset');
set(handles.axes3,'visible','off');
cla(handles.axes4,'reset');
set(handles.axes4,'visible','off');
axes(handles.axes2);
if normalized
    surf(t1,Omega,spectrogram*100);
else
    surf(t1,Omega,spectrogram);
end
view(2);
shading flat; 
axis tight;
colormap jet;
ylabel('Cyclic frequency','fontsize',10); 
xlabel('Time','fontsize',10); 
align([handles.axes1 handles.axes2],'VerticalAlignment','none', 'HorizontalAlignment','center');
set(handles.Wait,'String','Least-squares wavelet spectrogram');
pause(0.05);
handles=guidata(hObject);
handles.spec=spectrogram;
handles.specdB=zeros;
handles.stochastic=stoch_surf;
handles.t1=t1;
handles.window_size=window_size;
handles.chat1_cos=chat1_cos;
handles.chat1_sin=chat1_sin;
handles.cbar=0;
handles.cbar2=0;
handles.LSWS3D1=1; handles.LSWS3D2=1; handles.LSWS3D3=1;
handles.colorbar_VarordB=1;
handles.colorbar_check=1;
handles.color_bar=1;
guidata(hObject, handles);


% --------------------------------------------------------------------
function displays_Callback(hObject, eventdata, handles)
% hObject    handle to displays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tunit_Callback(hObject, eventdata, handles)
% hObject    handle to tunit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipanelforce,'visible','off');
set(handles.uipanellevel,'visible','off');
set(handles.uipanelwindow,'visible','off');
set(handles.uipanelaxes,'visible','on');
set(handles.uipanelarrow,'visible','off');
set(handles.uipanelfreq,'visible','off');
guidata(hObject, handles);

% --------------------------------------------------------------------
function dB_Callback(hObject, eventdata, handles)
% hObject    handle to dB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
if ~isfield(handles,'spec')
    cla(handles.axes2,'reset');
    set(handles.axes2,'visible','off');
    set(handles.Wait,'String','Please push LSWA first in the toolbar');
    pause(0.05);
    return;
end
set(handles.colorbarvariance,'state','off');
set(handles.colorbarstochastic,'state','off');
cla(handles.axes2,'reset');
set(handles.axes2,'visible','off');
cla(handles.axes3,'reset');
set(handles.axes3,'visible','off');
cla(handles.axes4,'reset');
set(handles.axes4,'visible','off');
normalized=get(handles.Normalizedcheckbox,'Value');
if ~normalized
      set(handles.Wait,'String','Please check Normalized mode first');
      set(handles.uipanellevel,'visible','on');
      pause(.5);
      return;
end
t1=handles.t1;
spec=handles.spec;
for i=1:length(spec(:,1))
    for j=1:length(spec(1,:))
        if spec(i,j)>0.999
            spec(i,j)=0.999;
        end
        if spec(i,j)<0.001
            spec(i,j)=0.001;
        end
    end
end
specdB=10*log10(spec./(1-spec)); % LSWA spectrogram in dB;
Omega=handles.Omega;
axes(handles.axes2);
LSWS3D3=handles.LSWS3D3;
if LSWS3D3==1
    mesh(t1,Omega,specdB);
    LSWS3D3=2;
else
    surf(t1,Omega,specdB);
    LSWS3D3=1;
end

shading flat; 
axis tight;
colormap jet;
ylabel('Cyclic frequency','fontsize',10); 
xlabel('Time','fontsize',10); 
align([handles.axes1 handles.axes2],'VerticalAlignment','none', 'HorizontalAlignment','center');
set(handles.Wait,'String','LSWS (power spectral density)');
pause(0.05);
handles.colorbar_VarordB=2;
handles.color_bar=1;
handles.colorbar_check=1;
handles.cbar=0;
handles.cbar2=0;
handles.specdB=specdB;
handles.LSWS3D3=LSWS3D3;
guidata(hObject,handles);

% --------------------------------------------------------------------
function PercentV_Callback(hObject, eventdata, handles)
% hObject    handle to PercentV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
if ~isfield(handles,'spec')
    cla(handles.axes2,'reset');
    set(handles.axes2,'visible','off');
    set(handles.Wait,'String','Please push LSWA first in the toolbar');
    pause(.5);
    return;
end
set(handles.colorbarvariance,'state','off');
set(handles.colorbarstochastic,'state','off');
cla(handles.axes2,'reset');
set(handles.axes2,'visible','off');
cla(handles.axes3,'reset');
set(handles.axes3,'visible','off');
cla(handles.axes4,'reset');
set(handles.axes4,'visible','off');
t1=handles.t1;
spec=handles.spec;
Omega=handles.Omega;
axes(handles.axes2);
LSWS3D1=handles.LSWS3D1;
normalized=get(handles.Normalizedcheckbox,'Value');
if normalized
    spectrum=spec*100;
    set(handles.Wait,'String','Least-squares wavelet spectrogram (percentage variance)');
    pause(0.05);
else
    spectrum=spec;
    set(handles.Wait,'String','Least-squares wavelet spectrogram (variance: signals)');
    pause(0.05);
end
if LSWS3D1==1
    mesh(t1,Omega,spectrum);
    LSWS3D1=2;
else
    surf(t1,Omega,spectrum);
    LSWS3D1=1;
end
shading flat;
axis tight;
colormap jet;
ylabel('Cyclic frequency','fontsize',10);
xlabel('Time','fontsize',10); 
align([handles.axes1 handles.axes2],'VerticalAlignment','none', 'HorizontalAlignment','center');
handles.colorbar_VarordB=1;
handles.colorbar_check=1;
handles.color_bar=1;
handles.cbar=0;
handles.LSWS3D1=LSWS3D1;
guidata(hObject,handles);



function X_panelaxes_Callback(hObject, eventdata, handles)
% hObject    handle to X_panelaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of X_panelaxes as text
%        str2double(get(hObject,'String')) returns contents of X_panelaxes as a double


% --- Executes during object creation, after setting all properties.
function X_panelaxes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to X_panelaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function texttimeunit_Callback(hObject, eventdata, handles)
% hObject    handle to texttimeunit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of texttimeunit as text
%        str2double(get(hObject,'String')) returns contents of texttimeunit as a double
TL=get(hObject,'String'); % gets the time label entered by the user
handles = guidata(hObject); 
if strcmp(TL,'s') || strcmp(TL,'S') || strcmp(TL,'Sec') || strcmp(TL,'sec') ...
        || strcmp(TL,'Secs') || strcmp(TL,'secs') || strcmp(TL,'Seconds') ...
        || strcmp(TL,'seconds') || strcmp(TL,'Second') || strcmp(TL,'second')
    axes(handles.axes1);
    xlabel(strcat('Time (s)'),'fontsize',10);
    axes(handles.axes2);
    ylabel('Frequency (Hz)','fontsize',10); 
    xlabel(strcat('Time (s)'),'fontsize',10);
    axes(handles.axes5);
    xlabel('Frequency (Hz)','fontsize',10);
elseif strcmp(TL,'min') || strcmp(TL,'h') || strcmp(TL,'d') || strcmp(TL,'y')
    axes(handles.axes1);
    xlabel(strcat('Time (', TL, ')'),'fontsize',10);
    axes(handles.axes2);
    ylabel(strcat('Cyclic frequency',' (c/', TL, ')'),'fontsize',10); 
    xlabel(strcat('Time (', TL, ')'),'fontsize',10);
    axes(handles.axes5);
    xlabel(strcat('Cyclic frequency',' (c/', TL, ')'),'fontsize',10);
elseif strcmp(TL,'mm') || strcmp(TL,'cm') || strcmp(TL,'in') || strcmp(TL,'ft')...
        || strcmp(TL,'yd') || strcmp(TL,'m') || strcmp(TL,'km') || strcmp(TL,'mi')
    axes(handles.axes1);
    xlabel(strcat('Distance (', TL, ')'),'fontsize',10);
    axes(handles.axes2);
    ylabel(strcat('Wavenumber',' (c/', TL, ')'),'fontsize',10); 
    xlabel(strcat('Distance (', TL, ')'),'fontsize',10);
    axes(handles.axes5);
    xlabel(strcat('Wavenumber',' (c/', TL, ')'),'fontsize',10);
else
    axes(handles.axes1);
    xlabel(strcat('(', TL, ')'),'fontsize',10);
    axes(handles.axes2);
    ylabel(strcat(' (c/', TL, ')'),'fontsize',10); 
    xlabel(strcat('(', TL, ')'),'fontsize',10);
    axes(handles.axes5);
    xlabel(strcat(' (c/', TL, ')'),'fontsize',10);
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function texttimeunit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to texttimeunit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Xpanelforce_Callback(hObject, eventdata, handles)
% hObject    handle to Xpanelforce (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xpanelforce as text
%        str2double(get(hObject,'String')) returns contents of Xpanelforce as a double
set(handles.uipanelforce,'visible','off');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Xpanelforce_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xpanelforce (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function X_paneltimeseries_Callback(hObject, eventdata, handles)
% hObject    handle to X_paneltimeseries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of X_paneltimeseries as text
%        str2double(get(hObject,'String')) returns contents of X_paneltimeseries as a double
set(handles.uipaneltimeseries,'visible','off');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function X_paneltimeseries_CreateFcn(hObject, eventdata, handles)
% hObject    handle to X_paneltimeseries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Xpanelwindow_Callback(hObject, eventdata, handles)
% hObject    handle to Xpanelwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xpanelwindow as text
%        str2double(get(hObject,'String')) returns contents of Xpanelwindow as a double
set(handles.uipanelwindow,'visible','off');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Xpanelwindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xpanelwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Xpanellevel_Callback(hObject, eventdata, handles)
% hObject    handle to Xpanellevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xpanellevel as text
%        str2double(get(hObject,'String')) returns contents of Xpanellevel as a double
set(handles.uipanellevel,'visible','off');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Xpanellevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xpanellevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Xpanelfrequencies_Callback(hObject, eventdata, handles)
% hObject    handle to Xpanelfrequencies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xpanelfrequencies as text
%        str2double(get(hObject,'String')) returns contents of Xpanelfrequencies as a double;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Xpanelfrequencies_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xpanelfrequencies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Xpaneltimeunit.
function Xpaneltimeunit_Callback(hObject, eventdata, handles)
% hObject    handle to Xpaneltimeunit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipanelaxes,'visible','off');
guidata(hObject, handles);

% --- Executes on button press in Xpanelforce.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to Xpanelforce (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Xpaneltimeseries.
function Xpaneltimeseries_Callback(hObject, eventdata, handles)
% hObject    handle to Xpaneltimeseries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipaneltimeseries,'visible','off');
guidata(hObject, handles);


% --------------------------------------------------------------------
function Tseriesproperties_Callback(hObject, eventdata, handles)
% hObject    handle to Tseriesproperties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Gridonoff_Callback(hObject, eventdata, handles)
% hObject    handle to Gridonoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Grid_ON_Callback(hObject, eventdata, handles)
% hObject    handle to Grid_ON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); 
axes(handles.axes1);
% grid on;
grid minor;
axes(handles.axes5);
% grid on;
grid minor;
guidata(hObject, handles);

% --------------------------------------------------------------------
function Grid_off_Callback(hObject, eventdata, handles)
% hObject    handle to Grid_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); 
axes(handles.axes1);
grid off;
axes(handles.axes5);
grid off;
guidata(hObject, handles);


% --------------------------------------------------------------------
function Time_Ser_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to Time_Ser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); 
if ~isfield(handles,'FirstTS')
    set(handles.Wait,'String','Please open a file containing time/data series');
    pause(.5);
    return;
end
cla(handles.axes1,'reset');
axes(handles.axes1);
T_S=handles.FirstTS;
set(handles.Wait,'String','Time/data series');
if contains(T_S,'.xlsx') || contains(T_S,'.xls')
    Data=xlsread(T_S);
else
    Data=load(T_S);
end
plot(Data(:,1)-Data(1,1),Data(:,2),'b');
set(gca,'ycolor','b')
axis tight;
hold off;
ylabel('Values','fontsize',10);
xlabel('Time','fontsize',10);
align([handles.axes1 handles.axes2],'VerticalAlignment','none', 'HorizontalAlignment','center');
guidata(hObject, handles);


function Wait_Callback(hObject, eventdata, handles)
% hObject    handle to Wait (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Wait as text
%        str2double(get(hObject,'String')) returns contents of Wait as a double


% --- Executes during object creation, after setting all properties.
function Wait_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Wait (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uipanelaxes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanelaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function colorbarvariance_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to colorbarvariance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
if ~isfield(handles,'colorbar_check')
    set(handles.colorbarvariance,'state','off');
    set(handles.axes2,'visible','off');
    set(handles.Wait,'String','Please push LSWA or LSCWA first in the toolbar');
    pause(.5);
    return;
end
colorbar_check=handles.colorbar_check;
normalized=get(handles.Normalizedcheckbox,'Value');
if colorbar_check==1 
    cbar=handles.cbar;
    if (cbar==1)
         cla(handles.axes3,'reset');
         set(handles.axes3,'visible','off');
         handles.cbar=0;
         guidata(hObject, handles);
         return;
    end
    colorbar_VarordB=handles.colorbar_VarordB;
    if colorbar_VarordB==1
       spec=handles.spec;
       axes(handles.axes3);
       cb=colorbar(handles.axes3);
       if normalized
           caxis([min(min(spec))*100 max(max(spec))*100]);
       else
           caxis([min(min(spec)) max(max(spec))]);
       end
       colormap(handles.axes3,jet);
       if normalized
           ylabel(cb, 'Percentage variance (%)','fontsize',10);
       else
           ylabel(cb, 'Variance (signals)','fontsize',10);
       end
    end
    % Colorbar for LSWS in PSD(dB)
    if colorbar_VarordB==2
       specdB=handles.specdB;
       axes(handles.axes3);
       cb=colorbar(handles.axes3);
       caxis([min(min(specdB)) max(max(specdB))]);
       colormap(handles.axes3,jet);
       ylabel(cb, 'Power spectral density (dB)','fontsize',10);
    end
end
%--------------------------------------------
% Colorbar for LSCWS
normalized=get(handles.Normalizedcheckbox,'Value');
if colorbar_check==2
    cbar=handles.cbar;
    if (cbar==1)
         cla(handles.axes3,'reset');
         set(handles.axes3,'visible','off');
         handles.cbar=0;
         guidata(hObject, handles);
         return;
    end
    xspec=handles.xspec;
    axes(handles.axes3);
    cb=colorbar(handles.axes3);
    if normalized
        caxis([min(min(xspec))*100 max(max(xspec))*100]);
    else
        caxis([min(min(xspec)) max(max(xspec))]);
    end
    colormap(handles.axes3,jet);
    if normalized
        ylabel(cb, 'Percentage variance (%)','fontsize',10);
    else
        ylabel(cb, 'Variance (signals)','fontsize',10); 
    end
  %  ylabel(cb, 'Ln (S1) +Ln (S2)');
end
%--------------------------------------------
% Colorbar for LSCWS
normalized=get(handles.Normalizedcheckbox,'Value');
if colorbar_check==3
    cbar=handles.cbar;
    if (cbar==1)
         cla(handles.axes3,'reset');
         set(handles.axes3,'visible','off');
         handles.cbar=0;
         guidata(hObject, handles);
         return;
    end
    xspec=handles.xspec;
    axes(handles.axes3);
    cb=colorbar(handles.axes3);
    mm=min(min(xspec));
    MM=max(max(xspec));     
    if normalized
        mm=mm*100;  MM=MM*100;
        ylabel(cb, 'Percentage variance (%)','fontsize',10);
    else
        ylabel(cb, 'Variance (signals)','fontsize',10);
    end
    colormap(handles.axes3,jet);
    if mm==MM
        caxis([mm-0.001 MM+0.001]);
    else
        caxis([mm MM]);
    end
  %  ylabel(cb, 'Ln (S1) +Ln (S2)');
end
handles.cbar=1;
guidata(hObject, handles);


% --------------------------------------------------------------------
function LSWA_Stochastic_Callback(hObject, eventdata, handles)
% hObject    handle to LSWA_Stochastic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
if ~isfield(handles,'spec')
    cla(handles.axes2,'reset');
    set(handles.axes2,'visible','off');
    set(handles.Wait,'String','Please push LSWA first in the toolbar');
    pause(0.05);
    return;
end
set(handles.colorbarvariance,'state','off');
set(handles.colorbarstochastic,'state','off');
cla(handles.axes2,'reset');
set(handles.axes2,'visible','off');
cla(handles.axes3,'reset');
set(handles.axes3,'visible','off');
cla(handles.axes4,'reset');
set(handles.axes4,'visible','off');

normalized=get(handles.Normalizedcheckbox,'Value');
spec=handles.spec;
stochastic=handles.stochastic;
t1=handles.t1;
Omega=handles.Omega;
LSWS3D2=handles.LSWS3D2;
axes(handles.axes2); 
if normalized
    spec=spec*100;
end
if LSWS3D2==1
    surf(t1,Omega,spec);
    colormap jet;
    if normalized
        hold on;
        freezeColors;
        mesh(t1,Omega,stochastic*100);
    else
        set(handles.Wait,'String','The level is defined for normalized spectrogram in here!');
        set(handles.uipanellevel,'visible','on');
    end 
    LSWS3D2=2;
elseif LSWS3D2==2
    mesh(t1,Omega,spec);
    colormap jet;
    if normalized
        hold on;
        freezeColors;
        mesh(t1,Omega,stochastic*100);
    else
        set(handles.Wait,'String','The level is defined for normalized spectrogram in here!');
        set(handles.uipanellevel,'visible','on');
    end
    LSWS3D2=3;
elseif LSWS3D2==3
    surf(t1,Omega,spec);
    colormap jet;
    if normalized
        hold on;
        freezeColors;
        surf(t1,Omega,stochastic*100);
    else
        set(handles.Wait,'String','The level is defined for normalized spectrogram in here!');
        set(handles.uipanellevel,'visible','on');
    end
    LSWS3D2=1;
end
shading flat; 
axis tight;
ylabel('Cyclic frequency','fontsize',10); 
xlabel('Time','fontsize',10);
if normalized
    colormap gray;
    [LOm,LT]=size(stochastic);
    LB=min(min(stochastic))*100;
    UB=stochastic(LOm, floor(LT/2))*100;
    if UB>LB
        caxis([LB UB+0.1*(UB-LB)]);
    end
    handles.stochSave=stochastic;
    set(handles.Wait,'String','The LSWS (PV) with its stochastic surface');
    pause(0.05);
    hold off
else
    set(handles.Wait,'String','Least-squares wavelet spectrogram (variance: signals)');
    pause(0.05);
end
align([handles.axes1 handles.axes2],'VerticalAlignment','none', 'HorizontalAlignment','center');
handles.cbar2=0;
handles.cbar=0;
handles.colorbar_VarordB=1;
handles.LSWS3D2=LSWS3D2;
handles.LSWSorLSCWSstoch=1;
handles.colorbar_check=1;
handles.color_bar=2;

guidata(hObject,handles);

% --------------------------------------------------------------------
function colorbarstochastic_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to colorbarstochastic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
if ~isfield(handles,'colorbar_check')
    set(handles.colorbarstochastic,'state','off');
    set(handles.axes2,'visible','off');
    set(handles.Wait,'String','Please push LSWA or LSCWA first in the toolbar');
      pause(.5);
      return;
end
color_bar=handles.color_bar;
if color_bar==1
   cla(handles.axes4,'reset');
   set(handles.axes4,'visible','off');
   set(handles.colorbarstochastic,'state','off');
   handles.cbar2=0;
   guidata(hObject, handles);
   return;
end
cbar2=handles.cbar2;
if (cbar2==1)
     cla(handles.axes4,'reset');
     set(handles.axes4,'visible','off');
     handles.cbar2=0;
     guidata(hObject, handles);
     return;
end

if ~isfield(handles,'LSWSorLSCWSstoch')
    set(handles.Wait,'String','Please go to Display in menu first');
    pause(.5);
    return;
end
LSWSorLSCWSstoch=handles.LSWSorLSCWSstoch;
if LSWSorLSCWSstoch==1
    stochastic=handles.stochastic;
    axes(handles.axes4);
    cb=colorbar(handles.axes4);
    [LOm,LT]=size(stochastic);
    LB=min(min(stochastic))*100;
    UB=stochastic(LOm, floor(LT/2))*100;
    if UB>LB
        caxis([LB UB+0.1*(UB-LB)]);
    end
    colormap(handles.axes4,gray);
    ylabel(cb, 'Critical percentage variance (%)','fontsize',10);
    handles.cbar2=1;
elseif LSWSorLSCWSstoch==2
    xstochastic=handles.xstochastic;
    axes(handles.axes4);
    cb=colorbar(handles.axes4);
    [LOm,LT]=size(xstochastic);
    LB=min(min(xstochastic));
    UB=xstochastic(LOm, floor(LT/2));
    if UB>LB
        caxis([LB UB+0.1*(UB-LB)]);
    end
    colormap(handles.axes4,gray);
    ylabel(cb, 'Critical percentage variance (%)','fontsize',10);
    handles.cbar2=1;
end
    
guidata(hObject, handles);


% --------------------------------------------------------------------
function LSSA_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to LSS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipanelforce,'visible','off');
set(handles.uipanellevel,'visible','off');
set(handles.uipanelwindow,'visible','off');
set(handles.uipanelarrow,'visible','off');
set(handles.uipanelfreq,'visible','off');
cla(handles.axes5,'reset');
set(handles.axes5,'visible','off');
if ~isfield(handles,'FirstTS')
    set(handles.Wait,'String','Please open a file containing time/data series');
      pause(.5);
      return;
end

T_S=handles.FirstTS;
set(handles.Wait,'String','Please wait (it may take a while)...');
pause(0.05);
guidata(hObject, handles);
Omega=handles.Omega;
if contains(T_S,'.xlsx') || contains(T_S,'.xls')
    Data=xlsread(T_S);
else
    Data=load(T_S);
end
tt = Data(:,1);
t1=tt-Data(1,1);
f = Data(:,2);
C_f=zeros;
[~,n_c]=size(Data);
if  (n_c>3)
     for i=1:length(t1)
        for j=1:length(t1)
            C_f(i,j)=Data(i,j+2);
        end
     end
     P_f=inv(C_f);
elseif n_c==3
     P_f =1./(Data(:,3).^2);
elseif n_c==2
     P_f =1;
end
%--------------------------------------------
content=get(handles.SignificanceLevel,'Value');
if content==1
    alpha=0.01;
elseif content==2
    alpha=0.05;
else
    alpha=0.1;
end
%--------------------------------------------
QQ(1)=get(handles.datum,'Value');
QQ(2)=get(handles.lineartrend,'Value');
QQ(3)=get(handles.quadratictrend,'Value');
QQ(4)=get(handles.cubictrend,'Value');
normalized=get(handles.Normalizedcheckbox,'Value');
force_sinusoids=get(handles.force,'String');
FS=str2num(force_sinusoids);
if isnan(FS)
    Q=QQ;
else
    Q=[QQ,FS];
end

IndexDatum1=get(handles.DatumShifts1,'String');
IndexDatum1=str2num(IndexDatum1);
if isnan(IndexDatum1)
   IndexDatum1=[];
end
Lt1=length(t1);
if any(floor(IndexDatum1)-IndexDatum1) || any(IndexDatum1>=Lt1) || any(IndexDatum1<1)
    set(handles.Wait,'String','Error! Datum shift indices must be natural numbers less than the size of series');
    set(handles.uipanelforce,'visible','on');
    return;
end
%-----------Run the LSSA---------------------
[spectrum, ~, lsslevel, r, xhat, COVcoeff, ~]=LSSA(t1, f, P_f, Omega, Q, IndexDatum1, alpha, normalized);
%--------------------------------------------
Lid=length(IndexDatum1);
LSSA_TrendsCoeff=zeros;
if Q(1)==1
    LSSA_CosSin=xhat(Lid+nnz(Q(1:4))+1:length(xhat));
    LSSA_TrendsCoeff=xhat(1:Lid+nnz(Q(1:4)));
elseif nnz(Q(2:4))>1
    LSSA_CosSin=xhat(nnz(Q(2:4))+1:length(xhat));
    LSSA_TrendsCoeff=xhat(1:nnz(Q(2:4)));
else
    LSSA_CosSin=xhat;
end

axes(handles.axes5);
if normalized
    spec=spectrum*100;
else
    spec=spectrum;
end
stem(Omega,spec, '.-b', 'MarkerSize',10,'linewidth',0.1);
hold on
plot(Omega,spec);
hold off
axis tight;
if normalized
    ylabel('Percentage variance (%)','fontsize',10);
else
    ylabel('Variance (signals)','fontsize',10);
end
xlabel('Cyclic frequency','fontsize',10); 
set(handles.Wait,'String','Least-squares spectrum');
pause(0.05);
handles=guidata(hObject);
handles.lsslevel=lsslevel;
handles.spectrum=spectrum;
handles.LSSA_Q = Q;
handles.LSSA_IndexDatum = IndexDatum1;
handles.LSSA_TrendsCoeff=LSSA_TrendsCoeff;
handles.LSSA_CosSin=LSSA_CosSin;
handles.t1=t1;
residual(:,1)=tt;
residual(:,2)=r;
if ~isempty(COVcoeff)
    handles.r_LSSA = residual;
    handles.LSSA_COVcoeff = COVcoeff;
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function LSspectrum_Callback(hObject, eventdata, handles)
% hObject    handle to LSspectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function LSwavelet_Callback(hObject, eventdata, handles)
% hObject    handle to LSwavelet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function LSS_PV_Callback(hObject, eventdata, handles)
% hObject    handle to LSS_PV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
if ~isfield(handles,'spectrum')
    cla(handles.axes5,'reset');
    set(handles.axes5,'visible','off');
    set(handles.Wait,'String','Please push LSSA first in the toolbar');
    pause(0.05);
    return;
end
cla(handles.axes5,'reset');
set(handles.axes5,'visible','off');
spectrum=handles.spectrum;
Omega=handles.Omega;
normalized=get(handles.Normalizedcheckbox,'Value');
axes(handles.axes5);
if normalized
    plot(Omega,spectrum*100); 
    ylabel('Percentage variance (%)','fontsize',10);
    set(handles.Wait,'String','Least-squares spectrum (percentage variance)');
    pause(0.05);
else
    plot(Omega,spectrum); 
    ylabel('Variance (signals)','fontsize',10);
    set(handles.Wait,'String','Least-squares spectrum (variance: signals)');
    pause(0.05);
end
axis tight;
xlabel('Cyclic frequency','fontsize',10); 
guidata(hObject,handles);

% --------------------------------------------------------------------
function LSS_SL_Callback(hObject, eventdata, handles)
% hObject    handle to LSS_SL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
if ~isfield(handles,'spectrum')
    cla(handles.axes5,'reset');
    set(handles.axes5,'visible','off');
    set(handles.Wait,'String','Please push LSSA first in the toolbar');
    pause(0.05);
    return;
end
cla(handles.axes5,'reset');
set(handles.axes5,'visible','off');
spectrum=handles.spectrum;
lsslevel=handles.lsslevel;
Omega=handles.Omega;
normalized=get(handles.Normalizedcheckbox,'Value');
axes(handles.axes5);
if normalized
    spec=spectrum*100;
else
    spec=spectrum;
end
stem(Omega,spec, '.-b', 'MarkerSize',10,'linewidth',0.1); 
hold on;
plot(Omega,spec, 'b');
if normalized
    ylabel('Percentage variance (%)','fontsize',10); 
    plot(Omega,lsslevel*100,'r');
    handles.lssCV=lsslevel(1);
    set(handles.Wait,'String','LSS (percentage variance) and its critical threshold');
    pause(0.05);
else
    ylabel('Variance (signals)','fontsize',10); 
    set(handles.Wait,'String','The level is defined for normalized spectrum in here!');
    set(handles.uipanellevel,'visible','on');
    pause(0.05);
end
hold off;
axis tight;
xlabel('Cyclic frequency','fontsize',10); 
guidata(hObject,handles);

% --------------------------------------------------------------------
function LSS_PSD_Callback(hObject, eventdata, handles)
% hObject    handle to LSS_PSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
if ~isfield(handles,'spectrum')
    cla(handles.axes5,'reset');
    set(handles.axes5,'visible','off');
    set(handles.Wait,'String','Please push LSSA first in the toolbar');
    pause(0.05);
    return;
end
normalized=get(handles.Normalizedcheckbox,'Value');
if ~normalized
      set(handles.Wait,'String','Please check Normalized mode first');
      set(handles.uipanellevel,'visible','on');
      pause(.5);
      return;
end
cla(handles.axes5,'reset');
set(handles.axes5,'visible','off');
spectrum=handles.spectrum;
for i=1:length(spectrum(:,1))
    for j=1:length(spectrum(1,:))
        if spectrum(i,j)>0.999
            spectrum(i,j)=0.999;
        end
        if spectrum(i,j)<0.001
            spectrum(i,j)=0.001;
        end
    end
end
spectrumdB=10*log10(spectrum./(1-spectrum)); % LSWA spectrogram in dB;
Omega=handles.Omega;
axes(handles.axes5);
plot(Omega,spectrumdB, '.-b', 'MarkerSize',10,'linewidth',0.1); 
axis tight;
ylabel('PSD (dB)','fontsize',10); 
xlabel('Cyclic frequency','fontsize',10);
set(handles.Wait,'String','LSS (power spectral density)');
pause(0.05);
%set(gca, 'FontSize' ,10) ;  %Set font size of x and y axis values
guidata(hObject,handles);


% --------------------------------------------------------------------
function TSerieswithouterrorbars_Callback(hObject, eventdata, handles)
% hObject    handle to TSerieswithouterrorbars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); 
if ~isfield(handles,'FirstTS')
    set(handles.Wait,'String','Please open a file containing 1st time/data series');
      pause(.5);
      return;
end
T_S=handles.FirstTS;

cla(handles.axes1,'reset');
axes(handles.axes1);
set(handles.Wait,'String','Time/data series');
if contains(T_S,'.xlsx') || contains(T_S,'.xls')
    Data=xlsread(T_S);
else
    Data=load(T_S);
end
Data0 = Data(:,1)-Data(1,1);
IndexDatum1=get(handles.DatumShifts1,'String');
IndexDatum1=str2num(IndexDatum1);
hold on
if ~isnan(IndexDatum1)
    IndexDatum1 = intersect(1:length(Data0),IndexDatum1);
    min_data=min(Data(:,2));
    max_data=max(Data(:,2));
    for k = 1:length(IndexDatum1)
        plot([Data0(IndexDatum1(k)) Data0(IndexDatum1(k))], [min_data max_data], '-b', 'LineWidth',1);
    end
end
plot(Data0, Data(:,2),'-db','MarkerFaceColor','b','LineWidth', 0.1,'MarkerSize',2);
set(gca,'ycolor','b')
axis tight;
ylabel('Values','fontsize',10);
xlabel('Time','fontsize',10);
hold off
align([handles.axes1 handles.axes2],'VerticalAlignment','none', 'HorizontalAlignment','center');
guidata(hObject, handles);

% --------------------------------------------------------------------
function TSerieswitherroebars_Callback(hObject, eventdata, handles)
% hObject    handle to TSerieswitherroebars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); 
if ~isfield(handles,'FirstTS')
    set(handles.Wait,'String','Please open a file containing 1st time/data series');
    pause(.5);
    return;
end
T_S=handles.FirstTS;
cla(handles.axes1,'reset');
axes(handles.axes1);
set(handles.Wait,'String','Please wait ...');
pause(0.05);
guidata(hObject, handles);
%--------------------------------------------
if contains(T_S,'.xlsx') || contains(T_S,'.xls')
    Data=xlsread(T_S);
else
    Data=load(T_S);
end
[~,n_c]=size(Data);
if  (n_c>3)
     std=zeros;
     for i=1:length(Data(:,1))
            std(i)=sqrt(Data(i,i+2));
     end
elseif n_c==3
     std=Data(:,3);
end
%--------------------------------------------
Data0 = Data(:,1)-Data(1,1);
IndexDatum1=get(handles.DatumShifts1,'String');
IndexDatum1=str2num(IndexDatum1);
hold on
if ~isnan(IndexDatum1)
    IndexDatum1 = intersect(1:length(Data0),IndexDatum1);
    if n_c>2
        min_data=min(Data(:,2)-std);
        max_data=max(Data(:,2)+std);
    else
        min_data=min(Data(:,2));
        max_data=max(Data(:,2));
    end
    for k = 1:length(IndexDatum1)
        plot([Data0(IndexDatum1(k)) Data0(IndexDatum1(k))], [min_data max_data], '-b','LineWidth',1);
    end
end
axis tight;
if n_c>2
    for i=1:length(Data(:,1))
       plot([Data0(i) Data0(i)],[Data(i,2)-std(i) Data(i,2)+std(i)],'-r',...
           'LineWidth',0.5);
       plot(Data0(i),Data(i,2),'db','MarkerFaceColor','b','MarkerSize',2);%,'LineWidth',0.01);
    end
else
    set(handles.Wait,'String','The standard deviations are not given!');
    plot(Data0,Data(:,2),'db','MarkerFaceColor','b','MarkerSize',2);
end
set(gca,'ycolor','b')
axis tight;
ylabel('Values','fontsize',10); 
xlabel('Time','fontsize',10);
align([handles.axes1 handles.axes2],'VerticalAlignment','none', 'HorizontalAlignment','center');
hold off
if n_c>2
    set(handles.Wait,'String','Time/data series');
end
guidata(hObject, handles);





% --------------------------------------------------------------------
function Help_1_Callback(hObject, eventdata, handles)
% hObject    handle to Help_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function Help_pdf_Callback(hObject, eventdata, handles)
% hObject    handle to Help_pdf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
winopen('User Manual.pdf')


% --- Executes during object creation, after setting all properties.
function TS2filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TS2filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function LSCWA_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to C_LSWA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set visibility of the panel "off"
handles = guidata(hObject); 
if ~isfield(handles,'FirstTS')
    set(handles.Wait,'String','Please open a file containing 1st time/data series');
      pause(.5);
      return;
end
if ~isfield(handles,'SecondTS')
    set(handles.Wait,'String','Please open a file containing 2nd time/data series');
      pause(.5);
      return;
end

set(handles.uipanelforce,'visible','off');
set(handles.uipanellevel,'visible','off');
set(handles.uipanelwindow,'visible','off');
set(handles.uipanelaxes,'visible','off');
set(handles.uipanelarrow,'visible','off');
set(handles.uipanelfreq,'visible','off');
T_S=handles.FirstTS;
T_S2=handles.SecondTS;
set(handles.Wait,'String','Please wait (it may take a while)...');
pause(0.05);
guidata(hObject, handles);
Omega=handles.Omega;
if contains(T_S,'.xlsx') || contains(T_S,'.xls')
    Data=xlsread(T_S);
else
    Data=load(T_S);
end
if contains(T_S2,'.xlsx') || contains(T_S2,'.xls')
    Data2=xlsread(T_S2);
else
    Data2=load(T_S2);
end
MinTim=min(Data(1,1),Data2(1,1));
tt = Data(:,1);
t1=tt-MinTim;
f1 = Data(:,2);
C_f1=zeros;
[~,n_c]=size(Data);
if  (n_c>3) 
    % formats the weight matrix (the inverse of the covariance matrix 
    % associated with the series)
     for i=1:length(t1)
        for j=1:length(t1)
            C_f1(i,j)=Data(i,j+2);
        end
     end
     P_f1=inv(C_f1);
elseif n_c==3
     P_f1 =1./(Data(:,3).^2);
elseif n_c==2
     P_f1 =ones(1,length(t1))';
end
%--------------------------------------------
tt2 = Data2(:,1);
t2=tt2-MinTim;
f2 = Data2(:,2);
C_f2=zeros;
[~,n_c2]=size(Data2);
if  (n_c2>3) 
    % formats the weight matrix (the inverse of the covariance matrix 
    % associated with the series)
     for i=1:length(t2)
        for j=1:length(t2)
            C_f2(i,j)=Data2(i,j+3);
        end
     end
     P_f2=inv(C_f2);
elseif n_c2==3
     P_f2 =1./(Data2(:,3).^2);
elseif n_c2==2
     P_f2 =1;
end
%--------------------------------------------
content=get(handles.SignificanceLevel,'Value');
if content==1
    alpha=0.01;
elseif content==2
    alpha=0.05;
else
    alpha=0.1;
end
%--------------------------------------------
M1=str2double(get(handles.SR,'String'));
if isnan(M1) || M1<1 % Checks sampling rate M for the 1st time series
    set(handles.SR,'String','');
    set(handles.Wait,'String','Error! Sampling rate must be a positive number');
    set(handles.uipanelwindow,'visible','on');
    return;
end
%--------------------------------------------
L11=str2double(get(handles.L1,'String'));
if isnan(L11) || L11<0
    set(handles.L1,'String','');
    set(handles.Wait,'String',...
        'Error! L1 must be a small nonnegative integer');
    set(handles.uipanelwindow,'visible','on');
    return;
end
%--------------------------------------------
L01=str2double(get(handles.L0,'String'));
if isnan(L01) || L01<0 || rem(L01,1)~=0
    set(handles.L0,'String','');
    set(handles.Wait,'String',...
        'Error! L0 must be a nonnegative integer');
    set(handles.uipanelwindow,'visible','on');
    return;
end
%--------------------------------------------
M2=str2double(get(handles.SR2,'String'));
if isnan(M2) || M2<1  % Checks sampling rate M for the 2nd time series
    set(handles.SR2,'String','');
    set(handles.Wait,'String','Error! Sampling rate must be a positive number');
    set(handles.uipanelwindow,'visible','on');
    return;
end
%--------------------------------------------
L12=str2double(get(handles.L_1,'String')); % L_1 for the second time series
if isnan(L12) || L12<0
    set(handles.L_1,'String','');
    set(handles.Wait,'String',...
        'Error! L1 must be a small nonnegative integer');
    set(handles.uipanelwindow,'visible','on');
    return;
end
%--------------------------------------------
L02=str2double(get(handles.L_0,'String')); % L_0 for the secind time series
if isnan(L02) || L02<0 || rem(L02,1)~=0
    set(handles.L_0,'String','');
    set(handles.Wait,'String',...
        'Error! L0 must be a nonnegative integer');
    set(handles.uipanelwindow,'visible','on');
    return;
end
%--------------------------------------------
% Adjust QQ1, the vector of constituents of known forms
QQ1(1)=get(handles.datum,'Value');
QQ1(2)=get(handles.lineartrend,'Value');
QQ1(3)=get(handles.quadratictrend,'Value');
QQ1(4)=get(handles.cubictrend,'Value');
force_sinusoids1=get(handles.force,'String');
FS1=str2num(force_sinusoids1);
if isnan(FS1)
    Q1=QQ1;
else
    Q1=[QQ1,FS1];
end
% Adjust QQ2, the vector of constituents of known forms
QQ2(1)=get(handles.datum2,'Value');
QQ2(2)=get(handles.lineartrend2,'Value');
QQ2(3)=get(handles.quadratictrend2,'Value');
QQ2(4)=get(handles.cubictrend2,'Value');
normalized=get(handles.Normalizedcheckbox,'Value');
force_sinusoids2=get(handles.force2,'String');
FS2=str2num(force_sinusoids2);
if isnan(FS2)
    Q2=QQ2;
else
    Q2=[QQ2,FS2];
end
%--------------------------------------------
IndexDatum1=get(handles.DatumShifts1,'String');
IndexDatum1=str2num(IndexDatum1);
if isnan(IndexDatum1)
   IndexDatum1=[];
end
Lt1=length(t1);
if any(floor(IndexDatum1)-IndexDatum1) || any(IndexDatum1>=Lt1) || any(IndexDatum1<1)
    set(handles.Wait,'String','Error! Datum shift indices must be natural numbers less than the size of series');
    set(handles.uipanelforce,'visible','on');
    return;
end
IndexDatum2=get(handles.DatumShifts2,'String');
IndexDatum2=str2num(IndexDatum2);
if isnan(IndexDatum2)
   IndexDatum2=[];
end
Lt2=length(t2);
if any(floor(IndexDatum2)-IndexDatum2) || any(IndexDatum2>=Lt2) || any(IndexDatum2<1)
    set(handles.Wait,'String','Error! Datum shift indices must be natural numbers less than the size of series');
    set(handles.uipanelforce,'visible','on');
    return;
end
%--------------------------------------------
MorlCoeff=get(handles.MorletCheck,'Value');
MorlCoeff=MorlCoeff*0.0125;
%-----------Run the LSCWA---------------------
[t_union, cross_spectrogram, cross_stoch_surf, phase_diff]=...
    LSCWA(t1, t2, f1, f2, P_f1, P_f2, Omega, M1, L11, L01, Q1, IndexDatum1, ...
    M2, L12, L02, Q2, IndexDatum2, alpha, MorlCoeff, normalized);
%---------------------------------------------
set(handles.colorbarvariance,'state','off');
set(handles.colorbarstochastic,'state','off');
cla(handles.axes2,'reset');
set(handles.axes2,'visible','off');
cla(handles.axes3,'reset');
set(handles.axes3,'visible','off');
cla(handles.axes4,'reset');
set(handles.axes4,'visible','off');
% set(gca,'xticklabel','');
% xlabel('');
axes(handles.axes2);
if normalized
    surf(t_union,Omega,cross_spectrogram*100); view(2);
else
    surf(t_union,Omega,cross_spectrogram); view(2);
end
shading flat; 
axis tight;
colormap jet;
ylabel('Cyclic frequency','fontsize',10);
xlabel('Time','fontsize',10); 
%set(gca, 'FontSize' ,12) ;  %Set font size of x and y axis values
align([handles.axes1 handles.axes2],'VerticalAlignment','none', 'HorizontalAlignment','center');
set(handles.Wait,'String','Least-squares cross-wavelet spectrogram');
pause(0.05);
handles=guidata(hObject);
handles.xspec=cross_spectrogram;
handles.xstochastic=cross_stoch_surf;
handles.t1=t1;
handles.t2=t2;
handles.t_union=t_union;
handles.cbar=0;
handles.cbar2=0;
handles.phase_diff=phase_diff;
handles.colorbar_VarordB=1;
handles.color_bar=1;
handles.colorbar_check=2;
handles.LSCWS3D1=1; handles.LSCWS3D2=1; handles.LSCWS3D3=1;
guidata(hObject, handles);


% --------------------------------------------------------------------
function TwoTime_Ser_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to TwoTime_Ser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); 
if ~isfield(handles,'FirstTS')
    set(handles.Wait,'String','Please open a file containing 1st time/data series');
      pause(.5);
      return;
end
if ~isfield(handles,'SecondTS')
    set(handles.Wait,'String','Please open a file containing 2nd time/data series');
      pause(.5);
      return;
end
T_S=handles.FirstTS;
T_S2=handles.SecondTS;
cla(handles.axes1,'reset');
axes(handles.axes1);
set(handles.Wait,'String','Both time/data series');
if contains(T_S,'.xlsx') || contains(T_S,'.xls')
    Data=xlsread(T_S);
else
    Data=load(T_S);
end
if contains(T_S2,'.xlsx') || contains(T_S2,'.xls')
    Data2=xlsread(T_S2);
else
    Data2=load(T_S2);
end
yyaxis left
minval=min(Data(1,1),Data2(1,1));
plot(Data(:,1)-minval,Data(:,2),'b');
set(gca,'ycolor','b')
hold on
ylabel('Values','fontsize',10);
yyaxis right
plot(Data2(:,1)-minval,Data2(:,2),'g');
set(gca,'ycolor','g')
ylabel('Values','fontsize',10);
xlabel('Time','fontsize',10);
axis tight
hold off
align([handles.axes1 handles.axes2],'VerticalAlignment','none', 'HorizontalAlignment','center');
guidata(hObject, handles);


% --- Executes on button press in datum2.
function datum2_Callback(hObject, eventdata, handles)
% hObject    handle to datum2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of datum2


% --- Executes on button press in lineartrend2.
function lineartrend2_Callback(hObject, eventdata, handles)
% hObject    handle to lineartrend2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lineartrend2


% --- Executes on button press in quadratictrend2.
function quadratictrend2_Callback(hObject, eventdata, handles)
% hObject    handle to quadratictrend2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of quadratictrend2


% --- Executes on button press in cubictrend2.
function cubictrend2_Callback(hObject, eventdata, handles)
% hObject    handle to cubictrend2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cubictrend2



function force2_Callback(hObject, eventdata, handles)
% hObject    handle to force2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of force2 as text
%        str2double(get(hObject,'String')) returns contents of force2 as a double


% --- Executes during object creation, after setting all properties.
function force2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to force2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function LSCSA_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to LSCS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'FirstTS')
    set(handles.Wait,'String','Please open a file containing 1st time/data series');
      pause(.5);
      return;
end
if ~isfield(handles,'SecondTS')
    set(handles.Wait,'String','Please open a file containing 2nd time/data series');
      pause(.5);
      return;
end
set(handles.uipanelforce,'visible','off');
set(handles.uipanellevel,'visible','off');
set(handles.uipanelwindow,'visible','off');
set(handles.uipanelaxes,'visible','off');
set(handles.uipanelarrow,'visible','off');
set(handles.uipanelfreq,'visible','off');
handles = guidata(hObject); 
cla(handles.axes5,'reset');
set(handles.axes5,'visible','off');
T_S=handles.FirstTS;
T_S2=handles.SecondTS;
set(handles.Wait,'String','Please wait (it may take a while)...');
pause(0.05);
guidata(hObject, handles);
Omega=handles.Omega;
if contains(T_S,'.xlsx') || contains(T_S,'.xls')
    Data=xlsread(T_S);
else
    Data=load(T_S);
end
if contains(T_S2,'.xlsx') || contains(T_S2,'.xls')
    Data2=xlsread(T_S2);
else
    Data2=load(T_S2);
end
MinTim=min(Data(1,1),Data2(1,1));
%--------------------------------------------
tt1 = Data(:,1);
t1=tt1-MinTim;
f1 = Data(:,2);
C_f1=zeros;
[~,n_c]=size(Data);
if  (n_c>3) 
    % formats the weight matrix (the inverse of the covariance matrix 
    % associated with the series)
     for i=1:length(t1)
        for j=1:length(t1)
            C_f1(i,j)=Data(i,j+2);
        end
     end
     P_f1=inv(C_f1);
elseif n_c==3
     P_f1 =1./(Data(:,3).^2);
elseif n_c==2
     P_f1 =1;
end
%--------------------------------------------
tt2 = Data2(:,1);
t2=tt2-MinTim;
f2 = Data2(:,2);
C_f2=zeros;
[~,n_c2]=size(Data2);
if  (n_c2>3) 
    % formats the weight matrix (the inverse of the covariance matrix 
    % associated with the second series)
     for i=1:length(t2)
        for j=1:length(t2)
            C_f2(i,j)=Data2(i,j+3);
        end
     end
     P_f2=inv(C_f2);
elseif n_c2==3
     P_f2 =1./(Data2(:,3).^2);
elseif n_c2==2
     P_f2 =1;
end
%--------------------------------------------
content=get(handles.SignificanceLevel,'Value');
if content==1
    alpha=0.01;
elseif content==2
    alpha=0.05;
else
    alpha=0.1;
end
%--------------------------------------------
QQ1(1)=get(handles.datum,'Value');
QQ1(2)=get(handles.lineartrend,'Value');
QQ1(3)=get(handles.quadratictrend,'Value');
QQ1(4)=get(handles.cubictrend,'Value');
force_sinusoids1=get(handles.force,'String');
FS1=str2num(force_sinusoids1);
if isnan(FS1)
    Q1=QQ1;
else
    Q1=[QQ1,FS1];
end
QQ2(1)=get(handles.datum2,'Value');
QQ2(2)=get(handles.lineartrend2,'Value');
QQ2(3)=get(handles.quadratictrend2,'Value');
QQ2(4)=get(handles.cubictrend2,'Value');
normalized=get(handles.Normalizedcheckbox,'Value');
force_sinusoids2=get(handles.force2,'String');
FS2=str2num(force_sinusoids2);
if isnan(FS2)
    Q2=QQ2;
else
    Q2=[QQ2,FS2];
end
IndexDatum1=get(handles.DatumShifts1,'String');
IndexDatum1=str2num(IndexDatum1);
if isnan(IndexDatum1)
   IndexDatum1=[];
end
Lt1=length(t1);
if any(floor(IndexDatum1)-IndexDatum1) || any(IndexDatum1>=Lt1) || any(IndexDatum1<1)
    set(handles.Wait,'String','Error! Datum shift indices must be natural numbers less than the size of series');
    set(handles.uipanelforce,'visible','on');
    return;
end
IndexDatum2=get(handles.DatumShifts2,'String');
IndexDatum2=str2num(IndexDatum2);
if isnan(IndexDatum2)
   IndexDatum2=[];
end
Lt2=length(t2);
if any(floor(IndexDatum2)-IndexDatum2) || any(IndexDatum2>=Lt2) || any(IndexDatum2<2)
    set(handles.Wait,'String','Error! Datum shift indices must be natural numbers less than the size of series');
    set(handles.uipanelforce,'visible','on');
    return;
end
%-----------Run the LSCSA--------------------
[cross_spectrum, CrossCritVal, phase_diff_cross_spectrum]=...
    LSCSA(t1, t2, f1, f2, P_f1, P_f2, Omega, Q1, Q2, IndexDatum1, IndexDatum2, alpha, normalized);
%--------------------------------------------
axes(handles.axes5);
if normalized
    xspec=cross_spectrum*100;
else
    xspec=cross_spectrum;
end
stem(Omega,xspec, '.-b', 'MarkerSize',10,'linewidth',0.1); 
hold on
plot(Omega,xspec); 
hold off
axis tight;
if normalized
    ylabel('Percentage variance (%)','fontsize',10);
else
    ylabel('Variance (signals)','fontsize',10);
end
xlabel('Cyclic frequency','fontsize',10); 
%set(gca, 'FontSize' ,10) ;  %Set font size of x and y axis values to 10
set(handles.Wait,'String','Least-squares cross-spectrum');
pause(0.05);
handles=guidata(hObject);
handles.LSCSlevel=CrossCritVal;
handles.phase_diff_cross_spectrum=phase_diff_cross_spectrum;
handles.xspectrum=cross_spectrum;
handles.t1=t1;
handles.t2=t2;
guidata(hObject, handles);


% --------------------------------------------------------------------
function TwoTimeSeriesWithout_Callback(hObject, eventdata, handles)
% hObject    handle to TwoTimeSeriesWithout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); 
if ~isfield(handles,'FirstTS')
    set(handles.Wait,'String','Please open a file containing 1st time/data series');
      pause(.5);
      return;
end
if ~isfield(handles,'SecondTS')
    set(handles.Wait,'String','Please open a file containing 2nd time/data series');
      pause(.5);
      return;
end
T_S=handles.FirstTS;
T_S2=handles.SecondTS;
cla(handles.axes1,'reset');
axes(handles.axes1);
set(handles.Wait,'String','Both time/data series');
if contains(T_S,'.xlsx') || contains(T_S,'.xls')
    Data=xlsread(T_S);
else
    Data=load(T_S);
end
if contains(T_S2,'.xlsx') || contains(T_S2,'.xls')
    Data2=xlsread(T_S2);
else
    Data2=load(T_S2);
end
yyaxis left;
minval=min(Data(1,1),Data2(1,1));
Data0 = Data(:,1)-minval;
IndexDatum1=get(handles.DatumShifts1,'String');
IndexDatum1=str2num(IndexDatum1);
hold on
if ~isnan(IndexDatum1)
    IndexDatum1 = intersect(1:length(Data0),IndexDatum1);
    min_data=min(Data(:,2));
    max_data=max(Data(:,2));
    for k = 1:length(IndexDatum1)
        plot([Data0(IndexDatum1(k)) Data0(IndexDatum1(k))], [min_data max_data], '-b','LineWidth',1);
    end
end
plot(Data0,Data(:,2),'-db','MarkerFaceColor','b','LineWidth', 0.01,'MarkerSize',1);
ylabel('Values','fontsize',10);
set(gca,'ycolor','b')
hold off
%--------------------------------------------
yyaxis right;
hold on
Data02 = Data2(:,1)-minval;
IndexDatum2=get(handles.DatumShifts2,'String');
IndexDatum2=str2num(IndexDatum2);
if ~isnan(IndexDatum2)
    IndexDatum2 = intersect(1:length(Data02),IndexDatum2);
    min_data2=min(Data2(:,2));
    max_data2=max(Data2(:,2));
    for k = 1:length(IndexDatum2)
        plot([Data02(IndexDatum2(k)) Data02(IndexDatum2(k))], [min_data2 max_data2], '-g','LineWidth',1);
    end
end
plot(Data02,Data2(:,2),'-sg','MarkerFaceColor','g','LineWidth', 0.01,'MarkerSize',2);
ylabel('Values','fontsize',10);
set(gca,'ycolor','g')
axis tight;
xlabel('Time','fontsize',10);
%set(gca,'fontsize',12);
align([handles.axes1 handles.axes2],'VerticalAlignment','none', 'HorizontalAlignment','center');
hold off
guidata(hObject, handles);

% --------------------------------------------------------------------
function TwoTSeriesWith_Callback(hObject, eventdata, handles)
% hObject    handle to TwoTSeriesWith (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); 
if ~isfield(handles,'FirstTS')
    set(handles.Wait,'String','Please open a file containing 1st time/data series');
      pause(.5);
      return;
end
if ~isfield(handles,'SecondTS')
    set(handles.Wait,'String','Please open a file containing 2nd time/data series');
      pause(.5);
      return;
end
T_S=handles.FirstTS;
T_S2=handles.SecondTS;
cla(handles.axes1,'reset');
axes(handles.axes1);

set(handles.Wait,'String','Please wait ...');
pause(0.05);
guidata(hObject, handles);

if contains(T_S,'.xlsx') || contains(T_S,'.xls')
    Data=xlsread(T_S);
else
    Data=load(T_S);
end
if contains(T_S2,'.xlsx') || contains(T_S2,'.xls')
    Data2=xlsread(T_S2);
else
    Data2=load(T_S2);
end
[~,n_c1]=size(Data);
if  (n_c1>3)
     std=zeros;
     for i=1:length(Data(:,1))
            std(i)=sqrt(Data(i,i+2));
     end
elseif n_c1==3
     std=Data(:,3);
end
yyaxis left
minval=min(Data(1,1),Data2(1,1));
Data0 = Data(:,1)-minval;
IndexDatum1=get(handles.DatumShifts1,'String');
IndexDatum1=str2num(IndexDatum1);
hold on
if ~isnan(IndexDatum1)
    IndexDatum1 = intersect(1:length(Data0),IndexDatum1);
    if n_c1>2
        min_data=min(Data(:,2)-std);
        max_data=max(Data(:,2)+std);
    else
        min_data=min(Data(:,2));
        max_data=max(Data(:,2));
    end
    for k = 1:length(IndexDatum1)
        plot([Data0(IndexDatum1(k)) Data0(IndexDatum1(k))], [min_data max_data], '-b','LineWidth',1);
    end
end
axis tight;
if n_c1>2
    for i=1:length(Data(:,1))
       plot([Data0(i) Data0(i)],[Data(i,2)-std(i) Data(i,2)+std(i)],'-r',...
           'LineWidth',0.5);
       plot(Data0(i),Data(i,2),'db','MarkerSize',1);%,'LineWidth',0.01);
    end
else
    set(handles.Wait,'String','The standard deviations are not given!');
    plot(Data0,Data(:,2),'db','MarkerFaceColor','b','MarkerSize',1);
end
set(gca,'ycolor','b')
axis tight;
ylabel('Values','fontsize',10); 
xlabel('Time','fontsize',10);
%--------------------------------------------
hold on
[~,n_c2]=size(Data2);
if  (n_c2>3)
     std2=zeros;
     for i=1:length(Data2(:,1))
            std2(i)=sqrt(Data2(i,i+3));
     end
elseif n_c2==3
     std2=Data2(:,3);
end
yyaxis right
Data02 = Data2(:,1)-minval;
IndexDatum2=get(handles.DatumShifts2,'String');
IndexDatum2=str2num(IndexDatum2);
if ~isnan(IndexDatum2)
    IndexDatum2 = intersect(1:length(Data02),IndexDatum2);
    if n_c2>2
        min_data2=min(Data2(:,2)-std2);
        max_data2=max(Data2(:,2)+std2);
    else
        min_data2=min(Data2(:,2));
        max_data2=max(Data2(:,2));
    end
    for k = 1:length(IndexDatum2)
        plot([Data02(IndexDatum2(k)) Data02(IndexDatum2(k))], [min_data2 max_data2], '-g','LineWidth',1);
    end
end
if n_c2>2
    for i=1:length(Data2(:,1))
       plot([Data02(i) Data02(i)],[Data2(i,2)-std2(i) Data2(i,2)+std2(i)],'-m',...
           'LineWidth',0.5);
       plot(Data02(i),Data2(i,2),'sg','MarkerFaceColor','g','MarkerSize',2);%,'LineWidth',0.01);
    end
else
    set(handles.Wait,'String','The standard deviations are not given!')
    plot(Data02,Data2(:,2),'sg','MarkerFaceColor','g','MarkerSize',2);
end
set(gca,'ycolor','g')
axis tight;
ylabel('Values','fontsize',10);
xlabel('Time','fontsize',10);
hold off
align([handles.axes1 handles.axes2],'VerticalAlignment','none', 'HorizontalAlignment','center');
if n_c1>2 && n_c2>2
    set(handles.Wait,'String','Both time/data series');
end

guidata(hObject, handles);


% --------------------------------------------------------------------
function Shading_Callback(hObject, eventdata, handles)
% hObject    handle to Shading (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function shading_faceted_Callback(hObject, eventdata, handles)
% hObject    handle to shading_faceted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); 
axes(handles.axes2);
shading faceted;
guidata(hObject, handles);

% --------------------------------------------------------------------
function Shading_Flat_Callback(hObject, eventdata, handles)
% hObject    handle to Shading_Flat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); 
axes(handles.axes2);
shading flat;
guidata(hObject, handles);


% --------------------------------------------------------------------
function LSCWS_Callback(hObject, eventdata, handles)
% hObject    handle to LSCWS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function LSCWS_PV_Callback(hObject, eventdata, handles)
% hObject    handle to LSCWS_PV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
if ~isfield(handles,'xspec')
    cla(handles.axes2,'reset');
    set(handles.axes2,'visible','off');
    set(handles.Wait,'String','Please push LSCWA first in the toolbar');
    pause(0.05);
    return;
end
set(handles.colorbarvariance,'state','off');
set(handles.colorbarstochastic,'state','off');
cla(handles.axes2,'reset');
set(handles.axes2,'visible','off');
cla(handles.axes3,'reset');
set(handles.axes3,'visible','off');
cla(handles.axes4,'reset');
set(handles.axes4,'visible','off');
t_union=handles.t_union;
xspec=handles.xspec;
Omega=handles.Omega;
LSCWS3D1=handles.LSCWS3D1;
axes(handles.axes2);
normalized=get(handles.Normalizedcheckbox,'Value');
if normalized
    xspec=xspec*100;
    set(handles.Wait,'String','LSCWS (percentage variance)');
    pause(0.05);
else
    set(handles.Wait,'String','LSCWS (variance: signals)');
    pause(0.05);
end
if LSCWS3D1==1
    mesh(t_union,Omega,xspec);
    LSCWS3D1=2;
elseif LSCWS3D1==2
    surf(t_union,Omega,xspec);
    LSCWS3D1=1;
end
axis tight;
shading flat;
colormap jet;
ylabel('Cyclic frequency','fontsize',10);
xlabel('Time','fontsize',10); 
%set(gca, 'FontSize' ,12) ;  %Set font size of x and y axis values
align([handles.axes1 handles.axes2],'VerticalAlignment','none', 'HorizontalAlignment','center');
handles.colorbar_VarordB=1;
handles.colorbar_check=2;
handles.color_bar=1;
handles.LSCWS3D1=LSCWS3D1;
handles.cbar=0;
guidata(hObject,handles);

% --------------------------------------------------------------------
function LSCWS_Stochastic_Callback(hObject, eventdata, handles)
% hObject    handle to LSCWS_Stochastic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
if ~isfield(handles,'xspec')
    cla(handles.axes2,'reset');
    set(handles.axes2,'visible','off');
    set(handles.Wait,'String','Please push LSCWA first in the toolbar');
    pause(0.05);
    return;
end
normalized=get(handles.Normalizedcheckbox,'Value');
set(handles.colorbarvariance,'state','off');
set(handles.colorbarstochastic,'state','off');
cla(handles.axes3,'reset');
set(handles.axes3,'visible','off');
cla(handles.axes4,'reset');
set(handles.axes4,'visible','off');
xspec=handles.xspec;
xstochastic=handles.xstochastic;
t_union=handles.t_union;
Omega=handles.Omega;
LSCWS3D2=handles.LSCWS3D2;

cla(handles.axes2,'reset');
set(handles.axes2,'visible','off');  
axes(handles.axes2); 
if normalized
    xsp=xspec*100;
else
    xsp=xspec;
end

if LSCWS3D2==1
    surf(t_union,Omega,xsp);
    colormap jet;
    if normalized
        hold on 
        freezeColors;
        mesh(t_union,Omega,xstochastic);
    else
        set(handles.Wait,'String','The level is defined for normalized cross-spectrogram in here!');
        set(handles.uipanellevel,'visible','on');
    end   
    LSCWS3D2=2;
elseif LSCWS3D2==2
    mesh(t_union,Omega,xsp);
    colormap jet;
    if normalized
        hold on 
        freezeColors;
        mesh(t_union,Omega,xstochastic);
    else
        set(handles.Wait,'String','The level is defined for normalized cross-spectrogram in here!');
        set(handles.uipanellevel,'visible','on');
    end
    LSCWS3D2=3;
elseif LSCWS3D2==3
    surf(t_union,Omega,xsp);
    colormap jet;
    if normalized
        hold on 
        freezeColors;
        surf(t_union,Omega,xstochastic);
    else
        set(handles.Wait,'String','The level is defined for normalized cross-spectrogram in here!');
        set(handles.uipanellevel,'visible','on');
    end
    LSCWS3D2=1;
end
if normalized
    [Lo,Lt]=size(xspec);
    LB=min(min(xstochastic));
    UB=xstochastic(Lo, floor(Lt/2));
    if UB>LB
        caxis([LB UB+0.1*(UB-LB)]);
    end
    colormap gray;
    handles.xstochSave=xstochastic;
    set(handles.Wait,'String','LSCWS (percentage variance) and its stochastic surface');
    pause(.05);
else
    set(handles.Wait,'String','LSCWS (variance: signals)');
    pause(0.05);
end
shading flat; 
axis tight;
ylabel('Cyclic frequency','fontsize',10); 
xlabel('Time','fontsize',10);
%hold off;
%set(gca, 'FontSize' ,12) ;  %Set font size of x and y axis values
align([handles.axes1 handles.axes2],'VerticalAlignment','none', 'HorizontalAlignment','center');
handles.cbar2=0;
handles.cbar=0;
handles.colorbar_VarordB=1;
handles.colorbar_check=2;
handles.color_bar=2;
handles.xspec=xspec;
handles.LSWSorLSCWSstoch=2;
handles.LSCWS3D2=LSCWS3D2;
guidata(hObject,handles);


% --------------------------------------------------------------------
function LSCS_Callback(hObject, eventdata, handles)
% hObject    handle to LSCS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function LSCS_PV_Callback(hObject, eventdata, handles)
% hObject    handle to LSCS_PV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
if ~isfield(handles,'xspectrum')
    cla(handles.axes5,'reset');
    set(handles.axes5,'visible','off');
    set(handles.Wait,'String','Please push LSCSA first in the toolbar');
    pause(.05);
    return;
end
cla(handles.axes5,'reset');
xspectrum=handles.xspectrum;
Omega=handles.Omega;
normalized=get(handles.Normalizedcheckbox,'Value');
axes(handles.axes5);
%figure;
if normalized
    xspec=xspectrum*100;
else
    xspec=xspectrum;
end
stem(Omega,xspec, '.-b', 'MarkerSize',10,'linewidth',0.1);
hold on
plot(Omega,xspec); 
hold off
axis tight;
if normalized
    ylabel('Percentage variance (%)','fontsize',10); 
    set(handles.Wait,'String','LSCS (percentage variance)');
    pause(0.05);
else
    ylabel('Variance (signals)','fontsize',10);
    set(handles.Wait,'String','LSCS (variance: signals)');
    pause(0.05);
end
xlabel('Cyclic frequency','fontsize',10); 
%set(gca, 'FontSize' ,10) ;  %Set font size of x and y axis values to 10
guidata(hObject,handles);

% --------------------------------------------------------------------
function LSCS_confidence_Callback(hObject, eventdata, handles)
% hObject    handle to LSCS_confidence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
if ~isfield(handles,'xspectrum')
    cla(handles.axes5,'reset');
    set(handles.axes5,'visible','off');
    set(handles.Wait,'String','Please push LSCSA first in the toolbar');
    pause(.05);
    return;
end
cla(handles.axes5,'reset');
xspectrum=handles.xspectrum;
LSCSlevel=handles.LSCSlevel;
Omega=handles.Omega;
normalized=get(handles.Normalizedcheckbox,'Value');
axes(handles.axes5);
if normalized
    xspec=xspectrum*100;
else
    xspec=xspectrum;
end
stem(Omega,xspec, '.-b', 'MarkerSize',10,'linewidth',0.1);
hold on;
plot(Omega,xspec); 
if normalized 
    ylabel('Percentage variance (%)','fontsize',10); 
    plot(Omega,LSCSlevel,'r');
    handles.LSCS_CV=LSCSlevel(1);
    set(handles.Wait,'String','LSCS (percentage variance) and its critical threshold');
    pause(0.05);
else
    ylabel('Variance (signals)','fontsize',10);
    set(handles.Wait,'String','The level is defined for normalized cross-spectrum in here!');
    set(handles.uipanellevel,'visible','on');
end
hold off;
axis tight;
xlabel('Cyclic frequency','fontsize',10); 
%set(gca, 'FontSize' ,10) ;  %Set font size of x and y axis values to 10
guidata(hObject,handles);


% --------------------------------------------------------------------
function Phases_Callback(hObject, eventdata, handles)
% hObject    handle to Phases (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
if ~isfield(handles,'xspec')
    cla(handles.axes2,'reset');
    set(handles.axes2,'visible','off');
    set(handles.Wait,'String','Please push LSCWA first in the toolbar');
    pause(0.05);
    return;
end
set(handles.colorbarvariance,'state','off');
set(handles.colorbarstochastic,'state','off');
cla(handles.axes2,'reset');
set(handles.axes2,'visible','off');
cla(handles.axes3,'reset');
set(handles.axes3,'visible','off');
cla(handles.axes4,'reset');
set(handles.axes4,'visible','off');
phase_diff=handles.phase_diff;
t_union=handles.t_union;
Omega=handles.Omega;
w=zeros;
Lt=length(t_union);
Lo=length(Omega);
xspec=handles.xspec;
xstochastic=handles.xstochastic;
%--------------------------------------------
% Re-scale for nicer display for the phase arrows only
x_t=t_union*(Omega(Lo)/t_union(Lt));

if ~isfield(handles,'X_arrow') || ~isfield(handles,'Y_arrow')
    set(handles.uipanelarrow,'visible','on');
    set(handles.Wait,'String','Please set up the density of the arrows and try again');
    pause(.5);
    return;
end

NumXArr=handles.X_arrow;
NumYArr=handles.Y_arrow;
u=zeros; v=zeros;

% Increment toward x-axis for phase arrows
t_inc=floor((Lt+1)/(NumXArr+1));
tt1=t_inc:t_inc:Lt-t_inc+1;
Ltt1=length(tt1);

% Increment toward y-axis for phase arrows
Om_inc=floor((Lo+1)/(NumYArr+1));
Om1=Om_inc:Om_inc:Lo-Om_inc+1;
LOm1=length(Om1);

xx=zeros(LOm1,Ltt1); 
yy=zeros(LOm1,Ltt1);
zz=zeros(LOm1,Ltt1); 

for i=1:LOm1
    xx(i,:)=x_t(tt1);
end
for i=1:LOm1
    yy(i,:)=repmat(Omega(Om1(i)),1,length(tt1));
end
for i=1:LOm1
    for j=1:Ltt1
       u(i,j)=cos(phase_diff(Om1(i),tt1(j)));
       v(i,j)=sin(phase_diff(Om1(i),tt1(j)));
    end
end

cla(handles.axes2,'reset');
axes(handles.axes2);
LSCWS3D3=handles.LSCWS3D3;
normalized=get(handles.Normalizedcheckbox,'Value');
if normalized
    xspec=xspec*100;
end
if LSCWS3D3==1
    surf(x_t,Omega,xspec);
    colormap jet;
    hold on
    if normalized 
        freezeColors;
        mesh(x_t,Omega,xstochastic);
    else
        set(handles.Wait,'String','The level is defined for normalized cross-spectrogram in here!');
        set(handles.uipanellevel,'visible','on');
    end 
    LSCWS3D3=2;
elseif LSCWS3D3==2
    mesh(x_t,Omega,xspec);
    colormap jet;
    hold on
    if normalized
        freezeColors;
        mesh(x_t,Omega,xstochastic);
    else
        set(handles.Wait,'String','The level is defined for normalized cross-spectrogram in here!');
        set(handles.uipanellevel,'visible','on');
    end 
    LSCWS3D3=3;
else
    surf(x_t,Omega,xspec);
    colormap jet;
    hold on
    if normalized
        freezeColors;
        surf(x_t,Omega,xstochastic);
    else
        set(handles.Wait,'String','The level is defined for normalized cross-spectrogram in here!');
        set(handles.uipanellevel,'visible','on');
    end 
    LSCWS3D3=1; 
end

shading flat; 
axis tight;
if normalized
    colormap gray;
    LB=min(min(xstochastic));
    UB=xstochastic(Lo, floor(Lt/2));
    if UB>LB
        caxis([LB UB+0.1*(UB-LB)]);
    end
    set(handles.Wait,'String','LSCWS (PV) and its stochastic surface and phase arrows');
    pause(.05);
    z=ones(LOm1,Ltt1)*100;
    handles.xstochSave=xstochastic;
else
    z=ones(LOm1,Ltt1)*max(max(xspec));
    set(handles.Wait,'String','LSCWS (variance: signals) and phase arrows');
    pause(.05);
end

xticks(x_t(tt1));
if t_union(tt1(Ltt1))-t_union(tt1(1))<1
    xticklabels(round(t_union(tt1),2));
    xtickangle(45);
else
    xticklabels(round(t_union(tt1),1));
end
axis tight;
ylabel('Cyclic frequency','fontsize',10); 
xlabel('Time','fontsize',10);

quiver3(xx,yy,z,u,v,zz,.4,'w','LineWidth',1.5);
axis tight;
hold off;

%set(gca, 'FontSize' ,12) ;  %Set font size of x and y axis values
align([handles.axes1 handles.axes2],'VerticalAlignment','none', 'HorizontalAlignment','center');
handles.cbar2=0;
handles.cbar=0;
handles.colorbar_VarordB=1;
handles.colorbar_check=2;
handles.color_bar=2;
handles.LSWSorLSCWSstoch=2;
handles.w=w;
handles.LSCWS3D3=LSCWS3D3;
guidata(hObject,handles);


% --------------------------------------------------------------------
function PhasesLSCS_Callback(hObject, eventdata, handles)
% hObject    handle to PhasesLSCS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
if ~isfield(handles,'xspectrum')
    cla(handles.axes5,'reset');
    set(handles.axes5,'visible','off');
    set(handles.Wait,'String','Please push LSCSA first in the toolbar');
    pause(0.05);
    return;
end

phase_diff_cross_spectrum=handles.phase_diff_cross_spectrum;
Omega=handles.Omega;

cla(handles.axes5,'reset');
set(handles.axes5,'visible','off');  
axes(handles.axes5); 

phase_diff_degree=phase_diff_cross_spectrum*180/pi;
plot(Omega,phase_diff_degree,'-db','MarkerFaceColor','g','MarkerSize',2);
axis tight;
xlabel('Cyclic frequency','fontsize',10); 
ylabel('Phase shifts in degree','fontsize',10);
handles.LSCSAPhases=phase_diff_degree;
guidata(hObject,handles);



function SR2_Callback(hObject, eventdata, handles)
% hObject    handle to SR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SR2 as text
%        str2double(get(hObject,'String')) returns contents of SR2 as a double


% --- Executes during object creation, after setting all properties.
function SR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function L_0_Callback(hObject, eventdata, handles)
% hObject    handle to L_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L_0 as text
%        str2double(get(hObject,'String')) returns contents of L_0 as a double


% --- Executes during object creation, after setting all properties.
function L_0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function L_1_Callback(hObject, eventdata, handles)
% hObject    handle to L_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L_1 as text
%        str2double(get(hObject,'String')) returns contents of L_1 as a double


% --- Executes during object creation, after setting all properties.
function L_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on datum and none of its controls.
function datum_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to datum (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in XPanelArrow.
function XPanelArrow_Callback(hObject, eventdata, handles)
% hObject    handle to XPanelArrow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipanelarrow,'visible','off');
guidata(hObject, handles);

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
content=get(hObject,'Value');
handles.X_arrow=content;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
content=get(hObject,'Value');
handles.Y_arrow=content;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Arrow_Density_Callback(hObject, eventdata, handles)
% hObject    handle to Arrow_Density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipanelforce,'visible','off');
set(handles.uipanellevel,'visible','off');
set(handles.uipanelwindow,'visible','off');
set(handles.uipanelaxes,'visible','off');
set(handles.uipanelarrow,'visible','on');
set(handles.uipanelfreq,'visible','off');
guidata(hObject, handles);


% --------------------------------------------------------------------
function ALLSSA_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ALLSSA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject); 

if ~isfield(handles,'FirstTS')
    set(handles.Wait,'String','Please open a file containing time/data series');
      pause(.05);
      return;
end

set(handles.uipanelforce,'visible','off');
set(handles.uipanellevel,'visible','off');
set(handles.uipanelwindow,'visible','off');
set(handles.uipanelarrow,'visible','off');
set(handles.uipanelfreq,'visible','off');
cla(handles.axes5,'reset');
set(handles.axes5,'visible','off');
T_S=handles.FirstTS;
set(handles.Wait,'String','Please wait (it may take a while)...');
pause(0.05);
guidata(hObject, handles);
Omega=handles.Omega;
if contains(T_S,'.xlsx') || contains(T_S,'.xls')
    Data=xlsread(T_S);
else
    Data=load(T_S);
end
tt = Data(:,1);
t1=tt-Data(1,1);
f = Data(:,2);
C_f=zeros;
[~,n_c]=size(Data);
if  (n_c>3)
     for i=1:length(t1)
        for j=1:length(t1)
            C_f(i,j)=Data(i,j+2);
        end
     end
     P_f=inv(C_f);
elseif n_c==3
     P_f =1./(Data(:,3).^2);
elseif n_c==2
     P_f =1;
end
%--------------------------------------------
content=get(handles.SignificanceLevel,'Value');
if content==1
    alpha=0.01;
elseif content==2
    alpha=0.05;
else
    alpha=0.1;
end
%--------------------------------------------
QQ(1)=get(handles.datum,'Value');
QQ(2)=get(handles.lineartrend,'Value');
QQ(3)=get(handles.quadratictrend,'Value');
QQ(4)=get(handles.cubictrend,'Value');
force_sinusoids=get(handles.force,'String');
FS=str2num(force_sinusoids);
if isnan(FS)
    Q=QQ;
else
    Q=[QQ,FS];
end

IndexDatum1=get(handles.DatumShifts1,'String');
IndexDatum1=str2num(IndexDatum1);
if isnan(IndexDatum1)
   IndexDatum1=[];
end
Lt1=length(t1);
if any(floor(IndexDatum1)-IndexDatum1) || any(IndexDatum1>=Lt1) || any(IndexDatum1<1)
    set(handles.Wait,'String','Error! Datum shift indices must be natural numbers less than the size of series');
    set(handles.uipanelforce,'visible','on');
    return;
end
%-----------Run the ALLSSA---------------------
[TrendsCoeff, SpecCosSin, r, Qnew, COVcoeff]=ALLSSA(t1, f, P_f, Omega, Q, IndexDatum1, alpha);
%----------------------------------------------
NewOmega=Qnew(5:length(Qnew));
ALLSS=zeros; u=1;
for k=1:2:length(SpecCosSin)
    chat=[SpecCosSin(k), SpecCosSin(k+1)];
    ALLSS(u)=sqrt(chat*chat');
    u=u+1;
end

axes(handles.axes5);
stem(NewOmega,ALLSS, '.-b', 'MarkerSize',10,'linewidth',0.1);

ylabel('Amplitude','fontsize',10);
xlabel('Cyclic frequency','fontsize',10); 
%set(gca, 'FontSize' ,10) ;  %Set font size of x and y axis values to 10
set(handles.Wait,'String','Antileakage least-squares spectrum');
pause(0.05);
handles=guidata(hObject);
handles.ALLSSA_Q = Q;
handles.ALLSSA_IndexDatum = IndexDatum1;
handles.SpecCosSin=SpecCosSin;
handles.TrendsCoeff=TrendsCoeff;
residual(:,1)=tt;
residual(:,2)=r;
if ~isempty(COVcoeff)
    handles.r_ALLSSA = residual;
    handles.ALLSSA_COVcoeff = COVcoeff;
end
handles.OptimalFrequencies=NewOmega;
handles.t1=t1;
guidata(hObject, handles);


% --- Executes on button press in MorletCheck.
function MorletCheck_Callback(hObject, eventdata, handles)
% hObject    handle to MorletCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MorletCheck
MorlCoeff=get(handles.MorletCheck,'Value');

if MorlCoeff~=0
    set(handles.L1,'String','6');
    set(handles.L_1,'String','6');
    set(handles.L0,'String','0');
    set(handles.L_0,'String','0');
else
    set(handles.L1,'String','2');
    set(handles.L0,'String','20');
    set(handles.L_1,'String','2');
    set(handles.L_0,'String','20');
end

guidata(hObject, handles);


% --- Executes on button press in Normalizedcheckbox.
function Normalizedcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to Normalizedcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
% Hint: get(hObject,'Value') returns toggle state of Normalizedcheckbox
handles = guidata(hObject);
if isfield(handles,'spectrum')
    handles = rmfield( handles, 'spectrum' );
end
if isfield(handles,'spec')
    handles = rmfield( handles, 'spec' );
end
if isfield(handles,'xspectrum')
    handles = rmfield( handles, 'xspectrum' );
end
if isfield(handles,'xspec')
    handles = rmfield( handles, 'xspec' );
end
guidata(hObject, handles);

% --------------------------------------------------------------------
function Save_work_Callback(hObject, eventdata, handles)
% hObject    handle to Save_work (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function SaveWorkspace_Callback(hObject, eventdata, handles)
% hObject    handle to SaveWorkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------
handles = guidata(hObject); 

%-----------LSSA save-----------------
if isfield(handles,'spectrum')
    LSSA_COVcoeff=(handles.LSSA_COVcoeff);
    spectrum=(handles.spectrum)';
    filter = {'*.mat';'*.xlsx';'*.dat';'*.txt'};
    [file, path] = uiputfile(filter,'File Selection','LSS');
    FullFileName = fullfile(path,file);
    LenFile=length(file);
    if LenFile>3
        if file(LenFile)=='x'
            xlswrite(FullFileName,spectrum);
        elseif file(LenFile-2:LenFile)=='mat'
            save(FullFileName, 'spectrum');   
        else
            save(FullFileName,'spectrum','-ascii');
        end
    end
end
LSSA_Est{1,1}='Intercept';
LSSA_Est{1,2}='Error';
LSSA_Est{1,3}='Slope';
LSSA_Est{1,4}='Error';
LSSA_Est{1,5}='Quad_Coeff';
LSSA_Est{1,6}='Error';
LSSA_Est{1,7}='Cub_Coeff';
LSSA_Est{1,8}='Error';
LSSA_Est{1,9}='Frequency';
LSSA_Est{1,10}='Cos_Coeff';
LSSA_Est{1,11}='Error';
LSSA_Est{1,12}='Sin_Coeff';
LSSA_Est{1,13}='Error';
LSSA_Est{1,14}='Amplitude';
LSSA_Est{1,15}='Error';
LSSA_Est{1,16}='Phase';
LSSA_Est{1,17}='Error';

if isfield(handles,'LSSA_TrendsCoeff')
    LSSA_TrendsCoeff=handles.LSSA_TrendsCoeff;
    LSSA_Q = handles.LSSA_Q;
    IndexDatum1 = handles.LSSA_IndexDatum;
    if isnan(IndexDatum1)
       IndexDatum1=[];
    end
    Lid=length(IndexDatum1);
    LSSA_Est{2,1}=0;LSSA_Est{2,2}=0;LSSA_Est{2,3}=0;LSSA_Est{2,4}=0;
    LSSA_Est{2,5}=0;LSSA_Est{2,6}=0;LSSA_Est{2,7}=0;LSSA_Est{2,8}=0;
    if LSSA_Q(1)==1
        for k=1:Lid+1
            LSSA_Est{k+1,1}=LSSA_TrendsCoeff(k);
            LSSA_Est{k+1,2}=sqrt(LSSA_COVcoeff(k,k));
        end
        k=k+1;
    else
        k=1;
    end
    
    if LSSA_Q(2)==1
        LSSA_Est{2,3}=LSSA_TrendsCoeff(k);
        LSSA_Est{2,4}=sqrt(LSSA_COVcoeff(k,k));
        k=k+1;
    end
    if LSSA_Q(3)==1
        LSSA_Est{2,5}=LSSA_TrendsCoeff(k);
        LSSA_Est{2,6}=sqrt(LSSA_COVcoeff(k,k));
        k=k+1;
    end
    if LSSA_Q(4)==1
        LSSA_Est{2,7}=LSSA_TrendsCoeff(k);
        LSSA_Est{2,8}=sqrt(LSSA_COVcoeff(k,k));
        k=k+1;
    end
end

if isfield(handles,'LSSA_CosSin')
    LSSA_CosSin=handles.LSSA_CosSin;
    force_sinusoids=get(handles.force,'String');
    FW=str2num(force_sinusoids);
    for j=1:length(FW)
        LSSA_Est{j+1,9}=FW(j);
        c=LSSA_CosSin(2*j-1); s=LSSA_CosSin(2*j); 
        ce=sqrt(LSSA_COVcoeff(k,k)); se=sqrt(LSSA_COVcoeff(k+1,k+1));
        cse=sqrt(LSSA_COVcoeff(k,k+1));
        LSSA_Est{j+1,10}=c;
        LSSA_Est{j+1,11}=ce;
        LSSA_Est{j+1,12}=s;
        LSSA_Est{j+1,13}=se;
        amp=sqrt(c^2+s^2); 
        ampE=sqrt((c*ce)^2+(s*se)^2 + 2*c*s*cse)/amp;
        LSSA_Est{j+1,14}= amp;
        LSSA_Est{j+1,15}= ampE;
        LSSA_Est{j+1,16}= 2*atan((amp-s)/c);
        if c<1e-7
            c=1e-7;
        end
        phaseE=2*sqrt(((s/c)*(amp-s))^2*ce^2+(s-amp)^2*se^2-2*(s/c)*(amp-s)^2*cse)/(abs(amp*c)*(1+(amp-s)^2/c^2));
        LSSA_Est{j+1,17}= phaseE;
        k=k+2;
    end
end


if isfield(handles,'LSSA_CosSin') || isfield(handles,'LSSA_TrendsCoeff')
    filter = {'*.mat';'*.xlsx';'*.dat';'*.txt'};
    [file, path] = uiputfile(filter,'File Selection','LSSA_Estimates');
    FullFileName = fullfile(path,file);
    LenFile=length(file);
    if LenFile>3
        if file(LenFile)=='x'
            xlswrite(FullFileName,LSSA_Est);
        elseif file(LenFile-2:LenFile)=='mat'
            save(FullFileName, 'LSSA_Est');   
        else
            fileID = fopen(FullFileName,'w');
            firstrow = '%9s %9s %12s %9s %16s %9s %16s %9s %17s %15s %9s %15s %9s %15s %9s %11s %9s \r\n';
            fprintf(fileID,firstrow,LSSA_Est{1,:});
            secondrow = '%9.3f %9.3f %12.3f %9.3f %16.3f %9.3f %16.3f %9.3f %17.3f %15.3f %9.3f %15.3f %9.3f %15.3f %9.3f %11.3f %9.3f\r\n';
            otherrows = '%9.3f %9.3f';
            SizeA=size(LSSA_Est); 
            fprintf(fileID,secondrow,LSSA_Est{2,:});
            for u=3:SizeA(1)
                if isempty(LSSA_Est{u,1}) 
                    fprintf(fileID,'%114.3f %15.3f %9.3f %15.3f %9.3f %15.3f %9.3f %11.3f %9.3f\r\n',LSSA_Est{u,9:17});
                else
                    fprintf(fileID,otherrows,LSSA_Est{u,1:2});
                    fprintf(fileID,'%95.3f %15.3f %9.3f %15.3f %9.3f %15.3f %9.3f %11.3f %9.3f\r\n',LSSA_Est{u,9:17});
                end
            end
            fclose(fileID);
        end
    end
end
if isfield(handles,'lssCV')
    lssCV=(handles.lssCV);
    filter = {'*.mat';'*.xlsx';'*.dat';'*.txt'};
    [file, path] = uiputfile(filter,'File Selection','LSSA_CV');
    FullFileName = fullfile(path,file);
    LenFile=length(file);
    if LenFile>3
        if file(LenFile)=='x'
            xlswrite(FullFileName,lssCV);
        elseif file(LenFile-2:LenFile)=='mat'
            save(FullFileName, 'lssCV');   
        else
            save(FullFileName,'lssCV','-ascii');
        end
    end
end

if isfield(handles,'LSSA_COVcoeff')
    r_LSSA=(handles.r_LSSA);
    filter = {'*.mat';'*.xlsx';'*.dat';'*.txt'};
    [file, path] = uiputfile(filter,'File Selection','LSSA_residual');
    FullFileName = fullfile(path,file);
    LenFile=length(file);
    if LenFile>3
        if file(LenFile)=='x'
            xlswrite(FullFileName,r_LSSA);
        elseif file(LenFile-2:LenFile)=='mat'
            save(FullFileName, 'r_LSSA');   
        else
            save(FullFileName,'r_LSSA','-ascii');
        end
    end
end
if isfield(handles,'LSSA_COVcoeff')
    LSSA_COVcoeff=(handles.LSSA_COVcoeff);
    filter = {'*.mat';'*.xlsx';'*.dat';'*.txt'};
    [file, path] = uiputfile(filter,'File Selection','LSSA_COVcoeff');
    FullFileName = fullfile(path,file);
    LenFile=length(file);
    if LenFile>3
        if file(LenFile)=='x'
            xlswrite(FullFileName,LSSA_COVcoeff);
        elseif file(LenFile-2:LenFile)=='mat'
            save(FullFileName, 'LSSA_COVcoeff');   
        else
            save(FullFileName,'LSSA_COVcoeff','-ascii');
        end
    end
end
%------------LSCSA save------------------------
if isfield(handles,'xspectrum')
    xspect=(handles.xspectrum);
    xspectrum=xspect(1);
    filter = {'*.mat';'*.xlsx';'*.dat';'*.txt'};
    [file, path] = uiputfile(filter,'File Selection','LSCS');
    FullFileName = fullfile(path,file);
    LenFile=length(file);
    if LenFile>3
        if file(LenFile)=='x'
            xlswrite(FullFileName,xspectrum);
        elseif file(LenFile-2:LenFile)=='mat'
            save(FullFileName, 'xspectrum');   
        else
            save(FullFileName,'xspectrum','-ascii');
        end
    end
end
if isfield(handles,'LSCS_CV')
    LSCS_CV=(handles.LSCS_CV);
    filter = {'*.mat';'*.xlsx';'*.dat';'*.txt'};
    [file, path] = uiputfile(filter,'File Selection','LSCSA_CV');
    FullFileName = fullfile(path,file);
    LenFile=length(file);
    if LenFile>3
        if file(LenFile)=='x'
            xlswrite(FullFileName,LSCS_CV);
        elseif file(LenFile-2:LenFile)=='mat'
            save(FullFileName, 'LSCS_CV');   
        else
            save(FullFileName,'LSCS_CV','-ascii');
        end
    end
end
if isfield(handles,'LSCSAPhases')
    LSCSAPhases=(handles.LSCSAPhases)';
    filter = {'*.mat';'*.xlsx';'*.dat';'*.txt'};
    [file, path] = uiputfile(filter,'File Selection','LSCSA_Phases');
    FullFileName = fullfile(path,file);
    LenFile=length(file);
    if LenFile>3
        if file(LenFile)=='x'
            xlswrite(FullFileName,LSCSAPhases);
        elseif file(LenFile-2:LenFile)=='mat'
            save(FullFileName, 'LSCSAPhases');   
        else
            save(FullFileName,'LSCSAPhases','-ascii');
        end
    end
end
%---------------LSWA save---------------------
if isfield(handles,'spec')
    spec=(handles.spec);
    filter = {'*.mat';'*.xlsx';'*.dat';'*.txt'};
    [file, path] = uiputfile(filter,'File Selection','LSWS');
    FullFileName = fullfile(path,file);
    LenFile=length(file);
    if LenFile>3
        if file(LenFile)=='x'
            xlswrite(FullFileName,spec);
        elseif file(LenFile-2:LenFile)=='mat'
            save(FullFileName, 'spec');   
        else
            save(FullFileName,'spec','-ascii');
        end
    end
end
if isfield(handles,'stochSave')
    stochastic=(handles.stochSave);
    filter = {'*.mat';'*.xlsx';'*.dat';'*.txt'};
    [file, path] = uiputfile(filter,'File Selection','LSWA_SS');
    FullFileName = fullfile(path,file);
    LenFile=length(file);
    if LenFile>3
        if file(LenFile)=='x'
            xlswrite(FullFileName,stochastic);
        elseif file(LenFile-2:LenFile)=='mat'
            save(FullFileName, 'stochastic');   
        else
            save(FullFileName,'stochastic','-ascii');
        end
    end
end

%---------------LSCWA save-------------------------
if isfield(handles,'xspec')
    xspec=(handles.xspec);
    filter = {'*.mat';'*.xlsx';'*.dat';'*.txt'};
    [file, path] = uiputfile(filter,'File Selection','LSCWS');
    FullFileName = fullfile(path,file);
    LenFile=length(file);
    if LenFile>3
        if file(LenFile)=='x'
            xlswrite(FullFileName,xspec);
        elseif file(LenFile-2:LenFile)=='mat'
            save(FullFileName, 'xspec');   
        else
            save(FullFileName,'xspec','-ascii');
        end
    end
end
if isfield(handles,'xstochSave')
    xstochastic=(handles.xstochSave)./100;
    filter = {'*.mat';'*.xlsx';'*.dat';'*.txt'};
    [file, path] = uiputfile(filter,'File Selection','LSCWA_SS');
    FullFileName = fullfile(path,file);
    LenFile=length(file);
    if LenFile>3
        if file(LenFile)=='x'
            xlswrite(FullFileName,xstochastic);
        elseif file(LenFile-2:LenFile)=='mat'
            save(FullFileName, 'xstochastic');   
        else
            save(FullFileName,'xstochastic','-ascii');
        end
    end
end
if isfield(handles,'phase_diff')
    phase_diff=(handles.phase_diff);
    filter = {'*.mat';'*.xlsx';'*.dat';'*.txt'};
    [file, path] = uiputfile(filter,'File Selection','LSCWA_Phases');
    FullFileName = fullfile(path,file);
    LenFile=length(file);
    if LenFile>3
        if file(LenFile)=='x'
            xlswrite(FullFileName,phase_diff);
        elseif file(LenFile-2:LenFile)=='mat'
            save(FullFileName, 'phase_diff');   
        else
            save(FullFileName,'phase_diff','-ascii');
        end
    end
end
%----------------ALLSSA save---------------------------
ALLSSA_Est{1,1}='Intercept';
ALLSSA_Est{1,2}='Error';
ALLSSA_Est{1,3}='Slope';
ALLSSA_Est{1,4}='Error';
ALLSSA_Est{1,5}='Quad_Coeff';
ALLSSA_Est{1,6}='Error';
ALLSSA_Est{1,7}='Cub_Coeff';
ALLSSA_Est{1,8}='Error';
ALLSSA_Est{1,9}='Frequency';
ALLSSA_Est{1,10}='Cos_Coeff';
ALLSSA_Est{1,11}='Error';
ALLSSA_Est{1,12}='Sin_Coeff';
ALLSSA_Est{1,13}='Error';
ALLSSA_Est{1,14}='Amplitude';
ALLSSA_Est{1,15}='Error';
ALLSSA_Est{1,16}='Phase';
ALLSSA_Est{1,17}='Error';

if isfield(handles,'TrendsCoeff')
    ALLSSA_COVcoeff=(handles.ALLSSA_COVcoeff);
    ALLSSA_TrendsCoeff=handles.TrendsCoeff;
    ALLSSA_Q = handles.ALLSSA_Q;
    IndexDatum1 = handles.ALLSSA_IndexDatum;
    if isnan(IndexDatum1)
       IndexDatum1=[];
    end
    Lid=length(IndexDatum1);
    ALLSSA_Est{2,1}=0;ALLSSA_Est{2,2}=0;ALLSSA_Est{2,3}=0;ALLSSA_Est{2,4}=0;
    ALLSSA_Est{2,5}=0;ALLSSA_Est{2,6}=0;ALLSSA_Est{2,7}=0;ALLSSA_Est{2,8}=0;
    if ALLSSA_Q(1)==1
        for k=1:Lid+1
            ALLSSA_Est{k+1,1}=ALLSSA_TrendsCoeff(k);
            ALLSSA_Est{k+1,2}=sqrt(ALLSSA_COVcoeff(k,k));
        end
        k=k+1;
    else
        k=1;
    end
    
    if ALLSSA_Q(2)==1
        ALLSSA_Est{2,3}=ALLSSA_TrendsCoeff(k);
        ALLSSA_Est{2,4}=sqrt(ALLSSA_COVcoeff(k,k));
        k=k+1;
    end
    if ALLSSA_Q(3)==1
        ALLSSA_Est{2,5}=ALLSSA_TrendsCoeff(k);
        ALLSSA_Est{2,6}=sqrt(ALLSSA_COVcoeff(k,k));
        k=k+1;
    end
    if ALLSSA_Q(4)==1
        ALLSSA_Est{2,7}=ALLSSA_TrendsCoeff(k);
        ALLSSA_Est{2,8}=sqrt(ALLSSA_COVcoeff(k,k));
        k=k+1;
    end
end

if isfield(handles,'SpecCosSin')
    ALLSSA_CosSin=handles.SpecCosSin;
    FW=handles.OptimalFrequencies;
    for j=1:length(FW)
        ALLSSA_Est{j+1,9}=FW(j);
        c=ALLSSA_CosSin(2*j-1); s=ALLSSA_CosSin(2*j);
        ce=sqrt(ALLSSA_COVcoeff(k,k)); se=sqrt(ALLSSA_COVcoeff(k+1,k+1));
        cse=sqrt(ALLSSA_COVcoeff(k,k+1));
        ALLSSA_Est{j+1,10}=c;
        ALLSSA_Est{j+1,11}=ce;
        ALLSSA_Est{j+1,12}=s;
        ALLSSA_Est{j+1,13}=se;
        amp=sqrt(c^2+s^2); 
        ampE=sqrt((c*ce)^2+(s*se)^2+ 2*c*s*cse)/amp;
        ALLSSA_Est{j+1,14}= amp;
        ALLSSA_Est{j+1,15}= ampE;
        ALLSSA_Est{j+1,16}= 2*atan((amp-s)/c);
        if c<1e-7
            c=1e-7;
        end
        phaseE=2*sqrt(((s/c)*(amp-s))^2*ce^2+(s-amp)^2*se^2-2*(s/c)*(amp-s)^2*cse)/(abs(amp*c)*(1+(amp-s)^2/c^2));
        ALLSSA_Est{j+1,17}= phaseE;
        k=k+2;
    end
end


if isfield(handles,'SpecCosSin') || isfield(handles,'TrendsCoeff')
    filter = {'*.mat';'*.xlsx';'*.dat';'*.txt'};
    [file, path] = uiputfile(filter,'File Selection','ALLSSA_Estimates');
    FullFileName = fullfile(path,file);
    LenFile=length(file);
    if LenFile>3
        if file(LenFile)=='x'
            xlswrite(FullFileName,ALLSSA_Est);
        elseif file(LenFile-2:LenFile)=='mat'
            save(FullFileName, 'ALLSSA_Est');   
        else
            fileID = fopen(FullFileName,'w');
            firstrow = '%9s %9s %12s %9s %16s %9s %16s %9s %17s %15s %9s %15s %9s %15s %9s %11s %9s \r\n';
            fprintf(fileID,firstrow,ALLSSA_Est{1,:});
            secondrow = '%9.3f %9.3f %12.3f %9.3f %16.3f %9.3f %16.3f %9.3f %17.3f %15.3f %9.3f %15.3f %9.3f %15.3f %9.3f %11.3f %9.3f\r\n';
            otherrows = '%9.3f %9.3f';
            SizeA=size(ALLSSA_Est); 
            fprintf(fileID,secondrow,ALLSSA_Est{2,:});
            for u=3:SizeA(1)
                if isempty(ALLSSA_Est{u,1}) 
                    fprintf(fileID,'%114.3f %15.3f %9.3f %15.3f %9.3f %15.3f %9.3f %11.3f %9.3f\r\n',ALLSSA_Est{u,9:17});
                else
                    fprintf(fileID,otherrows,ALLSSA_Est{u,1:2});
                    fprintf(fileID,'%95.3f %15.3f %9.3f %15.3f %9.3f %15.3f %9.3f %11.3f %9.3f\r\n',ALLSSA_Est{u,9:17});
                end
            end
            fclose(fileID);
        end
    end
end
if isfield(handles,'ALLSSA_COVcoeff')
    r_ALLSSA=(handles.r_ALLSSA);
    filter = {'*.mat';'*.xlsx';'*.dat';'*.txt'};
    [file, path] = uiputfile(filter,'File Selection','ALLSSA_residual');
    FullFileName = fullfile(path,file);
    LenFile=length(file);
    if LenFile>3
        if file(LenFile)=='x'
            xlswrite(FullFileName,r_ALLSSA);
        elseif file(LenFile-2:LenFile)=='mat'
            save(FullFileName, 'r_ALLSSA');   
        else
            save(FullFileName,'r_ALLSSA','-ascii');
        end
    end
end
if isfield(handles,'ALLSSA_COVcoeff')
    ALLSSA_COVcoeff=(handles.ALLSSA_COVcoeff);
    filter = {'*.mat';'*.xlsx';'*.dat';'*.txt'};
    [file, path] = uiputfile(filter,'File Selection','ALLSSA_COVcoeff');
    FullFileName = fullfile(path,file);
    LenFile=length(file);
    if LenFile>3
        if file(LenFile)=='x'
            xlswrite(FullFileName,ALLSSA_COVcoeff);
        elseif file(LenFile-2:LenFile)=='mat'
            save(FullFileName, 'ALLSSA_COVcoeff');   
        else
            save(FullFileName,'ALLSSA_COVcoeff','-ascii');
        end
    end
end



% --------------------------------------------------------------------
function SaveFigureAs_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFigureAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------

filter = {'*.jpeg';'*.png';'*.tiff'};
[file, path] = uiputfile(filter,'File Selection','Results_LSWA');
FullFileName = fullfile(path,file);
set(gcf, 'Color', 'w');
width=12;%8; %LSCS
height = 8;  % 540 pixels converted to inches
set(gcf, 'PaperUnits', 'inches','PaperPosition',[0 0 width height]);
set(gca,'LooseInset',get(gca,'TightInset'));
Lfile=length(file);
if Lfile>3
    if file(Lfile-2:Lfile)=='png'
        Fname=strcat('-d',file(Lfile-2:Lfile));
    else
        Fname=strcat('-d',file(Lfile-3:Lfile));
    end
    print( Fname,  FullFileName, '-r600');
end


% --------------------------------------------------------------------
function FirstTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to FirstTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile;
TS1 = fullfile(path,file);
if path 
    if contains(file,'.xlsx') || contains(file,'.xls')
        Data=xlsread(TS1);
    else
        Data=load(TS1);
    end
    t1=Data(:,1); Lt1=length(t1);
    M1=(Lt1-1)/(t1(Lt1)-t1(1));
    set(handles.SR,'String',num2str(round(M1,2)));
    LB = round(0.5/(t1(Lt1)-t1(1)),3);
    UB = round(M1/2,3);
    set(handles.FreqLB,'String',num2str(LB));
    set(handles.FreqUB,'String',num2str(UB));
    Num = str2double(get(handles.FreqQuantity,'String'));
    inc = (UB - LB)/Num;
    handles.Omega= LB:inc:UB-inc;
    handles.t1=t1;
    handles.FirstTS=TS1;
end
guidata(hObject, handles);

% --------------------------------------------------------------------
function SecondTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to SecondTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'t1')
    set(handles.Wait,'String','Please open a file containing 1st time/data series first!');
      pause(.5);
      return;
end

[file,path] = uigetfile;
TS2 = fullfile(path,file);
if path
    if contains(file,'.xlsx') || contains(file,'.xls')
        Data2=xlsread(TS2);
    else
        Data2=load(TS2);
    end
    t2=Data2(:,1); Lt2=length(t2);
    M2=(Lt2-1)/(t2(Lt2)-t2(1));
    set(handles.SR2,'String',num2str(round(M2,2)));
    handles.t2=t2;
    t1=handles.t1;
    Lt1=length(t1);
    LB=round(max(0.5/(t1(Lt1)-t1(1)),0.5/(t2(Lt2)-t2(1))),3);
    set(handles.FreqLB,'String',num2str(LB));
    M1 = str2double(get(handles.SR,'String'));
    UB=round(min(M1/2,M2/2),3);
    set(handles.FreqUB,'String',num2str(UB));
    Num = str2double(get(handles.FreqQuantity,'String'));
    inc = (UB - LB)/Num;
    handles.Omega= LB:inc:UB-inc;
    handles.SecondTS=TS2;
end
guidata(hObject, handles);

% --------------------------------------------------------------------
function SetFrequencies_Callback(hObject, eventdata, handles)
% hObject    handle to SetFrequencies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile;
Freq = fullfile(path,file);
A=fopen(Freq);
if A==-1
    set(handles.Wait,'String','The file containing frequencies is not opened! You may still continue...');
    LB = str2double(get(handles.FreqLB,'String'));
    UB = str2double(get(handles.FreqUB,'String'));
    Num = str2double(get(handles.FreqQuantity,'String'));
    inc = (UB - LB)/Num;
    handles.Omega= LB:inc:UB-inc;
    guidata(hObject, handles);
    return;
elseif contains(Freq,'.xlsx') || contains(Freq,'.xls')
    Omega=xlsread(Freq);
else
    Omega=load(Freq);
end
handles.Omega=Omega;
guidata(hObject, handles);


% --------------------------------------------------------------------
function FrequenciesAutomated_Callback(hObject, eventdata, handles)
% hObject    handle to FrequenciesAutomated (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipanelfreq,'visible','on');
set(handles.uipanelforce,'visible','off');
set(handles.uipanellevel,'visible','off');
set(handles.uipanelwindow,'visible','off');
set(handles.uipanelaxes,'visible','off');
set(handles.uipanelarrow,'visible','off');
guidata(hObject, handles);


function FreqLB_Callback(hObject, eventdata, handles)
% hObject    handle to FreqLB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FreqLB as text
%        str2double(get(hObject,'String')) returns contents of FreqLB as a double
LB = str2double(get(handles.FreqLB,'String'));
UB = str2double(get(handles.FreqUB,'String'));
Num = str2double(get(handles.FreqQuantity,'String'));
if 0<LB && LB<UB && 1<Num
    inc = (UB - LB)/Num;
    Omega = LB:inc:UB-inc;
    handles.Omega = Omega;
    set(handles.Wait,'String',strcat('OMEGA=[',num2str(Omega),']'));
else
    set(handles.Wait,'String','Invalid Lower Bound!');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function FreqLB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FreqLB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FreqQuantity_Callback(hObject, eventdata, handles)
% hObject    handle to FreqQuantity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FreqQuantity as text
%        str2double(get(hObject,'String')) returns contents of FreqQuantity as a double
LB = str2double(get(handles.FreqLB,'String'));
UB = str2double(get(handles.FreqUB,'String'));
Num = str2double(get(handles.FreqQuantity,'String'));
if 0<LB && LB<UB && 1<Num
    inc = (UB - LB)/Num;
    Omega = LB:inc:UB-inc;
    handles.Omega = Omega;
    set(handles.Wait,'String',strcat('OMEGA=[',num2str(Omega),']'));
else
    set(handles.Wait,'String','Invalid Quantity!');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function FreqQuantity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FreqQuantity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FreqUB_Callback(hObject, eventdata, handles)
% hObject    handle to FreqUB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FreqUB as text
%        str2double(get(hObject,'String')) returns contents of FreqUB as a double
LB = str2double(get(handles.FreqLB,'String'));
UB = str2double(get(handles.FreqUB,'String'));
Num = str2double(get(handles.FreqQuantity,'String'));
if 0<LB && LB<UB && 1<Num
    inc = (UB - LB)/Num;
    Omega = LB:inc:UB-inc;
    handles.Omega = Omega;
    set(handles.Wait,'String',strcat('OMEGA=[',num2str(Omega),']'));
else
    set(handles.Wait,'String','Invalid Upper Bound!');
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function FreqUB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FreqUB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% if get(handles.checkboxf1,'Value')
%     set(handles.checkboxf2,'Value',0);
% end
% guidata(hObject, handles);

% --- Executes on button press in XpanelFrequencies.
function XpanelFrequencies_Callback(hObject, eventdata, handles)
% hObject    handle to XpanelFrequencies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipanelfreq,'visible','off');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function uipanelfreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanelfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
