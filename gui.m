function varargout = gui(varargin)
% GUI MATLAB code for gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 01-May-2020 10:39:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
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


% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)

% Choose default command line output for gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(gcf,'name','基于二维L阵的高频地波雷达电离层回波方向估计'); 
% n=1:10;m=[12 34 45 67 76 78 65 45 53 52];
% plot(handles.img,n,m,'b-','LineWidth',1.5);
% colorbar;


global g_para
global g_signal
global g_array 
global g_axis_range;
global g_echos;
global dbf_mode;

g_signal.lamda = (3*10^8)/g_signal.freq;

testmode = 'dbf';
g_array.num = 16;
g_array.x_num = 8;
g_array.x_pos = g_array.span : g_array.span : (g_array.x_num)*g_array.span;
g_array.y_pos = 0 : g_array.span : (g_array.y_num-1)*g_array.span;

g_echos.theta.num = 45;
g_echos.phi.num = 45;
g_echos.theta.rad = g_echos.theta.num*g_para.rad;
g_echos.phi.rad = g_echos.phi.num*g_para.rad;
g_echos.num = 1;
g_echos.snr = 10;
g_echos.snapshot = 1000;
g_echos.t = [0:99]/1000;
g_echos.signal=sqrt(10^(g_echos.snr/10))*exp(j*2*pi*g_signal.freq*g_echos.t);  
% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global g_para g_signal g_array g_axis_range;
global g_echos;
global dbf_mode;
global ag_time;

switch dbf_mode
    case 'capon' 
        g_echos.theta.num = [15 15 15 15 15 15 15 15];
        g_echos.phi.num =   [65 55 45 45 35 25 15  5];
        g_echos.theta.rad = g_echos.theta.num*g_para.rad;
        g_echos.phi.rad = g_echos.phi.num*g_para.rad;
        g_echos.num = length(g_echos.theta.num);
        g_echos.snr = get(handles.edit1,'String');
        g_echos.snr = str2num(g_echos.snr);
        g_echos.snapshot = get(handles.edit2,'String'); 
        g_echos.snapshot = str2num(g_echos.snapshot);
%         g_echos.snr = 10;
%         g_echos.snapshot = 100; 
        g_echos.t = [0:99]/1000;
        g_echos.signal=rand(g_echos.num,g_echos.snapshot);
        abs_P=beamforming('capon');
        theta.rad=linspace(-pi/2,pi/2,181);
        theta.num = theta.rad/pi*180;
        phi.rad=linspace(-pi/2,pi/2,181);
        phi.num = phi.rad/pi*180;
        imagesc(handles.img,phi.num,theta.num,abs_P);
        hold on;
        title(handles.img,'二维Capon算法波束扫描功率谱');
        xlabel(handles.img,'方位角phi/degree');
        ylabel(handles.img,'俯仰角theta/degree');    
        colormap(handles.img,jet);
        colorbar;
        [M,~]=contour(phi.num,theta.num,abs_P,[-3.84,-3.84],':');
        [~,mynum]=size(M);
        Mx=M(1,2:mynum);
        My=M(2,2:mynum);
        xmax=max(Mx);
        xmin=min(Mx);
        ymax=max(My);
        ymin=min(My);
        xmin=round(xmin);
        xmax=round(xmax);
        ymin=round(ymin);
        ymax=round(ymax);
        xlabel(handles.img,{['方位角phi/degree'],...
                ['估计方位角：' num2str(xmin) '° ~' num2str(xmax) '°'],...
                ['估计俯仰角：' num2str(ymin) '° ~' num2str(ymax) '°']});
        rectangle(handles.img,'Position',[xmin ymin (xmax-xmin) (ymax-ymin)],'EdgeColor',[1 1 1],'LineWidth',2);
        set(handles.text5,'string',['核心代码运行时间为：',num2str(ag_time),'s']);

    case 'normal'
        g_echos.theta.num = [15 15 15 15 15 15 15 15];
        g_echos.phi.num =   [65 55 45 45 35 25 15  5];
        g_echos.theta.rad = g_echos.theta.num*g_para.rad;
        g_echos.phi.rad = g_echos.phi.num*g_para.rad;
        g_echos.num = length(g_echos.theta.num);
        g_echos.snr = get(handles.edit1,'String');
        g_echos.snr = str2num(g_echos.snr);
        g_echos.snapshot = get(handles.edit2,'String'); 
        g_echos.snapshot = str2num(g_echos.snapshot);
%         g_echos.snr = 10;
%         g_echos.snapshot = 100; 
        g_echos.t = [0:99]/1000;
        g_echos.signal=rand(g_echos.num,g_echos.snapshot);
        abs_P=beamforming('normal');
        theta.rad=linspace(-pi/2,pi/2,181);
        theta.num = theta.rad/pi*180;
        phi.rad=linspace(-pi/2,pi/2,181);
        phi.num = phi.rad/pi*180;
        imagesc(handles.img,phi.num,theta.num,abs_P);
        hold on;
        title(handles.img,'二维延迟相加算法波束扫描功率谱');
        xlabel(handles.img,'方位角phi/degree');
        ylabel(handles.img,'俯仰角theta/degree');    
        colormap(handles.img,jet);
        colorbar;
        [M,~]=contour(phi.num,theta.num,abs_P,[-3.84,-3.84],':');
        [~,mynum]=size(M);
        Mx=M(1,2:mynum);
        My=M(2,2:mynum);
        xmax=max(Mx);
        xmin=min(Mx);
        ymax=max(My);
        ymin=min(My);
        xmin=round(xmin);
        xmax=round(xmax);
        ymin=round(ymin);
        ymax=round(ymax);
        xlabel(handles.img,{['方位角phi/degree'],...
                ['估计方位角：' num2str(xmin) '° ~' num2str(xmax) '°'],...
                ['估计俯仰角：' num2str(ymin) '° ~' num2str(ymax) '°']});
        rectangle(handles.img,'Position',[xmin ymin (xmax-xmin) (ymax-ymin)],'EdgeColor',[1 1 1],'LineWidth',2);
        set(handles.text5,'string',['核心代码运行时间为：',num2str(ag_time),'s']);
end

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
global dbf_mode;
val=get(hObject,'Value');
switch val
case 1
    dbf_mode = 'normal';
case 2
    dbf_mode = 'capon';
end


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



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over popupmenu1.
function popupmenu1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


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
