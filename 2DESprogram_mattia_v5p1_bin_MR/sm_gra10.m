function varargout = sm_gra10(varargin)
% SM_GRA10 MATLAB code for sm_gra10.fig
%      SM_GRA10, by itself, creates a new SM_GRA10 or raises the existing
%      singleton*.
%
%      H = SM_GRA10 returns the handle to a new SM_GRA10 or the handle to
%      the existing singleton*.
%
%      SM_GRA10('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SM_GRA10.M with the given input arguments.
%
%      SM_GRA10('Property','Value',...) creates a new SM_GRA10 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sm_gra10_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sm_gra10_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sm_gra10

% Last Modified by GUIDE v2.5 27-Oct-2017 11:35:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sm_gra10_OpeningFcn, ...
                   'gui_OutputFcn',  @sm_gra10_OutputFcn, ...
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


% --- Executes just before sm_gra10 is made visible.
function sm_gra10_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sm_gra10 (see VARARGIN)

% Choose default command line output for sm_gra10
handles.output = hObject;
% panel_handles1=findobj(handles.figure1,'type','tpanel');
% set(panel_handles1,'parent',handles.figure1);

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes sm_gra10 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sm_gra10_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function path_Callback(hObject, eventdata, handles)
% hObject    handle to path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global a0 files
% global pat
% pat=get(handles.path,'String')
acquisition2(hObject, eventdata, handles)

% Hints: get(hObject,'String') returns contents of path as text
%        str2double(get(hObject,'String')) returns contents of path as a double


% --- Executes during object creation, after setting all properties.
function path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pathmat_Callback(hObject, eventdata, handles)
% hObject    handle to pathmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data T2 PhInt wpr mStep
set(hObject,'BackgroundColor','red');
pat=get(handles.pathmat,'String');

load(pat);

set(hObject,'BackgroundColor','white')
textLabel = sprintf('%d', length(T2));
set(handles.NT2, 'String', textLabel);
textLabel = [sprintf('%d', T2(1)) '-' sprintf('%d', T2(end))];
set(handles.textRangeT2, 'String', textLabel);
%% create global variables%%
setappdata(0, 'Data' ,Data);
setappdata(0, 'PhInt' ,PhInt);%storing matrix in GUI variables
setappdata(0, 'NumScans' ,NumScans);%storing matrix in GUI variables
setappdata(0, 'T2' ,T2);
setappdata(0, 'mStep' ,mStep);%storing matrix in GUI variables
setappdata(0, 'wpr' ,wpr);%storing matrix in GUI variables





% Hints: get(hObject,'String') returns contents of pathmat as text
%        str2double(get(hObject,'String')) returns contents of pathmat as a double


% --- Executes during object creation, after setting all properties.
function pathmat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pathmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pathcali_Callback(hObject, eventdata, handles)
% hObject    handle to pathcali (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'BackgroundColor','red');
pat=get(handles.pathcali,'String');
calibration(hObject,handles,pat)
    
% Hints: get(hObject,'String') returns contents of pathcali as text
%        str2double(get(hObject,'String')) returns contents of pathcali as a double


% --- Executes during object creation, after setting all properties.
function pathcali_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pathcali (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonCali.
% function pushbuttonCali_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbuttonCali (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% set(handles.uipanelCalibration,'visible','on')
% set(handles.text49,'visible','on')
% set(handles.uipanelMaps,'visible','off')
% set(handles.tpanel,'visible','off')


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function editBin0_Callback(hObject, eventdata, handles)
% hObject    handle to editBin0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% FT_mag_phase_unwrap.vi%%%%
Bin0function(handles)



% Hints: get(hObject,'String') returns contents of editBin0 as text
%        str2double(get(hObject,'String')) returns contents of editBin0 as a double


% --- Executes during object creation, after setting all properties.
function editBin0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBin0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject, 'Units', 'Normalized');
set(hObject, 'Position', [0 1 1 1]);



function editPPfile_Callback(hObject, eventdata, handles)
% hObject    handle to editPPfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
converter=1;
get(handles.checkboxNM,'Value')
if (get(handles.checkboxNM,'Value') == get(handles.checkboxNM,'Max'))
elseif (get(handles.checkboxEV,'Value') == get(handles.checkboxEV,'Max'))
    %converter=1/1239.84;
    converter=1129.125;
elseif (get(handles.checkboxCM,'Value') == get(handles.checkboxCM,'Max'))
    converter=1/10000000;
end
set(hObject,'BackgroundColor','red');
pat=get(handles.editPPfile,'String');
PumpProbe(hObject,handles,pat,converter)

% Hints: get(hObject,'String') returns contents of editPPfile as text
%        str2double(get(hObject,'String')) returns contents of editPPfile as a double


% --- Executes during object creation, after setting all properties.
function editPPfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPPfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonMaps.
% function pushbuttonMaps_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbuttonMaps (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% 
% set(handles.uipanelMaps,'visible','on')
% set(handles.tpanel,'visible','off')
% set(handles.uipanelCalibration,'visible','off')


function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editT2_Callback(hObject, eventdata, handles)
% hObject    handle to editT2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axesPP)
cla reset;
converter=1;
if (get(handles.checkboxNM,'Value') == get(handles.checkboxNM,'Max'))
elseif (get(handles.checkboxEV,'Value') == get(handles.checkboxEV,'Max'))
    converter=1239.84;
elseif (get(handles.checkboxCM,'Value') == get(handles.checkboxCM,'Max'))
    converter=10000000;
end
plotMapT2_bis(hObject,handles,converter)
set(hObject,'BackgroundColor','red');
pat=get(handles.editPPfile,'String');
PumpProbe(hObject,handles,pat,converter)
axes(handles.axesPP)
legend('2Dmap','Pump-Probe','Averaged Pump-Probe','Location','best')
%hold off
% Hints: get(hObject,'String') returns contents of editT2 as text
%        str2double(get(hObject,'String')) returns contents of editT2 as a double


% --- Executes during object creation, after setting all properties.
function editT2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editT2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMinProbe_Callback(hObject, eventdata, handles)
% hObject    handle to editMinProbe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinProbe as text
%        str2double(get(hObject,'String')) returns contents of editMinProbe as a double


% --- Executes during object creation, after setting all properties.
function editMinProbe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinProbe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMaxProbe_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxProbe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxProbe as text
%        str2double(get(hObject,'String')) returns contents of editMaxProbe as a double


% --- Executes during object creation, after setting all properties.
function editMaxProbe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxProbe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNumAverProbe_Callback(hObject, eventdata, handles)
% hObject    handle to editNumAverProbe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNumAverProbe as text
%        str2double(get(hObject,'String')) returns contents of editNumAverProbe as a double


% --- Executes during object creation, after setting all properties.
function editNumAverProbe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNumAverProbe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editZeroPadding_Callback(hObject, eventdata, handles)
% hObject    handle to editZeroPadding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editZeroPadding as text
%        str2double(get(hObject,'String')) returns contents of editZeroPadding as a double


% --- Executes during object creation, after setting all properties.
function editZeroPadding_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editZeroPadding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ediMinPump_Callback(hObject, eventdata, handles)
% hObject    handle to ediMinPump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ediMinPump as text
%        str2double(get(hObject,'String')) returns contents of ediMinPump as a double


% --- Executes during object creation, after setting all properties.
function ediMinPump_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ediMinPump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMaxPump_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxPump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxPump as text
%        str2double(get(hObject,'String')) returns contents of editMaxPump as a double


% --- Executes during object creation, after setting all properties.
function editMaxPump_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxPump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSpikeThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editSpikeThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSpikeThreshold as text
%        str2double(get(hObject,'String')) returns contents of editSpikeThreshold as a double


% --- Executes during object creation, after setting all properties.
function editSpikeThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSpikeThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editShiftTime_Callback(hObject, eventdata, handles)
% hObject    handle to editShiftTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editShiftTime as text
%        str2double(get(hObject,'String')) returns contents of editShiftTime as a double


% --- Executes during object creation, after setting all properties.
function editShiftTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editShiftTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonCursor.
function pushbuttonCursor_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCursor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2Dmap)
[Xcursor,Ycursor]=ginput(1);
textLabel = sprintf('%.3f', Xcursor);
set(handles.editXcursor, 'String', textLabel);
textLabel = sprintf('%.3f', Ycursor);
set(handles.editYcursor, 'String', textLabel);
Data_in=getappdata(0, 'DataAverageFree');
wprCutaver=getappdata(0, 'wprCutaver');
if (get(handles.checkboxNM,'Value') == get(handles.checkboxNM,'Max'))
IndexXcursor=find(wprCutaver>=Xcursor,1);
    else
        IndexXcursor=find(wprCutaver<=Xcursor,1);
end
Data_inFilter=getappdata(0, 'DataAverageFreeWindowed');   
axes(handles.axesXcursor)
if (get(handles.checkboxFS,'Value') == get(handles.checkboxFS,'Max'))
    x=linspace(0,length(Data_in(:,IndexXcursor)),length(Data_in(:,IndexXcursor)))*250/length(Data_in(:,IndexXcursor));
    plot(x,Data_in(:,IndexXcursor)*100,'LineWidth',1.5)
hold on
axes(handles.axesXcursor)
plot(x,Data_inFilter(:,IndexXcursor)*100,'r','LineWidth',1.5)
xlabel('t_1 [fs]');
ylabel('\Delta T/T [%]');
else
plot(Data_in(:,IndexXcursor)*100,'LineWidth',1.5)
hold on
axes(handles.axesXcursor)
plot(Data_inFilter(:,IndexXcursor)*100,'r','LineWidth',1.5)
xlabel('t_1 [Motor Steps]');
ylabel('\Delta T/T [%]');
end
legend('Signal average free','Signal average free and windowed')
set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
hold off



function editXcursor_Callback(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit32 as text
%        str2double(get(hObject,'String')) returns contents of edit32 as a double


% --- Executes during object creation, after setting all properties.
function editXcursor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editYcursor_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double


% --- Executes during object creation, after setting all properties.
function editYcursor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFWHMfilter_Callback(hObject, eventdata, handles)
% hObject    handle to editFWHMfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFWHMfilter as text
%        str2double(get(hObject,'String')) returns contents of editFWHMfilter as a double


% --- Executes during object creation, after setting all properties.
function editFWHMfilter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFWHMfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editThreshold as text
%        str2double(get(hObject,'String')) returns contents of editThreshold as a double


% --- Executes during object creation, after setting all properties.
function editThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit0Pad_Callback(hObject, eventdata, handles)
% hObject    handle to edit0Pad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit0Pad as text
%        str2double(get(hObject,'String')) returns contents of edit0Pad as a double


% --- Executes during object creation, after setting all properties.
function edit0Pad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit0Pad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNAverCal_Callback(hObject, eventdata, handles)
% hObject    handle to editNAverCal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNAverCal as text
%        str2double(get(hObject,'String')) returns contents of editNAverCal as a double


% --- Executes during object creation, after setting all properties.
function editNAverCal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNAverCal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editGaussCal_Callback(hObject, eventdata, handles)
% hObject    handle to editGaussCal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editGaussCal as text
%        str2double(get(hObject,'String')) returns contents of editGaussCal as a double


% --- Executes during object creation, after setting all properties.
function editGaussCal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editGaussCal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
converter=1;
if (get(handles.checkboxNM,'Value') == get(handles.checkboxNM,'Max'))
elseif (get(handles.checkboxEV,'Value') == get(handles.checkboxEV,'Max'))
    converter=1239.84;
elseif (get(handles.checkboxCM,'Value') == get(handles.checkboxCM,'Max'))
    converter=10000000;
end
SaveMap2D(hObject,handles,converter)


% --- Executes on button press in checkboxNM.
function checkboxNM_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxNM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxNM


% --- Executes on button press in checkboxEV.
function checkboxEV_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxEV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxEV


% --- Executes on button press in checkboxCM.
function checkboxCM_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxCM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxCM


% --- Executes on button press in checkboxFFT.
function checkboxFFT_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxFFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxFFT


% --- Executes on button press in checkboxDFT.
function checkboxDFT_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxDFT


% --- Executes on button press in checkboxNBP.
function checkboxNBP_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxNBP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxNBP



function editOrderG_Callback(hObject, eventdata, handles)
% hObject    handle to editOrderG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editOrderG as text
%        str2double(get(hObject,'String')) returns contents of editOrderG as a double


% --- Executes during object creation, after setting all properties.
function editOrderG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editOrderG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxSmoothing.
function checkboxSmoothing_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxSmoothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxSmoothing



function editNsmooth_Callback(hObject, eventdata, handles)
% hObject    handle to editNsmooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNsmooth as text
%        str2double(get(hObject,'String')) returns contents of editNsmooth as a double


% --- Executes during object creation, after setting all properties.
function editNsmooth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNsmooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNSmooth_Callback(hObject, eventdata, handles)
% hObject    handle to editNSmooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNSmooth as text
%        str2double(get(hObject,'String')) returns contents of editNSmooth as a double


% --- Executes during object creation, after setting all properties.
function editNSmooth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNSmooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxFS.
function checkboxFS_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxFS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxFS



function editCMmax_Callback(hObject, eventdata, handles)
% hObject    handle to editCMmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2Dmap)
cmin=str2double(get(handles.editCMmin,'String'));
cmax=str2double(get(handles.editCMmax,'String'));
caxis([cmin cmax])

% Hints: get(hObject,'String') returns contents of editCMmax as text
%        str2double(get(hObject,'String')) returns contents of editCMmax as a double


% --- Executes during object creation, after setting all properties.
function editCMmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCMmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCMmin_Callback(hObject, eventdata, handles)
% hObject    handle to editCMmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2Dmap)
cmin=str2double(get(handles.editCMmin,'String'));
cmax=str2double(get(handles.editCMmax,'String'));
caxis([cmin cmax])

% Hints: get(hObject,'String') returns contents of editCMmin as text
%        str2double(get(hObject,'String')) returns contents of editCMmin as a double


% --- Executes during object creation, after setting all properties.
function editCMmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCMmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxContour.
function checkboxContour_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxContour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxContour


% --- Executes on button press in tbutton.
% function tbutton_Callback(hObject, eventdata, handles)
% % hObject    handle to tbutton (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% %[ Map2DSave,T2In,wlpump,wprCutaver ]=dataset( hObject, eventdata, handles )
% 
% set(handles.uipanelCalibration,'visible','off')
% set(handles.uipanelMaps,'visible','off')
% set(handles.text49,'visible','off')
% set(handles.tpanel,'visible','on')

function edit31_Callback(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit31 as text
%        str2double(get(hObject,'String')) returns contents of edit31 as a double
st=get(handles.text46,'String');
 if strcmp(st,'insert T2 value')
T2In=getappdata(0, 'T2');
DataIn=getappdata(0, 'Data');
NumScans=getappdata(0, 'NumScans');
wlpump=getappdata(0, 'wlpump');
mStep=getappdata(0, 'mStep');
IndexPump=getappdata(0, 'IndexPump');
Bin0=str2double(get(handles.editBin0,'String'));
PhIntIn=getappdata(0, 'PhInt');
if (get(handles.checkboxNM,'Value') == get(handles.checkboxNM,'Max'))
wpr=getappdata(0, 'wpr');
else
    wpr=converter./getappdata(0, 'wpr');
end
ZeroPadding=str2double(get(handles.editZeroPadding,'String'));
    InterceptCal=str2double(get(handles.textCalibInter,'String'));
    SlopeCal=str2double(get(handles.CalibS,'String'));
    UnwrapPhase=1; %%user button from handles
    InvertPhase=0; %%%User button from handles
    padsize=[2^(nextpow2(length(mStep))+ZeroPadding)-length(mStep)];
    FTsize=length(mStep)+padsize;
    XX=linspace(1,FTsize,FTsize);
    XX=XX/FTsize;
    XX=fliplr(XX);
    if (get(handles.checkboxNM,'Value') == get(handles.checkboxNM,'Max'))
wlpump=(300000./(SlopeCal.*XX+InterceptCal));
else
   wlpump=converter./(300000./(SlopeCal.*XX+InterceptCal));
end
    
ProbeAverageNsquare=str2double(get(handles.editNumAverProbe,'String')); %from user
MinProbe=str2double(get(handles.editMinProbe,'String')); 
MaxProbe=str2double(get(handles.editMaxProbe,'String'));

if MinProbe<MaxProbe
    IndexProbe=find(wpr>=MinProbe & wpr<=MaxProbe);
else
    IndexProbe=find(wpr<=MinProbe & wpr>=MaxProbe);
end
 
wprCut=wpr(IndexProbe);
wprCutaver=zeros(floor(length(wprCut)/2^ProbeAverageNsquare),1);
DevStdInput=str2double(get(handles.editFWHMfilter,'String')); 

MaxPump=str2double(get(handles.editMaxPump,'String'));
MinPump=str2double(get(handles.ediMinPump,'String'));

if MaxPump>MinPump
IndexPump=find(wlpump>MinPump & wlpump<MaxPump);
else
    IndexPump=find(wlpump<MinPump & wlpump>MaxPump);
end

wlpump=wlpump(IndexPump);
%%%%%%%%%%%%%%%%%%%%%%shift application
speed=300000000;
sh=str2double(get(handles.edit36,'String'));
sh=sh*1000000000000;
f1=(speed./wlpump).*1000000000;
f2=f1+sh;
newpump=(speed./f2)*1000000000;
wlpump=newpump;
%%%%%%%%%%%%%%%%%%%%%%
Map2DSave=cell(length(T2In),1);
PhIntTrunc=cell(length(T2In),1);

    
for k=1:length(T2In)
%DataOut=DataIn{k}(Bin0+3:end,:); %%%selecting part of the matrix corrisponding to t1>0
DataOut=DataIn{k}(Bin0:end,:);
PhInt=PhIntIn{k}; %%%selected interferogram
spike_thrsh=str2double(get(handles.editSpikeThreshold,'String'));

%% FFT Unwrap.vi%%
    
    %PhIntPad= padarray(PhInt,padsize,'post');
    PhIntPad= padpost2(PhInt,padsize);
    PhIntPadRot=circshift(PhIntPad,-Bin0); 
    FFTPhIntPadRot=fft(PhIntPadRot, FTsize);
     if UnwrapPhase==1
        phase=unwrap(angle(FFTPhIntPadRot));
    end
    phase=angle(FFTPhIntPadRot);
    %phase=fliplr(phase);
    phase=flipud(phase);
    amplFFT=abs(FFTPhIntPadRot);
    if InvertPhase==1
        phase=(-1)*phase;
    end
    expIphase=exp(1i*phase);

    PhIntPadTrunc=PhIntPad(Bin0:end);
    FFTPhIntPadTrunc=fft(PhIntPadTrunc,FTsize).*expIphase;
    
    %% Probe Filter %%
Data=DataOut;
DataCut=Data(:,IndexProbe); %Data cut according to probe wavelength

%% Call function average_N.vi
DataCutaver=zeros(size(DataCut,1),floor(length(wprCut)/2^ProbeAverageNsquare));
for i=1:floor(length(wprCut)/2^ProbeAverageNsquare)
    wprCutaver(i)=mean(wprCut(i*2^ProbeAverageNsquare-2^ProbeAverageNsquare+1:i*2^ProbeAverageNsquare));
    DataCutaver(:,i)=mean(DataCut(:,i*2^ProbeAverageNsquare-2^ProbeAverageNsquare+1:i*2^ProbeAverageNsquare),2);
end
prwlFilt=wprCutaver;

%% Gaussian Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nOrder=str2double(get(handles.editOrderG,'String'));
% DevStdInput_rel=DevStdInput*size(DataCutaver,1); %the program works with indexes so std deviation must be used realtive to the employed vector
%     x=linspace(1,size(DataCutaver,1)*2,size(DataCutaver,1)*2);
%     gaussianF=(exp(-(x.^2./(2*DevStdInput_rel^2)).^nOrder))';
%     UniGaussWindow=ones(size(DataCutaver,1)*2,1).*gaussianF;
%     UniGaussWindow=UniGaussWindow(length(UniGaussWindow)/2:length(UniGaussWindow)-1);
if get(handles.gau,'Value') == get(handles.gau,'Max')
DevStdInput=str2double(get(handles.editFWHMfilter,'String')); 
nOrder=str2double(get(handles.editOrderG,'String'));    
DevStdInput_rel=DevStdInput*size(DataCutaver,1); %the program works with indexes so std deviation must be used realtive to the employed vector
    x=linspace(1,size(DataCutaver,1)*2,size(DataCutaver,1)*2);
   gaussianF=(exp(-(x.^2./(2*DevStdInput_rel^2)).^nOrder))';
    UniGaussWindow=ones(size(DataCutaver,1)*2,1).*gaussianF;
    UniGaussWindow=UniGaussWindow(length(UniGaussWindow)/2:length(UniGaussWindow)-1);
end 
   if get(handles.bh,'Value') == get(handles.bh,'Max')
        x=linspace(1,size(DataCutaver,1)*2,size(DataCutaver,1)*2);
        UniGaussWindow=blackmanharris(length(x));
        UniGaussWindow=UniGaussWindow(length(UniGaussWindow)/2:length(UniGaussWindow)-1);
   end

   if get(handles.han,'Value') == get(handles.han,'Max')
       DevStdInput=str2double(get(handles.editFWHMfilter,'String'));
        x=linspace(1,size(DataCutaver,1)*2,size(DataCutaver,1)*2);
        UniGaussWindow=tukeywin(length(x));
        tfun=bartlett(round(length(x)*DevStdInput));
        if DevStdInput<=1
        tfun=padarray(tfun,abs(length(x)-length(tfun)),'post');
        else
                tfun=tfun(1:length(x));
        end
         UniGaussWindow=UniGaussWindow.*tfun;
        UniGaussWindow=UniGaussWindow(length(UniGaussWindow)/2:length(UniGaussWindow)-1);
   end
    
    
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Phase_Spectra=1; %from user
    if Phase_Spectra==1
        phase=exp(1i*phase);
    end
   
NoPhase=1; %from user
if NoPhase==1;
    phase=ones(length(phase),1);
end;

%% Ft loop

Data_in=DataCutaver;
Data_in=Data_in-repmat(mean(Data_in,1),size(Data_in,1),1);
Filter=UniGaussWindow;
for i=1:size(Data_in,2)
SpikeInd=find(abs(Data_in(:,i))>spike_thrsh);
Data_in(SpikeInd,i)=0;
Data_inFilter(:,i)=Data_in(:,i).*Filter; %%multiply each column for the filter

end


%Data_inFilterPad= padarray(Data_inFilter,padsize,'post');
Data_inFilterPad= padpost2(Data_inFilter,padsize);
%%
lp=size(Data_inFilterPad,2);
support=expIphase;
phase_padding=expIphase;
for w=1:lp-1
    phase_padding=cat(2,phase_padding,support);
end
%phase_padding=flipud(phase_padding);
phase_padding=-1.*phase_padding;
%%
FFTData_inFilterPad=fft(Data_inFilterPad,FTsize);%.*phase;
Data_average_Free_windowed=Data_in;
Map2D=FFTData_inFilterPad;
%% Pump Filter
DataCut2D=Map2D(IndexPump,:); %Data cut according to probe wavelength
DataCut2D=DataCut2D.*phase_padding(IndexPump,:);

%% Call function average_N.vi  %%%Nel programma non media sul pump

Pump_spectrum=amplFFT(IndexPump);
square_root=0;
NormalizedBySpectrum=0;
if square_root==1
    DataCut2D=sqrt(DataCut2D);
end
if (get(handles.checkboxNBP,'Value') == get(handles.checkboxNBP,'Max'))
    for i=1:size(DataCut2D,2)
     DataCut2D(:,i)=DataCut2D(:,i)./(Pump_spectrum);
    end
end
%%%Application of Filtering %%%
if get(handles.buttondis,'Value') == get(handles.buttondis,'Max')
%DataCut2D_filt=sgolayfilt(imag(DataCut2D)',3,7);
DataCut2D_filt=imag(DataCut2D);
end
if get(handles.areal,'Value') == get(handles.areal,'Max')
%DataCut2D_filt=sgolayfilt(real(DataCut2D)',3,7);
DataCut2D_filt=real(DataCut2D);
end
if get(handles.buttonmag,'Value') == get(handles.buttonmag,'Max')
%DataCut2D_filt=sgolayfilt(abs(DataCut2D)',3,7);
DataCut2D_filt=abs(DataCut2D);
end
%DataCut2D_filt=sgolayfilt(real(DataCut2D)',3,7);
%DataCut2Dfilt=sgolayfilt(DataCut2D_filt',3,7);
DataCut2Dfilt=DataCut2D_filt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%change units
if get(handles.ucm,'Value') == get(handles.ucm,'Max')
  
    DataCut2Dfilt=fliplr(DataCut2Dfilt);
    DataCut2Dfilt=flipud(DataCut2Dfilt);
end
if get(handles.invax,'Value') == get(handles.invax,'Max')
     DataCut2Dfilt= DataCut2Dfilt';
end


%%%%%%%%%%%%%%%%%%%%%%%
Map2DSave{k}=DataCut2Dfilt;
PhIntTrunc{k}=PhIntPadTrunc;


end
if get(handles.ucm,'Value') == get(handles.ucm,'Max')
  wlpump=10000000./wlpump;
    %wlpump=fliplr(wlpump);
    wprCutaver=10000000./wprCutaver;
    %wprCutaver=flipud(wprCutaver);
end
% 
 handles.map=Map2DSave;
 handles.t=T2In;
 handles.pump=wlpump;
 handles.probe=wprCutaver;
 set(handles.text46,'String','T2 value');
 time1=num2str(T2In(1));
 time2=num2str(T2In(length(T2In)));
 set(handles.edit34,'String',time1);
 set(handles.edit35,'String',time2);
 end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Map2DSave=handles.map;
  T2In=handles.t;
  wlpump=handles.pump;
  wprCutaver=handles.probe;
T2value=str2double(get(handles.edit31,'String')); 
indexT2=find(T2In>=T2value,1);
B=cell2mat(Map2DSave(indexT2));
wpump=linspace(wlpump(1),wlpump(length(wlpump)),length(wprCutaver));
%%%%%%%%%%%%%%%%%%qui ce il plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.axes12)

if get(handles.invax,'Value') == get(handles.invax,'Max')

if get(handles.ucm,'Value') == get(handles.ucm,'Max')
    B=flipud(B);
    B=fliplr(B);
pcolor(wlpump,wprCutaver,B)
colormap jet
shading interp
hold on
[C,hfigc] = contour(wlpump, wprCutaver, B);
set(hfigc, ...
    'LineWidth',1.0, ...
    'Color', [1 1 1]);
xlabel('\nu_{1} [cm-1]');
ylabel('\nu_{3} [cm-1]'); 

else
pcolor(wlpump,wprCutaver,B)
colormap jet
shading interp
hold on
[C,hfigc] = contour(wlpump, wprCutaver, B);
set(hfigc, ...
    'LineWidth',1.0, ...
    'Color', [1 1 1]);
xlabel('\lambda_{pump} [nm]');
ylabel('\lambda_{probe} [nm]'); 
end

else

if get(handles.ucm,'Value') == get(handles.ucm,'Max')
    B=flipud(B);
    B=fliplr(B);
pcolor(wprCutaver,wlpump,B)
colormap jet
shading interp
hold on
[C,hfigc] = contour(wprCutaver, wlpump, B);
set(hfigc, ...
    'LineWidth',1.0, ...
    'Color', [1 1 1]);
xlabel('\nu_{3} [cm-1]');
ylabel('\nu_{1} [cm-1]'); 

else
pcolor(wprCutaver,wlpump,B)
colormap jet
shading interp
hold on
[C,hfigc] = contour(wprCutaver, wlpump, B);
set(hfigc, ...
    'LineWidth',1.0, ...
    'Color', [1 1 1]);
xlabel('\lambda_{probe} [nm]');
ylabel('\lambda_{pump} [nm]'); 
end



end



% pcolor(wprCutaver,wlpump,B)
% colormap jet
% shading flat
% hold on
% 
%  [C,hfigc] = contour(wprCutaver, wlpump, B);
% 
% 
% set(hfigc, ...
%     'LineWidth',1.0, ...
%     'Color', [1 1 1]);
% 
% 
% xlabel('\lambda_{probe} [nm]');
% ylabel('\lambda_{pump} [nm]'); 
% 
 hold on 
 line(wpump,wpump,'Color','k','LineStyle','--','LineWidth',3);
 
 hold off

 
 guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
  handles.map=cell(100,1);
  handles.t=0;
  handles.pump=0;
  handles.probe=0;
  handles.control=0;



function edit32_Callback(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit32 as text
%        str2double(get(hObject,'String')) returns contents of edit32 as a double


% --- Executes during object creation, after setting all properties.
function edit32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double


% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
 %hObject    handle to pushbutton10 (see GCBO)
 %eventdata  reserved - to be defined in a future version of MATLAB
 %handles    structure with handles and user data (see GUIDATA)
 Map2DSave=handles.map;
 T2In=handles.t;
 wlpump=handles.pump;
 wprCutaver=handles.probe;
 if (get(handles.radiobutton3,'Value') == get(handles.radiobutton3,'Max'))
axes(handles.axes12)
[Xcursor,Ycursor]=ginput(1);
textLabel = sprintf('%.3f', Xcursor);
set(handles.edit32, 'String', textLabel);
textLabel = sprintf('%.3f', Ycursor);
set(handles.edit33, 'String', textLabel);
%Data_in=getappdata(0, 'DataAverageFree');
%wprCutaver=getappdata(0, 'wprCutaver');
% if (get(handles.checkboxNM,'Value') == get(handles.checkboxNM,'Max'))
% IndexXcursor=find(wprCutaver>=Xcursor,1);
%     else
%         IndexXcursor=find(wprCutaver<=Xcursor,1);
% end

IndexXcursor=find(wprCutaver>=Xcursor,1);
IndexYcursor=find(wlpump>=Ycursor,1);
somma=0;
for i=1:length(T2In)
    A=cell2mat(Map2DSave(i));
   B(i)=A(IndexYcursor,IndexXcursor);
%    for raw=IndexYcursor-1:IndexYcursor+1
%        for col=IndexXcursor-1:IndexXcursor+1
%            somma=somma+A(raw,col);
%        end
%    end
%    somma=somma/9;
%    B(i)=somma;
%    somma=0;
end
time1=str2double(get(handles.edit34,'String'));
time2=str2double(get(handles.edit35,'String'));
        
 
%Data_inFilter=getappdata(0, 'DataAverageFreeWindowed');   
axes(handles.axes13)
%if (get(handles.checkboxFS,'Value') == get(handles.checkboxFS,'Max'))
    %x=linspace(0,length(Data_in(:,IndexXcursor)),length(Data_in(:,IndexXcursor)))*250/length(Data_in(:,IndexXcursor));
    %plot(x,Data_in(:,IndexXcursor)*100,'LineWidth',1.5)
    plot(T2In,B,'LineWidth',1.5)
    xlabel('T2[fs]');
ylabel('\Delta T/T');
xlim([time1 time2]);
hold off
N=length(T2In);
t=str2double(get(handles.edit31,'String'));
x1=str2double(get(handles.edit32,'String'));
y1=str2double(get(handles.edit33,'String'));
index=find(wprCutaver>=x1);
indy=find(wlpump>=y1);
indm=find(T2In>=t);
spec1=cell2mat(Map2DSave(indm(1)));
spec2=spec1(:,index(1));
spec3=spec1(indy(1),:);
axes(handles.axes15)
plot(wlpump,spec2,'LineWidth',1.5);
xlabel('\lambda_{pump} [nm]');
ylabel('\Delta T/T');
hold off
axes(handles.axes17)
plot(wprCutaver,spec3,'LineWidth',1.5);
xlabel('\lambda_{probe} [nm]');
ylabel('\Delta T/T');
hold off
N=length(T2In);
x1=str2double(get(handles.edit32,'String'));
y1=str2double(get(handles.edit33,'String'));
index=find(wprCutaver>=x1);
indy=find(wlpump>=y1);

for m=1:N
A=cell2mat(Map2DSave(m));
if m==1
puev=A(:,index(1));
prev=A(indy(1),:);
else
    puev=cat(2,puev,A(:,index(1)));
    prev=cat(1,prev,A(indy(1),:));
end
end
prev=prev';
axes(handles.axes16)
pcolor(T2In,wlpump,puev);
shading interp
hold on

  [C,hfigc] = contour(T2In,wlpump,puev);
 set(hfigc, ...
     'LineWidth',1.0, ...
     'Color', [1 1 1]);


xlabel('T2 [fs]');
ylabel('\lambda_{pump} [nm]'); 
xlim([time1 time2]);
%   set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
%   set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
%   set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)

 %colorbar;
hold off

axes(handles.axes18)
pcolor(T2In,wprCutaver,prev);
shading interp
hold on

  [C,hfigc] = contour(T2In,wprCutaver,prev);
 set(hfigc, ...
     'LineWidth',1.0, ...
     'Color', [1 1 1]);


xlabel('T2 [fs]');
 ylabel('\lambda_{probe} [nm]'); 
 xlim([time1 time2]);
%  set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
%  set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
%  set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
 %colorbar;
hold off
N=length(T2In);
for i=1:N
    A=cell2mat(Map2DSave(i));
    if i==1
    B=sum(A);
    else
       B=cat(1,B,sum(A)); 
    end
end
b=B';

axes(handles.axes14)
pcolor(T2In,wprCutaver,b);
colormap jet
shading interp
hold on

  [C,hfigc] = contour(T2In,wprCutaver,b);
 set(hfigc, ...
     'LineWidth',1.0, ...
     'Color', [1 1 1]);


xlabel('T2 [fs]');
ylabel('\lambda_{probe} [nm]'); 
xlim([time1 time2]);
%  set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
%  set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
%  set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
%  colorbar('EastOutside')
hold off
 end
 
 %%%%%%%%%%%%%%%%%%%%%%%PEAKS
 if (get(handles.radiobutton4,'Value') == get(handles.radiobutton4,'Max'))
     axes(handles.axes12)
[Xpoint,Ypoint]=ginput(1);
hold off
dvett=0;
avett=0;
j=str2double(get(handles.edit31,'String'));
jj=find(T2In>=j);
A=cell2mat(Map2DSave(jj(1)));
dim=size(A);
indy1=find(wlpump>=Ypoint);
indx1=find(wprCutaver>=Xpoint);
m1=length(wlpump)/length(wprCutaver);
m2=-1/m1;
vettprobe=1:length(wprCutaver);
vettpump1=round(m1.*(vettprobe-indx1(1))+indy1(1));
vettpump2=round(m2.*(vettprobe-indx1(1))+indy1(1));
for i=1:length(wprCutaver)
if vettpump1(i)>dim(1) || vettpump1(i)<=0
   dvett(i)=0;
else
    dvett(i)=A(vettpump1(i),i);
end

if vettpump2(i)>dim(1) || vettpump2(i)<=0
   avett(i)=NaN;
else
    avett(i)=A(vettpump2(i),i);
end
end
axes(handles.axes19)
plot(wprCutaver,dvett,'red','LineWidth',2);
hold on
plot(wprCutaver,avett,'LineWidth',2);
hold on
xlabel('\lambda_{probe}[nm]');
ylabel('amplitude'); 
hold on
legend('Diagonal','Off Diagonal','Location','NorthEastOutside');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%surface


N=length(T2In);
B=cell2mat(Map2DSave(1));
b=size(B);

for j=1:N
    A=cell2mat(Map2DSave(j));

for i=1:length(wprCutaver)
if vettpump1(i)>dim(1) || vettpump1(i)<=0
   Tdiag(i,j)=0;
else
    Tdiag(i,j)=A(vettpump1(i),i);
end

if vettpump2(i)>dim(1) || vettpump2(i)<=0
   Toff(i,j)=0;
else
    Toff(i,j)=A(vettpump2(i),i);
end
end


end
%Tmap=flipud(Tmap);
axes(handles.axes21)
contourf(T2In,wprCutaver,Tdiag,8);
xlabel('T2 [fs]');
ylabel('\lambda_{probe} [nm]'); 
zlabel('amplitude');
hold off
axes(handles.axes22)
%Toff=flipud(Toff);
contourf(T2In,wprCutaver,Toff,8);
xlabel('T2 [fs]');
ylabel('\lambda_{probe} [nm]'); 
zlabel('amplitude');
hold off
 end
%function edit32_Callback(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit32 as text
%        str2double(get(hObject,'String')) returns contents of edit32 as a double


% --- Executes during object creation, after setting all properties.
% function edit32_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit32 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end



%function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double


% --- Executes during object creation, after setting all properties.
% function edit33_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit33 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end



function edit34_Callback(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit34 as text
%        str2double(get(hObject,'String')) returns contents of edit34 as a double
 Map2DSave=handles.map;
 T2In=handles.t;
 wlpump=handles.pump;
 wprCutaver=handles.probe;
axes(handles.axes12)
Xcursor=str2double(get(handles.edit32,'String'));
Ycursor=str2double(get(handles.edit33,'String'));
%Data_in=getappdata(0, 'DataAverageFree');
%wprCutaver=getappdata(0, 'wprCutaver');
% if (get(handles.checkboxNM,'Value') == get(handles.checkboxNM,'Max'))
% IndexXcursor=find(wprCutaver>=Xcursor,1);
%     else
%         IndexXcursor=find(wprCutaver<=Xcursor,1);
% end

IndexXcursor=find(wprCutaver>=Xcursor,1);
IndexYcursor=find(wlpump>=Ycursor,1);
somma=0;
for i=1:length(T2In)
    A=cell2mat(Map2DSave(i));
   %B(i)=A(IndexYcursor,IndexXcursor);
   for raw=IndexYcursor-1:IndexYcursor+1
       for col=IndexXcursor-1:IndexXcursor+1
           somma=somma+A(raw,col);
       end
   end
   somma=somma/9;
   B(i)=somma;
   somma=0;
end
time1=str2double(get(handles.edit34,'String'));
time2=str2double(get(handles.edit35,'String'));
        
 
%Data_inFilter=getappdata(0, 'DataAverageFreeWindowed');   
axes(handles.axes13)
%if (get(handles.checkboxFS,'Value') == get(handles.checkboxFS,'Max'))
    %x=linspace(0,length(Data_in(:,IndexXcursor)),length(Data_in(:,IndexXcursor)))*250/length(Data_in(:,IndexXcursor));
    %plot(x,Data_in(:,IndexXcursor)*100,'LineWidth',1.5)
    plot(T2In,B,'LineWidth',1.5)
    xlabel('T2[fs]');
ylabel('\Delta T/T');
xlim([time1 time2]);
hold off
N=length(T2In);
x1=str2double(get(handles.edit32,'String'));
y1=str2double(get(handles.edit33,'String'));
index=find(wprCutaver>=x1);
indy=find(wlpump>=y1);

for m=1:N
A=cell2mat(Map2DSave(m));
if m==1
puev=A(:,index(1));
prev=A(indy(1),:);
else
    puev=cat(2,puev,A(:,index(1)));
    prev=cat(1,prev,A(indy(1),:));
end
end
prev=prev';
axes(handles.axes16)
pcolor(T2In,wlpump,puev);
shading interp
hold on

  [C,hfigc] = contour(T2In,wlpump,puev);
 set(hfigc, ...
     'LineWidth',1.0, ...
     'Color', [1 1 1]);


xlabel('T2 [fs]');
ylabel('\lambda_{pump} [nm]'); 
xlim([time1 time2]);
%   set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
%   set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
%   set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)

 %colorbar;
hold off

axes(handles.axes18)
pcolor(T2In,wprCutaver,prev);
shading interp
hold on

  [C,hfigc] = contour(T2In,wprCutaver,prev);
 set(hfigc, ...
     'LineWidth',1.0, ...
     'Color', [1 1 1]);


xlabel('T2 [fs]');
 ylabel('\lambda_{probe} [nm]'); 
 xlim([time1 time2]);
%  set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
%  set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
%  set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
 %colorbar;
hold off
N=length(T2In);
for i=1:N
    A=cell2mat(Map2DSave(i));
    if i==1
    B=sum(A);
    else
       B=cat(1,B,sum(A)); 
    end
end
b=B';

axes(handles.axes14)
pcolor(T2In,wprCutaver,b);
colormap jet
shading interp
hold on

  [C,hfigc] = contour(T2In,wprCutaver,b);
 set(hfigc, ...
     'LineWidth',1.0, ...
     'Color', [1 1 1]);


xlabel('T2 [fs]');
ylabel('\lambda_{probe} [nm]'); 
xlim([time1 time2]);
%  set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
%  set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
%  set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
%  colorbar('EastOutside')
hold off
% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double
 Map2DSave=handles.map;
 T2In=handles.t;
 wlpump=handles.pump;
 wprCutaver=handles.probe;
axes(handles.axes12)
Xcursor=str2double(get(handles.edit32,'String'));
Ycursor=str2double(get(handles.edit33,'String'));
%Data_in=getappdata(0, 'DataAverageFree');
%wprCutaver=getappdata(0, 'wprCutaver');
% if (get(handles.checkboxNM,'Value') == get(handles.checkboxNM,'Max'))
% IndexXcursor=find(wprCutaver>=Xcursor,1);
%     else
%         IndexXcursor=find(wprCutaver<=Xcursor,1);
% end

IndexXcursor=find(wprCutaver>=Xcursor,1);
IndexYcursor=find(wlpump>=Ycursor,1);
somma=0;
for i=1:length(T2In)
    A=cell2mat(Map2DSave(i));
   %B(i)=A(IndexYcursor,IndexXcursor);
   for raw=IndexYcursor-1:IndexYcursor+1
       for col=IndexXcursor-1:IndexXcursor+1
           somma=somma+A(raw,col);
       end
   end
   somma=somma/9;
   B(i)=somma;
   somma=0;
end
time1=str2double(get(handles.edit34,'String'));
time2=str2double(get(handles.edit35,'String'));
        
 
%Data_inFilter=getappdata(0, 'DataAverageFreeWindowed');   
axes(handles.axes13)
%if (get(handles.checkboxFS,'Value') == get(handles.checkboxFS,'Max'))
    %x=linspace(0,length(Data_in(:,IndexXcursor)),length(Data_in(:,IndexXcursor)))*250/length(Data_in(:,IndexXcursor));
    %plot(x,Data_in(:,IndexXcursor)*100,'LineWidth',1.5)
    plot(T2In,B,'LineWidth',1.5)
    xlabel('T2[fs]');
ylabel('\Delta T/T');
xlim([time1 time2]);
hold off
N=length(T2In);
x1=str2double(get(handles.edit32,'String'));
y1=str2double(get(handles.edit33,'String'));
index=find(wprCutaver>=x1);
indy=find(wlpump>=y1);

for m=1:N
A=cell2mat(Map2DSave(m));
if m==1
puev=A(:,index(1));
prev=A(indy(1),:);
else
    puev=cat(2,puev,A(:,index(1)));
    prev=cat(1,prev,A(indy(1),:));
end
end
prev=prev';
axes(handles.axes16)
pcolor(T2In,wlpump,puev);
shading interp
hold on

  [C,hfigc] = contour(T2In,wlpump,puev);
 set(hfigc, ...
     'LineWidth',1.0, ...
     'Color', [1 1 1]);


xlabel('T2 [fs]');
ylabel('\lambda_{pump} [nm]'); 
xlim([time1 time2]);
%   set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
%   set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
%   set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)

 %colorbar;
hold off

axes(handles.axes18)
pcolor(T2In,wprCutaver,prev);
shading interp
hold on

  [C,hfigc] = contour(T2In,wprCutaver,prev);
 set(hfigc, ...
     'LineWidth',1.0, ...
     'Color', [1 1 1]);


xlabel('T2 [fs]');
 ylabel('\lambda_{probe} [nm]'); 
 xlim([time1 time2]);
%  set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
%  set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
%  set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
 %colorbar;
hold off
N=length(T2In);
for i=1:N
    A=cell2mat(Map2DSave(i));
    if i==1
    B=sum(A);
    else
       B=cat(1,B,sum(A)); 
    end
end
b=B';

axes(handles.axes14)
pcolor(T2In,wprCutaver,b);
colormap jet
shading interp
hold on

  [C,hfigc] = contour(T2In,wprCutaver,b);
 set(hfigc, ...
     'LineWidth',1.0, ...
     'Color', [1 1 1]);


xlabel('T2 [fs]');
ylabel('\lambda_{probe} [nm]'); 
xlim([time1 time2]);
%  set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
%  set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
%  set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
%  colorbar('EastOutside')
hold off
% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
film(hObject, eventdata, handles)

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text58,'String','stop');

% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text58,'String','continue');
tempo2=get(handles.text57,'String');
[token,remain]=strtok(tempo2,'f');
j=str2double(token);
filmc( hObject, eventdata, handles,j )


% --- Executes when tpanel is resized.
function tpanel_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to tpanel (see GCBO)
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

s=get(hObject,'String');
val=get(hObject,'Value');

switch s{val}
    case 'CALIBRATION'
        set(handles.uipanelCalibration,'visible','on')
	set(handles.text49,'visible','on')
	set(handles.uipanelMaps,'visible','off')
	set(handles.tpanel,'visible','off')
    case 'MAPS'
       set(handles.uipanelMaps,'visible','on')
	set(handles.tpanel,'visible','off')
	set(handles.uipanelCalibration,'visible','off')
    case 'T2 ANALYSIS'
        set(handles.uipanelCalibration,'visible','off')
	set(handles.uipanelMaps,'visible','off')
	set(handles.text49,'visible','off')
	set(handles.tpanel,'visible','on')
    set(handles.text46,'String','insert T2 value');
    set(handles.pushbutton10,'String','Cursor');
  end
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in uv.
function uv_Callback(hObject, eventdata, handles)
% hObject    handle to uv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of uv
set(handles.vis,'Value',0);

% --- Executes on button press in vis.
function vis_Callback(hObject, eventdata, handles)
% hObject    handle to vis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uv,'Value',0);
% Hint: get(hObject,'Value') returns toggle state of vis



function edit36_Callback(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit36 as text
%        str2double(get(hObject,'String')) returns contents of edit36 as a double


% --- Executes during object creation, after setting all properties.
function edit36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in areal.
function areal_Callback(hObject, eventdata, handles)
% hObject    handle to areal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of areal
set(handles.buttonmag,'Value',0);
set(handles.buttondis,'Value',0);

% --- Executes on button press in buttonmag.
function buttonmag_Callback(hObject, eventdata, handles)
% hObject    handle to buttonmag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of buttonmag
set(handles.areal,'Value',0);
set(handles.buttondis,'Value',0);

% --- Executes on button press in buttondis.
function buttondis_Callback(hObject, eventdata, handles)
% hObject    handle to buttondis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of buttondis
set(handles.buttonmag,'Value',0);
set(handles.areal,'Value',0);


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
set(handles.radiobutton4,'Value',0);
set(handles.diagpanel,'Visible','off')

% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.radiobutton3,'Value',0);
set(handles.diagpanel,'Visible','on')
% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in gau.
function gau_Callback(hObject, eventdata, handles)
% hObject    handle to gau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gau
set(handles.bh,'Value',0);
set(handles.han,'Value',0);

% --- Executes on button press in bh.
function bh_Callback(hObject, eventdata, handles)
% hObject    handle to bh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.gau,'Value',0);
set(handles.han,'Value',0);
% Hint: get(hObject,'Value') returns toggle state of bh


% --- Executes on button press in han.
function han_Callback(hObject, eventdata, handles)
% hObject    handle to han (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of han
set(handles.gau,'Value',0);
set(handles.bh,'Value',0);


% --- Executes on button press in ucm.
function ucm_Callback(hObject, eventdata, handles)
% hObject    handle to ucm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ucm


% --- Executes on button press in invax.
function invax_Callback(hObject, eventdata, handles)
% hObject    handle to invax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of invax


% --- Executes on button press in sm.
function sm_Callback(hObject, eventdata, handles)
% hObject    handle to sm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sm
