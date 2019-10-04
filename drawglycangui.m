function varargout = drawglycangui(varargin)
% DRAWGLYCANGUI MATLAB code for drawglycangui.fig
%      DRAWGLYCANGUI, by itself, creates a new DRAWGLYCANGUI or raises the existing
%      singleton*.
%
%      H = DRAWGLYCANGUI returns the handle to a new DRAWGLYCANGUI or the handle to
%      the existing singleton*.
%
%      DRAWGLYCANGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DRAWGLYCANGUI.M with the given input arguments.
%
%      DRAWGLYCANGUI('Property','Value',...) creates a new DRAWGLYCANGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before drawglycangui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to drawglycangui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2016, Research Foundation for State University of New York. All rights reserved
%

% Edit the above text to modify the response to help drawglycangui

% Last Modified by GUIDE v2.5 06-Jun-2019 19:12:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @drawglycangui_OpeningFcn, ...
                   'gui_OutputFcn',  @drawglycangui_OutputFcn, ...
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


% --- Executes just before drawglycangui is made visible.
function drawglycangui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to drawglycangui (see VARARGIN)

% Choose default command line output for drawglycangui
handles.output = hObject;
defoptname = {'orientation','monosacsize','msperiwidth',...
    'bondwidth','perpendicularmonosac',...
    'linkinfodist','linkinfotheta','bondbreaksiglength',...
    'showlink','fontsize','workingmode',...
    'structspacing','aaspacing','fileout',...
    'visible','inputformat','displaybrackettext',...
    'specialoptions','linkfontsize','sortbranches',...
    'figurehandle'};

defoptvalue = {'left',.5,1,...
    2,{'Fuc','Xyl','Deoxyhexose','Pentose'},...
    .7,-30,.5,...
    'yes',12,'',...
    1,.75,'',...
    'on','iupac','no',...
    {},12,true,...
    []};
handles.optionstruct = usg('gen',defoptname,defoptvalue);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes drawglycangui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = drawglycangui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on key press with focus on iupacstring and none of its controls.
function iupacstring_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to iupacstring (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
key = get(gcf,'CurrentKey');
if(strcmp (key , 'return'))
    run_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    iupacstr = get(handles.iupacstring,'String');
    [gly,pep,~,~] = distggp(iupacstr,'IUPAC');
    if isempty(gly) && isempty(pep)
        errordlg('Looks like there is no glycan, peptide or glycopeptide here?')
        return
    end
    drawglycan(iupacstr,handles.optionstruct);
catch
    errordlg('DrawGlycan cannot recognize your input.')
end
%% temp


% --- Executes on selection change in poporientation.
function poporientation_Callback(hObject, eventdata, handles)
% hObject    handle to poporientation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns poporientation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from poporientation
contents = cellstr(get(hObject,'String'));
orientation = contents{get(hObject,'Value')};
handles.optionstruct.orientation = lower(orientation);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function poporientation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poporientation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popmonosacsize.
function popmonosacsize_Callback(hObject, eventdata, handles)
% hObject    handle to popmonosacsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popmonosacsize contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popmonosacsize
contents = cellstr(get(hObject,'String'));
monosacsize = contents{get(hObject,'Value')};
handles.optionstruct.monosacsize = str2num(monosacsize);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function popmonosacsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popmonosacsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showlinkage.
function showlinkage_Callback(hObject, eventdata, handles)
% hObject    handle to showlinkage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showlinkage
showlink = get(hObject,'Value');
switch showlink
    case 0
        handles.optionstruct.showlink = 'no';
    case 1
        handles.optionstruct.showlink = 'yes';
end
guidata(hObject,handles)

% --- Executes on selection change in poplinkfontsize.
function poplinkfontsize_Callback(hObject, eventdata, handles)
% hObject    handle to poplinkfontsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns poplinkfontsize contents as cell array
%        contents{get(hObject,'Value')} returns selected item from poplinkfontsize
contents = cellstr(get(hObject,'String'));
fontsize = contents{get(hObject,'Value')};
handles.optionstruct.linkfontsize = str2num(fontsize);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function poplinkfontsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poplinkfontsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in poplinktextspacing.
function poplinktextspacing_Callback(hObject, eventdata, handles)
% hObject    handle to poplinktextspacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns poplinktextspacing contents as cell array
%        contents{get(hObject,'Value')} returns selected item from poplinktextspacing
contents = cellstr(get(hObject,'String'));
linkinfodist = contents{get(hObject,'Value')};
handles.optionstruct.linkinfodist = str2num(linkinfodist);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function poplinktextspacing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poplinktextspacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iupacstring_Callback(hObject, eventdata, handles)
% hObject    handle to iupacstring (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iupacstring as text
%        str2double(get(hObject,'String')) returns contents of iupacstring as a double


% --- Executes during object creation, after setting all properties.
function iupacstring_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iupacstring (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function icon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to icon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate icon
if ispc
    [pathstr,~,~] = fileparts(which('drawglycangui.m'));
    imshow([pathstr,'\DrawGlycanLogo.jpg'])
else
    imshow([pathstr,'/DrawGlycanLogo.jpg'])
end
