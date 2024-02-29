
function varargout = main_window(varargin)
% MAIN_WINDOW MATLAB code for main_window.fig
%      MAIN_WINDOW, by itself, creates a new MAIN_WINDOW or raises the existing
%      singleton*.
%
%      H = MAIN_WINDOW returns the handle to a new MAIN_WINDOW or the handle to
%      the existing singleton*.
%
%      MAIN_WINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_WINDOW.M with the given input arguments.
%
%      MAIN_WINDOW('Property','Value',...) creates a new MAIN_WINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_window_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_window_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES.


% Edit the above text to modify the response to help main_window

% Last Modified by GUIDE v2.5 05-May-2022 21:50:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_window_OpeningFcn, ...
                   'gui_OutputFcn',  @main_window_OutputFcn, ...
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


% --- Executes just before main_window is made visible.
function main_window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main_window (see VARARGIN)

% Choose default command line output for main_window
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main_window wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_window_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'Enable','off')
error_eps = eval(get(handles.error_field, 'String'));
energy = eval(get(handles.energy_field, 'String'));
sample_num = eval(get(handles.sample_field, 'String'));
k_value = eval(get(handles.k_value, 'String'));
global node_energy_constants;
    results = iot_main(node_energy_constants, length(node_energy_constants), energy, error_eps, sample_num, k_value);
if (k_value == 0)
    set(handles.onehop, 'string', mat2str(results(1, :)));
    set(handles.twohop, 'string', mat2str(results(2, :)));
else
    set(handles.khop, 'string', mat2str(results(1, :)));
end
set(hObject,'Enable','on')


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global node_energy_constants;
global BS;

global start_energy;
start_energy = 0;
global node_paths;
node_paths = [];
global node_coordinates;
node_count = eval(get(handles.point_num, 'String'));
BS = node_count;
node_coordinates = rand(node_count,2);

alpha = 2;
teta = 1;
sigma_power_z = 100;
D=squareform(pdist(node_coordinates));
node_energy_constants = -1 * (D.^alpha)*teta*sigma_power_z;
showNodesOnAxes(handles.node_place);


% --- Executes on button press in city_button.
function city_button_Callback(hObject, eventdata, handles)
% hObject    handle to city_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
columns = [0.2 0.3 0.6 0.8];
rows = [0.1 0.3 0.5 0.7 0.9 0.95];
global start_energy;
start_energy = 0;
global node_paths;
node_paths = [];
global node_coordinates;
node_coordinates = [];
for col=columns
    for row=rows
        node_coordinates(end + 1, :) = [col row];
    end
end

global node_energy_constants
alpha = 2;
teta = 1;
sigma_power_z = 100;
D=squareform(pdist(node_coordinates));
node_energy_constants = -1 * (D.^alpha)*teta*sigma_power_z;
for i=1:length(node_coordinates)
    for j=1:length(node_coordinates)
        node_a = node_coordinates(i, :);
        node_b = node_coordinates(j, :);
        if node_a(1) ~= node_b(1) && node_a(2) ~= node_b(2)
            node_energy_constants(i,j) = node_energy_constants(i,j) * 1000;
        end
    end
end
showNodesOnAxes(handles.node_place);


function showNodesOnAxes(figure_h)
global node_coordinates;
global BS;
global node_energy_vector;
global start_energy;
global node_paths;
%scatter_colors = zeros(length(node_coordinates), 3);
%scatter_colors(SOURCE) = 0;
%if start_energy ~= 0
    %brightness = ((start_energy - node_energy_vector') ./ start_energy) .* 0.8;
    %scatter_colors = brightness * ones(1, 3);
%end
%scatter_colors(BS, :) = [1 0 0];
%scatter_colors(1, :) = [0 0 1];
axes(figure_h);
cla reset;
%scatter(node_coordinates(:,1), node_coordinates(:,2), 10, scatter_colors, "filled")
iot_device = imread('iot_device.png');
iot_server = imread('iot_server.png');
hold on
for i = 1:length(node_coordinates)
    coordinate = node_coordinates(i, :);
    if i == BS
        image([coordinate(1) + 0.02, coordinate(1) - 0.02], [coordinate(2) + 0.02, coordinate(2) - 0.02], iot_server);  
    else
        %brightness = 1 - ((start_energy - node_energy_vector(i)) ./ start_energy) .* 0.8;
        image([coordinate(1) + 0.02, coordinate(1) - 0.02], [coordinate(2) + 0.02, coordinate(2) - 0.02], iot_device);  
    end
end
for node_pair = node_paths'
    node_a = node_coordinates(node_pair(1), :);
    node_b = node_coordinates(node_pair(2), :);
    dp = node_b - node_a;
    quiver(node_a(1), node_a(2), dp(1), dp(2), 0);
end
hold off
xlim([0, 1]);
ylim([0, 1]);

function point_num_Callback(hObject, eventdata, handles)
% hObject    handle to point_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of point_num as text
%        str2double(get(hObject,'String')) returns contents of point_num as a double


% --- Executes during object creation, after setting all properties.
function point_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to point_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function energy_field_Callback(hObject, eventdata, handles)
% hObject    handle to energy_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of energy_field as text
%        str2double(get(hObject,'String')) returns contents of energy_field as a double


% --- Executes during object creation, after setting all properties.
function energy_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to energy_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sample_field_Callback(hObject, eventdata, handles)
% hObject    handle to sample_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sample_field as text
%        str2double(get(hObject,'String')) returns contents of sample_field as a double


% --- Executes during object creation, after setting all properties.
function sample_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sample_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function error_field_Callback(hObject, eventdata, handles)
% hObject    handle to error_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of error_field as text
%        str2double(get(hObject,'String')) returns contents of error_field as a double


% --- Executes during object creation, after setting all properties.
function error_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to error_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function compression_Callback(hObject, eventdata, handles)
% hObject    handle to compression (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of compression as text
%        str2double(get(hObject,'String')) returns contents of compression as a double


% --- Executes during object creation, after setting all properties.
function compression_CreateFcn(hObject, eventdata, handles)
% hObject    handle to compression (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_leach.
function run_leach_Callback(hObject, eventdata, handles)
% hObject    handle to run_leach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'Enable','off')
error_eps = eval(get(handles.error_field, 'String'));
energy = eval(get(handles.energy_field, 'String'));
sample_num = eval(get(handles.sample_field, 'String'));
p_node_election = eval(get(handles.node_election, 'String'));
compression = eval(get(handles.compression, 'String'));
global node_energy_constants;
    results = leach_main(node_energy_constants, length(node_energy_constants), energy, error_eps, sample_num, p_node_election, 1-compression);
    set(handles.leach, 'string', mat2str(results));
set(hObject,'Enable','on')


% --- Executes on button press in pushbutton_pegasis.
function pushbutton_pegasis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pegasis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'Enable','off')
error_eps = eval(get(handles.error_field, 'String'));
energy = eval(get(handles.energy_field, 'String'));
sample_num = eval(get(handles.sample_field, 'String'));
set(handles.text_pegasis, 'string', '');
global node_energy_constants;
    results = pegasis_main(node_energy_constants, length(node_energy_constants), energy, error_eps, sample_num);
    set(handles.text_pegasis, 'string', mat2str(results));
set(hObject,'Enable','on')



function node_election_Callback(hObject, eventdata, handles)
% hObject    handle to node_election (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of node_election as text
%        str2double(get(hObject,'String')) returns contents of node_election as a double


% --- Executes during object creation, after setting all properties.
function node_election_CreateFcn(hObject, eventdata, handles)
% hObject    handle to node_election (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k_value_Callback(hObject, eventdata, handles)
% hObject    handle to k_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k_value as text
%        str2double(get(hObject,'String')) returns contents of k_value as a double


% --- Executes during object creation, after setting all properties.
function k_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over run_leach.
function run_leach_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to run_leach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton6.
function step_leach_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global node_energy_vector;
global node_paths;
global BS;
global node_energy_constants;

error_eps = eval(get(handles.error_field, 'String'));
p_node_election = eval(get(handles.node_election, 'String'));
compression = eval(get(handles.compression, 'String'));
[paths, energy] = leach_step(node_energy_constants, BS, node_energy_vector, 1 - error_eps, p_node_election, compression);
node_paths = paths;
node_energy_vector = energy;
showNodesOnAxes(handles.node_place);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global node_energy_vector;
global node_energy_constants;
global start_energy;
energy = eval(get(handles.energy_field, 'String'));
start_energy = energy;
node_energy_vector = ones(1, length(node_energy_constants)) * energy;


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global node_energy_vector;
global node_paths;
global BS;
global node_energy_constants;

error_eps = eval(get(handles.error_field, 'String'));
k_value = eval(get(handles.k_value, 'String'));
[paths, energy] = k_hop_step(node_energy_constants, BS, node_energy_vector, 1 - error_eps, k_value);
node_paths = paths;
node_energy_vector = energy;
showNodesOnAxes(handles.node_place);
