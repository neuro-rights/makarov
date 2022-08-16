function varargout = struct2d(varargin)
%   SYNTAX
%   struct2d
%   DESCRIPTION
%   This GUI script creates planar meshes using PDE toolbox. It allows the user
%   to define the planar geometry of the transmission line to be studied.
%   The planar geometry consists of the projection in the x-y plane of all
%   borders between different domains.  To define the geometry, the user
%   may delimit up to a total of ten rectangles and ellipses. After
%   delimiting the shapes to be included in the geometry, the user presses
%   the "View mesh" button to view a triangular mesh.  If the mesh is
%   satisfactory, the user presses the button labeled "Accept mesh" to save
%   the data and then the "Close" button to exit. 
%   When the user presses the "Accept mesh" button, STRUCT2D saves the
%   planar mesh to the file struct2d.mat.  This file contains two matrices:
%
%   P - 3-by-M array of nodes; P(1,I) is the x-coordinate of the Ith node,
%       P(2,I) is the y-coordinate, and P(3,I) is zero and represents the
%       z-coordinate of the node.
%   t - 4-by-N array of triangles; t(1:3,J) are the three indices in P of
%       the vertices of the Jth triangle.  In addition, t(4,J) is the
%       subdomain number of the Jth triangle.
%   
%   Authors: A. Marut, S.N. Makarov
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @struct2d_OpeningFcn, ...
                   'gui_OutputFcn',  @struct2d_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before struct2d is made visible.
function struct2d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to struct2d (see VARARGIN)

global STRUCT2D_BASE_DIRECTORY;
STRUCT2D_BASE_DIRECTORY = fileparts(mfilename('fullpath'));
addpath(STRUCT2D_BASE_DIRECTORY);
addpath(fullfile(STRUCT2D_BASE_DIRECTORY, 'codes'));

% Choose default command line output for struct2d
handles.output = hObject;

if ispc
    editbox_bgcolor = 'white';
else
    editbox_bgcolor = get(0,'defaultUicontrolBackgroundColor');
end

handles.        NAME = 1;
handles.          XC = 2;
handles.          YC = 3;
handles.      SQUARE = 4;
handles.      CIRCLE = 5;
handles.       WIDTH = 6;
handles.      HEIGHT = 7;
handles.      DIEL = 8;
handles.DIEL_R = 9;

handles.NUM_RECTS = 10;
handles.NUM_POLY_VERTICES = 8;

HORIZONTAL_ALIGNMENT = 'center';

handles.rect_inputs = zeros(handles.NUM_RECTS, 9);
handles.poly_inputs = zeros(2, handles.NUM_POLY_VERTICES);

re_tag_prefixes = {'include', 'xc',   'yc',   'square',   'circle',      'width', 'height', 'diel',   'diel_r'};
re_styles =       {'edit',    'edit', 'edit', 'checkbox', 'radiobutton', 'edit',  'edit',   'checkbox', 'edit'};

for row = 1:handles.NUM_RECTS;
    for col = 1:length(re_tag_prefixes)
        handles.rect_inputs(row, col) = uicontrol( ...
            'Style',    re_styles{col}, ...
            'FontSize', 8.0, ...
            'Tag',      strcat(re_tag_prefixes{col}, num2str(row)), ...
            'String',   '', ...
            'Callback', '', ...
            'HorizontalAlignment', HORIZONTAL_ALIGNMENT);
        if strcmp(re_styles{col}, 'edit')
            set(handles.rect_inputs(row, col), 'BackgroundColor', editbox_bgcolor);
        end
    end
    set(handles.rect_inputs(row, handles.NAME), 'String', ['S'  num2str(row)], 'Enable', 'off');
    set(handles.rect_inputs(row, handles.SQUARE), 'Value', 1.0, 'Callback', strcat('struct2d(''shape_square_Callback'',gcbo,[],guidata(gcbo),', num2str(row), ')' ));
    set(handles.rect_inputs(row, handles.CIRCLE), 'Value', 0.0, 'Callback', strcat('struct2d(''shape_circle_Callback'',gcbo,[],guidata(gcbo),', num2str(row), ')' ));
end


for col=1:handles.NUM_POLY_VERTICES
    x = 63*col - 15;
    handles.poly_lbl_x(col)    = uicontrol( 'Style', 'text', 'FontSize', 8.0, 'String', 'Name', 'Callback', '');
    handles.poly_inputs(1,col) = uicontrol( 'Style', 'edit', 'FontSize', 8.0, 'Tag', ['x' num2str(col)], 'String', '', 'Callback', '', 'HorizontalAlignment', HORIZONTAL_ALIGNMENT, 'BackgroundColor', editbox_bgcolor);
    handles.poly_lbl_y(col)    = uicontrol( 'Style', 'text', 'FontSize', 8.0, 'String', 'Value(s)', 'Callback', '');
    handles.poly_inputs(2,col) = uicontrol( 'Style', 'edit', 'FontSize', 8.0, 'Tag', ['y' num2str(col)], 'String', '', 'Callback', '', 'HorizontalAlignment', HORIZONTAL_ALIGNMENT, 'BackgroundColor', editbox_bgcolor);
end

set(handles.edit_triangle_size, 'HorizontalAlignment', HORIZONTAL_ALIGNMENT);

if exist('dimensions.m','file')
    load_prev(handles);
end

% Update handles structure
guidata(hObject, handles);

set_position(gcf, [10 100 700 1], 'pixels'); % the height (4th element) gets set to the minimum value
p = get_position(gcf, 'normalized');
set_position(gcf, [(1-min([p(3:4); 1 1]))/2  p(3:4)], 'normalized');
figure1_ResizeFcn(gcf, [], handles);

% UIWAIT makes struct2d wait for user response (see UIRESUME)
uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = struct2d_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global STRUCT2D_BASE_DIRECTORY;
rmpath(STRUCT2D_BASE_DIRECTORY);
rmpath(fullfile(STRUCT2D_BASE_DIRECTORY, 'codes'));
clear MAT05_STRUCT2D_BASE_DIRECTORY;

% Get default command line output from handles structure
if isempty(handles)
    % user pressed X button
else
    varargout{1} = handles.output;
    close gcf;
end



% --- Executes on button press in a square.
function shape_square_Callback(hObject, eventdata, handles, source)
% hObject    handle to the square (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% source     an integer specifying in which row this event originated

square = handles.rect_inputs(source, handles.SQUARE);
circle = handles.rect_inputs(source, handles.CIRCLE);
set(square, 'Value', get(square, 'Max'));
set(circle, 'Value', get(circle, 'Min'));
guidata(hObject, handles);



% --- Executes on button press in a circle.
function shape_circle_Callback(hObject, eventdata, handles, source)
% hObject    handle to the circle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% source     an integer specifying in which row this event originated

circle = handles.rect_inputs(source, handles.CIRCLE);
square = handles.rect_inputs(source, handles.SQUARE);
set(circle, 'Value', get(circle, 'Max'));
set(square, 'Value', get(square, 'Min'));
guidata(hObject, handles);



% --- Executes on button press in cb_include_polygon.
function cb_include_polygon_Callback(hObject, eventdata, handles)
% hObject    handle to cb_include_polygon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This function does nothing.



% --- Executes during object creation, after setting all properties.
function edit_triangle_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_triangle_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit_triangle_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_triangle_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This function does nothing.



% --- Executes on button press in btn_view_mesh.
function btn_view_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to btn_view_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[p e t mid] = make_mesh(handles);
f = figure;
pdemesh(p{mid},e{mid},t{mid});
axis equal;



% --- Calculates the mesh from the data entered by the user
function [p,e,t,mid] = make_mesh(handles)
% First extract the parameterization information from the GUI
param_names  = {};
param_values = {};
param_dim_vector = [];
mid = []; % index of design variation for which parameters have the middle values
for col=1:8
    current_param_name   = get(handles.poly_inputs(1,col),'String');
    current_param_values = get(handles.poly_inputs(2,col),'String');
    if isempty(current_param_name) || isempty(current_param_values)
        continue
    end
    param_names{end+1}  = current_param_name;
    param_values{end+1} = str2num(current_param_values);
    param_dim_vector(end+1) = length(param_values{end});
    mid(end+1) = ceil(param_dim_vector(end)/2);
    temp = param_values{end};
    mid(end) = temp(mid(end));
end

% Recursively enumerate all design variations
param_matrix = zeros(0,1);
for i = 1:length(param_dim_vector)
    current_param_values = param_values{i};
    param_matrix = [repmat(param_matrix, 1, length(current_param_values)); reshape(repmat(current_param_values, size(param_matrix, 2), 1),1,[])];
    % Each column of param_matrix has one combination of parameter values.
end
setappdata(handles.figure1, 'Parameter_List',   param_names);
setappdata(handles.figure1, 'Parameter_Matrix', param_matrix);
if isempty(param_matrix)
    mid = 1;
else
    mid = find(all(param_matrix == repmat(mid', 1, size(param_matrix,2)),1));
end

p = cell(1,size(param_matrix, 2));
e = cell(1,size(param_matrix, 2));
t = cell(1,size(param_matrix, 2));

sf = get(handles.edit_set_formula, 'String');
if isempty(regexpi(sf, '\S')) % if sf contains only whitespace
    sf = generate_set_formula(handles);
    set(handles.edit_set_formula, 'String', sf);
end

included_shapes = zeros(1,handles.NUM_RECTS);
for shape_index = 1:handles.NUM_RECTS;
    included_shapes(shape_index) = ~isempty(regexpi(sf, ['S' num2str(shape_index) '(\D|$)']));
end

for param_matrix_col = 1:size(param_matrix,2)
    % Initialize the parameter values
    for param_matrix_row = 1:length(param_names)
        eval([param_names{param_matrix_row}  ' = ' num2str(param_matrix(param_matrix_row, param_matrix_col)) ';']);
    end
    gd = [];        % geometry description matrix (see DECSG)
    ns_cell = {};   % Cell array precursor of Name Space matrix
    POLYGON = 2;
    RECTANGLE = 3;
    ELLIPSE = 4;

    for rect=1:handles.NUM_RECTS
        inputs = handles.rect_inputs(rect,:);
        if included_shapes(rect)
            xc_vec = eval(get(inputs(handles.XC),     'String'));
            yc_vec = eval(get(inputs(handles.YC),     'String'));
            width  = eval(get(inputs(handles.WIDTH),  'String'));
            height = eval(get(inputs(handles.HEIGHT), 'String'));
            if get(inputs(handles.DIEL), 'Value') == 1.0
                diel_r = eval(get(inputs(handles.DIEL_R), 'String'));
            end
            shape_name_base = strcat('S', num2str(rect));
            for xi = 1:length(xc_vec)
                xc = xc_vec(xi);
                if length(xc_vec)>1
                    shape_name_x_suffix = strcat('_',num2str(xi));
                else
                    shape_name_x_suffix = '';
                end
                for yi = 1:length(yc_vec)
                    yc = yc_vec(yi);
                    if length(yc_vec)>1
                        shape_name_y_suffix = strcat('_',num2str(yi));
                    else
                        shape_name_y_suffix = '';
                    end
                    if get(inputs(handles.SQUARE), 'Value') == 1.0
                        % draw a rectangle
                        xmin = xc - width/2;
                        xmax = xc + width/2;
                        ymin = yc - height/2;
                        ymax = yc + height/2;
                        gd   = [gd [RECTANGLE; 4; xmin; xmax; xmax; xmin; ymax; ymax; ymin; ymin]];
                        ns_cell{end+1} = strcat(shape_name_base, shape_name_x_suffix, shape_name_y_suffix);                       
                    else
                        % draw an ellipse
                        gd   = [gd [ELLIPSE; xc; yc; width/2; height/2; 0; 0; 0; 0; 0]]; 
                        ns_cell{end+1} = strcat(shape_name_base, shape_name_x_suffix, shape_name_y_suffix);                       
                    end
                end
            end
        end
    end

    ns = char(ns_cell)';

    trisize = str2num(get(handles.edit_triangle_size, 'String'));
    if isempty(trisize)
        trisize = inf;
    end

    dl = decsg(gd, sf, ns);
    [curr_p curr_e curr_t] = initmesh(dl, 'Hmax', trisize, 'Hgrad', 1.75);  % do not change Hgrad
    curr_p = curr_p/1000; % mm
    p{param_matrix_col} = curr_p;
    e{param_matrix_col} = curr_e;
    t{param_matrix_col} = curr_t;
end



% --- Generates set formula corresponding to union of all defined shapes.
% Currently returns the empty string.
function sf = generate_set_formula(handles)
sf = '';
for shape_index = 1:handles.NUM_RECTS
    if         ~isempty(regexpi(get(handles.rect_inputs(shape_index, handles.XC    ), 'String'), '\S')) ...
            && ~isempty(regexpi(get(handles.rect_inputs(shape_index, handles.YC    ), 'String'), '\S')) ...
            && ~isempty(regexpi(get(handles.rect_inputs(shape_index, handles.WIDTH ), 'String'), '\S')) ...
            && ~isempty(regexpi(get(handles.rect_inputs(shape_index, handles.HEIGHT), 'String'), '\S'))
    % The above treats a field that contains only whitespace as empty.
        sf = [sf  'S'  num2str(shape_index)  ' + '];
    end
end
if ~isempty(sf)
    sf = sf(1:end-3); % Delete the trailing ' + '
end
    
% --- Executes on button press in btn_accept_mesh.
function btn_accept_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to btn_accept_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Save the user data
global MAT05_STRUCT2D_BASE_DIRECTORY;
dimensions_m = fopen(fullfile(MAT05_STRUCT2D_BASE_DIRECTORY,'dimensions.m'), 'wt');
fprintf(dimensions_m, '%%DIMENSIONS.M - Script used by STRUCT2D.M\n');
for i=1:handles.NUM_RECTS
    fprintf(dimensions_m, '  re_name{%d} = ''%s'';\n', i, get(handles.rect_inputs(i,handles.NAME), 'String'));
    fprintf(dimensions_m, '    re_xc{%d} = ''%s'';\n', i, get(handles.rect_inputs(i,handles.XC),   'String'));
    fprintf(dimensions_m, '    re_yc{%d} = ''%s'';\n', i, get(handles.rect_inputs(i,handles.YC),   'String'));
    if get(handles.rect_inputs(i,handles.SQUARE),'Value')
        fprintf(dimensions_m, ' re_shape{%d} = ''rectangle'';\n', i);
    else
        fprintf(dimensions_m, ' re_shape{%d} = ''ellipse'';\n', i);
    end
    fprintf(dimensions_m, ' re_width{%d} = ''%s'';\n', i, get(handles.rect_inputs(i,handles.WIDTH),'String'));
    fprintf(dimensions_m, 're_height{%d} = ''%s'';\n', i, get(handles.rect_inputs(i,handles.HEIGHT),'String'));
    fprintf(dimensions_m, 're_diel(%d) = %s;\n', i, mat2str(logical(get(handles.rect_inputs(i,handles.DIEL),'Value'))));
    fprintf(dimensions_m, 're_diel_r{%d} = ''%s'';\n\n', i, get(handles.rect_inputs(i,handles.DIEL_R),'String'));
end

sf = get(handles.edit_set_formula, 'String');
if isempty(regexpi(sf, '\S')) % if sf contains only whitespace
    sf = generate_set_formula(handles);
    set(handles.edit_set_formula, 'String', sf);
end
fprintf(dimensions_m, '\nset_formula = ''%s'';\n', sf);

fprintf(dimensions_m, '\n\ninclude_poly = %s;\n\n', mat2str(logical(get(handles.cb_include_polygon,'Value'))));
for i=1:handles.NUM_POLY_VERTICES
    fprintf(dimensions_m, 'poly_x{%d} = ''%s'';\n', i, get(handles.poly_inputs(1,i),'String'));
    fprintf(dimensions_m, 'poly_y{%d} = ''%s'';\n', i, get(handles.poly_inputs(2,i),'String'));
end
fprintf(dimensions_m, '\n\ntriangle_size = ''%s'';\n', get(handles.edit_triangle_size,'String'));
fclose(dimensions_m);

[p e tri mid] = make_mesh(handles);

param_names  = getappdata(handles.figure1, 'Parameter_List');
param_matrix = getappdata(handles.figure1, 'Parameter_Matrix');

variations_txt = fopen(fullfile(MAT05_STRUCT2D_BASE_DIRECTORY, 'variations.txt'), 'wt');
fprintf(variations_txt, 'Design variations and corresponding parameter values:\n');

for index=1:length(p)
    P             =   p{index};
    P(3,:)        =   0;
    t             =   tri{index};
    t(1:3,:)      =   sort(t(1:3,:),1);
    %--------------------------------------------------------------
    % Temporary solution
    N = max(t(4, :));
    t_ = t;
    leg1 = P(:,t(1,:))-P(:,t(2,:)); % vector from second vertex to first for every triangle
    leg2 = P(:,t(3,:))-P(:,t(2,:)); % vector from second vertex to third for every triangle
    areas = sqrt(sum(cross(leg1, leg2).^2))/2; % vector containing area of every triangle
    area_size = zeros(1,N);
    for n = 1:N
        area_size(n) = sum(areas(t(4, :)==n));
        % area_size(n) is the total area of the triangles with subdomain number n
    end
    [dummy, index1] = sort(area_size);
    % index1(n) is the nth-smallest subdomain
    for n = 1:N
        in = find(t(4, :)==index1(n));
        % indices in t of those triangles in the nth-smallest subdomain
        t_(4, in) = n;
    end   
    t = t_;
    t(1:3,:)      =   sort(t(1:3,:),1);
    %--------------------------------------------------------------
    save_data     =   struct('P', P, 't', t, 'e', e{index});
    mat_filename  = ['struct2d_var'  num2str(index)];
    fprintf(variations_txt, ['\n' mat_filename '\n']);
    
    for param_index = 1:length(param_names)
        save_data = setfield(save_data, param_names{param_index}, param_matrix(param_index, index));
        fprintf(variations_txt, [param_names{param_index} ' = ' num2str(param_matrix(param_index, index)) '\n']);
    end
    save(mat_filename, '-struct', 'save_data');
    if index == mid
        figure;
        viewer(P, t);
        view(0, 90);
    end
end

fclose(variations_txt);



% --- Loads the contents of the fields from the last time struct2d was run
function load_prev(handles)
clear dimensions;
dimensions;
for i=1:min([handles.NUM_RECTS  length(re_xc)])
    if exist('re_name', 'var')
        set(handles.rect_inputs(i,handles.NAME),'String',re_name{i});
    end
    set(handles.rect_inputs(i,handles.XC),'String',re_xc{i});
    set(handles.rect_inputs(i,handles.YC),'String',re_yc{i});
    set(handles.rect_inputs(i,handles.SQUARE),'Value',0);
    set(handles.rect_inputs(i,handles.CIRCLE),'Value',0);
    if strcmpi(re_shape{i},'rectangle')
        set(handles.rect_inputs(i,handles.SQUARE),'Value',1);
    else
        set(handles.rect_inputs(i,handles.CIRCLE),'Value',1);
    end
    set(handles.rect_inputs(i,handles.WIDTH),'String',re_width{i});
    set(handles.rect_inputs(i,handles.HEIGHT),'String',re_height{i});
    set(handles.rect_inputs(i,handles.DIEL),'Value',re_diel(i));
    set(handles.rect_inputs(i,handles.DIEL_R),'String',re_diel_r{i});
end
if ~exist('set_formula', 'var') % backwards compatibility
    set_formula = generate_default_set_formula(find(re_inc));
end
set(handles.edit_set_formula,'String',set_formula);
set(handles.cb_include_polygon,'Value',include_poly);
for i=1:handles.NUM_POLY_VERTICES
    set(handles.poly_inputs(1,i),'String',poly_x{i});
    set(handles.poly_inputs(2,i),'String',poly_y{i});
end
set(handles.edit_triangle_size,'String',triangle_size);



% --- Generates a set formula as the union of all included shapes.
function set_formula = generate_default_set_formula(included_shape_indices)
% included_shape_indices  list of the indices of rectangles and ellipses
%                         included in the design
set_formula = '';
for i=included_shape_indices;
    set_formula = strcat(set_formula, ' + S', num2str(i));
end
if ~isempty(set_formula)
    set_formula = set_formula(4:end);   % Get rid of trailing ' + '
end



% --- Executes on button press in btn_close.
function btn_close_Callback(hObject, eventdata, handles)
% hObject    handle to btn_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);



% --- Executes on resizing the figure window.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to callback object (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

window_width = get_width(gcf, 'pixels');
window_height = get_height(gcf, 'pixels');

if isempty(handles)
    return;
end

if window_width < 489
    window_width = 489;
    set_width(gcf, window_width, 'pixels');
end

if window_height < 27*handles.NUM_RECTS+231
    window_height = 27*handles.NUM_RECTS+231;
    set_height(gcf, window_height, 'pixels');
    p = get_position(gcf, 'normalized');
    if p(2)+p(4) > 1    % the top of the window is above the top of the screen
        set_y(gcf, 1-p(4), 'normalized');
        set_y(gcf, get_y(gcf, 'pixels')-16, 'pixels'); % make the title bar visible
    end
end

CHSP = 10;                      % Component (or Control) Horizontal SPacing
CVSP = (window_height-21*handles.NUM_RECTS-159)/(handles.NUM_RECTS+12);  % Component (or Control) Vertical   SPacing
WIDTH_NAME_LABEL = 35;
WIDTH_SHAPE_LABEL = 32;
SPACING_SHAPE_COLUMN = 1;
WIDTH_DIEL_LABEL = 34;
WIDTH_DIEL_R_LABEL = 66;  % Assumption: WIDTH_DIEL_R_LABEL <= WIDTH_TEXT_FIELD
WIDTH_CHECKBOX = 18;
WIDTH_RADIO_BUTTON = 18;
WIDTH_TEXT_FIELD = (window_width + 3 - (11*CHSP ...
    + max([WIDTH_NAME_LABEL  WIDTH_CHECKBOX])  ...
    + max([WIDTH_SHAPE_LABEL    WIDTH_CHECKBOX+SPACING_SHAPE_COLUMN+WIDTH_RADIO_BUTTON]) ...
    + max([WIDTH_DIEL_LABEL   WIDTH_CHECKBOX]) + 2))/5;

TFW = WIDTH_TEXT_FIELD;  % Text Field Width


set_position(handles.lbl_title,        [   0     11*CVSP+(21+CVSP)*handles.NUM_RECTS+140    window_width                      26                  ], 'pixels');
set_position(handles.lbl_units,        [   0     10*CVSP+(21+CVSP)*handles.NUM_RECTS+127    window_width                      18                  ], 'pixels');

% Rectangle/ellipse frame:

set_position(handles.frame1,           [ CHSP-1                 7*CVSP+119                 window_width-15   2*CVSP+(21+CVSP)*handles.NUM_RECTS+13], 'pixels');

% Contents of the frame:

set_position(handles.lbl_rect_ellipse, [                 2*CHSP-  2  9*CVSP+(21+CVSP)*handles.NUM_RECTS+123   108                       15], 'pixels');
set_position(handles.lbl_shape_name,   [                 2*CHSP      8*CVSP+(21+CVSP)*handles.NUM_RECTS+115   WIDTH_NAME_LABEL          15], 'pixels');
set_position(handles.lbl_xc,           [0.1*window_width+3*CHSP+  5  8*CVSP+(21+CVSP)*handles.NUM_RECTS+115   12                        15], 'pixels');
set_position(handles.lbl_yc,           [0.3*window_width+4*CHSP- 36  8*CVSP+(21+CVSP)*handles.NUM_RECTS+115   12                        15], 'pixels');
set_position(handles.lbl_shape,        [0.4*window_width+5*CHSP- 48  8*CVSP+(21+CVSP)*handles.NUM_RECTS+115   WIDTH_SHAPE_LABEL         15], 'pixels');
set_position(handles.lbl_width,        [0.5*window_width+6*CHSP- 49  8*CVSP+(21+CVSP)*handles.NUM_RECTS+115   28                        15], 'pixels');
set_position(handles.lbl_height,       [0.7*window_width+7*CHSP- 93  8*CVSP+(21+CVSP)*handles.NUM_RECTS+115   31                        15], 'pixels');
set_position(handles.lbl_diel,       [0.8*window_width+8*CHSP- 98  8*CVSP+(21+CVSP)*handles.NUM_RECTS+115   WIDTH_DIEL_LABEL        15], 'pixels');
set_position(handles.lbl_diel_r, [0.9*window_width+9*CHSP-122  8*CVSP+(21+CVSP)*handles.NUM_RECTS+115   WIDTH_DIEL_R_LABEL  15], 'pixels');

for row = 1:handles.NUM_RECTS;
    y = 8*CVSP + (21+CVSP)*(handles.NUM_RECTS - row) + 120;
    
    x = 2 * CHSP; % the x-coordinate of the label 'Name'
    x = x + 0;  % center the text field under the label
    
    set_position(handles.rect_inputs(row, handles.NAME), [x y  WIDTH_NAME_LABEL 21], 'pixels');
    
    x = x + WIDTH_NAME_LABEL + CHSP;

    set_position(handles.rect_inputs(row, handles.XC), [x y TFW 21], 'pixels');
    
    x = x + TFW + CHSP;    
        
    set_position(handles.rect_inputs(row, handles.YC), [x y TFW 21], 'pixels');
    
    x = x + TFW + CHSP;
    x = x + max([0  (WIDTH_SHAPE_LABEL - (WIDTH_CHECKBOX+SPACING_SHAPE_COLUMN+WIDTH_RADIO_BUTTON))/2]);    
        
    set_position(handles.rect_inputs(row, handles.SQUARE), [x y  WIDTH_CHECKBOX 21], 'pixels');
    
    x = x + WIDTH_CHECKBOX + SPACING_SHAPE_COLUMN;
    
    set_position(handles.rect_inputs(row, handles.CIRCLE), [x y  WIDTH_RADIO_BUTTON 21], 'pixels');
    
    x = x + WIDTH_RADIO_BUTTON + max([0  (WIDTH_SHAPE_LABEL - (WIDTH_CHECKBOX+SPACING_SHAPE_COLUMN+WIDTH_RADIO_BUTTON))/2]) + CHSP;
    
    set_position(handles.rect_inputs(row, handles.WIDTH), [x y TFW 21], 'pixels');
    
    x = x + TFW + CHSP;    
        
    set_position(handles.rect_inputs(row, handles.HEIGHT), [x y TFW 21], 'pixels');
    
    x = x + TFW + CHSP + max([0  (WIDTH_DIEL_LABEL-WIDTH_CHECKBOX)/2]);

    set_position(handles.rect_inputs(row, handles.DIEL), [x y  WIDTH_CHECKBOX 21], 'pixels');
    
    x = x + WIDTH_CHECKBOX + max([0  (WIDTH_DIEL_LABEL-WIDTH_CHECKBOX)/2]) + CHSP;

    set_position(handles.rect_inputs(row, handles.DIEL_R), [x y TFW 21], 'pixels');
    
end

set_position(handles.lbl_set_formula,  [  CHSP- 1   6*CVSP+101                       66   15], 'pixels');
set_position(handles.edit_set_formula, [2*CHSP+65   6*CVSP+ 98   window_width-3*CHSP-61   21], 'pixels');

% Polygon frame:

TFW = (window_width + 2 - (4 + handles.NUM_POLY_VERTICES) * CHSP - max([WIDTH_CHECKBOX  WIDTH_NAME_LABEL]))/handles.NUM_POLY_VERTICES;

set_position(handles.frame2,             [  CHSP-1  2*CVSP+25  window_width-2*CHSP+5  3*CVSP+69], 'pixels');
set_position(handles.lbl_polygons,       [2*CHSP-1  5*CVSP+85  57 15], 'pixels');
set_position(handles.lbl_include_poly,   [2*CHSP    4*CVSP+77  WIDTH_NAME_LABEL 15], 'pixels');
set_position(handles.cb_include_polygon, [2*CHSP+9  4*CVSP+61  WIDTH_CHECKBOX 15], 'pixels');

for col=1:handles.NUM_POLY_VERTICES
    x = (TFW + CHSP)*(col - 1) + 3 * CHSP + max([WIDTH_CHECKBOX  WIDTH_NAME_LABEL]);
    set_position(handles.poly_lbl_x(col),    [x  4*CVSP+77  TFW  15], 'pixels');
    set_position(handles.poly_inputs(1,col), [x  4*CVSP+58  TFW  21], 'pixels');
    set_position(handles.poly_lbl_y(col),    [x  3*CVSP+45  TFW  15], 'pixels');
    set_position(handles.poly_inputs(2,col), [x  3*CVSP+26  TFW  21], 'pixels');
end

% Things at the bottom of the window:
set_position(handles.lbl_tri_size,       [               CHSP-  4  CVSP    92                       22], 'pixels');
set_position(handles.edit_triangle_size, [             2*CHSP+ 86  CVSP-1  window_width-6*CHSP-358  26], 'pixels');
set_position(handles.btn_view_mesh,      [window_width-3*CHSP-272  CVSP-1  92                       26], 'pixels');
set_position(handles.btn_accept_mesh,    [window_width-2*CHSP-180  CVSP-1  92                       26], 'pixels');
set_position(handles.btn_close,          [window_width-  CHSP- 88  CVSP-1  92                       26], 'pixels');



% Returns the position [x y width height] of the UI object represented by a
% handle using the requested units of measurement
function p = get_position(h, units)
old_units = get(h, 'Units');
set(h, 'Units', units);
p = get(h, 'Position');
set(h, 'Units', old_units);



% Changes the position [x y width height] of the UI object represented by a
% handle using the requested units of measurement
function set_position(h, position, units);
old_units = get(h, 'Units');
set(h, 'Units', units);
set(h, 'Position', position);
set(h, 'Units', old_units);



% Returns the x-coordinate of the UI object represented by a
% handle using the requested units of measurement
function x = get_x(h, units)
p = get_position(h, units);
x = p(1);



% Changes the x-coordinate of the UI object represented by a
% handle using the requested units of measurement
function set_x(h, x, units);
p = get_position(h, units);
p(1) = x;
set_position(h, p, units);



% Returns the y-coordinate of the UI object represented by a
% handle using the requested units of measurement
function y = get_y(h, units)
p = get_position(h, units);
y = p(2);



% Changes the y-coordinate of the UI object represented by a
% handle using the requested units of measurement
function set_y(h, y, units);
p = get_position(h, units);
p(2) = y;
set_position(h, p, units);



% Returns the width of the UI object represented by a
% handle using the requested units of measurement
function width = get_width(h, units)
p = get_position(h, units);
width = p(3);



% Changes the width of the UI object represented by a
% handle using the requested units of measurement
function set_width(h, width, units);
p = get_position(h, units);
p(3) = width;
set_position(h, p, units);



% Returns the height of the UI object represented by a
% handle using the requested units of measurement
function height = get_height(h, units)
p = get_position(h, units);
height = p(4);



% Changes the x-coordinate of the UI object represented by a
% handle using the requested units of measurement
function set_height(h, height, units);
p = get_position(h, units);
p(4) = height;
set_position(h, p, units);



function edit_set_formula_Callback(hObject, eventdata, handles)
% hObject    handle to edit_set_formula (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_set_formula as text
%        str2double(get(hObject,'String')) returns contents of edit_set_formula as a double



% --- Executes during object creation, after setting all properties.
function edit_set_formula_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_set_formula (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object deletion, before destroying properties.
function lbl_width_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to lbl_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


