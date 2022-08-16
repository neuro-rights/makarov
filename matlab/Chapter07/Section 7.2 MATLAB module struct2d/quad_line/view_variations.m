function view_variations
%   SYNTAX
%   view_variations
%   DESCRIPTION
%   This function displays the variations of a design. Call VIEW_VARIATIONS
%   with no arguments.  Use the slider to iterate through the variations.
%   
%   Authors: A. Marut, S.N. Makarov
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

dir_contents = dir;
num_variations = 0;
for i=1:length(dir_contents)
    if regexpi(dir_contents(i).name, '^struct2d_var')
        num_variations = num_variations + 1;
    end
end
fig = figure;
set(fig, 'Color', [0.95 0.95 0.95], 'Toolbar', 'figure', 'Name', '2D geometry sweep');
ax = axes;
set(ax, 'Units', 'normalized', 'OuterPosition', [0 0.05 1 0.95]);

if num_variations > 1
    slider = uicontrol('Style', 'slider', ...
        'Units', 'normalized', 'Position', [0 0 1 0.05], ...
        'BackgroundColor', [192 203 225]/255, ...
        'Min', 1, 'Max', num_variations, ...
        'Value', 1, ...
        'Callback', @slider_callback, ...
        'SliderStep', [1/(num_variations-1)  1/(sqrt(num_variations-1))]);
    slider_callback(slider, []);
elseif num_variations == 1
    load struct2d_var1;
    cd codes
    viewer(P,t);
    cd ..
    view(0,90);
    xlabel('x [m]');
    ylabel('y [m]');
    grid on;
end

function slider_callback(hObject, eventdata)
val = get(hObject, 'Value');
val = round(val);
set(hObject, 'Value', val);
matfilename = ['struct2d_var'  num2str(val)];
param_names = who('-file', matfilename);
title_string = '';
load(matfilename);
for i = 1:length(param_names);
    name = param_names{i};
    if strcmp(name, 'P') || strcmp(name, 't')
        continue;
    end
    if ~isempty(title_string)
        title_string = [title_string ', '];
    end
    title_string = [title_string  param_names{i} ' = ' num2str(eval(param_names{i}))];
end
cd codes
viewer(P,t);
cd ..
view(0,90);
xlabel('x [m]');
ylabel('y [m]');
grid on;
title(title_string);

