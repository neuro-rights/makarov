function E73
%   SYNTAX
%   E73
%   DESCRIPTION
%   This module provides an upper estimate for eddy currents generated in a
%   heterogeneous conducting medium. 
%   G.M. Noetscher, S.N. Makarov, F. Scire-Scappuzzo, A. Pascual-Leone, "A simple absolute
%   estimate of peak eddy currents induced by TMS using the GR model," IEEE
%   Trans Magn., 2013;49 (9):4999–5003.
%    
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.
%   Module E72

    %   Two elliptical ring axes (major axis and minor axis)
    global a b; a = 0.035; b = 0.035; 
    %   Coil position (in the xy-plane)
    center = [0 0 0];
    %   Uniform conductivity, S/m
    global sigma; sigma = 1.65;
    %   Size of the xy plane 
    global W; W = 0.5;
    %   Separation from the coil
    global Z; Z = 0.75; 
    %   Center frequency, Hz
    global freq; freq = 5000;
    %   Peak current, A
    global I; I = 1800; 
    %   Quantity to plot 
    global type; type = 1; 
    %   Body type
    global body; body = 'person_pregnant';

    %% Step a
    % Create the GUI figure
    f = figure('Units', 'normalized', 'Position', [1/4, 1/4, 1/2, 1/2]);
    % Assign the GUI a name to appear in the window title
    set(f, 'Name', 'E73 - Upper estimate for eddy current density generated in a heterogeneous conducting medium');
    set(f, 'NumberTitle', 'off');

    %% Step b
    % Create the axes window within the figure
    ha = axes('Units', 'normalized', 'Position', [1/12, 1/6, 2/3, 2/3]);

    %% Step c
    % Initialize the GUI: create a default plot in the axes
    staticsW(a, b, center, sigma, W, Z, freq, I, type, body);

    %% Step d1
    % Add an editable field to change the conductivity value
    SIGMA = uicontrol('Style', 'edit', 'Units', 'normalized',...
    'Position', [5/6, 0.86, 0.1, 0.07], 'String', '1.65', ...
    'FontUnits', 'normalized', 'FontSize', 0.4,...
    'BackgroundColor', [1 1 1],...
    'Callback', {@SIGMA_Callback});    

    %% Step d2
    % Add an editable field to change the size of the xy-plane
    Width = uicontrol('Style', 'edit', 'Units', 'normalized',...
    'Position', [5/6, 0.72, 0.1, 0.07], 'String', '0.5', ...
    'FontUnits', 'normalized', 'FontSize', 0.4,...
    'BackgroundColor', [1 1 1],...
    'Callback', {@Width_Callback});    

    %% Step d3
    % Add an editable field to change the plane separation from the coil
    ZC = uicontrol('Style', 'edit', 'Units', 'normalized',...
    'Position', [5/6, 0.58, 0.1, 0.07], 'String', '0.75', ...
    'FontUnits', 'normalized', 'FontSize', 0.4,...
    'BackgroundColor', [1 1 1],...
    'Callback', {@ZC_Callback});    

    %% Step d4
    % Add an editable field to change the signal frequency
    TAU = uicontrol('Style', 'edit', 'Units', 'normalized',...
    'Position', [5/6, 0.44, 0.1, 0.07], 'String', '5000', ...
    'FontUnits', 'normalized', 'FontSize', 0.4,...
    'BackgroundColor', [1 1 1],...
    'Callback', {@FREQ_Callback});

    %% Step d5
    % Add an editable field to change the peak current
    Current = uicontrol('Style', 'edit', 'Units', 'normalized',...
    'Position', [5/6, 0.30, 0.1, 0.07], 'String', '1800', ...
    'FontUnits', 'normalized', 'FontSize', 0.4,...
    'BackgroundColor', [1 1 1],...
    'Callback', {@Current_Callback});

    %% Step d6
    % Add an editable field to change the loop radius a
    RADA = uicontrol('Style', 'edit', 'Units', 'normalized',...
    'Position', [5/6-0.03, 0.16, 0.07, 0.07], 'String', '0.035', ...
    'FontUnits', 'normalized', 'FontSize', 0.4,...
    'BackgroundColor', [1 1 1],...
    'Callback', {@RADA_Callback});

    %% Step d7
    % Add an editable field to change the loop radius b
    RADB = uicontrol('Style', 'edit', 'Units', 'normalized',...
    'Position', [5/6+0.05, 0.16, 0.07, 0.07], 'String', '0.035', ...
    'FontUnits', 'normalized', 'FontSize', 0.4,...
    'BackgroundColor', [1 1 1],...
    'Callback', {@RADB_Callback});
    %%  Step d8 
    % Add a pop-up menu to change quantity to be plotted 
    TYPE = uicontrol('Style', 'popup', 'Units', 'normalized',...
    'Position', [0.20, 0.90, 0.1, 0.08], 'String', 'mag(J)|Jx|Jy|mag(B)|Bx|By|Bz',...
    'FontUnits', 'normalized', 'FontSize', 0.4,...
    'BackgroundColor', [1 1 1],...
    'Callback', {@TYPE_Callback});

    %%  Step d9 
    % Add a pop-up menu to change color palette
    COLOR = uicontrol('Style', 'popup', 'Units', 'normalized',...
    'Position', [0.40, 0.90, 0.1, 0.08], 'String', 'jet|cool|gray',...
    'FontUnits', 'normalized', 'FontSize', 0.4,...
    'BackgroundColor', [1 1 1],...
    'Callback', {@COLOR_Callback});

   %%  Step d10 
    % Add a pop-up menu to change body shape
    BODY = uicontrol('Style', 'popup', 'Units', 'normalized',...
    'Position', [0.60, 0.90, 0.14, 0.08], 'String', 'person_pregnant|person1_arms_up|none',...
    'FontUnits', 'normalized', 'FontSize', 0.4,...
    'BackgroundColor', [1 1 1],...
    'Callback', {@BODY_Callback});
    %% Step e0
    % Add a static text field "hit Enter"
    E1 = uicontrol('Style', 'text', 'Units', 'normalized',...
    'Position', [5/6-0.05, 0.92, 0.2, 0.08], 'String', 'Update value and then hit Enter:',...   
    'FontUnits', 'normalized', 'FontSize', 0.3,...
    'BackgroundColor', [.8 .8 .8]);

    %% Step e1
    % Add a static text field for conductivity
    E2 = uicontrol('Style', 'text', 'Units', 'normalized',...
    'Position', [5/6-0.05, 0.78, 0.2, 0.08], 'String', 'Uniform conductivity, S/m',...   
    'FontUnits', 'normalized', 'FontSize', 0.3,...
    'BackgroundColor', [.8 .8 .8]);

    %% Step e2
    % Add a static text field for the domain size
    E3 = uicontrol('Style', 'text', 'Units', 'normalized',...
    'Position', [5/6-0.05, 0.64, 0.2, 0.08], 'String', 'Size of the xy-plane, m',...   
    'FontUnits', 'normalized', 'FontSize', 0.3,...
    'BackgroundColor', [.8 .8 .8]);

    %% Step e3
    % Add a static text field for the separation
    E4 = uicontrol('Style', 'text', 'Units', 'normalized',...
    'Position', [5/6-0.05, 0.50, 0.2, 0.08], 'String', 'Plane separation from the coil, m',...   
    'FontUnits', 'normalized', 'FontSize', 0.3,...
    'BackgroundColor', [.8 .8 .8]);

    %% Step e4
    % Add a static text field for the frequency
    E5 = uicontrol('Style', 'text', 'Units', 'normalized',...
    'Position', [5/6-0.05, 0.36, 0.2, 0.08], 'String', 'Center frequency, Hz',...
    'FontUnits', 'normalized', 'FontSize', 0.3,...
    'BackgroundColor', [.8 .8 .8]);

   %% Step e5
    % Add a static text field for the peak current
    E6 = uicontrol('Style', 'text', 'Units', 'normalized',...
    'Position', [5/6-0.05, 0.22, 0.2, 0.08], 'String', 'Peak current, A',...
    'FontUnits', 'normalized', 'FontSize', 0.3,...
    'BackgroundColor', [.8 .8 .8]);

   %% Step e6
    % Add a static text field for both loop radii
    E6 = uicontrol('Style', 'text', 'Units', 'normalized',...
    'Position', [5/6-0.05, 0.08, 0.2, 0.08], 'String', 'Semi-major loop axes a, b, m',...
    'FontUnits', 'normalized', 'FontSize', 0.3,...
    'BackgroundColor', [.8 .8 .8]);

    %% Step f1
    % Define the callback function that changes the uniform conductivity
    function SIGMA_Callback(hObject, handles)
        user_entry = str2double(get(hObject, 'string'));
        if isnan(user_entry)
            errordlg('You must enter a numeric value','Bad Input','modal')
            uicontrol(hObject); return;
        end
        % Proceed with callback...
        sigma = user_entry;
        staticsW(a, b, center, sigma, W, Z, freq, I, type, body);        
    end

    %% Step f2
    % Define the callback function that changes the plane size
    function Width_Callback(hObject, handles)
        user_entry = str2double(get(hObject, 'string'));
        if isnan(user_entry)
            errordlg('You must enter a numeric value','Bad Input','modal')
            uicontrol(hObject); return;
        end
        % Proceed with callback...
        W = user_entry;        
        staticsW(a, b, center, sigma, W, Z, freq, I, type, body);        
    end

    %% Step f3
    % Define the callback function that changes the plane separation
    function ZC_Callback(hObject, ~)
        user_entry = str2double(get(hObject, 'string'));
        if isnan(user_entry)
            errordlg('You must enter a numeric value','Bad Input','modal')
            uicontrol(hObject); return;
        end
        % Proceed with callback...
        Z = user_entry;
        staticsW(a, b, center, sigma, W, Z, freq, I, type, body);
    end

    %% Step f4
    % Define the callback function that changes the pulse frequency
    function FREQ_Callback(hObject, handles)
        user_entry = str2double(get(hObject, 'string'));
        if isnan(user_entry)
            errordlg('You must enter a numeric value','Bad Input','modal')
            uicontrol(hObject); return;
        end
        % Proceed with callback...
        freq = user_entry;
        staticsW(a, b, center, sigma, W, Z, freq, I, type, body);
    end

    %% Step f5
    % Define the callback function that changes peak pulse current
    function Current_Callback(hObject, handles)
        user_entry = str2double(get(hObject, 'string'));
        if isnan(user_entry)
            errordlg('You must enter a numeric value','Bad Input','modal')
            uicontrol(hObject); return;
        end
        % Proceed with callback...
        I = user_entry;    
        staticsW(a, b, center, sigma, W, Z, freq, I, type, body);        
    end

    %% Step f6
    % Define the callback function that the semi-major loop axis a
    function RADA_Callback(hObject, handles)
        user_entry = str2double(get(hObject, 'string'));
        if isnan(user_entry)
            errordlg('You must enter a numeric value','Bad Input','modal')
            uicontrol(hObject); return;
        end
        % Proceed with callback...
        a = user_entry;    
        staticsW(a, b, center, sigma, W, Z, freq, I, type, body);        
    end

   %% Step f7
    % Define the callback function that the semi-major loop axis b
    function RADB_Callback(hObject, handles)
        user_entry = str2double(get(hObject, 'string'));
        if isnan(user_entry)
            errordlg('You must enter a numeric value','Bad Input','modal')
            uicontrol(hObject); return;
        end
        % Proceed with callback...
        b = user_entry;    
        staticsW(a, b, center, sigma, W, Z, freq, I, type, body);        
    end

    %% Step f8
    % Define the callback function that changes field quantity to be
    % plotted
    function TYPE_Callback(hObject, handles)        
        type = get(TYPE, 'Value');        
        staticsW(a, b, center, sigma, W, Z, freq, I, type, body);  
    end

    %% Step f9
    % Define the callback function that changes field quantity to be
    % plotted
    function COLOR_Callback(hObject, handles)        
        val = get(COLOR, 'Value');
        if val == 1
            colormap(jet)
        elseif val == 2
            colormap(cool)
        elseif val == 3
            colormap(gray)
        end
    end

  %% Step f10
    % Define the callback function that changes the body shape
    function BODY_Callback(hObject, handles)        
        val = get(BODY, 'Value');
        if val == 1
            body = 'person_pregnant';
        elseif val == 2
            body = 'person1_arms_up';
        elseif val == 3
            body = 'none';
        end        
        staticsW(a, b, center, sigma, W, Z, freq, I, type, body);  
    end
end

    %%  Define the fields
    function staticsW(a, b, center, sigma, W, Z, freq, I, type, body);
    %   STATICSW plots the fields in a plane      
        x = linspace(-W/2, W/2, 25);            %   x-grid
        y = x;                                  %   y-grid
        [X0, Y0] = ndgrid(x, y);
        r0x = reshape(X0, size(X0, 1)*size(X0, 2), 1);
        r0y = reshape(Y0, size(Y0, 1)*size(Y0, 2), 1);
        r0z = -Z*ones(size(r0x));
        r0 = [r0x r0y r0z];
        %   Find the fields
        [A1 B1, Contour1] = staticsH(a, b, center-[a 0 0], Z, r0);
        [A2 B2, Contour2] = staticsH(a, b, center+[a 0 0], Z, r0);
        A = A1 - A2;
        B = B1 - B2;
        Ax = reshape(A(:,1), size(X0));
        Ay = reshape(A(:,2), size(X0));
        Az = reshape(A(:,3), size(X0));
        Bx = reshape(B(:,1), size(X0));
        By = reshape(B(:,2), size(X0));
        Bz = reshape(B(:,3), size(X0));
        Jx = sigma*I*2*pi*freq*Ax;
        Jy = sigma*I*2*pi*freq*Ay;
        Jz = sigma*I*2*pi*freq*Az;
        %   Find quantity to be plotted
        if type == 1
            PLOT = sqrt(Jx.^2 + Jy.^2 + Jz.^2);
        elseif type == 2
            PLOT = Jx;
        elseif type == 3
            PLOT = Jy;
        elseif type == 4
            PLOT = sqrt(Bx.^2 + By.^2 + Bz.^2);         
        elseif type == 5
            PLOT = Bx;
        elseif type == 6
            PLOT = By;
        elseif type == 7
            PLOT = Bz;                
        end
        MAX = max(max(abs(PLOT)));    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Figure
        rotate3d on; view (45, 30);  
        [C, h] = contourf(X0, Y0, PLOT, 15,'LineWidth', 0.1);     
        line(Contour1(1, :), Contour1(2, :), Contour1(3, :), 'Color', 'm', 'LineWidth', 3);
        line(Contour2(1, :), Contour2(2, :), Contour2(3, :), 'Color', 'm', 'LineWidth', 3);  
        if strcmp(body, 'person_pregnant') load person_pregnant; end;
        if strcmp(body, 'person1_arms_up') load person1_arms_up; end;
        if strcmp(body, 'person1_bent_over') load person1_bent_over; end;
        if strcmp(body, 'none') end;
        if exist('P')
            P(:, 3) = P(:, 3)  - max(P(:, 3)) + Z;
            patch('vertices', P, 'faces', t, 'EdgeColor',[0.75 0.75 0.75], 'LineWidth', 1, 'FaceAlpha', 0.01, 'Clipping', 'off');
        end
        axis equal; axis tight; xlabel ('x, m'); ylabel('y, m'); colorbar;    
        axis([-W/2 W/2 -W/2 W/2 0 Z]); grid on;
        if type <=3
            title('Eddy current density distribution, A/m^2 - upper limit');
        else
            title('Magnetic flux distribution, T - upper limit');
        end
        text('Units', 'normalized', 'Position',[0.7, 0.9, 0],'String', strcat('Max. abs. value=', num2str(MAX)));    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    function [A, B, Contour] = staticsH(a, b, center, Z, r0)
    %   STATICSH finds magnetic field generated by an ellipse of a uniform
    %   current 
    %   By Sheila Werth, Kaung Myat Win, and S. Makarov - March 2012
    %   Ellipse is centered in the xy-plane
    %   a - major axis, m
    %   b - minor axis, m
    %   loop current of 1A is assumed
    %   center - lop center (a 1 by 3 vector)
    %   r0 - observation point for the B-field, m (a N by 3 vector) 
    %   A - magnetic vector potential (N, 3), T*m
    %   B - magnetic flux density (N, 3), T*m
        %   EM DATA
        mu0 = 1.25663706e-006;    %    magnetic permeability of vacuum(~air)

        %   DISCRETIZATION
        t = linspace(0, 2*pi, 100);  %   parameteric ellipse form
        dt = t(2) - t(1); t = t + dt/2; t(end) = [];
        %   INTEGRATION
        B = zeros(size(r0, 1), 3);  %   magnetic flux
        A = zeros(size(r0, 1), 3);  %   magnetic vector potential
        rshift0 = r0 - repmat(center, [size(r0, 1) 1]);

        for m = 1:length(t)
            x  = [+a*cos(t(m)) +b*sin(t(m)) 0];       %   parameterization
            dl = [-a*sin(t(m)) +b*cos(t(m)) 0]*dt;    %   parameterization
            dL =      repmat(dl, [size(r0, 1) 1]);
            r  = rshift0 - repmat(x, [size(r0, 1) 1]); 
            R  = dot(r, r, 2);
            R1 = R.^0.5;
            R3 = R.^1.5;
            R1 = repmat(R1, [1, 3]);
            R3 = repmat(R3, [1, 3]);        
            B  = B + cross(dL, r, 2)./R3;
            A  = A + dL./R1;
        end
        B = mu0*B/(4*pi);
        A = mu0*A/(4*pi);    
        Contour(1, :) = a*cos(t)-center(1);
        Contour(2, :) = b*sin(t)-center(2);
        Contour(3, :) = Z*ones(size(t))-center(3);
        Contour(:, end+1) = Contour(:, 1);
    end