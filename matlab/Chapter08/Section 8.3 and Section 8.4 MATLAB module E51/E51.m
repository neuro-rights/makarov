function fig_hdl = E51
%   SYNTAX
%   E51
%   DESCRIPTION
%   This module is an accurate MoM solution for steady-state current flow
%   in simple shapes of finite conductivity with attached electrodes. An
%   arbitrary number of electrodes may be used, with arbitrary voltages.
%   The module also computes surface charge distributions, electric current
%   distribution in a plane, and electric current at a point.
%   
%   Mesh generator for planar shapes: Copyright DISTMESH 2004-2012 Per-Olof
%   Persson. 
%   GUI: Mr. Xingchi Dai for NEVA EM
%   Algorithm and its implementation: Sergey N Makarov
%    
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

%%  EM constants
eps0        = 8.85418782e-012;                  %   dielectric permittivity of vacuum(~air)
mu0         = 1.25663706e-006;                  %   magnetic permeability of vacuum(~air)
c0          = 1/sqrt(eps0*mu0);                 %   speed of light in vacuum(~air)
eta0        = sqrt(mu0/eps0);                   %   vacuum/air impedance, ~377 ohms

%%  GUI window - input parameters
global P;
global P0;
global PP;
global t;
global t0;
global tt;
global cond;
global S;
global c;
global Area;
global Center;
global normals;
global normals0;
global nnormals;
global Size;
global Points;
global epsr;
global contrasts;
global COLOR;   
global TRANSPARENCY;
global ProjectFile; ProjectFile = []; 
global FirstNontrivial;
     
%   Electrodes (12 max)                           
strge.NumberOfElectrodes     = 2;             
strge.TypeOfElectrodes       = ['V'; 'V'; 'V'; 'V'; 'V'; 'V'; 'V'; 'V'; 'V'; 'V'; 'V'; 'V'];
strge.VoltageOfElectrodes    = [1 -1 1 -1 1 -1 1 -1 1 -1 1 -1];                                       
strge.PositionOfElectrodes   = [[0.0 -0.025 0.0]; [0.0 0.025 0.0]; [0.0 -0.1 0.0]; [0.0 -0.1 0.0]; [0.0 -0.1 0.0]; [0.0 -0.1 0.0]; [0.0 -0.1 0.0]; [0.0 -0.1 0.0]; [0.0 -0.1 0.0]];  
strge.RadiusOfElectrodes     = [0.005 0.005 0.005 0.005 0.005 0.005 0.005 0.005 0.005 0.005 0.005 0.005]; 
strge.NumberOfTriangles      = [75 75 75 75 75 75 75 75 75 75 75 75]; 
strge.ParameterOfElectrodes  =  [1 1 1 1 1 1 1 1 1 1 1 1];
%  Bodies (4 max)

strge.color                  = [[0 0.5 0]; [0 0 1]; [1 0 0]; [0 0 1]; [0 1 0]; [0 0 1]; [0 0 1]; [0 0 1]; [0 0 1]; [0 0 1]; [0 0 1]; [0 0 1]]; 
strge.transparency           = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
strge.Conductivity           = [0.5 0.5 0.1 0.7 0.6 0.1 0.6 0.5 0.1 0.1 0.1 0.1];
strge.include{1} = 'yes';
strge.include{2} = 'no';
strge.include{3} = 'no';
strge.include{4} = 'no';
strge.edges{1}    = 'no';
strge.edges{2}    = 'no';
strge.edges{3}    = 'no';
strge.edges{4}    = 'no';
%   Mesh
strge.Indicator     = 0;
strge.OuterSize     = 0; 
strge.Tracker       = [];

strge.FileName{1}      = 'cylinder1_05cm_1cm.mat';
strge.FileName{2}      = 'cylinder2_15cm_7cm.mat';
strge.FileName{3}      = 'sphere1_4cm.mat';
strge.FileName{4}      = 'sphere2_1cm.mat';

strge.DielConst     = ones(4, 1);
strge.PositionX     = zeros(4, 1);
strge.PositionY     = zeros(4, 1);
strge.PositionZ     = zeros(4, 1);
strge.Scaling       = ones(4, 1);

%%  Output graphics parameters
%   Surface charge and general visualization (uitable)
%   Parameters to scale charge distribution for better visual inspection
strsc.positive = 0.2;         %   positive charge densities higher than this number times the average
%   positive charge density are assigned the same value
strsc.negative = 0.2;         %   negative charge densities smaller than this number times the average
%   negative charge density are assigned the same value

%   Visualization of the electric current in the observation plane (uitable)
strop.yes = 'no';           %   include (yes) or not (no) plot of the current
strop.potential = 'yes';    %   include (yes) or not (no) plot of the electric potential
strop.planetype   = 'xy';   %   xy, xz, or yz plane
strop.planex = 0.0;         %   plane center: x in m
strop.planey = 0.0;         %   plane center: y in m
strop.planez = 0.0;         %   plane center: z in m
strop.planesizex   = 0.05;  %   plane length in m
strop.planesizey   = 0.07;  %   plane width in m
strop.divisionsx   = 40;    %   divisions with respect to length
strop.divisionsy   = 40;    %   divisions with respect to width
strop.arrow = 3;            %   relative arrow size versus default size

%   Parameters for the observation point(s) (used to obtain exact values of the field) (uitable)
stroc.yes         = 'no';     %   'yes' - present; 'no' - absent
stroc.x           = 0.0;       %   x position in m
stroc.y           = 0.0;       %   y position in m
stroc.z           = 0.0;       %   z position in m
stroc.size        = 1.5;       %   relative marker size versus default size

%%   Numerical parameters (uitable)
global R; R   = 3;        %    dimensionless radius of an enclosing sphere for precise integration
%    R=0 - only self integrals are calculated precisely
%    R=10- integrals for all neighbor triangles whose
%    center-to-center distances from the observation triangle
%    are less than ten times the effective triangle size are calculated
%    precisely
global gauss; gauss = 4;      %    Number of integration points in the Gaussian quadrature
%    Numbers 1, 4, 7, 13, and 25 are permitted
global precision; precision = 'double'; 

global solver; solver = 'exact';

%%  End of GUI window - input parameters

%   Input Parameters--------
%%  Global parameters of output results
strout.phi          = 0;                        %    Electric potential at the observation point (if any)
strout.PatchesM     = 0;                        %    Number of triangular patches in the metal mesh
strout.PatchesD     = 0;                        %    Number of triangular patches in the dielectric mesh
strout.quality      = 0;                        %    Minimum triangle quality
strout.time1        = 0;                        %    CPU time in sec for filling the MoM matrix
strout.time2        = 0;                        %    CPU time in sec for solving the system of MoM eqs.
strout.ChargeCond  = 0;                         %    Charge of electrodes, C
strout.ChargeDielectric  = 0;                   %    Charge of dielectric object, C
strout.ChargeTotal  = 0;                        %    Sum of charges for the entire structure, C
strout.CurrentElectrodesmA      = [0 0 0 0 0 0 0 0 0 0 0 0];  %    Electrode current
strout.VoltageElectrodes      = [0 0 0 0 0 0 0 0 0 0 0 0];    %    Electrode current
strout.CurrentObsDensity   = [0 0 0];                  %    Current at the obs. point
strout.AreaElectrodes      = [0 0 0 0 0 0 0 0 0 0 0 0];
%    Current at the obs. point

%%  End of GUI window - input parameters
% Initialize handles structure
handles = struct();
handles.geometry = 0;   %   The real handles
handles.simulate = 0;   %   Only an indicator that the simulations are complete
handles.object = 2;
% Create all UI controls
build_gui();
geometry(0);

%%  Functions
  function outputgraphics(par)
        %%  Definition of nodal points in the observation plane (array Points) and initializaton of the E-field/Electric current  
        if (par==0)|(par==1)
            if strcmp(strop.planetype, 'xy') nx = 0; ny = 0; nz = 1; end
            if strcmp(strop.planetype, 'xz') nx = 0; ny = 1; nz = 0; end
            if strcmp(strop.planetype, 'yz') nx = 1; ny = 0; nz = 0; end
            [Points, dummy] = plate0(strop.planesizex, strop.planesizey, strop.divisionsx, strop.divisionsy, 0, nx, ny, nz);
            Points(:, 1) = Points(:, 1) +  strop.planex;
            Points(:, 2) = Points(:, 2) +  strop.planey;
            Points(:, 3) = Points(:, 3) +  strop.planez;
            strop.Points = Points;        
            %   Conductivity/ map            
            Conductivity        = zeros(size(Points, 1), 1);
            inside              = logical(check(P0, t0, normals0, Points));                                         %   column of size Points
            Conductivity(inside)= cond(1);
         
            for count = 1:strge.Indicator
                inside_sub{count}                = logical(check(PP{count}, tt{count}, nnormals{count}, Points));   %   column of size Points
                Conductivity(inside_sub{count})  = cond(count+1);
                inside = inside&(~inside_sub{count});                                                               %   inside the main structure but outside the others
            end
          
            strop.Conductivity  = Conductivity;
            strop.E              = zeros(size(Points));
            strop.Potential      = zeros(size(Points, 1), 1);     
            
            %%  Find the E-field/Electric current in a plane                    
            if strcmp(strop.yes, 'yes')
                msg = 1;
                strop.E(inside, :)  = efield(Points(inside, :), c, P, t, Center, Area, normals, Size, R, eps0, msg);
                for count = 1:strge.Indicator
                    msg = count + 1;
                    strop.E(inside_sub{count}, :)  = efield(Points(inside_sub{count}, :), c, P, t, Center, Area, normals, Size, R, eps0, msg);
                end
                strop.Current      = strop.E.*repmat(strop.Conductivity, 1, 3);
            end        
            %%  Find the potential in a plane
            if strcmp(strop.potential, 'yes')
                msg = 1;
                strop.Potential = potential(Points, c, P, t, Center, Area, normals, Size, R, eps0, msg);
            end       
        end
        if (par==0)|(par==2)
            %%  Find the E-field/Electric current at an observation point(s)            
            strout.phi        = 'none selected';
            if strcmp(stroc.yes, 'yes')
                msg = 0;
                points = [stroc.x stroc.y stroc.z];
                %   Conductivity map
                Conductivity        = zeros(size(points, 1), 1);
                inside              = logical(check(P0, t0, normals0, points));
                Conductivity(inside)= cond(1);
                for count = 1:strge.Indicator
                    temp    = logical(check(PP{count}, tt{count}, nnormals{count}, points));
                    Conductivity(temp)= cond(count+1);
                end
                strout.ConductivityObsPoint  = Conductivity;            
                E = efield(points, c, P, t, Center, Area, normals, Size, R, eps0, msg);
                strout.CurrentObsDensity   = E.*repmat(strout.ConductivityObsPoint, 1, 3);
                strout.phi = potential(points, c, P, t, Center, Area, normals, Size, R, eps0, msg);                
            end
        end
        %%   Output graphics
        set(0, 'CurrentFigure', handles.figure1);
        if handles.geometry ~= 0; delete(handles.geometry); end
        io = 1;
        handles.geometry = subplot(1,1,1,'Parent',handles.uipanel2);
        graphics_E51(io, strop, stroc, P, t, c, Area, strsc, strge, COLOR, TRANSPARENCY);    
    end

    function strout = clearoutput(strout)
        %   Clear results anyways
        strout.phi          = 0;        %   Electric potential at the observation point (if any)
        strout.time1        = 0;        %   CPU time in sec for filling the MoM matrix
        strout.time2        = 0;        %   CPU time in sec for solving the system of MoM eqs.
        strout.ChargeCond   = 0;        %   Total charge o
        strout.ChargeDielectric = 0;    %   Total charge of the dielectric object
        strout.ChargeTotal  = 0;        %   Sum of charges for the entire structure
        strout.CurrentElectrodesmA = [0 0 0 0 0 0 0 0 0 0];    %    Electrode current
        strout.CurrentObsDensity   = [0 0 0];             %    Current at the obs. point
    end

    function geometry(par)        
        count = 0; strge.Tracker = [];
        if ~par
            P = []; t = []; normals = [];     
        end
        FirstNontrivial = [];        
        for m = 1:4
            if strcmp(strge.include{m}, 'yes')                
                S  = load(strge.FileName{m}, '-mat');                
                cond(count+1) = strge.Conductivity(m);
                if count > 0
                    S.normals = -S.normals;   
                end                
                S.P = S.P*strge.Scaling(m);              
                S.P(:, 1) = S.P(:, 1) + strge.PositionX(m);
                S.P(:, 2) = S.P(:, 2) + strge.PositionY(m);
                S.P(:, 3) = S.P(:, 3) + strge.PositionZ(m);  
                temp = S.t + size(P, 1);
                temp(:, 4) = 0;
                if ~par|count>0
                    t = [t; temp];
                    P = [P; S.P];
                    normals   = [normals; S.normals];            
                end                                                                  
                if count == 0   %   outer shape
                    FirstNontrivial = m;
                    t0 = t; normals0 = normals; P0 = P; strge.OuterSize = size(t0, 1);
                    if ~par
                        [P, t, normals, No, si] = electrodes(P, t, normals, strge.NumberOfElectrodes, strge.PositionOfElectrodes, strge.RadiusOfElectrodes, strge.NumberOfTriangles, strge.ParameterOfElectrodes);
                    end                     
                    contrasts     = ones(size(t, 1), 1);          %    always the same for the single body in contact with air      
                    COLOR         = repmat(strge.color(m, :), size(t, 1), 1);         % 
                    TRANSPARENCY  = repmat(strge.transparency(m), size(t, 1), 1);     %    
                    strge.Tracker = [strge.Tracker; m*ones(size(t, 1), 1)];
                else
                    PP{count} = S.P;                   
                    tt{count} = S.t;
                    nnormals{count} = -S.normals;              
                    contrasts = [contrasts; (strge.Conductivity(FirstNontrivial)-strge.Conductivity(m))/...
                                            (strge.Conductivity(FirstNontrivial)+strge.Conductivity(m))*ones(size(S.t, 1), 1)];
                    COLOR     = [COLOR; repmat(strge.color(m, :), size(S.t, 1), 1)]; 
                    TRANSPARENCY    = [TRANSPARENCY; repmat(strge.transparency(m), size(S.t, 1), 1)];                              
                    strge.Tracker = [strge.Tracker; m*ones(size(S.t, 1), 1)];
                end  
                count = count + 1;                 
            end    
        end
        strge.Indicator = count-1;    %   Counted total number of inner meshes excluding the outer shape
       
        NM = length(find(t(:, 4)>0));
        ND = size(t, 1) - NM;
        strout.PatchesD = ND;
        strout.PatchesM = NM;
        c = zeros(size(t, 1), 1);           
        %%  Definition of nodal points in the observation plane (array Points) and initializaton of the E-field
        if strcmp(strop.planetype, 'xy') nx = 0; ny = 0; nz = 1; end
        if strcmp(strop.planetype, 'xz') nx = 0; ny = 1; nz = 0; end
        if strcmp(strop.planetype, 'yz') nx = 1; ny = 0; nz = 0; end
        [Points, dummy] = plate0(strop.planesizex, strop.planesizey, strop.divisionsx, strop.divisionsy, 0, nx, ny, nz);       
        Points(:, 1) = Points(:, 1) +  strop.planex;
        Points(:, 2) = Points(:, 2) +  strop.planey;
        Points(:, 3) = Points(:, 3) +  strop.planez;
        strop.Points = Points;
        strop.E      = zeros(size(Points));
        
        %%  Input graphics
        strout.quality  = min(simpqual(P, t));
        set(0, 'CurrentFigure', handles.figure1);
        if handles.geometry ~= 0; delete(handles.geometry); handles.geometry = 0; end
        handles.geometry = subplot(1,1,1,'Parent',handles.uipanel2);
        io = 0;  
        graphics_E51(io, strop, stroc, P, t, c, Area, strsc, strge, COLOR, TRANSPARENCY);
    end

    function simulate()
        h    = waitbar(0, 'Please wait - filling the MoM matrix');
        time1 = cputime;
        %%  Parameter initialization for the combined mesh
        Metal = [];
        for m = 1:strge.NumberOfElectrodes
            temp = find(t(:, 4)==m);
            Metal = [Metal; temp];
        end        
        NonMetal    = find(t(:, 4)==0);
        t           = [t(Metal, :); t(NonMetal, :)];    %   metal patches up front   
        normals     = [normals(Metal, :); normals(NonMetal, :)];    %   metal patches up front 
        NM          = length(Metal);
        ND          = size(t, 1)- NM;
        %   Matrix storage - up front
        ZD      = zeros(NM+ND, NM+ND, precision);    %   full matrix - covers everything    
        ZM      = zeros(NM, NM+ND, precision); 
        
        %   Other arrays
        Center = zeros(length(t), 3);   %   face center
        Area   = zeros(length(t), 1);   %   face area        
        Size   = zeros(length(t), 1);   %   face size defined as distance from center to furthest vertex        
        
        %%   Find base parameters for all faces (metal or dielectric)  
        %   Metal and dielectric
        
        for m = 1:NM+ND
            Vertexes        = P(t(m, 1:3)', :)';
            r1              = Vertexes(:, 1);
            r2              = Vertexes(:, 2);
            r3              = Vertexes(:, 3);
            tempv           = cross(r2-r1, r3-r1);  %   definition (*)
            temps           = sqrt(tempv(1)^2 + tempv(2)^2 + tempv(3)^2);
            normalcheck     = tempv'/temps;
            if sum(normalcheck.*normals(m, :))<0;   %   rearrange vertices to have exactly the outer normal
                t(m, 2:3)   = t(m, 3:-1:2);         %   by definition (*)
            end
            Area(m)         = temps/2;
            tempc           = (r1+r2+r3)/3;
            Center(m, :)    = tempc';
            Size(m)         = sqrt(max([sum((tempc-r1).^2) sum((tempc-r2).^2) sum((tempc-r3).^2)]));
        end
        
        %%  Accurate calculation of potential integrals and MoM matrix
        [coeffS, weightsS, IndexS]  = tri(7, 5);
        if gauss == 1;  [coeffS, weightsS, IndexS]  = tri(1, 1); end;
        %   if gauss == 3;  [coeffS, weightsS, IndexS]  = tri(3, 2); end;
        if gauss == 4;  [coeffS, weightsS, IndexS]  = tri(4, 3); end;
        %   if gauss == 6;  [coeffS, weightsS, IndexS]  = tri(6, 3); end;
        if gauss == 7;  [coeffS, weightsS, IndexS]  = tri(7, 5); end;
        %   if gauss == 9;  [coeffS, weightsS, IndexS]  = tri(9, 5); end;
        if gauss == 13; [coeffS, weightsS, IndexS]  = tri(13, 7); end;
        if gauss == 25; [coeffS, weightsS, IndexS]  = tri(25, 10); end;
        
        if ~sum(gauss==[1, 4, 7, 13, 25])
            h = errordlg('Number of nodes in the Gaussian quadrature must be 1, 4, 7, 13, or 25');
            return;
        end
        Size                        = sqrt(Size);
        ObsPoint                    = zeros(IndexS*length(t), 3);
        for p =1:IndexS
            ObsPoint(p:IndexS:end, :) = coeffS(1, p)*P(t(:, 1), :) +  coeffS(2, p)*P(t(:, 2), :) +  coeffS(3, p)*P(t(:, 3), :);
            ObsIndex = repmat([1:IndexS], 1, size(t, 1))'; 
            WeightsS = repmat(weightsS,  size(t, 1), 1);
        end
        
        %%  Filling ZM/ZD/Z
        %   Prepare distance matrix DIST for ZM
        %   First row - distances from the center of face #1 to all other face centers
        %   Second row - distances from the center of face #2 to all other face centers, etc.               
        %   MoM matrix ZM initialization (center point integration)    
        
        T = size(t, 1);        
        for m = 1:NM
            temp = Center- repmat(Center(m, :), T, 1);        
            ZM(m, :) = 1./sqrt(sum(temp.*temp, 2));
        end             
        %   MoM matrix ZD initialization (center point integration)            
        for m = 1:NM+ND    %   rowwise here!
            temp        = Center- repmat(Center(m, :), T, 1);    %   this is rn - rm     
            ZD(m, :)    = +(sum(repmat(normals(m, :), T, 1).*temp, 2)./sum(temp.*temp, 2).^1.5)';
        end
        for n = 1:NM+ND
            ZD(:, n) = ZD(:, n)*Area(n);
            ZM(:, n) = ZM(:, n)*Area(n);
        end
        
        %   Loop over columns of impedance matrix ZM/ZD!
        %   n is the number of the inner triangle (every column has the only the n-th inner triangle)
        for n =1:NM+ND
            normdummy    = normals(n, :);
            temp = Center- repmat(Center(n, :), T, 1);
            index   = find(sqrt(sum(temp.*temp, 2))./(Size*Size(n))<=R+1e-16);   % index is local
            index1  = index(index<=NM);                   
            r1      = P(t(n, 1), :);    %   row
            r2      = P(t(n, 2), :);    %   row
            r3      = P(t(n, 3), :);    %   row                   
            m1 = repmat(index1*IndexS-IndexS, [1 IndexS])';            
            m3 = m1(:) + ObsIndex(1:IndexS*length(index1));            
            [I, IRho] = potint(r1, r2, r3, normdummy, ObsPoint(m3, :)); %   I was calculated with the area Area(n)
            ZM(index1, n) = sum(WeightsS(1:length(index1), :).*reshape(I, IndexS, length(index1))', 2);      
            
            m1 = repmat(index*IndexS-IndexS, [1 IndexS])';
            m3 = m1(:) + ObsIndex(1:IndexS*length(index));
            I = potint2(r1, r2, r3, normdummy, ObsPoint(m3, :)); %   I was calculated with the area Area(n)
            J = zeros(length(index)*IndexS, 1);
            for p = 1:IndexS
                J(p:IndexS:end) = sum(I(p:IndexS:end, :).*normals(index, :), 2);
            end
            ZD(index, n) = sum(WeightsS(1:length(index), :).*reshape(J, IndexS, length(index))', 2);                        
            waitbar(n/(NM+ND));       
        end   
        %   Full matrix ZM
        ZM = +ZM/(4*pi*eps0);
        
        %%   Add coefficients to ZD MoM matrix
        for m = 1:NM+ND
            ZD(m, :) = ZD(m, :)*contrasts(m)/(4*pi*eps0);       %   this is (sig1-sig2)/(sig1+sig2)*n*int
            ZD(m, m) = 1/(2*eps0);                              %   this is only the diagonal term (sig/(2*eps))
        end     
        close(h);
        clear ObsPoint ObsIndex WeightsS 
        
        %%  RHS identification and Z-matrix adjustment
        %%  Find partial electrode currents (this is why do we need the full ZD)
        ElectrodeCurrent    = zeros(strge.NumberOfElectrodes, NM+ND);
        ElectrodeVoltage    = zeros(strge.NumberOfElectrodes, NM+ND);
        RHS     =  zeros(size(t, 1), 1);
        for m = 1:strge.NumberOfElectrodes                          
            index = find(t(:, 4)==m);
            ElectrodeCurrent(m, :) = -sum(ZD((index), :).*repmat(Area(index), 1, NM+ND), 1)*...
                                      strge.Conductivity(FirstNontrivial);   %   Total electrode current (after dot with c)
            ElectrodeVoltage(m, :) = +mean(ZM((index), :).*repmat(Area(index), 1, NM+ND), 1)/mean(Area(index));          
                                                                             %   Average "electrode" voltage (after dot with c)
            if strcmp(strge.TypeOfElectrodes(m), 'V')   
                RHS(index)   = strge.VoltageOfElectrodes(m);
                ZD(index, :) = ZM(index, :);      
            else
                %   Voltage of electrodes is now the total current
                CurrentDensity  = strge.VoltageOfElectrodes(m)/sum(Area(index));
                RHS(index)      = CurrentDensity;  
                for n = 1:length(index)
                    ZD(index(n), :)    = -strge.Conductivity(FirstNontrivial)*ZD(index(n), :);
                end
            end         
        end
        LastRow  = sum(ElectrodeCurrent, 1);     %   Total current of all electrodes for current conservation
        clear ZM; 
        
        %% Current conservation law (critical here)
        %   enforce the current conservation law
        for m = 1:NM+ND-1
            RHS(m)    = RHS(m)  - RHS(NM+ND);
            ZD(m, :) = ZD(m,: ) - ZD(NM+ND, :);
        end
        % Last row
        ZD(NM+ND, :)         = LastRow;       
        RHS(NM+ND) = 0;               
        time1 = cputime - time1;
        
%         %% Charge conservation law (critical here)
%         %   enforce the current conservation law
%         for m = 1:NM+ND-1
%             RHS(m)    = RHS(m)  - RHS(NM+ND);
%             ZD(m, :) = ZD(m,: ) - ZD(NM+ND, :);
%         end
%         % Last row
%         ZD(NM+ND, :)         = Area';       
%         RHS(NM+ND) = 0;               
%         time1 = cputime - time1;
        
        %%   Solve MoM equations
        if (length(t)>1000)
            h    = waitbar(0, 'Please wait - solving the MoM equations');
            waitbar(0.5);
        end
        time2 = cputime;
        
        %   Jacobi (diagonal) preconditioner
        DIAG = diag(ZD);
        for m = 1:size(ZD, 1)    %   rowwise
            ZD(m, :) = ZD(m, :)/DIAG(m);
        end
        RHS = RHS./DIAG;               
         
        if strcmp(solver, 'iterative')
            [c, flag] = gmres(ZD, RHS, 10, 1e-3);
        else
             c = ZD\RHS;                 %   solution for charge coefficients
             c = double(c); 
        end                    
        clear ZD;                        %   clear impedance matrix
        time2 = cputime - time2;
        if (length(t)>1000)
            close(h);
        end
        Charge              = c'.*Area';                  %   charges for every face
        ChargeCond          = sum(Charge(t(:, 4)>0));     %   total charge of the conductors, C        
        ChargeDielectric    = sum(Charge(t(:, 4)==0));    %   total charge of dielectric object, C
        ChargeTotal         = ChargeCond + ChargeDielectric;           
        %% Find electrode current             
        for m = 1:strge.NumberOfElectrodes
            strout.AreaElectrodes(m)      = sum(Area((t(:, 4)==m)))*10000;                  %   in cm^2
            Correction = 0.0001*strout.AreaElectrodes(m)/(pi*strge.RadiusOfElectrodes(m)^2);
            strout.CurrentElectrodesmA(m) = Correction*sum(ElectrodeCurrent(m, :).'.*c)*1000;    %current, mA
            % total electrode current
            strout.VoltageElectrodes(m) = sum(ElectrodeVoltage(m, :).'.*c);
            % average voltage at "electrode" surface  
        end       
      
        %%  Fields and output graphics
        handles.simulate = 1;
        outputgraphics(0);
        
        %%  GUI window - output parameters
        strout.PatchesM = NM;                       %   Number of triangular patches in the metal mesh
        strout.PatchesD = ND;                       %   Number of triangular patches in the dielectric mesh
        strout.quality  = min(simpqual(P, t));
        strout.time1    = time1;                    %   CPU time in sec for filling the MoM matrix
        strout.time2    = time2;                    %   CPU time in sec for solving the system of MoM eqs.
        strout.ChargeCond = ChargeCond;             %   Total charge of the conductors        
        strout.ChargeDielectric =ChargeDielectric;  %   Total charge of the dielectric object
        strout.ChargeTotal   = ChargeTotal;         %   Sum of charges for the entire structure
    end

%% ---------------------------------------------------------------------------
    function build_gui()
        % Creation of all uicontrols
        %---Basic Position parameters ---
        
        %% global Position Parameters
        %global Position_local;  Position_local  = [0.2, 0.27, 2.3/12, 6/12];
        global Position_b1;     Position_b1     = [1/12, 0.5/12, 4/12, 1/12];
        global Position_b2;     Position_b2     = [7/12, 0.5/12, 4/12, 1/12];
        global Position_b3;     Position_b3     = [0.1/12, 10.5/12,11.9/12, 1.5/12];
        global Position_table;  Position_table  = [0.1/12, 2/12, 1, 8/12];
        Position_units = 'pixels';
        
        % --- FIGURE -------------------------------------
        handles.figure1 = figure( ...
            'Tag', 'figure1', ...
            'Units', 'characters', ...
            'Position', [103.8 23.0769230769231 135.6 38.4615384615385], ...
            'Name', 'E51- Conductors with direct-current electrodes', ...
            'MenuBar', 'none', ...
            'NumberTitle', 'off', ...
            'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
            'Resize', 'on',...
            'Toolbar','figure');
        set( handles.figure1, ...
            'Units', 'pixels' );
        
        %get your display size
        screenSize = get(0, 'ScreenSize');
        
        %calculate the center of the display
        position = get( handles.figure1, ...
            'Position' );
        position(2) = (screenSize(4)-position(4))/2;
        position(1) = (screenSize(3)-position(3))/2;
        %center the window
        set( handles.figure1, ...
            'Position', position );
        
        %table position
        Position_local = get(handles.figure1,'position');
        Position_local(1) = 0.5*Position_local(1) ;
        Position_local(2) = 0.5*Position_local(2) ;
        Position_local(3) = 0.48*Position_local(3);
        Position_local(4) = Position_local(4);
        
        
        hToolLegend1 = findall(handles.figure1,'tag','Standard.NewFigure');
        hToolLegend2 = findall(handles.figure1,'tag','Standard.FileOpen');
        hToolLegend3 = findall(handles.figure1,'tag','Standard.SaveFigure');
        hToolLegend4 = findall(handles.figure1,'tag','Standard.EditPlot');
        hToolLegend5 = findall(handles.figure1,'tag','Exploration.Pan');
        hToolLegend6 = findall(handles.figure1,'tag','DataManager.Linking');
        hToolLegend7 = findall(handles.figure1,'tag','Annotation.InsertLegend');
        hToolLegend8 = findall(handles.figure1,'tag','Exploration.Brushing');
        hToolLegend9 = findall(handles.figure1,'tag','Standard.PrintFigure');
        hToolLegend10 = findall(handles.figure1,'tag','Exploration.ZoomIn');
        hToolLegend11 = findall(handles.figure1,'tag','Exploration.ZoomOut');
        
        
        delete(hToolLegend1);
        delete(hToolLegend2);
        delete(hToolLegend4);
        delete(hToolLegend5);
        delete(hToolLegend6);
        delete(hToolLegend7);
        delete(hToolLegend8);
        delete(hToolLegend10);
        delete(hToolLegend11);
        
        % Make a Print-view
        hToolbar = findall(gcf,'tag','FigureToolBar');
        hPrintButton = findall(hToolbar,'tag','Standard.PrintFigure');
        set(hPrintButton, 'ClickedCallback','printpreview(gcbf)', ...
            'TooltipString','Print Preview');
        
        
        % --- Menus --------------------------------------
        
        handles.fFileMenu      =   uimenu(...       % File menu
            'Parent', handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Project');
        
%         uimenu(...
%             'Parent', handles.fFileMenu,...
%             'HandleVisibility','callback', ...
%             'Label','New Project',...
%             'Callback', {@New_Callback});
%         
%         
%         function New_Callback(~, ~)
%             E51;
%         end
        
        uimenu(...
            'Parent', handles.fFileMenu,...
            'HandleVisibility','callback', ...
            'Label','Open',...
            'Callback', {@Open_Callback});
        
        function Open_Callback(~, ~)            
            file = uigetfile('*.mat', 'Select project file to open');
            ProjectFile = file;          
            if ~isequal(file, 0)                
                M   = load(file);
                strsc               = M.strsc;
                strge               = M.strge;
                strop               = M.strop;    
                strout              = M.strout;    
                stroc               = M.stroc;
                R                   = M.R;
                gauss               = M.gauss;   
                precision           = M.precision; 
                solver              = M.solver;            
                epsr                = M.epsr;                       
                set(0, 'CurrentFigure', handles.figure1);
                geometry(0);
            end
        end
        
        uimenu(...
            'Parent', handles.fFileMenu,...
            'HandleVisibility','callback', ...
            'Label','Save',...
            'Callback', {@Save_Callback});
        
        function Save_Callback(~, ~)
            if ~isempty(ProjectFile)
                save(ProjectFile);
            else
                save('E51project');
            end
        end
        
        uimenu(...
            'Parent', handles.fFileMenu,...
            'HandleVisibility','callback', ...
            'Label','Save As',...
            'Callback', {@SaveAs_Callback});
        
        function SaveAs_Callback(~,~)
            file = uiputfile('*.mat','Save project file as');
            if ~isequal(file, 0)
                save(file);
            end
        end
         
        handles.fFigureMenu      =   uimenu(...       % File menu
            'Parent', handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Save Figure',...
            'callback',@SaveFigureCallback);
        
        function SaveFigureCallback(~,~) 
           file = uiputfile({'*.fig';'*.png';'*.jpeg'},'Save figure as');
            if ~isequal(file, 0)
                if handles.geometry ~= 0; saveas(handles.geometry,file);end
                if handles.simulate ~= 0; saveas(handles.simulate,file);end
            end
        end
                
        % --- Electrode Parameters Menu ---------------------------
        handles.fFirstMenu   =   uimenu(...             % File menu
            'Parent', handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Electr. Param.');
        
        handles.fParametersMenu      =   uimenu(...       % File menu
            'Parent', handles.fFirstMenu,...
            'HandleVisibility','callback', ...
            'Label','Electrode parameters',...
            'Callback', @fParametersMenuCallback);
        
        function fParametersMenuCallback(~, ~)
            strge0 = strge;
            % Callback function run when_______________________________
            if strge.NumberOfElectrodes>0
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                local = figure('Units', Position_units, 'Position', Position_local);
                set(local, 'Name', 'Electrode setup', 'NumberTitle', 'off', 'MenuBar', 'none');
                b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
                b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
                b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'This is the panel to control electrode parameters';'Press ENTER after changing each parameter'},'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for m = 1:strge.NumberOfElectrodes     
                    if strcmp(strge.TypeOfElectrodes(m), 'I')
                        temp(m) = 1e3*strge.VoltageOfElectrodes(m);
                    else
                        temp(m) = strge.VoltageOfElectrodes(m);
                    end
                    data(m, :) =  {strge.TypeOfElectrodes(m), temp(m), strge.PositionOfElectrodes(m, 1), strge.PositionOfElectrodes(m, 2), strge.PositionOfElectrodes(m, 3),...
                                                                 strge.RadiusOfElectrodes(m), strge.NumberOfTriangles(m), strge.ParameterOfElectrodes(m)};                  
                end
                columnname =   {'Type (V or I)' 'V(V) or I(mA)', 'Pos. x, m', 'Pos. y, m', 'Pos. z, m', 'Rad., m', 'No of triangles', 't-Size: boundary vs. center'};
                columnformat = {'char', 'char', 'char', 'char', 'char', 'char', 'char', 'char'};
                columneditable =  [true true true true true true true true];

                uitable('Units', 'normalized', 'Position', Position_table,...
                    'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                    'ColumnName', columnname, 'ColumnFormat', columnformat, ...
                    'ColumnEditable', columneditable, 'ToolTipString', 'Change electrodes',...
                    'ColumnWidth',{85,85,85,85,85},...
                    'CellEditCallback', {@edit_callback3});
            else
                h = msgbox('The number of electrodes should be greater than 1', 'The number of layers is less than 2', 'warn');
            end
            
            function edit_callback3(hObject, ~)
                data = get(hObject,'Data');                
                for m = 1:strge.NumberOfElectrodes
                    strge.TypeOfElectrodes(m)           = data{m, 1};
                    strge.VoltageOfElectrodes(m)        = data{m, 2};
                    strge.PositionOfElectrodes(m, 1)    = data{m, 3};                
                    strge.PositionOfElectrodes(m, 2)    = data{m, 4};
                    strge.PositionOfElectrodes(m, 3)    = data{m, 5};
                    strge.RadiusOfElectrodes(m)         = data{m, 6};             
                    strge.RadiusOfElectrodes(m)         = data{m, 6};     
                    strge.NumberOfTriangles(m)          = data{m, 7};
                    strge.ParameterOfElectrodes(m)      = data{m, 8};
                    if strcmp(strge.TypeOfElectrodes(m), 'I')
                        strge.VoltageOfElectrodes(m) = 1e-3*strge.VoltageOfElectrodes(m);
                    end
                end                
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
%                 ind1 = sum(sum(abs(strge.PositionOfElectrodes - strge0.PositionOfElectrodes)));
%                 ind2 = sum(abs(strge.RadiusOfElectrodes - strge0.RadiusOfElectrodes));
%                 ind3 = sum(abs(strge.NumberOfTriangles - strge0.NumberOfTriangles));
%                 ind4 = sum(~strcmp(strge.TypeOfElectrodes, strge0.TypeOfElectrodes));
%                 if ind1>0|ind2>0|ind3|ind4
%                     geometry(0);              
%                 end
                geometry(0);
                strout = clearoutput(strout);
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strge = strge0;
                delete(local);
            end
            
        end
        
        %----- The Menu of human body --------------       
        handles.fDielectricMenu   =   uimenu(...            %    File menu
            'Parent', handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Cond. Setup');
     
         handles.fDielParametersMenu      =   uimenu(...       % File menu
        'Parent', handles.fDielectricMenu,...
        'HandleVisibility','callback', ...
        'Label','Conducting objects and their parameters',...
        'Callback', @fDielParametersMenuCallback);
    
         function fDielParametersMenuCallback(~, ~)
            strge0 = strge;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setup', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'This is the panel for cond. shapes (may be replaced by any other shapes)';'Press ENTER after changing each parameter'},'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
            for m = 1:4               
                data(m, :) =  {strge.include{m}, strge.Conductivity(m), strge.Scaling(m), strge.PositionX(m), strge.PositionY(m), strge.PositionZ(m), strge.color(m, 1), strge.color(m, 2), strge.color(m, 3), strge.transparency(m), strge.edges{m}};                  
            end
            columnname      =   {'Include', 'Cond., S/m', 'Scale', 'Pos. X, m', 'Pos. Y, m', 'Pos. Z, m', 'Color (r)', 'Color (g)', 'Color (b)', 'Inv. Transpar.', 'Show mesh'};
            columnformat    =   {'char', 'char', 'char', 'char', 'char', 'char', 'char', 'char', 'char', 'char', 'char'};
            columneditable  =  [true];  
            rnames{1, 1}    = 'Cylinder (5x1cm)';
            rnames{1, 2}    = 'Cylinder (15x7cm)';
            rnames{1, 3}    = 'Sphere (4cm)';
            rnames{1, 4}    = 'Sphere (1cm)';            

            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'ColumnName', columnname, 'ColumnFormat', columnformat,  'RowName', rnames,...
                'ColumnEditable', columneditable, 'ToolTipString', 'Change body composition',...
                'ColumnWidth',{60, 75, 75, 75, 75, 80, 80},...
                'CellEditCallback', {@edit_callback4});
                    
            function edit_callback4(hObject, ~)
                data = get(hObject,'Data');                
                for m = 1:4  
                    strge.include{m}            = data{m, 1};
                    strge.Conductivity(m)       = data{m, 2};  
                    strge.Scaling(m)            = data{m, 3};
                    strge.PositionX(m)          = data{m, 4};
                    strge.PositionY(m)          = data{m, 5};
                    strge.PositionZ(m)          = data{m, 6};
                    strge.color(m, 1)           = data{m, 7};
                    strge.color(m, 2)           = data{m, 8};
                    strge.color(m, 3)           = data{m, 9};    
                    strge.transparency(m)       = data{m, 10};
                    strge.edges{m}              = data{m, 11};
                end   
            end
            
            function B1_Callback(~, ~, ~)
                geometry(0);
                strout = clearoutput(strout);
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strge = strge0;
                delete(local);
            end
            
        end         
                  
        % --- Output ----
        handles.fOutputSetupMenu    =   uimenu(...       % File menu
            'Parent',handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Output Setup',...
            'callback',@fOutputSetupCallback);
        
        function fOutputSetupCallback(~, ~)
            strop0 = strop;
            strsc0 = strsc;
            stroc0 = stroc;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Output Setup', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'This is the panel to setup the parameters for the output graphics (what to show)';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Plot current density in a plane',strop.yes;...
                'Plot potential (16 levels) in a plane',strop.potential;...
                'Plane type (xy,xz,yz)',strop.planetype;...
                'Plane center x,m',strop.planex;...
                'Plane center y,m',strop.planey;...
                'Plane center z,m',strop.planez;...
                'Plane length in m',strop.planesizex;...
                'Plane width in m',strop.planesizey;...
                'Relative arrow size versus default size',strop.arrow;...                
                'Scale charge, positive bound',strsc.positive;...
                'Scale charge, negative bound',strsc.negative;...
                'Observation point(s) display',stroc.yes;...
                'x position of observation point, m',stroc.x;...
                'y position of observation point, m',stroc.y;...
                'z position of observation point, m',stroc.z;...
                'Relative marker size versus default size',stroc.size;...
                'Plane-divisions with respect to length',strop.divisionsx;...
                'Plane-divisions with respect to width',strop.divisionsy};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'Change Upper Plate shape',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback1});
            function edit_callback1(hObject, ~)
                data = get(hObject,'Data');
                if (data{28}>1)|(data{29}>1)
                    h = errordlg('Charge scaling factors must be between 0 and 1');
                    return;
                end                         
                strop.yes               = data{19};
                strop.potential         = data{20};
                strop.planetype         = data{21};
                strop.planex            = data{22};
                strop.planey            = data{23};
                strop.planez            = data{24};
                strop.planesizex        = data{25};
                strop.planesizey        = data{26};
                strop.arrow             = data{27};
                strsc.positive          = data{28};
                strsc.negative          = data{29};
                stroc.yes               = data{30};
                stroc.x                 = data{31};
                stroc.y                 = data{32};
                stroc.z                 = data{33};
                stroc.size              = data{34};
                strop.divisionsx        = data{35};
                strop.divisionsy        = data{36};                      
            end
          function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);   
                strOutGraphics.Point = (stroc.x~=stroc0.x)|(stroc.y~=stroc0.y)|(stroc.z~=stroc0.z)|(~strcmp(stroc.yes, stroc0.yes)&strcmp(stroc.yes, 'yes'));
                strOutGraphics.Point = strOutGraphics.Point&strcmp(stroc.yes, 'yes');
                strOutGraphics.Plane = (strop.planex~=strop0.planex)|(strop.planey~=strop0.planey)|(strop.planez~=strop0.planez)...
                                       |(~strcmp(strop.planetype, strop0.planetype))|(strop.divisionsx~=strop0.divisionsx)|(strop.divisionsy~=strop0.divisionsy)...
                                       |(strop.divisionsx~=strop0.divisionsx)|(strop.divisionsy~=strop0.divisionsy)|(~strcmp(strop.yes, strop0.yes)&strcmp(strop.yes, 'yes')); 
                strOutGraphics.Plane = strOutGraphics.Plane&strcmp(strop.yes, 'yes');
                if handles.simulate
                    if strOutGraphics.Point&strOutGraphics.Plane
                        outputgraphics(0);
                    end
                    if ~strOutGraphics.Point&strOutGraphics.Plane
                        outputgraphics(1);
                    end
                    if  strOutGraphics.Point&~strOutGraphics.Plane
                        outputgraphics(2);
                    end
                    if ~strOutGraphics.Point&~strOutGraphics.Plane
                        outputgraphics(3);
                    end                    
                else
                    geometry(0);
                end
                delete(local);
            end
            
            
            function B2_Callback(~, ~, ~)
                strop = strop0;
                strsc = strsc0;
                stroc = stroc0;
                delete(local);
            end
            
        end
        
        % -- Modeling Setup ----
        
        handles.fModelingSetupMenu    =   uimenu(...       % File menu
            'Parent',handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Modeling Setup',...
            'callback',@fModelingSetupCallback);
        
        function fModelingSetupCallback(~, ~)
            R_temp         = R;
            gauss_temp     = gauss;
            precision_temp = precision; 
            solver_temp = solver; 
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Modeling Setup', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'This is the panel to setup the numerical modeling parameters';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Dimensionless radius of integration sphere',R;...
                '# of integration points',gauss;...
                'MoM matrix precision (double/single)', precision;...
                'Solver (exact or iterative-faster)', solver};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'Setup Modeling parameters',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback2});
            
            function edit_callback2(hObject, ~)
                data = get(hObject,'Data');
                R                       = data{5};
                gauss                   = data{6};
                precision               = data{7};
                solver                  = data{8};
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry(1);
                strout = clearoutput(strout);
                delete(local);
            end            
            
            function B2_Callback(~, ~, ~)
                R = R_temp;
                gauss = gauss_temp;
                precision = precision_temp; 
                solver = solver_temp;
                delete(local);
            end
            
        end
                
        % -- Output Results -----
        handles.fOutputDataMenu    =   uimenu(...       % File menu
            'Parent',handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Output Data',...
            'callback',@fOutputDataCallback);
        
        function fOutputDataCallback(~, ~)
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Output Data', 'NumberTitle', 'off', 'MenuBar', 'none');
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Exit', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', 'Output data and results. When current electrodes(sources) are used, the average voltage at the source surface is reported ', 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ColumnName =   {'Name',   'Value'};
            columneditable  =  [true];
            data = {...             
                'Total charge of the electrodes, C',strout.ChargeCond;...           
                'Total charge of the body surface, C',strout.ChargeDielectric;...
                'Sum of charges for entire structure, C',strout.ChargeTotal;... 
                'Current x-comp. at obs. point, A/m^2',strout.CurrentObsDensity(1);...
                'Current y-comp. at obs. point, A/m^2',strout.CurrentObsDensity(2);...
                'Current z-comp. at obs. point, A/m^2',strout.CurrentObsDensity(3);...
                'Electric potential at obs. point, V',strout.phi;...
                'Number of triangular patches in the metal mesh',strout.PatchesM;...
                'Number of triangular patches in the dielectric mesh',strout.PatchesD;...
                'Minimum triangle quality in the mesh',strout.quality;...
                'CPU time in sec for filling the MoM matrix',strout.time1;...
                'CPU time in sec for solving the system of MoM eqs.',strout.time2};            
            
            if (strcmp(stroc.yes,'no'))== 1; data{14+strge.NumberOfElectrodes} = 'N.A.'; data{15+strge.NumberOfElectrodes} = 'N.A.'; data{16+strge.NumberOfElectrodes} = 'N.A.'; data{17+strge.NumberOfElectrodes} = 'N.A.'; end;
            
            data1{1, 1} = 'empty';
            data1{1, 2} = 'empty';
            data2{1, 1} = 'empty';
            data2{1, 2} = 'empty';
            data3{1, 1} = 'empty';
            data3{1, 2} = 'empty';            
            for m = 1:strge.NumberOfElectrodes
                string = strcat('Current, mA at electrode #', num2str(m));
                data1{m, 1} = string;
                data1{m, 2} = strout.CurrentElectrodesmA(m);
                string = strcat('Voltage, V at electrode #', num2str(m));
                data2{m, 1} = string;
                data2{m, 2} = strout.VoltageElectrodes(m);  
                string = strcat('Electrode area, cm^2 #', num2str(m));
                data3{m, 1} = string;
                data3{m, 2} = strout.AreaElectrodes(m);
            end          
            data = [data1; data2; data3; data];                          

            handles.outputtable = uitable('Units', 'normalized', 'Position', Position_table,...
                'Data',data,'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'Change Upper Plate shape',...
                'ColumnWidth',{200,100});
            
            function B2_Callback(~, ~, ~)
                delete(local);
            end
            
        end
       % --- View Direction ---
        handles.fViewMenu    =   uimenu(...       % File menu
            'Parent',handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','View',...
            'callback',@fViewCallback);
        
        function fViewCallback(~,~)
            
        % --- RADIO BUTTONS -------------------------------------
       handles.figure3 = figure( ...
                'Tag', 'figure1', ...
                'Units', 'characters', ...
                'Position', [103.8 38.6153846153846 28.4 20], ...
                'Name', 'View direction', ...
                'MenuBar', 'none', ...
                'NumberTitle', 'off', ...
                'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
                'Resize', 'on');
            set( handles.figure3, ...
                'Units', 'pixels' );
            
            %calculate the center of the display
            position = get(handles.figure1,'position');
            position(1) = position(1)*0.75;
            position(2) = position(2)*1.3;
            position(3) = position(3)*0.25;
            position(4) = position(4)*0.5;
            %center the window
            set( handles.figure3, ...
                'Position', position ); 
            
            handles.uipanel3 = uibuttongroup( ...
                'Parent', handles.figure3, ...
                'Tag', 'uipanel1', ...
                'Units', 'normalized', ...
                'Position', [0.0141823161189358 0.0334928229665072 0.965277777777778 0.961722488038277], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'Title', 'View direction');
            
            set(handles.uipanel3,'SelectionChangeFcn',@uibuttongroup1_SelectionChangeFcn);
            
            
            function uibuttongroup1_SelectionChangeFcn(~,eventdata)
                switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
                    case 'Current'
                        set(0, 'CurrentFigure', handles.figure1);                 
                        view(46, 41);
                    case 'Front'
                        set(0, 'CurrentFigure', handles.figure1);
                        view(90,0);
                    case 'Rear'
                        set(0, 'CurrentFigure', handles.figure1);
                        view(270,0);
                    case 'Top'
                        set(0, 'CurrentFigure', handles.figure1);
                        view(90,90);
                    case 'Bottom'
                        set(0, 'CurrentFigure', handles.figure1);
                        view(90,-90);
                    case 'Left'
                        set(0, 'CurrentFigure', handles.figure1);
                        view(0,0);
                    case 'Right'
                        set(0, 'CurrentFigure', handles.figure1);
                        view(180,0);
                    otherwise
                        errordlg('Fatal error! Contact us!');
                end
            end
                   
            handles.radiobutton7 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'Current', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.86 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Default');
            
            
            handles.radiobutton8 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'Front', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.72 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Front');
                        
            handles.radiobutton9 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'Rear', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.58 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Rear');
                        
            handles.radiobutton10 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'Top', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.44 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Top');
            
            handles.radiobutton11 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'Bottom', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.28 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Bottom');
            
            handles.radiobutton12 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'Left', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.14 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Left');
                
            handles.radiobutton13 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'Right', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.00 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Right');
        end
                
        %--- Help Menu ----
        handles.fHelpMenu    =   uimenu(...       % File menu
            'Parent',handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Help',...
            'callback',@fHelpSetupCallback);
        
       function fHelpSetupCallback(~, ~)
            HelpPath = 'index_s51.html';
            web(HelpPath);
        end
        
        
        % --- PANELS -------------------------------------
        handles.uipanel1 = uipanel( ...
            'Parent', handles.figure1, ...
            'Tag', 'uipanel1', ...
            'Units', 'normalized', ...
            'Position', [0.0162241887905605 0.822 0.967551622418879 0.156], ...
            'FontSize', 7.6666666666666, ...
            'FontUnits', 'pixels', ...
            'Title', 'Abstract');
        
        handles.uipanel2 = uipanel( ...
            'Parent', handles.figure1, ...
            'Tag', 'uipanel2', ...
            'Units', 'normalized', ...
            'Position', [0.0162241887905605 0.104 0.967551622418879 0.718], ...
            'FontSize', 7.6666666666666, ...
            'FontUnits', 'pixels', ...
            'Title', 'Figure');
        % --- ABSTRACT PANEL----------------------------------
        handles.abstract = uicontrol( ...
            'Parent', handles.uipanel1, ...
            'Tag', 'Abstract Description', ...
            'HorizontalAlignment', 'left',...
            'Style', 'text', ...
            'Units', 'normalized', ...
            'Position', [0 0.2 1 0.7], ...
            'FontSize', 8, ...
            'FontUnits', 'pixels', ...
            'String', 'This module is an accurate MoM solution for steady-state current flow in simple shapes of finite conductivity with attached electrodes. An arbitrary number of electrodes may be used, with arbitrary voltages/currents. The shape may contain another shape of a different conductivity. The module also computes surface charge distributions, electric current distribution in a plane and at a point.');

        % --- PUSHBUTTONS -------------------------------------
        handles.pushbutton2 = uicontrol( ...
            'Parent', handles.figure1, ...
            'Tag', 'pushbutton2', ...
            'Style', 'pushbutton', ...
            'Units', 'normalized', ...
            'Position', [0.70 0.0271317829457364 0.105413105413105 0.0542635658914729], ...
            'FontSize', 7.6666666666666, ...
            'FontUnits', 'pixels', ...
            'String', 'Reload',...
            'callback',@GeometryCallback);
        
        function GeometryCallback(~,~)            
            handles.simulate = 0;
            strout = clearoutput(strout);
            geometry(0);     
        end
        
        handles.pushbutton3 = uicontrol( ...
            'Parent', handles.figure1, ...
            'Tag', 'pushbutton3', ...
            'Style', 'pushbutton', ...
            'Units', 'normalized', ...
            'Position', [0.85 0.0271317829457364 0.106837606837607 0.0542635658914729], ...
            'FontSize', 7.6666666666666, ...
            'FontUnits', 'pixels', ...
            'String', 'Exit',...
            'Callback',@ExitCallback);
        
        function ExitCallback(~,~)            
            button = questdlg('Save project before closing?');            
            switch button
                case 'Yes',
                   save('E51project');
                    close all;
                case'No',
                    close all;
                case 'Cancel',
                    quit cancel;
            end            
        end
            
        handles.pushbutton1 = uicontrol( ...
            'Parent', handles.figure1, ...
            'Tag', 'pushbutton1', ...
            'Style', 'pushbutton', ...
            'Units', 'normalized', ...
            'Position', [0.55 0.0271317829457364 0.105413105413105 0.0542635658914729], ...
            'FontSize', 7.6666666666666, ...
            'FontUnits', 'pixels', ...
            'String', 'Simulate',...
            'callback',@SimulateCallback);
        
        function SimulateCallback(~,~)
            if strge.NumberOfElectrodes<2
                errordlg('You must introduce electrodes','Bad Input','modal');
                return;
            end
            simulate();
        end
        
        handles.pushbutton4 = uicontrol( ...
            'Parent', handles.figure1, ...
            'Tag', 'pushbutton3', ...
            'Style', 'pushbutton', ...
            'Units', 'normalized', ...
            'Position', [0.15 0.0271317829457364 0.15 0.0542635658914729], ...
            'FontSize', 7.6666666666666, ...
            'FontUnits', 'pixels', ...
            'String', 'Separate window',...
            'Callback',@OpenZoomCallback);
        
        function OpenZoomCallback(~,~)
             handles.figure2 = figure( ...
                'Tag', 'figure2', ...
                'Units', 'characters', ...
                'Position', [103.8 23.0769230769231 50 20], ...
                'Name', 'Figure window', ...
                'MenuBar', 'figure', ...
                'NumberTitle', 'off', ...
                'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
                'Resize', 'on',...
                'Toolbar','figure');
            set( handles.figure2, ...
                'Units', 'pixels' );
            %calculate the center of the display
            position = get(handles.figure1,'position');
            position(1) = position(1)*0.75;
            position(2) = position(2)*0.65;        
            position(3) = position(3)*0.75;
            position(4) = position(4)*0.75;
            %center the window
            set( handles.figure2, 'Position', position ); 
            if handles.geometry ~= 0 ; copyobj(handles.geometry,handles.figure2);end        
        end        
        
        handles.pushbutton5 = uicontrol( ...
            'Parent', handles.figure1, ...
            'Tag', 'pushbutton3', ...
            'Style', 'pushbutton', ...
            'Units', 'normalized', ...
            'Position', [0.35 0.0271317829457364 0.15 0.0542635658914729], ...
            'FontSize', 7.6666666666666, ...
            'FontUnits', 'pixels', ...
            'String', 'Set electrodes',...
            'Callback',@ElectrodeCallback);
        
        function ElectrodeCallback(~,~)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
            prompt = {'Required number of electrodes'};
            dlg_title = 'Required number of electrodes';
            num_lines = 1;
            def = {num2str(2)};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            user_entry = str2double(answer);
            if isnan(user_entry)
                errordlg('You must enter a numeric value','Bad Input','modal');
                return;
            end
            if user_entry<=1
                errordlg('Minimum number of electrodes is two','Bad Input','modal');
                return;
            end
            if ~isempty(user_entry)
                if user_entry>strge.NumberOfElectrodes
                    for m = 1:user_entry
                        strge.RadiusOfElectrodes(m) = strge.RadiusOfElectrodes(1);  
                    end                                            
                end
                strge.NumberOfElectrodes = user_entry; 
            end       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            handles.figure3  = figure( ...
                'Tag', 'figure2', ...
                'Units', 'characters', ...
                'Position', [103.8 23.0769230769231 50 20], ...
                'Name', 'Figure window', ...
                'MenuBar', 'none', ...
                'NumberTitle', 'off', ...
                'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
                'Resize', 'on',...
                'Toolbar','figure');
            set(handles.figure3, 'Units', 'pixels' );
            %calculate the center of the display
            position = get(handles.figure1,'position');
            position(1) = position(1)*0.75;
            position(2) = position(2)*0.65;        
            position(3) = position(3)*0.75;
            position(4) = position(4)*0.75;
            %center the window
            set(handles.figure3, 'Position', position );       
            %   Clear history and arrays P, t, normals (only load the outer shape)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            P           = P0; 
            t           = t0; 
            t(:, 4)     = 0;
            normals     = normals0;            
            P = P*strge.Scaling(1);
            P(:, 1) = P(:, 1) + strge.PositionX(1);
            P(:, 2) = P(:, 2) + strge.PositionY(1);
            P(:, 3) = P(:, 3) + strge.PositionZ(1);                                   
            t(:, 4) = 0;                    
            %   Identify electrode nodes
            for m = 1:strge.NumberOfElectrodes
                title('Mouse click and hit Enter');
                patch('Faces', t(:, 1:3), 'Vertices', P, 'FaceColor', [1 1 0], 'EdgeColor', 'k', 'FaceAlpha', 1.0);
                axis 'equal';    axis 'tight', set(gca, 'YDir','normal');
                xlabel('x, m'); ylabel('y, m'); view(46, 41); grid on;
                dcm_obj = datacursormode(handles.figure3);
                set(dcm_obj,'DisplayStyle', 'datatip', 'SnapToDataVertex', 'on', 'Enable','on');
                pause  
                c_info = getCursorInfo(dcm_obj);
                position = c_info.Position;
                strge.PositionOfElectrodes(m, :) = position;             
                clf;
                figure(handles.figure3);
            end                               
            [P, t, normals, No, si] = electrodes(P, t, normals, strge.NumberOfElectrodes, strge.PositionOfElectrodes, strge.RadiusOfElectrodes, strge.NumberOfTriangles, strge.ParameterOfElectrodes);
            
            N       = size(t, 1);
            V       = zeros(N, 1);
            Color   = repmat([1 1 0], N, 1);
            for m = 1:strge.NumberOfElectrodes    
                %   Update color array
                index = find(t(:, 4)==m);
                for n = 1:length(index)
                    Color(index(n), :) = [strge.VoltageOfElectrodes(m)-min(strge.VoltageOfElectrodes) 0.3 max(strge.VoltageOfElectrodes)-strge.VoltageOfElectrodes(m)];
                end   
            end
            title('Hit Enter');
            patch('Faces', t(:, 1:3), 'Vertices', P, 'FaceColor', 'flat', 'FaceVertexCData', Color);                
            axis 'equal';    axis 'tight', set(gca, 'YDir','normal');
            xlabel('x, m'); ylabel('y, m'); view(41, 46); grid on;            
            pause
            delete(handles.figure3);
            geometry(1);    
            figure(handles.figure1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end                
        %   End of electrode callback   
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Pnew] = rotatex(P, xangle)
    %%   Rotate the mesh about the x-axis
    %   The local coordinate system is used with the origin at the center of gravity
    %   SNM Summer 2012
    anglex = xangle/180*pi;
    %   Rotation about the x-axis  
    Pnew(:, 1) = +P(:, 1);
    Pnew(:, 2) = +P(:, 2)*cos(anglex) - P(:, 3)*sin(anglex);
    Pnew(:, 3) = +P(:, 2)*sin(anglex) + P(:, 3)*cos(anglex);
    Pnew(:, 2) = Pnew(:, 2);
    Pnew(:, 3) = Pnew(:, 3);
end

function [Pnew] = rotatey(P, yangle)
    %%   Rotate the mesh about the y-axis
    %   The local coordinate system is used with the origin at the center of gravity
    %   SNM Summer 2012
    angley = yangle/180*pi;
    %   Rotation about the y-axis    
    Pnew(:, 1) = +P(:, 1)*cos(angley) + P(:, 3)*sin(angley);
    Pnew(:, 2) = +P(:, 2);
    Pnew(:, 3) = -P(:, 1)*sin(angley) + P(:, 3)*cos(angley);
    Pnew(:, 1) = Pnew(:, 1);
    Pnew(:, 3) = Pnew(:, 3);
end

function [Pnew] = rotatez(P, zangle)
    %%   Rotate the mesh about the z-axis
    %   The local coordinate system is used with the origin at the center of gravity
    %   SNM Summer 2012
    anglez = zangle/180*pi;
    %   Rotation about the z-axis         
    Pnew(:, 1)  = +P(:, 1)*cos(anglez) - P(:, 2)*sin(anglez);
    Pnew(:, 2)  = +P(:, 1)*sin(anglez) + P(:, 2)*cos(anglez);
    Pnew(:, 3)  = +P(:, 3); 
    Pnew(:, 1) = Pnew(:, 1);
    Pnew(:, 2) = Pnew(:, 2);
end

function [P, t] = plate0(L, W, Nx, Ny, indicator, nx, ny, nz)
    %   Create structured mesh for a base plate on the size L by W with the normal vector nx, ny, nz  
    %   Use indicator=1 for a non-uniform mesh
    %   SNM Summer 2012
    %   Create mesh in the xy-plane
    if indicator == 0
        x = [0:Nx]/Nx*L;    %   uniform grid
        y = [0:Ny]/Ny*W;    %   uniform grid
    else
        x = L/2*(1 - cos(pi*[0:Nx]/Nx));    %   non-uniform grid
        y = W/2*(1 - cos(pi*[0:Ny]/Ny));    %   non-uniform grid
    end
    P(:, 1) = repmat(x, 1, length(y))';
    P(:, 2) = reshape(repmat(y, length(x), 1), 1, length(x)*length(y))';
    P(:, 1) = P(:, 1) - mean(P(:, 1));
    P(:, 2) = P(:, 2) - mean(P(:, 2));
    P(:, 3) = 0;
    %   Define triangles 
    t = []; temp = [];
    for n = 1:Ny
        for m = 1:Nx
            t1 = [m m+1    m+Nx+2]' + (n-1)*(Nx+1);
            t2 = [m m+Nx+1 m+Nx+2]' + (n-1)*(Nx+1);
            temp = [temp t1 t2];
        end
    end
    t = temp'; 
    %   Rotate mesh
    theta   = acos(nz)*sign(ny + eps);              %   rotation about the x-axis following the right-hand rule
    phi     = asin(nx/sqrt(nx^2+ny^2+eps));         %   rotation about the z-axis (the same way) 
    theta   = theta/pi*180;
    phi     = phi/pi*180; 
    P       = rotatex(P, theta);
    P       = rotatez(P, phi);
end

function [P] = circle0(R, Tr, par)    
    %   Create a uniform or nonuniform (ratio 1:5) mesh for a base circle 
    %   of radius R, normal vector nx, ny, nz, and with approximately M 
    %   triangles
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson
    %   Adopted and simplified: SNM Summer 2012 
    if Tr<20
        %   Create a simple structured triangular mesh
        M = Tr-1;
        x = R*[cos(2*pi*[0:M]/(M+1)) 0];
        y = R*[sin(2*pi*[0:M]/(M+1)) 0];
        P(:, 1) = x';
        P(:, 2) = y';
        P(:, 3) = 0;
        t(:, 1) = [1:M+1  ]';
        t(:, 2) = [2:M+1 1]';
        t(:, 3) = (M+2);
        return;
    end
    
    h0    = 2.5*R/sqrt(Tr);                             %   Approximate desired edge length
    xmin  = -R; xmax = +R; ymin = -R; ymax = +R;        %   Bounding box

    dptol   = 0.0005;                                   %   Termination criterion
    ttol = 0.1;                                         %   Criterion for retriangulation                  
    Fscale  = 1.2;                                      %   Inner "internal pressure" (increased by 0.01)
    deltat  = 0.2;                                      %   Ratio between force and movement
    geps    = 0.001*h0;                                 %   Relative boundary tolerance
    deps    = sqrt(eps)*h0;                             %   Accuracy of gradient calculation

    %   Create initial distribution in bounding box 
    %   Equilateral triangles - faster convergence
    [x, y] = meshgrid(xmin:h0:xmax, ymin:h0*sqrt(3)/2:ymax);
    x(2:2:end, :) = x(2:2:end, :) + h0/2;                 %   Shift even rows
    P = [x(:), y(:)];                                     %   List of vertices

    %   Remove vertices outside the region
    dist = dcircle(P, R);                               %   A column of signed distances
    P    = P(dist<geps, :);                             %   Keep only d<0 (geps) points
    N    = size(P, 1);                                  %   Number of vertices N
    
    Pold = inf;                                         %   For the first iteration
    count = 0;
    
    while 1        
        %   Retriangulation by the Delaunay algorithm
        if max(sqrt(sum((P-Pold).^2,2))/h0)>ttol        %   Any large movement?
            Pold    = P;                                %   Save current positions
            dt      = DelaunayTri(P);                   %   2D Delaunay triangulation
            t       = dt.Triangulation;                 %   List of triangles
            ic      = incenters(dt);                    %   Centroids of triangles
            dist    = dcircle(ic, R);                   %   A column of signed distances
            t       = t(dist<-geps, :);                 %   Keep interior trianges
            bars    = edges(TriRep(t, P));              %   All sorted mesh edges [M, 2]
            M       = size(bars, 1);                    %   Length of the array of edges       
        end
        count   = count + 1;

        %   Move mesh points based on edge lengths L1 and forces F
        barvec  = P(bars(:, 1),:)-P(bars(:, 2),:);            %   List of edge vectors [M, 2]
        L1       = sqrt(sum(barvec.^2,2));                    %   Edge lengths [M, 1]
        barsc   = 0.5*(P(bars(:, 1), :) + P(bars(:, 2), :));  %   Edge centers [M, 2]
        hbars   = hcircle(barsc, R, par);                     %   Element sizes for edge centers [M, 1]
        hbars   = hbars/sqrt(sum(hbars.^2)/M);                %   Element sizes normalized by its average [M, 1] 
        L0      = Fscale*sqrt(sum(L1.^2)/M)*hbars;            %   Desired (non-uniform) edge lengths [M, 1] 
        F       = max(L0-L1, 0);                              %   Edge forces [M, 1]
        Fvec    = repmat(F./L1, 1, 2).*barvec;                %   Edge forces (x,y components) [M, 2]

        Ftot = full(sparse(bars(:, [1,1,2,2]), ones(size(F))*[1,2,1,2], [Fvec, -Fvec], N, 2));
        P = P + deltat*Ftot;                                %   Update node positions

        %   Bring outside vertices back to the boundary
        dist = dcircle(P, R);                           %   Find distances
        ind   = dist>0;                                 %   Find vertices outside (d>0)
        dgradx = (dcircle([P(ind, 1)+deps, P(ind, 2)], R)-dist(ind))/deps;   %   Numerical
        dgrady = (dcircle([P(ind, 1), P(ind, 2)+deps], R)-dist(ind))/deps;   %   gradient
        dgrad2 = dgradx.^2 + dgrady.^2;
        P(ind, :) = P(ind, :) - ...
            [dist(ind).*dgradx./dgrad2, dist(ind).*dgrady./dgrad2];     %   Projection

        %   Termination criterion: All interior nodes move less than dptol (scaled)
        factor = max(sqrt(sum(deltat*Ftot(dist<-geps,:).^2, 2))/h0); 
        if factor<dptol, break; end
        if count>2^13 break; end;
    end
    P(:, 3) = 0;
end
function h = hcircle(P, R, par)
    %   Compute element size function for a circle with radius R 
    %   0<par<1
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson (PhD thesis, pp. 38-39)
    h = 1/par - (1/par-1)*sqrt(sum(P.^2, 2))/R;     %   Relative edge length
end
function d = dcircle(P, R)
    %   Compute signed distance function for a circle with radius R 
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson (PhD thesis, pp. 38-39)

    d = sqrt(P(:, 1).^2 + P(:, 2).^2) - R;
end


function [] = box(L, W, H, XC, YC, ZC, Color, Transparency, LineWidth, EdgeColor)
%%  Draws a rectangular object in 3D (a plane, a brick, or a line)
%   Usage: box(2, 2, 2, 0, 0, 0, [1 1 0], 0, 1, [ 0 0 0]);
%   SNM Summer 2012
%   __________L__________
%   |         |y        |
%   |W        |         |
%   |         *(XC,YC,ZC)   
%   |         |         |
%   |_________|_________|x
%       

    if (L>0) & (W>0)
        hr = patch([-L/2 -L/2 +L/2 +L/2]+XC, [-W/2 +W/2 +W/2 -W/2]+YC, [-H/2 -H/2 -H/2 -H/2]+ZC, ...
            Color, 'FaceAlpha', Transparency, 'LineWidth', LineWidth, 'EdgeColor', EdgeColor);   %   bottom
        hr = patch([-L/2 -L/2 +L/2 +L/2]+XC, [-W/2 +W/2 +W/2 -W/2]+YC, [+H/2 +H/2 +H/2 +H/2]+ZC, ...
            Color, 'FaceAlpha', Transparency, 'LineWidth', LineWidth, 'EdgeColor', EdgeColor);   %   top
    end
    if (W>0) & (H>0)    
        hr = patch([+L/2 +L/2 +L/2 +L/2]+XC, [-W/2 -W/2 +W/2 +W/2]+YC, [-H/2 +H/2 +H/2 -H/2]+ZC, ...
            Color, 'FaceAlpha', Transparency, 'LineWidth', LineWidth, 'EdgeColor', EdgeColor);   %   right   
        hr = patch([-L/2 -L/2 -L/2 -L/2]+XC, [-W/2 -W/2 +W/2 +W/2]+YC, [-H/2 +H/2 +H/2 -H/2]+ZC, ...
            Color, 'FaceAlpha', Transparency, 'LineWidth', LineWidth, 'EdgeColor', EdgeColor);   %   left
    end
    if (L>0) & (H>0)
        hr = patch([-L/2 -L/2 +L/2 +L/2]+XC, [+W/2 +W/2 +W/2 +W/2]+YC, [-H/2 +H/2 +H/2 -H/2]+ZC, ...
            Color, 'FaceAlpha', Transparency, 'LineWidth', LineWidth, 'EdgeColor', EdgeColor);   %   back   
        hr = patch([-L/2 -L/2 +L/2 +L/2]+XC, [-W/2 -W/2 -W/2 -W/2]+YC, [-H/2 +H/2 +H/2 -H/2]+ZC, ...
            Color, 'FaceAlpha', Transparency, 'LineWidth', LineWidth, 'EdgeColor', EdgeColor);   %   front
    end
    if (W==0)&(H==0)    %  x-line
        hrl = line([-L/2+XC L/2+XC], [YC YC], [ZC ZC], 'LineWidth', 3, 'Color', Color);  
    end
    if (L==0)&(H==0)    %  y-line
        hrl = line([XC XC], [-W/2+YC W/2+YC], [ZC ZC], 'LineWidth', 3, 'Color', Color);  
    end
    if (L==0)&(W==0)    %  z-line
        hrl = line([XC XC], [YC YC], [-H/2+ZC H/2+ZC], 'LineWidth', 3, 'Color', Color);  
    end
    h = 1; 

end
 
function [P, t] = combine(varargin)
    %%   Combine multiple meshes into one mesh 
    %   Usage: [P, t] = combine(P1, P2, P3, t1, t2, t3);
    %   SNM Summer 2012    
    MN = size(varargin, 2); %   total number of arguments
    
    %   Combine arrays of vertices
    Ptemp = [];
    for m = 1:MN/2
        Ptemp = [Ptemp varargin{m}'];
    end
    P = Ptemp';
   
    %   Combine arrays of triangles
    temp = 0;
    ttemp = [];
    for m = MN/2+1:MN        
        ttemp = [ttemp (varargin{m}'+temp)];
        temp = temp + size(varargin{m-MN/2}, 1);
    end
    t = ttemp';    
end

function [fv] = cone(X, Y, Z, E1, E2, E3, M, R, H, Enormalize);
%   My own cone plot
%   SNM Summer 2012
%   Test with
    %     M = 24;             %   number of divisions
    %     R = 0.2;            %   radius, m  
    %     H = 1;              %   height ,m 
    %     X = [0 1]';
    %     Y = [0 1]';
    %     Z = [0 0];
    %     E1 = [1 1]';
    %     E2 = [0.3 0.5]';
    %     E3 = [1 0]';
    %%  Cone
    %   Create triangular mesh for the bottom
    x = R*[cos(2*pi*[0:M]/(M+1)) 0];
    y = R*[sin(2*pi*[0:M]/(M+1)) 0];
    P1(:, 1) = x'; 
    P1(:, 2) = y';
    P1(:, 3) = 0;
    t1(:, 1) = [1:M+1  ]';
    t1(:, 2) = [2:M+1 1]';
    t1(:, 3) = (M+2);
    %   Create triangular mesh for the top
    P2 = P1;
    P2(end, :) = [0 0 H];
    t2 = t1 + M + 2;

    %%  Combine into one mesh
    P = [P1' P2']'; 
    P(:, 3) = P(:, 3) - (max(P(:, 3)) + min(P(:, 3)))/2;   %   Observation point is at the center
    P = P/(max(P(:, 3)) - min(P(:, 3)));                   %   The arrow has the total height of 1
    t = [t1' t2']';

    %%  Prepare the field data
    temp = sqrt(E1.*E1 + E2.*E2 + E3.*E3);
    NX   = E1./temp;        %   unit vector
    NY   = E2./temp;        %   unit vector
    NZ   = E3./temp;        %   unit vector
    % sizearrow = temp/max(temp);     %   normalized size of the arrow
    sizearrow = temp/Enormalize;      %   absolute size of the arrow
    unitsize = max([abs(X(2)-X(1)) abs(Y(2)-Y(1)) abs(Z(2)-Z(1))]);     %   unit cell size
    
    %%  Clone the mesh, rotate copiles, and move copies    
    tconj = [];
    Pconj = [];
    N     = size(X, 1);
    K0    = repmat([0 0 0]', [1 size(P, 1)])';
    for m = 1:N
        PP = P; 
        %   scaling
        PP = PP*sizearrow(m)*unitsize;
        %   rotation (Rodrigues' rotation formula)
        theta = acos(NZ(m));
        k     = [-NY(m) +NX(m) 0]';
        K     = K0;
        if dot(k, k)>1e-6
            K     = repmat(k, [1 size(PP, 1)])'/sqrt(dot(k, k));
        end
        PP = PP*cos(theta) + cross(K, PP, 2)*sin(theta) + K.*repmat(dot(K, PP, 2), [1 3])*(1-cos(theta));
        %   movement
        PP(:, 1) = PP(:, 1) + X(m);
        PP(:, 2) = PP(:, 2) + Y(m);
        PP(:, 3) = PP(:, 3) + Z(m);        
        %   combining
        Pconj = [Pconj PP'];
        tconj = [tconj t'+length(PP)*(m-1)];
    end
    P = Pconj';
    t = tconj';
    fv.vertices = P;
    fv.faces = t;
end
 
function [] = contourca(strop, Potential)
    %   map of 16 colors
    N = 16;
    map =   [...
             0         0    0.7500;...
             0         0    1.0000;...
             0    0.2500    1.0000;...
             0    0.5000    1.0000;...
             0    0.7500    1.0000;...
             0    1.0000    1.0000;...
        0.2500    1.0000    0.7500;...
        0.5000    1.0000    0.5000;...
        0.7500    1.0000    0.2500;...
        1.0000    1.0000         0;...
        1.0000    0.7500         0;...
        1.0000    0.5000         0;...
        1.0000    0.2500         0;...
        1.0000         0         0;...
        0.7500         0         0;...
        0.5000         0         0];...

    scale     = max(Potential) - min(Potential);
    levels    = linspace(min(Potential)+scale/(N+1), max(Potential)-scale/(N+1), size(map, 1));
    Potential = reshape(Potential, strop.divisionsx+1, strop.divisionsy+1)';

    if strcmp(strop.planetype, 'xy') 
        x = strop.Points(1:strop.divisionsx+1,     1);
        y = strop.Points(1:strop.divisionsx+1:end, 2);
        z =  strop.planez;
        [D, color] = contource(x, y, Potential, levels, map); 
        for m =1:length(D)
            D{m}(3, :) =  strop.planez;
        end
    end
    if strcmp(strop.planetype, 'xz') 
        x = strop.Points(1:strop.divisionsx+1,     1);
        z = strop.Points(1:strop.divisionsx+1:end, 3);
        y =  strop.planey;
        [D, color] = contource(x, z, Potential, levels, map); 
        for m =1:length(D)
            D{m}(3, :) = D{m}(2, :);
            D{m}(2, :) =  strop.planey;
        end
    end
    if strcmp(strop.planetype, 'yz') 
        x = strop.planex;
        y = strop.Points(1:strop.divisionsx+1,     2);
        z = strop.Points(1:strop.divisionsx+1:end, 3);
        [D, color] = contource(y, z, Potential, levels, map); 
        for m =1:length(D)
            D{m}(3, :) = D{m}(2, :);
            D{m}(2, :) = D{m}(1, :);
            D{m}(1, :) =  strop.planex;
        end
    end

    if scale>1e-6
        for m = 1:length(D)
            line(D{m}(1, :), D{m}(2, :),  D{m}(3, :), 'LineWidth', 1.5, 'Color', color(m, :));
        end
        xlabel('x'); ylabel('y'); zlabel('z'); grid on; axis equal; view(46, 41);
    else
        disp('Electric potential is nearly constant in this plane - contour plot should not be used');       
        xlabel('x'); ylabel('y'); zlabel('z'); grid on; axis equal; view(46, 41);
    end
        
end

function [D, color] = contource(x, y, Potential, levels, map)
    %   Extracts curve data from the contourc matrix into cell array D
    %   SNM Summer 2012

    C = contourc(x, y, Potential, levels);
    %   Read C- matrix into a cell array D
    count = 1; count1 =1; S = size(C, 2);
    while 1
        number      = C(2, count);
        level       = C(1, count);
        color(count1, :) = map(levels==level, :);
        D{count1}   = C(:, count+1:count+number);
        if count+number==S break; end;
        count       = count + number + 1;
        count1      = count1 + 1;
    end
 
end

function E = efield(Points, c, P, t, centers, areas, normals, sizes, R, eps, msg)
%   Electric field from free/polarization charges
%   Vectorized for an arbitrary number of observation points
%   SNM Summer 2012
    if msg >0
        string = strcat('Computing current densities in a plane for domain #', num2str(msg));
        h    = waitbar(0, string);
    end
    if msg == 0
        h    = waitbar(0, 'Please wait - computing current density at the obs. point');
    end
    
%%  Find the field at a number of points    
    E       = zeros(size(Points));      %   standard format 
    %   contribution of triangle m to all points is sought
    for m =1:size(centers, 1)
        temp        = repmat(centers(m, :), size(Points, 1), 1) - Points;   %   these are distances to the observation point
        DIST        = sqrt(dot(temp, temp, 2));                            %   single column
        index   = find(DIST./(sizes(m)^2)<=R+1e-16);   % index into Points
        r1      = P(t(m, 1), :);    %   row
        r2      = P(t(m, 2), :);    %   row
        r3      = P(t(m, 3), :);    %   row   
        I           = areas(m)*temp./repmat(DIST.^3, 1, 3);              %   integral, standard format
        I(index, :) = potint2(r1, r2, r3, normals(m, :), Points(index, :));%   I was calculated with the area 
        E = E + (- c(m)*I/(4*pi*eps));
        waitbar(m/size(centers, 1));
    end
    close(h);
end

function [fv] = kreuz(X, Y, Z, M, R, H);
%   Our own cross plot
%   SNM/Yana Summer 2012
P            = [ 0.05 0.05 -0.05;
                 0.05 0.05 0.05;
                 0.5  0.05 -0.05;
                 0.5  0.05  0.05;
                 0.5  -0.05 -0.05;
                 0.5  -0.05  0.05
                 0.05 -0.05  0.05
                 -0.05 0.05 -0.05;
                 -0.05 0.05 0.05;
                 -0.5  0.05 -0.05;
                 -0.5  0.05  0.05;
                 -0.5  -0.05 -0.05;
                 -0.5  -0.05  0.05
                 -0.05 -0.05  0.05
                 -0.05  0.5  0.05
                  0.05  0.5  0.05
                 -0.05 0.5 -0.05
                  0.05 0.5 -0.05
                  -0.05 0.05 0.5
                  0.05 0.05 0.5
                  0.05 -0.05 0.5
                  -0.05 -0.05 0.5
                  -0.05 0.05 -0.5
                  0.05 0.05 -0.5
                  0.05 -0.05 -0.5
                  -0.05 -0.05 -0.5
                  -0.05 -0.5  0.05
                  0.05 -0.5  0.05
                 -0.05 -0.5 -0.05
                  0.05 -0.5 -0.05
                  0.05 -0.05  -0.05
                  -0.05 -0.05  -0.05];
t              = [ 1 2 4 3                  %   rectangles here
                 4 3 5 6
                 2 4 6 7 
                 8 10 11 9
                 10 11 13 12
                 11 9 14 13
                 9 2 16 15
                 15 17 18 16
                 2 16 18 1
                 9 19 20 2
                 19 20 21 22
                 23 24 25 26
                 27 28 30 29
                 1 3 5 31
                 7 31 5 6
                 10 8 32 12
                 13 12 32 14
                 15 17 8 9
                 17 18 1 8
                 20 2 7 21
                 9 14 22 19
                 22 21 7 14
                 7 31 30 28
                 32 14 27 29
                 27 28 7 14
                 29 30 31 32
                 1 31 25 24
                 8 32 26 23
                 8 1 24 23
                 31 32 26 25];
                          
    P = P*H;
    P(:, 1) = P(:, 1) + X;
    P(:, 2) = P(:, 2) + Y;
    P(:, 3) = P(:, 3) + Z;        

    fv.vertices = P;
    fv.faces = t;
    
end

function [Pnew] = move(P, x, y, z)
    %%   Move the center of the mesh by x y z
    %   SNM Summer 2012
    Pnew(:, 1) = +P(:, 1) + x;
    Pnew(:, 2) = +P(:, 2) + y;
    Pnew(:, 3) = +P(:, 3) + z;
end

function [centers, normals] = normcenters(P, t);
    %%   Find centers and outer normal vectors for convex bodies
    tr              = TriRep(t(:, 1:3), P);
    centers         = incenters(tr);
    normals         = faceNormals(tr);  %   normals from Tri class
    temp            = sign(dot(normals, centers, 2));
    temp(temp==0)   = 1;
    normals     = normals.*repmat(temp, 1, 3); 
end

function Potential = potential(Points, c, P, t, centers, areas, normals, sizes, R, eps, msg)
%   Electric field from free/polarization charges
%   Vectorized for an arbitrary number of observation points
%   SNM Summer 2012
    if msg ==1
        h    = waitbar(0,  'Please wait - computing electric potentials in a plane');
    else
        h    = waitbar(0, 'Please wait - computing electric potential at the obs. point');
    end
%%  Find the potential at a number of points    
    Potential       = zeros(size(Points, 1), 1);                                %   Standard format 
    %   contribution of triangle m to all points is sought
    for m =1:size(centers, 1)
        temp        = repmat(centers(m, :), size(Points, 1), 1) - Points;       %   These are distances to the observation point
        DIST        = sqrt(dot(temp, temp, 2));                                 %   Single column
        index   = find(DIST./(sizes(m)^2)<=R+1e-16);   % index into Points
        r1      = P(t(m, 1), :);    %   row
        r2      = P(t(m, 2), :);    %   row
        r3      = P(t(m, 3), :);    %   row   
        I           = areas(m)./DIST;                                           %   Integral, standard format
        [I(index), IRho] = potint(r1, r2, r3, normals(m, :), Points(index, :)); %   I was calculated with the area 
        Potential = Potential + c(m)*I/(4*pi*eps);
        waitbar(m/size(centers, 1));
    end
    close(h);
end

function [I, IRho] = potint(r1, r2, r3, norm, ObsPoint)
%   Potential integrals 1/r and vec(r)/r for a single triangle
%   The integrals are not divided by the area
%   Vectorized for an arbitrary number of observation points
%   Comment: new version; no alpha; uses the last formula in Eq. (5) of Wilton et al. 1984, p. 279
%   SNM June 2012

    N = size(ObsPoint, 1);
    
    Z = zeros(N, 9);
    S = zeros(N, 3);

    %   Create r+ and r- coordinates
    temp = [r2 r1 r3 r1 r3 r2];
    r    = repmat(temp, N, 1);

    %   Create l coordinates
    normabsl1 = sqrt(sum((r2-r1).^2));
    normabsl2 = sqrt(sum((r3-r1).^2));
    normabsl3 = sqrt(sum((r3-r2).^2));
    temp0 = [(r2-r1)./normabsl1 (r3-r1)./normabsl2 (r3-r2)./normabsl3];
    l    = repmat(temp0, N, 1);

    %   Create unit normal to the edges of the triangle
    u(1:3) = +cross((r2-r1)./normabsl1, norm);
    u(4:6) = -cross((r3-r1)./normabsl2, norm);
    u(7:9) = +cross((r3-r2)./normabsl3, norm);
    u      = repmat(u, N, 1);

    %   Create projection vector of the observation point
    NORM = repmat(norm, N, 1);
    ndot = sum(ObsPoint.*NORM, 2);
    p    = ObsPoint - repmat(ndot, 1, 3).*NORM;
    
    % Calculation of the analytic formula
    count = 0;
    for c1 = 0:2    
        %   Distance of observation point perpendicular to the plane with triangle
        temporary   = ObsPoint - r(:, [1:3] + 3*(count + 1));
        distanceobs = abs(sum(NORM.*temporary, 2));     

        %   Calculate p+ and p-
        d1 = sum(NORM.*r(:, [1:3] + 3*count), 2);
        d2 = sum(NORM.*r(:, [1:3] + 3*(count + 1)), 2);

        pplus  = r(:, [1:3] + 3*count)       - NORM.*repmat(d1, 1, 3);
        pminus = r(:, [1:3] + 3*(count + 1)) - NORM.*repmat(d2, 1, 3);

        %   Calculate l+ and l-
        lplus  = sum(l(:, [1:3] + 3*c1).*(pplus - p), 2);
        lminus = sum(l(:, [1:3] + 3*c1).*(pminus - p), 2);

        %   Perpendicular distance from projection vector to edge
        P0 = abs( sum(u(:, [1:3] + 3*c1).*(pminus-p), 2) ); 

        %   Distances to l+ and l- from projection vector
        PPLUS  = sqrt(P0.*P0 + lplus.*lplus);
        PMINUS = sqrt(P0.*P0 + lminus.*lminus);

        %   Vector containing line P0 measured
        PHAT = (pminus - p - repmat(lminus, 1, 3).*l(:, [1:3] + 3*c1))./repmat(P0, 1, 3);       

        %   Distances to l+ and l- from observation point
        RPLUS  = sqrt(PPLUS.*PPLUS + distanceobs.*distanceobs);
        RMINUS = sqrt(PMINUS.*PMINUS + distanceobs.*distanceobs);
        R0	   = sqrt(P0.*P0 + distanceobs.*distanceobs);

        %   A value of one term of the analytic sum 1/R
        d1      = P0.*log((RPLUS+lplus)./(RMINUS+lminus));        
        d2      = distanceobs.*(atan((P0.*lplus)./(R0.*R0 + distanceobs.*RPLUS)));        
        d3      = distanceobs.*(atan((P0.*lminus)./(R0.*R0 + distanceobs.*RMINUS)));
        
        dotPHATu= sum(PHAT.*u(:, [1:3] + 3*c1), 2);
        S(:, c1+1)   = (d1 - d2 + d3).*dotPHATu;

        %   A value of one term of the analytic sum RHO /R
        d1 = R0.*R0.*log((RPLUS + lplus)./(RMINUS + lminus));
        d2 = lplus.*RPLUS - lminus.*RMINUS;
        d3 = d1 + d2;
        Z(:, [1:3] + c1*3) = repmat(d3, 1, 3).*u(:, [1:3] + c1*3);        
        count = count + 2;
    end
    
    %   Contribution is zero when the projection point is on the edge (Wilton et al. 1984, p. 279)
    S(isnan(S)) = 0;  S(isinf(S)) = 0;
    I       = sum(S, 2);    
    IRho    = 0.5*(Z(:, 1:3) + Z(:, 4:6) + Z(:, 7:9));   
end

function [Int] = potint2(r1, r2, r3, normal, ObsPoint)
%   Potential integrals grad(1/r) for a single triangle
%   The integrals are not divided by the area
%   Vectorized for an arbitrary number of observation points
%   SNM June 2012
%   Test (comparison with Wang et al, 2003):
% %     clear all;
% %     format long;
% %     r1 = [62.5 25.0 0];
% %     r2 = [62.5 25.0 2];
% %     r3 = [62.5 37.5 0];
% %     ObsPoint = [62.5 0.0 0.0];
% %     tempv           = cross(r2-r1, r3-r1);  %   correct normal sign!
% %     temps           = sqrt(tempv(1)^2 + tempv(2)^2 + tempv(3)^2);
% %     normal          = tempv/temps;
% %     Area            = temps/2;   
% %     %   Analytical integral
% %      [coeff, weights, IndexF] = tri(25, 10);
% %      Int = [0 0 0];
% %      for m = 1:length(coeff)     
% %          Point = coeff(1, m)*r1 + coeff(2, m)*r2 + coeff(3, m)*r3;
% %          R = ObsPoint - Point;
% %          Int = Int - Area*weights(m)*R/(sqrt(sum(R.*R)))^3;
% %      end
% %      Int
% %      Int = potint2(r1, r2, r3, normal, ObsPoint)
 
    N           = size(ObsPoint, 1);
    I           = zeros(N, 3);
    S           = zeros(N, 9);
    Int         = zeros(N, 3);
    BetaTerm    = zeros(N, 3);

    %   Create r+ and r- coordinates
    temp = [r2 r1 r3 r1 r3 r2];
    r    = repmat(temp, N, 1);

    %   Create l coordinates
    normabsl1 = sqrt(sum((r2-r1).^2));
    normabsl2 = sqrt(sum((r3-r1).^2));
    normabsl3 = sqrt(sum((r3-r2).^2));
    temp = [(r2-r1)./normabsl1 (r3-r1)./normabsl2 (r3-r2)./normabsl3];
    l    = repmat(temp, N, 1);

    %   Create unit normal to the edges of the triangle
    u(1:3) = +cross((r2-r1)./normabsl1, normal);
    u(4:6) = -cross((r3-r1)./normabsl2, normal);
    u(7:9) = +cross((r3-r2)./normabsl3, normal);
    u      = repmat(u, N, 1);

    %   Create projection vector of the observation point
    NORM = repmat(normal, N, 1);
    ndot = sum(ObsPoint.*NORM, 2);
    p    = ObsPoint - repmat(ndot, 1, 3).*NORM;

    %   Is the projection point on the triangle edge or its continuation?
    %   Calculate midpoints of the edges of a triangle
    temp     = 0.5*[r1+r2 r1+r3 r2+r3];
    midpoint = repmat(temp, N, 1);
    %   Calculate vectors from observation point project to midpoints
    vector = midpoint - [p p p];
      %   Move the observation  point (projection) from the edge
    factor = 1e-6;
    factor1 = factor*normabsl1;
    factor2 = factor*normabsl2;
    factor3 = factor*normabsl3;
    check1 = sum(vector(:, 1:3).*u(:, 1:3), 2);
    check2 = sum(vector(:, 4:6).*u(:, 4:6), 2);
    check3 = sum(vector(:, 7:9).*u(:, 7:9), 2);	
    index  = (abs(check1)<factor1); ObsPoint(index, :) =  ObsPoint(index, :) - factor1*u(index, 1:3);    
    index  = (abs(check2)<factor2); ObsPoint(index, :) =  ObsPoint(index, :) - factor2*u(index, 4:6);
    index  = (abs(check3)<factor3); ObsPoint(index, :) =  ObsPoint(index, :) - factor3*u(index, 7:9);
    ndot   = sum(ObsPoint.*NORM, 2);
    p      = ObsPoint - repmat(ndot, 1, 3).*NORM;
    
    % Calculation of the analytic formula
    count = 0;
    for c1 = 0:2    
        %   Distance of observation point perpendicular to the plane with triangle
        temporary   = ObsPoint - r(:, [1:3] + 3*(count + 1));
         %%  Changes compared to potint.m   
        distanceobs = sum(NORM.*temporary, 2);     
         %%  End of changes compared to potint.m   
         
        %   Calculate p+ and p-
        d1 = sum(NORM.*r(:, [1:3] + 3*count), 2);
        d2 = sum(NORM.*r(:, [1:3] + 3*(count + 1)), 2);

        pplus  = r(:, [1:3] + 3*count)       - NORM.*repmat(d1, 1, 3);
        pminus = r(:, [1:3] + 3*(count + 1)) - NORM.*repmat(d2, 1, 3);

        %   Calculate l+ and l-
        lplus  = sum(l(:, [1:3] + 3*c1).*(pplus - p), 2);
        lminus = sum(l(:, [1:3] + 3*c1).*(pminus - p), 2);

        %   Perpendicular distance from projection vector to edge
        P0 = abs( sum(u(:, [1:3] + 3*c1).*(pminus-p), 2) );

        %   Distances to l+ and l- from projection vector
        PPLUS  = sqrt(P0.*P0 + lplus.*lplus);
        PMINUS = sqrt(P0.*P0 + lminus.*lminus);

        %   Vector containing line P0 measured
        PHAT = (pminus - p - repmat(lminus, 1, 3).*l(:, [1:3] + 3*c1))./repmat(P0, 1, 3);

        %   Distances to l+ and l- from observation point
        RPLUS  = sqrt(PPLUS.*PPLUS + distanceobs.*distanceobs);
        RMINUS = sqrt(PMINUS.*PMINUS + distanceobs.*distanceobs);
        R0	   = sqrt(P0.*P0 + distanceobs.*distanceobs);
        
        %%  Changes compared to potint.m
        %   A value of one term of the analytic sum 1/R         
        d1      = atan(P0.*lplus./(R0.*R0 + abs(distanceobs).*RPLUS));
        d2      = atan(P0.*lminus./(R0.*R0 + abs(distanceobs).*RMINUS));
        d3      = log((RPLUS + lplus)./(RMINUS + lminus));
        
        %% End of changes compared to potint.m
        
        dotPHATu= sum(PHAT.*u(:, [1:3] + 3*c1), 2);        
        
        %%  Changes compared to potint.m        
        BetaTerm(:, c1+1) = dotPHATu.*(d1 - d2);
        S(:, 1 + 3*c1) = d3.*u(:, 1 + 3*c1);
        S(:, 2 + 3*c1) = d3.*u(:, 2 + 3*c1);
        S(:, 3 + 3*c1) = d3.*u(:, 3 + 3*c1);
        %% End of changes compared to potint.m
        count = count + 2;      
    end
    
    
    %%  Changes compared to potint.m         
    Beta = sum(BetaTerm, 2);
    
    I(:, 1) = S(:, 1) + S(:, 4) + S(:, 7);
    I(:, 2) = S(:, 2) + S(:, 5) + S(:, 8);
    I(:, 3) = S(:, 3) + S(:, 6) + S(:, 9);
    
    %   find sign of distanceobs
    Sign = sign(distanceobs); 
    
    %   value of integral for 1/R
    Int(:, 1) = -normal(1)*Sign.*Beta - I(:, 1);
    Int(:, 2) = -normal(2)*Sign.*Beta - I(:, 2);
    Int(:, 3) = -normal(3)*Sign.*Beta - I(:, 3);
    
    %   Contribution is zero when the projection point is on the edge (Wilton et al. 1984, p. 279)    
    Int(isnan(Int)) = 0;  Int(isinf(Int)) = 0;    
    %% End of changes compared to potint.m
end

function q = simpqual(P, t)
%   Facet quality - radius ratio 
%   Original code: DISTMESH 2004-2012 Per-Olof Persson

    a = sqrt(sum((P(t(:, 2), :) - P(t(:, 1), :)).^2,2));
    b = sqrt(sum((P(t(:, 3), :) - P(t(:, 1), :)).^2,2));
    c = sqrt(sum((P(t(:, 3), :) - P(t(:, 2), :)).^2,2));
    r = 1/2*sqrt((b+c-a).*(c+a-b).*(a+b-c)./(a+b+c));
    R = a.*b.*c./sqrt((a+b+c).*(b+c-a).*(c+a-b).*(a+b-c));
    q = 2*r./R;
    
end

function [coeff, weights, IndexF] = tri(arg1, arg2);
    % Creates integration points and weights for triangles 
    % Syntax:
    % [coeff, weights, IndexF] = tri(3) -   will do barycentric subdivision with
    %                                       3*3 = 9 integration points
    % [coeff, weights, IndexF] = tri(7,5) - will use Gaussian quadrature of fifth 
    %                                       order with seven integration points
    %
    % Description:
    %
    % Uses Gaussian quadratures (with author's permission) from 
    % http://www.cs.kuleuven.ac.be/~nines/research/ecf/ecf.html
    % Formulas of third (4 integration points), fifth (7 integration points)
    % seventh (13 integration points), and tenth (25 integration points) order
    % of accuracy may be created.
    %
    % Uses the "edge" method for barycentric subdivision of arbitrary order,
    % where the edges of smaller triangles (similar to the original one) are
    % equally subdivided. This gives the desired barycentric points.
    %
    % Antenna Lab, ECE, WPI 2002
    nargs = nargin; 

    if  nargs < 1
        error('Requires at least 1 input');
    end

    if nargs == 1
    %   Barycentric triangle subdivision - coefficients for vertexes only
    %   M - subdivision order (number of subtriangles is M*M)
        M = arg1;	
        if M < 2 
            coeff(:,1) = [1/3 1/3 1/3]';
            weights(1) =   1;
        end

        coeff   = zeros(3, M*M);
        k       = 1;
        scale   = 1;
        eps     = 2 + 1e-9;   % scaling

        if mod(M,3) == 0
            N   =   M/3*2;
        elseif mod(M,3) == 2
            N   =   M/3*2 - 1/3;
        else
            N   =   (M-1)/3*2;
        end    

        %   Border loop - starts with the outer border of integration points
        %   and then goes inside - "triangle" by "triangle" 
        for m   = 1 : N
            div         = M - m - floor(m/eps);          %   integer - edge is divided into div (XX-jump-XX-jump)
            scale       = div/M;                         %   real - relative
            alpha       = (1 + 2*scale)/3;
            beta        = (1 - scale)/3;    
            coeff1      = [alpha beta beta]';    %   p1 new 
            coeff2      = [beta alpha beta]';    %   p2 new
            coeff3      = [beta beta alpha]';    %   p3 new
            %   first edge
            for n = 1 : div
                vector      = coeff1*(div - n + 1)/div + coeff2*(n - 1)/div;
                coeff(:,k)  = vector;
                k = k + 1;
            end
            %   second edge
            for n = 1 : div
                vector      = coeff2*(div - n + 1)/div + coeff3*(n - 1)/div;
                coeff(:,k)  = vector;
                k = k + 1;
            end
            %   third edge
            for n = 1 : div
                vector      = coeff3*(div - n + 1)/div + coeff1*(n - 1)/div;
                coeff(:,k)  = vector;
                k = k + 1;
            end  
        end

        %   Center point 
        if (3*floor(M/3) ~= M)
            coeff(:,k) = [1/3 1/3 1/3]';
        end	
        weights = 1/size(coeff,2)*ones(1,size(coeff,2)); 
    end

    if nargs == 2
    %   Gaussian quadrature formulae
    %   arg - number of integration points 
        if arg1 == 1    % first order (center)               
            coeff(:,1)  = [1/3 1/3 1/3]';
            weights(1)  =   1;
        end

        if (arg1 == 3) & (arg2 == 2)    %   second order (sides)                
            coeff(:,1)  = [1/2 1/2 0]';
            coeff(:,2)  = [0 1/2 1/2]';
            coeff(:,3)  = [1/2 0 1/2]';
            weights(1)  =   1/3;
            weights(2)  =   1/3;
            weights(3)  =   1/3;
        end
        if (arg1 == 4) & (arg2 == 3)    % third order                
            a1          =   0.6;
            b1          =   0.2;       
            coeff(:,1)  = [1/3 1/3 1/3]';
            coeff(:,2)  = [a1  b1  b1]';
            coeff(:,3)  = [b1  a1  b1]';
            coeff(:,4)  = [b1  b1  a1]';       
            weights(1)  =   -27/48;
            weights(2)  =    25/48;
            weights(3)  =    25/48;
            weights(4)  =    25/48;    	
        end
        if (arg1 == 6) & (arg2 == 3)   % third order (sides)     
            coeff(:,1)  = [1/2 1/2 0]';
            coeff(:,2)  = [0 1/2 1/2]';
            coeff(:,3)  = [1/2 0 1/2]';
            weights(1)  =   0.016666666;
            weights(2)  =   0.016666666;
            weights(3)  =   0.016666666;        
            a1          =   0.666666666;
            b1          =   0.166666666;       
            coeff(:,4)  = [a1  b1  b1]';
            coeff(:,5)  = [b1  a1  b1]';
            coeff(:,6)  = [b1  b1  a1]';       
            weights(4)  =    0.15;
            weights(5)  =    0.15;
            weights(6)  =    0.15;    	
            weights     = weights*2;
        end    

        if (arg1 == 7) & (arg2 == 5)    % fifth order                
            a1          =   0.797426985353087;
            b1          =   0.101286507323456;       
            a2          =   0.059715871789770;
            b2          =   0.470142064105115;	
            coeff(:,1)  = [1/3 1/3 1/3]';
            coeff(:,2)  = [a1  b1  b1]';
            coeff(:,3)  = [b1  a1  b1]';
            coeff(:,4)  = [b1  b1  a1]';       
            coeff(:,5)  = [a2  b2  b2]';
            coeff(:,6)  = [b2  a2  b2]';
            coeff(:,7)  = [b2  b2  a2]';    
            weights(1)    =   0.2250000;
            weights(2)    =   0.1259392;
            weights(3)    =   0.1259392;
            weights(4)    =   0.1259392;    
            weights(5)    =   0.1323942;
            weights(6)    =   0.1323942;
            weights(7)    =   0.1323942;        
       end
       if (arg1 == 9) & (arg2 == 5)    % fifth order (sides)
            coeff(:,1)  = [1 0 0]';
            coeff(:,2)  = [0 1 0]';
            coeff(:,3)  = [0 0 1]';
            coeff(:,4)  = [1/2 1/2 0]';
            coeff(:,5)  = [0 1/2 1/2]';
            coeff(:,6)  = [1/2 0 1/2]';
            a1          =   0.62283903060711;
            b1          =   0.18858048469644;       
            coeff(:,7)  = [a1  b1  b1]';
            coeff(:,8)  = [b1  a1  b1]';
            coeff(:,9)  = [b1  b1  a1]';
            weights(1)  =   0.01027006767296;
            weights(2)  =   0.01027006767296;
            weights(3)  =   0.01027006767296;        
            weights(4)  =   0.03098774943413;
            weights(5)  =   0.03098774943413;
            weights(6)  =   0.03098774943413;                          
            weights(7)  =   0.12540884955956;  
            weights(8)  =   0.12540884955956;
            weights(9)  =   0.12540884955956; 
            weights = weights*2;
       end
       if (arg1 == 13) & (arg2 == 7)    % seventh order
           a1 = 0.4793080678;
           b1 = 0.2603459660;
           a2 = 0.8697397941;
           b2 = 0.0651301029;
           a3 = 0.6384441885;			
           b3 = 0.3128654960;
           c3 = 0.0486903154;
           coeff(:,1)  = [1/3 1/3 1/3]';
           coeff(:,2)  = [a1  b1  b1]';
           coeff(:,3)  = [b1  a1  b1]';
           coeff(:,4)  = [b1  b1  a1]';       
           coeff(:,5)  = [a2  b2  b2]';
           coeff(:,6)  = [b2  a2  b2]';
           coeff(:,7)  = [b2  b2  a2]';    			       
           coeff(:,8) = [a3  b3  c3]'; 
           coeff(:,9) = [a3  c3  b3]'; 
           coeff(:,10) = [b3  a3  c3]'; 
           coeff(:,11) = [b3  c3  a3]'; 
           coeff(:,12) = [c3  a3  b3]'; 
           coeff(:,13) = [c3  b3  a3]'; 
           weights(1)  =    -0.14957004; 
           weights(2)  =	 0.1756152574; 
           weights(3)  =	 0.1756152574; 
           weights(4)  =	 0.1756152574; 
           weights(5)  =	 0.0533472356; 
           weights(6)  =	 0.0533472356; 
           weights(7)  =	 0.0533472356; 
           weights(8)  =	 0.0771137608; 
           weights(9)  =	 0.0771137608; 
           weights(10) =	 0.0771137608; 
           weights(11) =	 0.0771137608; 
           weights(12) =	 0.0771137608; 
           weights(13) =	 0.0771137608;       
       end
       if (arg1 == 25) & (arg2 == 10)    % tenth order       
            a1   =   0.1498275788;
            b1   =   0.4250862106;		
            a2   =   0.9533822650;
            b2   =   0.0233088675;		
            a3   =   0.6283074002; 
            b3   =   0.2237669736; 
            c3   =   0.1479256262; 		
            a4   =   0.6113138262; 
            b4   =   0.3587401419; 
            c4   =   0.0299460319;		
            a5   =   0.8210720699; 
            b5   =   0.1432953704; 
            c5   =   0.0356325597;    
            coeff(:,1)  = [1/3 1/3 1/3]';
            coeff(:,2)  = [a1  b1  b1]';
            coeff(:,3)  = [b1  a1  b1]';
            coeff(:,4)  = [b1  b1  a1]';       
            coeff(:,5)  = [a2  b2  b2]';
            coeff(:,6)  = [b2  a2  b2]';
            coeff(:,7)  = [b2  b2  a2]';    			       		
            coeff(:,8)  = [a3  b3  c3]'; 
            coeff(:,9)  = [a3  c3  b3]'; 
            coeff(:,10) = [b3  a3  c3]'; 
            coeff(:,11) = [b3  c3  a3]'; 
            coeff(:,12) = [c3  a3  b3]'; 
            coeff(:,13) = [c3  b3  a3]';        		
            coeff(:,14) = [a4  b4  c4]'; 
            coeff(:,15) = [a4  c4  b4]'; 
            coeff(:,16) = [b4  a4  c4]'; 
            coeff(:,17) = [b4  c4  a4]'; 
            coeff(:,18) = [c4  a4  b4]'; 
            coeff(:,19) = [c4  b4  a4]';        		
            coeff(:,20) = [a5  b5  c5]'; 
            coeff(:,21) = [a5  c5  b5]'; 
            coeff(:,22) = [b5  a5  c5]'; 
            coeff(:,23) = [b5  c5  a5]'; 
            coeff(:,24) = [c5  a5  b5]'; 
            coeff(:,25) = [c5  b5  a5]';        		
            weights(1) =     0.03994725237;
            weights(2) =     0.03556190112;
            weights(3) =     0.03556190112;
            weights(4) =     0.03556190112;        
            weights(5) =     0.00411190935;
            weights(6) =     0.00411190935;            
            weights(7) =     0.00411190935;                 
            weights(8) =     0.02271529614;
            weights(9) =     0.02271529614;
            weights(10) =    0.02271529614;                
            weights(11) =    0.02271529614;
            weights(12) =    0.02271529614;
            weights(13) =    0.02271529614;        
            weights(14) =    0.01867992812;
            weights(15) =    0.01867992812;
            weights(16) =    0.01867992812;                
            weights(17) =    0.01867992812;
            weights(18) =    0.01867992812;
            weights(19) =    0.01867992812;        
            weights(20) =    0.01544332844;
            weights(21) =    0.01544332844;
            weights(22) =    0.01544332844;
            weights(23) =    0.01544332844;
            weights(24) =    0.01544332844;
            weights(25) =    0.01544332844;
            weights     = weights*2;
        end
    end

    IndexF = size(coeff,2);
end
function [ ] = graphics_E51(io, strop, stroc, P, t, c, Area, strsc, strge, COLOR, TRANSPARENCY);
    %%   Program I/O graphics
    %    variable io is zero for the input graphics and one for the output graphics 
    warning off;
    %%   Draw the observation plane (always shown when included)
     if strcmp(strop.yes, 'yes')|strcmp(strop.potential, 'yes')
        if strcmp(strop.planetype, 'xy')  
            box(strop.planesizex, strop.planesizey, 0, strop.planex, strop.planey, strop.planez, 'y', 0.1, 0.5, 'k');
        end
        if strcmp(strop.planetype, 'xz')   
            box(strop.planesizex, 0, strop.planesizey, strop.planex, strop.planey, strop.planez, 'y', 0.1, 0.5, 'k');
        end
        if strcmp(strop.planetype, 'yz')   
           box(0, strop.planesizex, strop.planesizey,  strop.planex, strop.planey, strop.planez, 'y', 0.1, 0.5, 'k');
        end
    end
    
    %%  Draw the observation point(s)
    if strcmp(stroc.yes, 'yes')
        scaling = 0.05*max(strop.planesizex, strop.planesizey)*stroc.size;
        fv      = kreuz(stroc.x, stroc.y, stroc.z, 4, 0.1*scaling, 1*scaling);
        patch(fv, 'FaceColor', 'g', 'EdgeColor', 'k'); 
    end
    hold on
    %%   Input graphics (within the main GUI window)    
    if io == 0            
        Color = COLOR;
        Transparency = TRANSPARENCY;
        SCALE = strge.VoltageOfElectrodes(1:strge.NumberOfElectrodes);
        SCALE = (SCALE-min(SCALE))/(max(SCALE)-min(SCALE));
        for m = 1:strge.NumberOfElectrodes    
            %   Update color array
            index = find(t(:, 4)==m);            
            for n = 1:length(index)
                Color(index(n), :) = [SCALE(m) 0 1-SCALE(m)];
                Transparency(index(n)) = 1;
            end   
        end 
        %   Main plot
        patch('Faces', t(:, 1:3), 'Vertices', P, 'FaceColor', 'flat', 'FaceVertexCData', Color, 'EdgeColor', 'none', ...
            'FaceAlpha', 'flat', 'AlphaDataMapping', 'none', 'FaceVertexAlphaData', Transparency);       
        %   Plot edges
        for m = 1:4                    
            if strcmp(strge.edges{m}, 'yes');
                index = (strge.Tracker == m);
                patch('Faces', t(index, 1:3), 'Vertices', P, 'FaceColor', 'none', 'EdgeColor', 'y'); 
            end
        end
        xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); 
        title('Problem geometry (with electrodes)');      
    end
    %%   Output graphics 
    if io == 1
       %%  Prepare common variables
        tM = t(t(:, 4)>0, :);
        tD = t(t(:, 4)==0, :);
        cM = c(t(:, 4)>0, :);
        cD = c(t(:, 4)==0, :);

        Ptemp = P'; ttemp = t'; ind = size(t, 1);  

        ctempM = cM'; ttempM = tM'; 
        XM = reshape(Ptemp(1, ttempM(1:3, :)),[3, size(ttempM, 2)]);
        YM = reshape(Ptemp(2, ttempM(1:3, :)),[3, size(ttempM, 2)]);
        ZM = reshape(Ptemp(3, ttempM(1:3, :)),[3, size(ttempM, 2)]);  

        ctempD = cD'; ttempD = tD'; 
        XD = reshape(Ptemp(1, ttempD(1:3, :)),[3, size(ttempD, 2)]);
        YD = reshape(Ptemp(2, ttempD(1:3, :)),[3, size(ttempD, 2)]);
        ZD = reshape(Ptemp(3, ttempD(1:3, :)),[3, size(ttempD, 2)]);          
         %%  Charge scaling (scaled variable is "charge")
        cplus       = max(c);
        cminus      = min(c);
        charge      = c;
        charge(charge>=+cplus*strsc.positive)  = cplus*strsc.positive;
        charge(charge<=+cminus*strsc.negative) = cminus*strsc.negative;
        
        %%   Interpolate charge density for vertexes - global   
        ctemp = charge';
        cv    = zeros(1, length(Ptemp));
        for m = 1:size(Ptemp, 2)
            [p, q] = find(ttemp(1:3, :)-m==0);
            if (length(q))
                cv(m) = sum(ctemp(q), 2)/length(q);
            end    
        end        
        CD = cv(ttempD(1:3, :));
        CM = cv(ttempM(1:3, :));      
%      
%         if (~strcmp(strop.yes, 'yes'))&(~strcmp(strop.potential, 'yes'))   
%             if strge.Indicator==0
%                 patch(XD, YD, ZD, CD, 'FaceAlpha', 1.0, 'EdgeColor', 'k', 'FaceLighting', 'none');
%             else
%                 patch(XD, YD, ZD, CD, 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'FaceLighting', 'none');
%             end
%             patch(XM, YM, ZM, [0.5 0.5 0.5], 'FaceAlpha', 1.0, 'EdgeColor', 'k', 'FaceLighting', 'none');              
%             xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); colorbar;
%             title('Solution: Surface charge density in C/m^2')
%         end        
%         if strcmp(strop.yes, 'yes')|strcmp(strop.potential, 'yes')       
             patch(XD, YD, ZD, CD, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceLighting', 'none');            
%             %patch(XM, YM, ZM, CM, 'FaceAlpha', 1.0, 'EdgeColor', 'none', 'FaceLighting', 'none');
             patch('Faces', tM(:, 1:3), 'Vertices', P, 'FaceColor', [0.5 0.5 0.5]);            
%             xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); colorbar;
%             if strcmp(strop.yes, 'yes')
%                 title('Solution: Surface charge density in C/m^2 and electric current in the observation plane'); 
%             else
%                  title('Solution: Surface charge density in C/m^2 and electric potential in the observation plane'); 
%             end
%         end
        if strcmp(strop.yes, 'yes')      
            %%   Draw electric field/electric current (final) in the observation plane
            strop.Current = strop.Current + 1e-12*max(max(strop.Current))*rand(size(strop.Current)); %   randomize a bit            
            normalize = sqrt(sum(strop.Current.*strop.Current, 2));
            normalize = max(normalize);
            normalize = normalize/strop.arrow;                    
            fv = cone(strop.Points(:, 1), strop.Points(:, 2), strop.Points(:, 3), strop.Current(:, 1), strop.Current(:, 2), strop.Current(:, 3), 36, 0.25, 1, normalize);
            h = patch(fv, 'FaceColor', 'y', 'EdgeColor', 'none'); 
            camproj('orthographic');
            light('Position',[1 3 2]); light('Position',[-3 -1 3]); material shiny;     
        end
        if strcmp(strop.potential, 'yes')
            contourca(strop, strop.Potential);
        end             
    end
    
    %%   General settings 
    axis('equal'); axis('tight'); grid on; 
    view(46, 41);
end

function [P, t, normals, No, si] = electrodes(P, t, normals, NumberOfElectrodes, PositionOfElectrodes, RadiusOfElectrodes, NumberOfTriangles, ParameterOfElectrodes);
    tolerance = 1e-3; 
    h    = waitbar(0.5, 'Please wait - (re)computing the structure');
    %   Return new arrays P, t, normals
    %   Returns node numbers, arrays of neighbor triangles si, 
    %   and indexes in t into particular electrodes
    No = [];
    %   Triangles attached to every vertex (cell array)
    si = cell(size(P, 1), 1);
    for m = 1:size(P, 1)
        temp  = (t(:,1)==m)|(t(:,2)==m)|(t(:,3)==m);
        si{m} = find(temp>0);
    end    
    if NumberOfElectrodes==0
        close(h);    
        return;
    end    
    %   Find nearest vertex number  
    for m = 1:NumberOfElectrodes
        position = PositionOfElectrodes(m, :);      
        temp1 = P - repmat(position, size(P, 1), 1);
        temp2 = sum(temp1.*temp1, 2);        
        [dummy, index]    = min(temp2); 
        No(m) = index(1);       
    end      
    P(No, :) = PositionOfElectrodes(1:NumberOfElectrodes, :);
    
    %   Add/remove neighbor nodes
    Remove = [];
    Add    = [];
    for m = 1:NumberOfElectrodes          
        %   Find all edges attached to the vertex
        edges  = [];
        for n = 1:length(si{No(m)})
            nodes = t(si{No(m)}(n), 1:3);
            nodes = setdiff(nodes, No(m));
            edges = [edges; [No(m) nodes(1)]; [No(m) nodes(2)]];
        end
        edges           = unique(edges, 'rows');
        NodeNumber = edges(:, 2);   %    a vector
        DIST = (P(NodeNumber, :) - repmat(P(No(m), :), length(NodeNumber), 1));
        dist = sqrt(sum(DIST.*DIST, 2));
        NodeNumber = NodeNumber(dist<RadiusOfElectrodes(m)*(1+tolerance)); 
        Remove          = [Remove; NodeNumber];  
        
        %   Find all nodes within the circle
        position =  PositionOfElectrodes(m, :);
        temp  = ((P(:, 1)-position(1)).^2 + (P(:, 2)-position(2)).^2 + (P(:, 3)-position(3)).^2)<RadiusOfElectrodes(m)^2*(1+2*tolerance);
        Remove          = [Remove; find(temp)]; 
        if NumberOfTriangles(m)>12 
            Remove          = [Remove; No(m)]; 
        end
        Remove = unique(Remove);
        
        MEAN            = P(No(m), :);
        NORMAL          = sum(normals(si{No(m)}, :), 1)/size(edges, 1);
        NORMAL          = NORMAL/norm(NORMAL);       
        PP              = circle0(RadiusOfElectrodes(m), NumberOfTriangles(m), ParameterOfElectrodes(m));       
        %   Rotation (Rodrigues' rotation formula)
        theta = acos(NORMAL(3));
        k     = [-NORMAL(2) +NORMAL(1) 0]';
        K     = repmat([0 0 0]', [1 size(PP, 1)])';
        if dot(k, k)>1e-6
            K     = repmat(k, [1 size(PP, 1)])'/sqrt(dot(k, k));
        end
        PP = PP*cos(theta) + cross(K, PP, 2)*sin(theta) + K.*repmat(dot(K, PP, 2), [1 3])*(1-cos(theta));
        %   movement
        PP(:, 1) = PP(:, 1) + MEAN(1);
        PP(:, 2) = PP(:, 2) + MEAN(2);
        PP(:, 3) = PP(:, 3) + MEAN(3);        
        Add      = [Add; PP];    
    end

    P(Remove, :) = [];
    P = [P; Add];
    for m =  1:NumberOfElectrodes
        if NumberOfTriangles(m)<12
            P = [P; PositionOfElectrodes(m, :)];
        end
    end
    P = unique(P, 'rows');
    [t, normals]    = MyRobustCrust(P);
    si = cell(size(P, 1), 1);
    for m = 1:size(P, 1)
        temp  = (t(:,1)==m)|(t(:,2)==m)|(t(:,3)==m);
        si{m} = find(temp>0);
    end    
    Center = 1/3*(P(t(:, 1),:) + P(t(:, 2),:) + P(t(:, 3),:));

    for m =  1:NumberOfElectrodes
        if NumberOfTriangles(m)>12      
            %   Identify new electrode triangles     
            position =  PositionOfElectrodes(m, :);
            %   Find vertex number
            temp = (Center(:, 1)-position(1)).^2 + (Center(:, 2)-position(2)).^2 + (Center(:, 3)-position(3)).^2;
            [dummy, index] = min(temp);
            Normal = normals(index, :);
            tempn  = sum(normals.*repmat(Normal, size(t, 1), 1), 2);
            temp   = (temp<RadiusOfElectrodes(m)^2)&(tempn>(1-tolerance));
            t(temp, 4) = m;     
        end    
        if NumberOfTriangles(m)<=12
            %   Triangles attached to every new vertex (cell array)        
            position =  PositionOfElectrodes(m, :);
            %   Find vertex number
            temp  = (abs(P(:, 1)-position(1))<eps)&(abs(P(:, 2)-position(2))<eps)&(abs(P(:, 3)-position(3))<eps);
            No(m)    = find(temp);     
            t(si{No(m)}, 4) = m;        
        end 
    end
    close(h);    
end

function [inside] = check(P0, t0, normals0, Points)
    %   inside +1; outside +0    
    inside = zeros(size(Points, 1), 1);     %   column
    Center0 = 1/3*(P0(t0(:, 1),:) + P0(t0(:, 2),:) + P0(t0(:, 3),:));
    N = size(t0, 1);

    for m = 1: size(Points, 1)
        temp1 = repmat(Points(m, :), N, 1) - Center0;
        temp2 = sum(temp1.*temp1, 2);
        [dummy, index] = min(temp2);
        vector = Points(m, :) - Center0(index, :);
        inside(m) = sum(vector.*normals0(index, :))<0;
    end
end

function [t,tnorm]=MyRobustCrust(p)
        %error check

        if nargin>1
            error('The only input must be the Nx3 array of points');
        end

        [n]=size(p,2);
        if n ~=3
               error('Input 3D points must be stored in a Nx3 array');
        end 
        clear  n
    
        %%   Main
        starttime=clock;

        %add points to the given ones, this is usefull
        %to create outside tetraedroms
        %tic
        [p,nshield]=AddShield(p);
        %fprintf('Added Shield: %4.4f s\n',toc)

        %tic
        tetr=delaunayn(p);%creating tedraedron
        tetr=int32(tetr);%use integer to save memory
        %fprintf('Delaunay Triangulation Time: %4.4f s\n',toc)



        %connectivity data
        %find triangles to tetraedrom and tetraedrom to triangles connectivity data
        %tic
        [t2tetr,tetr2t,t]=Connectivity(tetr);
        %fprintf('Connectivity Time: %4.4f s\n',toc)


        %tic
        [cc,r]=CC(p,tetr);%Circumcenters of tetraedroms
        %fprintf('Circumcenters Time: %4.4f s\n',toc)
        clear n


        %tic
        [tbound,Ifact]=Marking(p,tetr,tetr2t,t2tetr,cc,r,nshield);%Flagging tetraedroms as inside or outside
        %fprintf('Walking Time: %4.4f s\n',toc)

        %recostructed raw surface
        t=t(tbound,:);
        % Ifact=Ifact(tbound);
        clear tetr tetr2t t2tetr

        %m
        %tic
        [t,tnorm]=ManifoldExtraction(t,p);
        %fprintf('Manifold extraction Time: %4.4f s\n',toc)



        time=etime(clock,starttime);
        %fprintf('Total Time: %4.4f s\n',time)
        
end

%% Circumcenters
function [cc,r]=CC(p,tetr)
%finds circumcenters from a set of tetraedroms


%batch the code to save memory
ntetr=size(tetr,1);
cutsize=25000;
i1=1;i2=cutsize;
r=zeros(ntetr,1);
cc=zeros(ntetr,1);

if i2>ntetr
    i2=ntetr;%to trigeer the terminate criterion
end


while 1




%points of tetraedrom
p1=(p(tetr(i1:i2,1),:));
p2=(p(tetr(i1:i2,2),:));
p3=(p(tetr(i1:i2,3),:));
p4=(p(tetr(i1:i2,4),:));

%vectors of tetraedrom edges
v21=p1-p2;
v31=p3-p1;
v41=p4-p1;




%Solve the system using cramer method
d1=sum(v41.*(p1+p4)*.5,2);
d2=sum(v21.*(p1+p2)*.5,2);
d3=sum(v31.*(p1+p3)*.5,2);

det23=(v21(:,2).*v31(:,3))-(v21(:,3).*v31(:,2));
det13=(v21(:,3).*v31(:,1))-(v21(:,1).*v31(:,3));
det12=(v21(:,1).*v31(:,2))-(v21(:,2).*v31(:,1));

Det=v41(:,1).*det23+v41(:,2).*det13+v41(:,3).*det12;


detx=d1.*det23+...
    v41(:,2).*(-(d2.*v31(:,3))+(v21(:,3).*d3))+...
    v41(:,3).*((d2.*v31(:,2))-(v21(:,2).*d3));

dety=v41(:,1).*((d2.*v31(:,3))-(v21(:,3).*d3))+...
    d1.*det13+...
    v41(:,3).*((d3.*v21(:,1))-(v31(:,1).*d2));

detz=v41(:,1).*((v21(:,2).*d3)-(d2.*v31(:,2)))...
    +v41(:,2).*(d2.*v31(:,1)-v21(:,1).*d3)...
    +d1.*(det12);

%Circumcenters
cc(i1:i2,1)=detx./Det;
cc(i1:i2,2)=dety./Det;
cc(i1:i2,3)=detz./Det;

%Circumradius
r(i1:i2)=realsqrt((sum((p2-cc(i1:i2,:)).^2,2)));%ciecum radius

if i2==ntetr
    break;%terminate criterion
end

i1=i1+cutsize;
i2=i2+cutsize;

if i2>ntetr
    i2=ntetr;%to trigeer the terminate criterion
end


end

end



%% Connectivity

function [t2tetr,tetr2t,t]=Connectivity(tetr)

%Gets conectivity relantionships among tetraedroms

numt = size(tetr,1);
vect = 1:numt;
t = [tetr(:,[1,2,3]); tetr(:,[2,3,4]); tetr(:,[1,3,4]);tetr(:,[1,2,4])];%triangles not unique
[t,j,j] = unique(sort(t,2),'rows');%triangles
t2tetr = [j(vect), j(vect+numt), j(vect+2*numt),j(vect+3*numt)];%each tetraedrom has 4 triangles


% triang-to-tetr connectivity

nume = size(t,1);
tetr2t  = zeros(nume,2,'int32');
count= ones(nume,1,'int8');
for k = 1:numt

    for j=1:4
        ce = t2tetr(k,j);
        tetr2t(ce,count(ce)) = k;
        count(ce)=count(ce)+1;
    end

end


end      % connectivity()



%% Marking
function [tbound,Ifact]=Marking(p,tetr,tetr2t,t2tetr,cc,r,nshield)
%The more important routine to flag tetredroms as outside or inside

%costants for the algorithm

TOLLDIFF=.01;%tollerance decrease at each iteration
% (the higher the value the more robust but slower is the algorithm. It is also required
% a higher MAXLEVEL value to rach the end of iterations. );


INITTOLL=.99;%starting tollerance

MAXLEVEL=10/TOLLDIFF;%maximum  reachable level 
BRUTELEVEL=MAXLEVEL-50;%level to start  brute continuation



%preallocation
np=size(p,1)-nshield;%nshield = number of shield points put at the end of array
numtetr=size(tetr,1);
nt=size(tetr2t,1);
% deleted=true(numtetr,1);%deleted tetraedroms
% checked=false(numtetr,1);%checked tetraedroms
onfront=false(nt,1);%tetraedroms that need to be checked
% countchecked=0;%counter of checked tetraedroms


%First flag as outside tetraedroms with Shield points

%unvectorized
% for i=1:numtetr
%     for j=1:4
%         if tetr(i,j)>np;
%             deleted(i)=true;
%             checked(i)=true;
%             onfront(t2tetr(i,:))=true;
%             countchecked=countchecked+1;
%             break
%         end
%     end
% end

%vectorized
deleted=any(tetr>np,2);%deleted tetraedroms
checked=deleted;%checked tetraedroms
onfront(t2tetr(checked,:))=true;
countchecked=sum(checked);%counter of checked tetraedroms


%tollerances to mark as in or out
toll=zeros(nt,1)+INITTOLL;
level=0;

%intersection factor
%it is computed from radius of the tetraedroms circumscribed sphere
% and the distance between their center
Ifact=IntersectionFactor(tetr2t,cc,r);
clear cc r




%         Now we scan all tetraedroms. When one is scanned puts on front is
%         neighbor. This means that now  the neighobor can be checked too.
%         At the begining only tetraedroms with shield points are on front,
%         because we are sure the are out. Tetraedrom with high
%         intersection factor  will be marked as equal else different. When
%         I say high i mean under a set tollerance that becames lower as
%         the algorithm progresses. This Aims to avoid errors propagation
%         when a tetraedrom is wrong marked.
%
ids=1:nt;
queue=ids(onfront);
nt=length(queue);
while countchecked<numtetr && level<MAXLEVEL
    level=level+1;%level of scan reached

    for i=1:nt%loop trough triangles <-----better is check only unchecked

        id=queue(i);

        tetr1=tetr2t(id,1);tetr2=tetr2t(id,2);%tetraedroms linked to triangle under analysis
        if  tetr2==0 %do not check boundary triangles
            onfront(id)=false;
            continue

        elseif (checked(tetr1) && checked(tetr2)) %tetraedroms are already checked
            onfront(id)=false;
            continue

        end

        if Ifact(id)>=toll(id) %flag as equal
            if checked(tetr1)%find the checked one between the two
                deleted(tetr2)=deleted(tetr1) ;%flag as equal
                checked(tetr2)=true;%check
                countchecked=countchecked+1;
                onfront(t2tetr(tetr2,:))=true;%put on front all tetreadrom triangles
            else
                deleted(tetr1)=deleted(tetr2) ;%flag as equal
                checked(tetr1)=true;%check
                countchecked=countchecked+1;
                onfront(t2tetr(tetr1,:))=true;%put on front all tetreadrom triangles
            end
            onfront(id)=false;%remove from front


        elseif Ifact(id)<-toll(id)%flag as different
            if checked(tetr1)%find the checked one between the two
                deleted(tetr2)=~(deleted(tetr1)) ;%flag as different
                checked(tetr2)=true;%check
                countchecked=countchecked+1;
                onfront(t2tetr(tetr2,:))=true;%put on front all tetreadrom triangles
            else
                deleted(tetr1)=~(deleted(tetr2)) ;%flag as different
                checked(tetr1)=true;%check
                countchecked=countchecked+1;
                onfront(t2tetr(tetr1,:))=true;%put on front all tetreadrom triangles
            end
            onfront(id)=false;%remove from front


        else
            toll(id)=toll(id)-TOLLDIFF;%tolleraces were too high next time will be lower

        end



    end

    if level==BRUTELEVEL %brute continuation(this may appens when there are almost null volume tetraedroms)
        beep
        warning('Brute continuation necessary')
        onfront(t2tetr(~(checked),:))=true;%force onfront collocation
    end

    %update the queue
    queue=ids(onfront);
    nt=length(queue);

end










%extract boundary triangles
 tbound=BoundTriangles(tetr2t,deleted);


% this is the raw surface and needsimprovements to be used in CAD systems.
% Maybe in my next revision I will add surface post treatments. Anyway for
% grafical purpose this should be good.



%Output Data
numchecked=countchecked/numtetr;
if level==MAXLEVEL
    %warning([num2str(level),' th level was reached\n'])
else
    %fprintf('%4.0f th level was reached\n',level)
end
%fprintf('%4.4f %% of Tetraedroms were checked\n',numchecked*100)



end






%% AddShield
function [pnew,nshield]=AddShield(p)

%adds outside points to the given cloud forming outside tetraedroms

%shield points are very good in detectinf outside tetraedroms. Unfortunatly
%delunany triangulation with these points can be even of 50% slower.

%find the bounding box
maxx=max(p(:,1));
maxy=max(p(:,2));
maxz=max(p(:,3));
minx=min(p(:,1));
miny=min(p(:,2));
minz=min(p(:,3));

%give offset to the bounding box
step=max(abs([maxx-minx,maxy-miny,maxz-minz]));

maxx=maxx+step;
maxy=maxy+step;
maxz=maxz+step;
minx=minx-step;
miny=miny-step;
minz=minz-step;

N=10;%number of points of the shield edge

step=step/(N*N);%decrease step, avoids not unique points



nshield=N*N*6;

%creating a grid lying on the bounding box
vx=linspace(minx,maxx,N);
vy=linspace(miny,maxy,N);
vz=linspace(minz,maxz,N);




[x,y]=meshgrid(vx,vy);
facez1=[x(:),y(:),ones(N*N,1)*maxz];
facez2=[x(:),y(:),ones(N*N,1)*minz];
[x,y]=meshgrid(vy,vz-step);
facex1=[ones(N*N,1)*maxx,x(:),y(:)];
facex2=[ones(N*N,1)*minx,x(:),y(:)];
[x,y]=meshgrid(vx-step,vz);
facey1=[x(:),ones(N*N,1)*maxy,y(:)];
facey2=[x(:),ones(N*N,1)*miny,y(:)];

%add points to the p array
pnew=[p;
    facex1;
    facex2;
    facey1;
    facey2;
    facez1;
    facez2];

% figure(4)
% plot3(pnew(:,1),pnew(:,2),pnew(:,3),'.g')

end



%% BoundTriangles
function tbound=BoundTriangles(tetr2t,deleted)
%extracts boundary triangles from a set tetr2t connectivity and form the
%deleted vector which tells tetraedroms that are marked as out

nt=size(tetr2t,1);%number of totals triangles

tbound=true(nt,2);%inizilize to keep shape in next operation

ind=tetr2t>0;%avoid null index
tbound(ind)=deleted(tetr2t(ind));%mark 1 for deleted 0 for kept tetraedroms

tbound=sum(tbound,2)==1;%bounary triangles only have one tetraedrom

end


%% Intersection factor
function Ifact=IntersectionFactor(tetr2t,cc,r)
nt=size(tetr2t,1);
Ifact=zeros(nt,1);%intersection factor
%it is computed from radius of the tetraedroms circumscribed sphere
% and the distance between their center
i=tetr2t(:,2)>0;
distcc=sum((cc(tetr2t(i,1),:)-cc(tetr2t(i,2),:)).^2,2);%distance between circumcenters
Ifact(i)=(-distcc+r(tetr2t(i,1)).^2+r(tetr2t(i,2)).^2)./(2*r(tetr2t(i,1)).*r(tetr2t(i,2)));

%unvectorized
% for i=1:nt
%     if tetr2t(i,2)>0 %jump boundary tetraedrom
%         distcc=sum((cc(tetr2t(i,1),:)-cc(tetr2t(i,2),:)).^2,2);%distance between circumcenters
%         %intersection factor
%         Ifact(i)=(-distcc+r(tetr2t(i,1))^2+r(tetr2t(i,2))^2)/(2*r(tetr2t(i,1))*r(tetr2t(i,2)));
%     end
% end
end




%% Manifold Extraction

function [t,tnorm]=ManifoldExtraction(t,p)
%Given a set of trianlges,
%Buils a manifolds surface with the ball pivoting method.



% building the etmap

numt = size(t,1);
vect = 1:numt;                                                             % Triangle indices
e = [t(:,[1,2]); t(:,[2,3]); t(:,[3,1])];                                  % Edges - not unique
[e,j,j] = unique(sort(e,2),'rows');                                        % Unique edges
te = [j(vect), j(vect+numt), j(vect+2*numt)];
nume = size(e,1);
e2t  = zeros(nume,2,'int32');

clear vect j
ne=size(e,1);
np=size(p,1);


count=zeros(ne,1,'int32');%numero di triangoli candidati per edge
etmapc=zeros(ne,4,'int32');
for i=1:numt

    i1=te(i,1);
    i2=te(i,2);
    i3=te(i,3);



    etmapc(i1,1+count(i1))=i;
    etmapc(i2,1+count(i2))=i;
    etmapc(i3,1+count(i3))=i;


    count(i1)=count(i1)+1;
    count(i2)=count(i2)+1;
    count(i3)=count(i3)+1;
end

etmap=cell(ne,1);
for i=1:ne

    etmap{i,1}=etmapc(i,1:count(i));

end
clear  etmapc

tkeep=false(numt,1);%all'inizio nessun trinagolo selezionato


%Start the front

%building the queue to store edges on front that need to be studied
efront=zeros(nume,1,'int32');%exstimate length of the queue

%Intilize the front


         tnorm=Tnorm(p,t);%get traingles normals
         
         %find the highest triangle
         [foo,t1]=max( (p(t(:,1),3)+p(t(:,2),3)+p(t(:,3),3))/3);

         if tnorm(t1,3)<0
             tnorm(t1,:)=-tnorm(t1,:);%punta verso l'alto
         end
         
         %aggiungere il ray tracing per verificare se il triangolo punta
         %veramente in alto.
         %Gli altri triangoli possono essere trovati sapendo che se un
         %triangolo ha il baricentro pi? alto sicuramente contiene il punto
         %pi? alto. Vanno analizzati tutto i traingoli contenenti questo
         %punto
         
         
            tkeep(t1)=true;%primo triangolo selezionato
            efront(1:3)=te(t1,1:3);
            e2t(te(t1,1:3),1)=t1;
            nf=3;%efront iterato
      

while nf>0


    k=efront(nf);%id edge on front

    if e2t(k,2)>0 || e2t(k,1)<1 || count(k)<2 %edge is no more on front or it has no candidates triangles

        nf=nf-1;
        continue %skip
    end
  
   
      %candidate triangles
    idtcandidate=etmap{k,1};

    
     t1=e2t(k,1);%triangle we come from
    
   
        
    %get data structure
%        p1
%       / | \
%  t1 p3  e1  p4 t2(idt)
%       \ | /  
%        p2
         alphamin=inf;%inizilizza
          ttemp=t(t1,:);
                etemp=e(k,:);
                p1=etemp(1);
                p2=etemp(2);
                p3=ttemp(ttemp~=p1 & ttemp~=p2);%terzo id punto
        
                
         %plot for debug purpose
%          close all
%          figure(1)
%          axis equal
%          hold on
%          
%          fs=100;
%         
%          cc1=(p(t(t1,1),:)+p(t(t1,2),:)+p(t(t1,3),:))/3;
%          
%          trisurf(t(t1,:),p(:,1),p(:,2),p(:,3))
%          quiver3(cc1(1),cc1(2),cc1(3),tnorm(t1,1)/fs,tnorm(t1,2)/fs,tnorm(t1,3)/fs,'b');
%                 
       for i=1:length(idtcandidate)
               t2=idtcandidate(i);
               if t2==t1;continue;end;
                
               %debug
%                cc2=(p(t(t2,1),:)+p(t(t2,2),:)+p(t(t2,3),:))/3;
%          
%                 trisurf(t(t2,:),p(:,1),p(:,2),p(:,3))
%                 quiver3(cc2(1),cc2(2),cc2(3),tnorm(t2,1)/fs,tnorm(t2,2)/fs,tnorm(t2,3)/fs,'r');
%                
%                

               
                ttemp=t(t2,:);
                p4=ttemp(ttemp~=p1 & ttemp~=p2);%terzo id punto
        
   
                %calcola l'angolo fra i triangoli e prendi il minimo
              
                
                [alpha,tnorm2]=TriAngle(p(p1,:),p(p2,:),p(p3,:),p(p4,:),tnorm(t1,:));
                
                if alpha<alphamin
                    
                    alphamin=alpha;
                    idt=t2;  
                    tnorm(t2,:)=tnorm2;%ripristina orientazione   
                     
                    %debug
%                      quiver3(cc2(1),cc2(2),cc2(3),tnorm(t2,1)/fs,tnorm(t2,2)/fs,tnorm(t2,3)/fs,'c');
                    
                end
                %in futuro considerare di scartare i trianoli con angoli troppi bassi che
                %possono essere degeneri
                
       end

   %update front according to idttriangle
          tkeep(idt)=true;
        for j=1:3
            ide=te(idt,j);
           
            if e2t(ide,1)<1% %Is it the first triangle for the current edge?
                efront(nf)=ide;
                nf=nf+1;
                e2t(ide,1)=idt;
            else                     %no, it is the second one
                efront(nf)=ide;
                nf=nf+1;
                e2t(ide,2)=idt;
            end
        end
        
     
        

         nf=nf-1;%per evitare di scappare avanti nella coda e trovare uno zero
end

t=t(tkeep,:);
tnorm=tnorm(tkeep,:);

end


%% TriAngle
function  [alpha,tnorm2]=TriAngle(p1,p2,p3,p4,planenorm)

%per prima cosa vediamo se il p4 sta sopra o sotto il piano identificato
%dalla normale planenorm e il punto p3

test=sum(planenorm.*p4-planenorm.*p3);



%Computes angle between two triangles
v21=p1-p2;
v31=p3-p1;

tnorm1(1)=v21(2)*v31(3)-v21(3)*v31(2);%normali ai triangoli
tnorm1(2)=v21(3)*v31(1)-v21(1)*v31(3);
tnorm1(3)=v21(1)*v31(2)-v21(2)*v31(1);
tnorm1=tnorm1./norm(tnorm1);



v41=p4-p1;
tnorm2(1)=v21(2)*v41(3)-v21(3)*v41(2);%normali ai triangoli
tnorm2(2)=v21(3)*v41(1)-v21(1)*v41(3);
tnorm2(3)=v21(1)*v41(2)-v21(2)*v41(1);
tnorm2=tnorm2./norm(tnorm2);
alpha=tnorm1*tnorm2';%coseno dell'angolo
%il coseno considera l'angolo fra i sempipiani e non i traigoli, ci dice
%che i piani sono a 180 se alpha=-1 sono concordi se alpha=1, a 90?

alpha=acos(alpha);%trova l'angolo

%Se p4 sta sopra il piano l'angolo ? quello giusto altrimenti va maggiorato
%di 2*(180-alpha);

if test<0%p4 sta sotto maggioriamo
   alpha=alpha+2*(pi-alpha);
end

%         fs=100;
%          cc2=(p1+p2+p3)/3;
%        quiver3(cc2(1),cc2(2),cc2(3),tnorm1(1)/fs,tnorm1(2)/fs,tnorm1(3)/fs,'m');
%        cc2=(p1+p2+p4)/3;
%               quiver3(cc2(1),cc2(2),cc2(3),tnorm2(1)/fs,tnorm2(2)/fs,tnorm2(3)/fs,'m');

%vediamo se dobbiamo cambiare l'orientazione del secondo triangolo
%per come le abbiamo calcolate ora tnorm1 t tnorm2 non rispettano
%l'orientamento
testor=sum(planenorm.*tnorm1);
if testor>0 
    tnorm2=-tnorm2;
end

end


%% Tnorm

function tnorm1=Tnorm(p,t)
    %Computes normalized normals of triangles


    v21=p(t(:,1),:)-p(t(:,2),:);
    v31=p(t(:,3),:)-p(t(:,1),:);

    tnorm1(:,1)=v21(:,2).*v31(:,3)-v21(:,3).*v31(:,2);%normali ai triangoli
    tnorm1(:,2)=v21(:,3).*v31(:,1)-v21(:,1).*v31(:,3);
    tnorm1(:,3)=v21(:,1).*v31(:,2)-v21(:,2).*v31(:,1);

    L=sqrt(sum(tnorm1.^2,2));

    tnorm1(:,1)=tnorm1(:,1)./L;
    tnorm1(:,2)=tnorm1(:,2)./L;
    tnorm1(:,3)=tnorm1(:,3)./L;
end
