function fig_hdl = E34
%   SYNTAX
%   E34
%   DESCRIPTION
% This module is an accurate MoM solution for modeling a capacitive 2D
% touchpad with a lattice of rows/columns ("mutual capacitance" method). The
% default scheme: All electrodes are assigned voltage V1=1V. Ground plane
% and the human finger phantom are assigned voltage V2=0V.
% Induced charge Q on every electrode is reported in pC which gives “mutual capacitance” in the form
% C=Q(V2-V1).
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
global t;
global c;
global Area;
global Center;
global normals;
global Size;
global Points;
global epsr;
global contrasts;

%%   Lattice
global NumberOfPads;
global PadSize;     
global PadSpacing;  
global PadGap;   
global Triangles;   
global Row;

NumberOfPads            = 5;
PadSize                 = 0.0033;
PadSpacing              = 0.005;
PadGap                  = 0.0005;
Triangles               = 1500; 
Row                     = 2;

input.NumberOfPads  = NumberOfPads;
input.PadSize       = PadSize;
input.PadSpacing    = PadSpacing;
input.PadGap        = PadGap;
input.Triangles     = Triangles;
input.Row           = Row;          %   Driven row

%% Dielectric brick
strbrick.L          = 0.0;            %    brick length, m
strbrick.W          = 0.0;            %    brick width, m
strbrick.H          = 0.001;          %    brick height, m
strbrick.X          = 0.0;            %    center position x, m
strbrick.Y          = 0.0;            %    center position y, m
strbrick.Z          = 0;              %    center position z, m        
strbrick.epsr       = 3;              %    relative dielectric constant of the brick
strbrick.par  =1.0;                   %    uniform (0) or non-uniform (1) grid
strbrick.Tr  = 400;                 

%%  Finger
strcylinder.yes        = 'yes';       %     yes or no  
strcylinder.R          = 0.003;       %     cylinder radius, m
strcylinder.H          = 0.02;        %     cylinder height, m
strcylinder.X          = +0.0025;     %     center position x, m
strcylinder.Y          = -0.0050;     %     center position y, m
strcylinder.Z          = 0.0106;      %     center position z, m
strcylinder.ax         = 0;           %     rotate about the x-axis, deg
strcylinder.ay         = 0;           %     rotate about the y-axis, deg
strcylinder.az         = 0;           %     rotate about the z-axis, deg
strcylinder.par        = 1;           %    uniform (0) or non-uniform (1) grid
strcylinder.Tr         = 40;          %    approximate number of triangular patches

strei.epsr = 1;                             %    dielectric constant of the external medium

%   Data for the plate  - second conductor('plate' is always selected) (uitable)
strge.L(2)              = 0.01;         %   plate length, m
strge.W(2)              = 0.01;         %   plate width, m
strge.X(2)              = 0;           %   center position x, m
strge.Y(2)              = 0;           %   center position y, m
strge.Z(2)              = 0.0;        %   center position z, m
strge.ax(2)             = 0;           %   rotate about the x-axis, deg
strge.ay(2)             = 0;           %   rotate about the y-axis, deg
strge.az(2)             = 0;           %   rotate about the z-axis, deg
strge.par(2)            = 0.3;           %   uniform (0) or non-uniform (1) grid
strge.Tr(2)             = 200;         %   approximate number of triangular patches

%%  Geometry parameters - dielectric object parameters
%   Parameters for a dielectric object (uimenu)
objecttype = 'brick';  %   'none' or 'sphere' or 'brick' or 'cylinder' is allowed

%%  Output graphics parameters
%   Surface charge and general visualization (uitable)
%   Parameters to scale charge distribution for better visual inspection
strsc.positive = 0.25;         %   positive charge densities higher than this number times the average
%   positive charge density are assigned the same value
strsc.negative = 0.25;         %   negative charge densities smaller than this number times the average
%   negative charge density are assigned the same value

%   Visualization of the E-field in the observation plane (uitable)
strop.yes = 'no';           %   include (yes) or not (no) plot of the E-field
strop.potential = 'no';     %   include (yes) or not (no) plot of the electric potential
strop.planetype   = 'xy';   %   xy, xz, or yz plane
strop.planex = 0.0;         %   plane center: x in m
strop.planey = -0.0025;         %   plane center: y in m
strop.planez = 0.005;       %   plane center: z in m
strop.planesizex   = 0.025; %   plane length in m
strop.planesizey   = 0.025; %   plane width in m
strop.divisionsx   = 16;    %   divisions with respect to length
strop.divisionsy   = 16;    %   divisions with respect to width
strop.arrow = 1.5;          %   relative arrow size versus default size

%   Parameters for the observation point(s) (used to obtain exact values of the field) (uitable)
stroc.yes         = 'no';     %   'yes' - present; 'no' - absent
stroc.x           = 0.0;       %   x position in m
stroc.y           = 0.0;       %   y position in m
stroc.z           = 0.0;       %   z position in m
stroc.size        = 1.5;       %   relative marker size versus default size

%%   Numerical parameters (uitable)
global R; R   = 5;        %    dimensionless radius of an enclosing sphere for precise integration
%    R=0 - only self integrals are calculated precisely
%    R=10- integrals for all neighbor triangles whose
%    center-to-center distances from the observation triangle
%    are less than ten times the effective triangle size are calculated
%    precisely
global gauss; gauss = 7;      %    Number of integration points in the Gaussian quadrature
%    Numbers 1, 4, 7, 13, and 25 are permitted

%%  End of GUI window - input parameters

%   Input Parameters--------
%%  Global parameters of output results
strout.E(1,1)       = 0;
strout.E(1,2)       = 0; 
strout.E(1,3)       = 0;        %    E-field at the observation point (if any)
strout.phi          = 0;        %    Electric potential at the observation point (if any)
strout.PatchesM     = 0;        %    Number of triangular patches in the metal mesh
strout.PatchesF     = 0;        %    Number of triangular patches in the finger
strout.PatchesD     = 0;        %    Number of triangular patches in the dielectric mesh
strout.quality      = 0;        %    Minimum triangle quality
strout.time1        = 0;        %    CPU time in sec for filling the MoM matrix
strout.time2        = 0;        %    CPU time in sec for solving the system of MoM eqs.
strout.SelfCapacitance  = zeros(input.NumberOfPads, input.NumberOfPads); 
strout.ChargeTotal  = 0;

%% End of GUI window - input parameters

% Initialize handles structure
handles = struct();
handles.geometry = 0;   %   The real handles
handles.simulate = 0;   %   Only an indicator that the simulations are complete
handles.object = 2;
% Create all UI controls
build_gui();
geometry();

%%  Functions
   function outputgraphics
        %%  Definition of nodal points in the observation plane (array Points) and initializaton of the E-field
        if strcmp(strop.planetype, 'xy') nx = 0; ny = 0; nz = 1; end
        if strcmp(strop.planetype, 'xz') nx = 0; ny = 1; nz = 0; end
        if strcmp(strop.planetype, 'yz') nx = 1; ny = 0; nz = 0; end
        [Points, dummy] = plate0(strop.planesizex, strop.planesizey, strop.divisionsx, strop.divisionsy, 0, nx, ny, nz);
        Points(:, 1) = Points(:, 1) +  strop.planex;
        Points(:, 2) = Points(:, 2) +  strop.planey;
        Points(:, 3) = Points(:, 3) +  strop.planez;
        strop.Points = Points;
        strop.E              = zeros(size(Points));
        strop.Potential      = zeros(size(Points, 1), 1);
        %%  Find the E-field in a plane
        if strcmp(strop.yes, 'yes')
            msg = 1;
            strop.E = efield(Points, c, P, t, Center, Area, normals, Size, R, eps0, msg);
        end
         %%  Find the potential in a plane
        if strcmp(strop.potential, 'yes')
            msg = 1;
            strop.Potential = potential(Points, c, P, t, Center, Area, normals, Size, R, eps0, msg);
        end
        %%  Find the E-field at an observation point(s)
        strout.E          = 'none selected';
        strout.phi        = 'none selected';
        if strcmp(stroc.yes, 'yes')
            msg = 0;
            points = [stroc.x stroc.y stroc.z];
            strout.E = efield(points, c, P, t, Center, Area, normals, Size, R, eps0, msg);
            strout.phi = potential(points, c, P, t, Center, Area, normals, Size, R, eps0, msg);
        end      
        %%   Output graphics
        set(0, 'CurrentFigure', handles.figure1);
        if handles.geometry ~= 0; delete(handles.geometry); end
        io = 1;
        handles.geometry = subplot(1,1,1,'Parent',handles.uipanel2);
        graphics_S32(io, strop, stroc, strei, P, t, c, Area, strsc, strbrick);    
   end

    function cleaning
        %   cleaning figure without re-plotting
        strout = clearoutput(strout); io = 0;
        if handles.geometry ~= 0; delete(handles.geometry); handles.geometry = 0;end
        handles.geometry = subplot(1,1,1,'Parent',handles.uipanel2);
        graphics_S32(io, strop, stroc, strei, P, t, c, Area, strsc, strbrick);
    end

    function strout = clearoutput(strout)
        %   Clear results anyways
        strout.Cnum         = 0;        %    Numerical capacitance in pF       
        strout.E(1,1)       = 0;
        strout.E(1,2)       = 0; 
        strout.E(1,3)       = 0;        %    E-field at the observation point (if any)
        strout.phi          = 0;        %    Electric potential at the observation point (if any)
        strout.time1        = 0;        %    CPU time in sec for filling the MoM matrix
        strout.time2        = 0;        %    CPU time in sec for solving the system of MoM eqs.
        strout.ChargeLower  = 0;        %   Total charge of the first conductor
        strout.ChargeUpper  = 0;        %   Total charge of the second conductor
        strout.ChargeDielectric = 0;    %   Total charge of the dielectric object
        strout.ChargeTotal  = 0;        %   Sum of charges for the entire structure
    end

    function geometry()        
        %%  Complete geometry - lattice       
        [PD{1}, tD{1}] = lattice(input);
        index = (tD{1}(:, 4) ==0);
        tD{1}(index, 4) = 0.5; 
        PD{1}(:, 3) = 0; 
        [centersD{1}, normalsD{1}] = normcenters(PD{1}, tD{1});        
        %   Data for the dielectric brick is based on the lattice data
        strbrick.L          = max(PD{1}(:, 1))-min(PD{1}(:, 1));            %    brick length, m
        strbrick.W          = max(PD{1}(:, 2))-min(PD{1}(:, 2));            %    brick width, m       
        strbrick.X          = (max(PD{1}(:, 1))+min(PD{1}(:, 1)))/2;        %    center position x, m
        strbrick.Y          = (max(PD{1}(:, 2))+min(PD{1}(:, 2)))/2;        %    center position y, m    
        BrickTop = strbrick.H/2 + strbrick.Z; 
        
        %%  Complete geometry - finger
        if strcmp(strcylinder.yes, 'yes')
            [PD{2}, tD{2}]  = cylinder(strcylinder.R, strcylinder.H, strcylinder.Tr, strcylinder.par);
            [PD{2}] = rotatex(PD{2}, strcylinder.ax);
            [PD{2}] = rotatey(PD{2}, strcylinder.ay);
            [PD{2}] = rotatez(PD{2}, strcylinder.az);
            [centersD{2}, normalsD{2}] = normcenters(PD{2}, tD{2});
            %   movement
            PD{2}(:, 1) = PD{2}(:, 1) + strcylinder.X;
            PD{2}(:, 2) = PD{2}(:, 2) + strcylinder.Y;
            PD{2}(:, 3) = PD{2}(:, 3) + strcylinder.Z;
            centersD{2}(:, 1) = centersD{2}(:, 1) + strcylinder.X;
            centersD{2}(:, 2) = centersD{2}(:, 2) + strcylinder.Y;
            centersD{2}(:, 3) = centersD{2}(:, 3) + strcylinder.Z;
            tD{2}(:, 4) = 1000; cD{2} = zeros(size(tD{1}, 1), 1);    %   finger will have a fourth index of 1000      
            NF = size(tD{2}, 1);
            strout.PatchesF = NF;
            FingerBottom = min(PD{2}(:, 3));
            if FingerBottom<BrickTop
                PD{2}(:, 3)= PD{2}(:, 3) + (BrickTop-FingerBottom) + 1e-4;
                centersD{2}(:, 3) = centersD{2}(:, 3) + (BrickTop-FingerBottom) + 1e-4;
                strcylinder.az = strcylinder.az + (BrickTop-FingerBottom) + 1e-4;
            end         
        else
            strout.PatchesF = 0;
        end        
        
        %%  Combining meshes together (while keeping the fourth index)
        P = [];
        t = [];
        normals = [];
        M = 1;
        if strcmp(strcylinder.yes, 'yes')
            M = 2;
        end
        for m = 1:M
            tD{m}(:, 1:3)   = tD{m}(:, 1:3) + size(P, 1);
            P               = [P' PD{m}']';
            t               = [t' tD{m}']';
            normals         = [normals' normalsD{m}']';     %   normals
        end
        
        %%   Complete geometry - dielectric mesh
        tD = [];
        PD = [];               
        [PD{1}, tD{1}] = brick(strbrick.L, strbrick.W, strbrick.H, strbrick.Tr, strbrick.par);     
        [centersD{1}, normalsD{1}] = normcenters(PD{1}, tD{1});
        PD{1}(:, 1) = PD{1}(:, 1) + strbrick.X;
        PD{1}(:, 2) = PD{1}(:, 2) + strbrick.Y;
        PD{1}(:, 3) = PD{1}(:, 3) + strbrick.Z;
        centersD{1}(:, 1) = centersD{1}(:, 1) + strbrick.X;
        centersD{1}(:, 2) = centersD{1}(:, 2) + strbrick.Y;
        centersD{1}(:, 3) = centersD{1}(:, 3) + strbrick.Z;
        epsr = strbrick.epsr;
        contrast{1} = (epsr - strei.epsr)/(epsr + strei.epsr);
        tD{1}(:, 4) = 0; cD{1} = zeros(size(tD{1}, 1), 1); %   all dielectrics will have a fourth index zero!!!        
        
        %%  Combining metal and dielectric mesh together
        NM = size(t, 1);
        contrasts = [];    %   no dielectric contrast for metal paches
        for m = 1:length(tD)
            tD{m}(:, 1:3)   = tD{m}(:, 1:3) + size(P, 1);
            P               = [P' PD{m}']';
            t               = [t' tD{m}']';
            normals         = [normals' normalsD{m}']';     %   dielectric normals
            contrasts       = [contrasts' contrast{m}*ones(1, size(tD{m}, 1))]'; % starts with diel.
        end
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
        io = 0;  strei.normalize = 1;
        graphics_S32(io, strop, stroc, strei, P, t, c, Area, strsc, strbrick);        
    end

    function simulate()
        h    = waitbar(0, 'Please wait - filling the MoM matrix');
        time1 = cputime;
        %%  Parameter initialization for the combined mesh
        ND = strout.PatchesD;
        NM = strout.PatchesM;
        NF = strout.PatchesF;
        Center = zeros(length(t), 3);   %   face center
        Area   = zeros(length(t), 1);   %   face area
        Normal = zeros(length(t), 3);   %   face normal
        Size   = zeros(length(t), 1);   %   face size defined as distance from center to furthest vertex
        
        %%   Find base parameters for all faces (metal or dielectric)                
        %   Metal
        for m = 1:NM
            Vertexes        = P(t(m, 1:3)', :)';
            r1              = Vertexes(:, 1);
            r2              = Vertexes(:, 2);
            r3              = Vertexes(:, 3);
            tempv           = cross(r2-r1, r3-r1);
            temps           = sqrt(tempv(1)^2 + tempv(2)^2 + tempv(3)^2);
            normals(m, :)   = tempv'/temps; %   these normals are consistent with pot2
            Area(m)         = temps/2;
            tempc           = (r1+r2+r3)/3;
            Center(m, :)    = tempc';
            Size(m)         = sqrt(max([sum((tempc-r1).^2) sum((tempc-r2).^2) sum((tempc-r3).^2)]));
        end
        %   Dielectric
        for m = NM+1:NM+ND
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
        
        %%  Filling ZM
        %   Prepare distance matrix DIST for ZM
        %   First row - distances from the center of face #1 to all other face centers
        %   Second row - distances from the center of face #2 to all other face centers, etc.
        DIST = zeros(NM, NM + ND);
        for m = 1:NM
            temp = Center'- repmat(Center(m, :)', 1, length(t));
            DIST(m, :) = sqrt(dot(temp, temp));
        end
        %   MoM matrix ZM initialization (center point integration)
        ZM = 1./DIST;
        for n = 1:NM+ND
            ZM(:, n) = ZM(:, n)*Area(n);
        end
        %   Loop over columns of impedance matrix ZM!
        %   n is the number of the inner triangle (every column has the only the n-th inner triangle)
        for n =1:NM+ND
            normdummy    = normals(n, :);
            index   = find(DIST(1:NM, n)'./(Size(1:NM)'*Size(n))<=R+1e-16);   % index is local
            r1      = P(t(n, 1), :);    %   row
            r2      = P(t(n, 2), :);    %   row
            r3      = P(t(n, 3), :);    %   row
       
            m1 = repmat(index'*IndexS-IndexS, [1 IndexS])';
            m3 = m1(:) + ObsIndex(1:IndexS*length(index));
            [I, IRho] = potint(r1, r2, r3, normdummy, ObsPoint(m3, :)); %   I was calculated with the area Area(n)
            ZM(index, n) = sum(WeightsS(1:length(index), :).*reshape(I, IndexS, length(index))', 2);     
            
            waitbar(NM/(NM+ND)*n/length(t));
        end
        %   Full matrix ZM
        ZM = +ZM/(4*pi*eps0);        
        
        %%  Filling ZD
        %   Prepare distance matrix DIST for ZD
        %   First row - distances from the center of face #1 to all other face centers
        %   Second row - distances from the center of face #2 to all other face centers, etc.
        %   m shifts by NM, but not for DIST
        ZD      = zeros(ND, NM + ND);
        DIST    = zeros(ND, NM + ND);
        for m = 1:ND    %   rowwise here!
            mm = m + NM;
            temp        = Center'- repmat(Center(mm, :)', 1, length(t));    %   this is rn - rm
            DIST(m, :)  = sqrt(dot(temp, temp));
            ZD(m, :)    = +sum(repmat(normals(mm, :)', 1, length(t)).*temp, 1)./(DIST(m, :)).^3;
        end
        for n = 1:NM+ND
            ZD(:, n) = ZD(:, n)*Area(n);
        end
        %   Loop over columns of impedance matrix ZD!
        %   n is the number of the inner triangle (every column has the only the n-th inner triangle)
        %   over which the exact integration is done (NM - metal + ND - dielectric)
        rownumber = NM+1:NM+ND; %   global row number
        for n = 1:NM+ND
            index   = find(DIST(:, n)'./(Size(rownumber)'*Size(n))<=R+1e-16);   % index is local (1:ND)
            r1      = P(t(n, 1), :);    %   row
            r2      = P(t(n, 2), :);    %   row
            r3      = P(t(n, 3), :);    %   row
            
            m1 = repmat(index'*IndexS-IndexS, [1 IndexS])';
            m3 = m1(:) + ObsIndex(1:IndexS*length(index)) + NM*IndexS;
            I = potint2(r1, r2, r3, normals(n, :), ObsPoint(m3, :)); %   I was calculated with the area Area(n)
            J = zeros(length(index)*IndexS, 1);
            for p = 1:IndexS
                J(p:IndexS:end) = sum(I(p:IndexS:end, :).*normals(index+NM, :), 2);
            end
            ZD(index, n) = sum(WeightsS(1:length(index), :).*reshape(J, IndexS, length(index))', 2);       
            waitbar(NM/(NM+ND) + ND/(NM+ND)*n/length(t));
        end
        
        %%   Add coefficients to ZD MoM matrix
        for m = 1:ND
            ZD(m, :) = ZD(m, :)*contrasts(m)/(4*pi*eps0*strei.epsr);
        end
        for m = 1:ND
            mm = m + NM;
            ZD(m, mm) = 1/(2*eps0*strei.epsr);
        end
        close(h);
        
        %%  Construct the full MoM matrix
        Z(1:NM, :)          = ZM;
        Z(NM+1:NM+ND, :)    = ZD;
        clear ZM; clear ZD;
        
        %%   The full voltage vector
        V1 = 1;
        V                       =  zeros(length(t), 1);
        V(t(:, 4)>0.5&t(:, 4)<1000)       =  V1;
        
        %% Charge conservation law (critical here)
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %  Enforce the charge conservation law for the row-column pair        
%         index1 = find(t(:, 4)==input.Row);
%         index2 = find(t(:, 4)==input.Column+input.NumberOfPads);
%         index = sort([index1; index2]);
%         for m = 1:length(index)-1
%             n = index(m);
%             V(n)    = V(n)    - V(index(end));
%             Z(n, :) = Z(n,: ) - Z(index(end), :);
%         end
%         % Last row
%         Z(index(end), :) = 0;
%         Z(index(end), index) = Area(index)';        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         %   enforce charge consevation law for the entire conductor structure
%         for m = 1:NM-1
%             V(m)    = V(m) - V(end);
%             Z(m, :) = Z(m,: ) - Z(NM, :);
%         end
%         % Last metal row
%         Z(NM, :) = 0;
%         Z(NM, 1:NM) = Area(1:NM)';
%         V(end) = 0;        
            
%         %  Enforce the charge conservation law for the finger only
%         for m = NM-NF+1:NM          
%             Z(m, :) = Z(m,: ) - Z(NM, :);
%         end
%         if NF>0
%             % Last row
%             Z(NM, :) = 0;
%             index = (t(:, 4) == 1000);
%             Z(NM, index) = Area(index)';        
%             V(NM) = 0;
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Enforce the charge conservation law (resulting in zero charge) for the dielectric object
%         for m = NM+1:NM+ND-1
%             V(m)    = V(m)    - V(NM+ND);
%             Z(m, :) = Z(m,: ) - Z(NM+ND, :);
%         end
%         % Last row
%         Z(NM+ND, NM+1:NM+ND) = ones(1, NM+ND-NM).*Area(NM+1:NM+ND)';
%         Z(NM+ND, 1:NM) = 0;
%         V(NM+ND) = 0;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        time1 = cputime - time1;
        
        %%   Solve MoM equations
        if (length(t)>1000)
            h    = waitbar(0, 'Please wait - solving the MoM equations');
            waitbar(0.5);
        end
        time2 = cputime;
        
        %   Jacobi (diagonal) preconditioner
        DIAG = diag(Z);
        for m = 1:size(Z, 1)    %   rowwise
            Z(m, :) = Z(m, :)/DIAG(m);
        end
        V = V./DIAG;
        
        c       = Z\V;                  %   solution for charge coefficients
        clear Z;                        %   clear impedance matrix
        time2 = cputime - time2;
        if (length(t)>1000)
            close(h);
        end
        
        Charge              = c'.*Area';                  %   charges for every face
        ChargeTotal         = sum(Charge);
        
        for m = 1:input.NumberOfPads
            for n = 1:input.NumberOfPads
                index = n + (m-1)*input.NumberOfPads;
                strout.SelfCapacitance(m, n) = strbrick.epsr*sum(Charge(t(:, 4)==index))/V1; 
            end
        end
        
    
        %%  Fields and output graphics
        handles.simulate = 1;
        outputgraphics;
        
        %%  GUI window - output parameters      
        strout.ChargeTotal      = ChargeTotal;
        strout.E;                                           %   E-field at the observation point (if any)
        strout.PatchesM;                                    %   Number of triangular patches in the metal mesh
        strout.PatchesD;                                    %   Number of triangular patches in the dielectric mesh
        strout.PatchesF;                                    %   Number of triangular patches in the finger mesh
        strout.quality  = min(simpqual(P, t));
        strout.time1    = time1;                            %   CPU time in sec for filling the MoM matrix
        strout.time2    = time2;                            %   CPU time in sec for solving the system of MoM eqs.     
        strout
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
            'Name', 'E34-Self-capacitance touchpad', ...
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
        
        uimenu(...
            'Parent', handles.fFileMenu,...
            'HandleVisibility','callback', ...
            'Label','New Project',...
            'Callback', {@New_Callback});
                
        function New_Callback(~, ~)
            E34;
        end
        
        uimenu(...
            'Parent', handles.fFileMenu,...
            'HandleVisibility','callback', ...
            'Label','Open',...
            'Callback', {@Open_Callback});
        
        function Open_Callback(~, ~)
            
            file = uigetfile('*.mat', 'Select project file to open');
            if ~isequal(file, 0)
                %% Data for the first plate
                M   = load(file);
                input               = M.input;
                strsc               = M.strsc;
                strge               = M.strge;
                strop               = M.strop;
                strei               = M.strei;
                objecttype          = M.objecttype;
                stroc               = M.stroc;
                R                   = M.R;
                gauss               = M.gauss;                
                strout              = M.strout;
                epsr                = M.epsr;
                strbrick            = M.strbrick;
                strcylinder         = M.strcylinder;
                switch objecttype
                    case 'sphere'
                        handles.object = 1;
                    case 'brick'
                        handles.object = 2;
                    case 'cylinder'
                        handles.object = 3;
                end               
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
            end
        end
        
        uimenu(...
            'Parent', handles.fFileMenu,...
            'HandleVisibility','callback', ...
            'Label','Save',...
            'Callback', {@Save_Callback});
        
        function Save_Callback(~, ~)
            save('E32project');
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
        
        % --- Lattice menu ---------------------------
        handles.fFirstMenu   =   uimenu(...       % File menu
            'Parent', handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Electrode Pattern');
        
        handles.fPadsMenu      =   uimenu(...       % File menu
            'Parent', handles.fFirstMenu,...
            'HandleVisibility','callback', ...
            'Label','Number of pads',...
            'Callback', @fPadsCallback);
        
        function fPadsCallback(~, ~)
            % Callback function run when the Potential menu item is selected
            prompt = {'Number of pads (number of rows/columns'};
            dlg_title = 'Number of pads';
            num_lines = 1;
            def = {num2str( NumberOfPads )};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            user_entry = str2double(answer);
            if isnan(user_entry)
                errordlg('You must enter a numeric value','Bad Input','modal');
                return;
            end
            if ~isempty(user_entry)
                NumberOfPads = user_entry;
                input.NumberOfPads = NumberOfPads;
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);              
            end
        end
                
        handles.fParametersMenu      =   uimenu(...       % File menu
            'Parent', handles.fFirstMenu,...
            'HandleVisibility','callback', ...
            'Label','Parameters',...
            'Callback', @fFirstParametersMenuCallback);
        
        function fFirstParametersMenuCallback(~, ~)
            strge0 = strge;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Lattice setup', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'This is the panel to control the lattice parameters';'Press ENTER after changing each parameter'},'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            data = {...
                'Circle diameter, m',                       input.PadSize;...
                'Circle spacing, m',                        input.PadSpacing;...
                'Circle gap, m',                            input.PadGap;...
                'Number of triangles in the pad',           input.Triangles};          
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'Change Upper Plate shape',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback});
            function edit_callback(hObject, ~)
                data = get(hObject,'Data');                
                input.PadSide       = data{5};
                input.PadSpacing    = data{6};
                input.PadGap        = data{7};                
                input.Triangles     = data{8};             
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry()
                strout = clearoutput(strout);
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strge = strge0;
                delete(local);
            end
            
        end
        
        % --- Finger Menu ----------------------------        
        handles.fSecondMenu    =   uimenu(...       % File menu
            'Parent',handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Finger');
                
         handles.fObjectTypeMenu      =   uimenu(...       % File menu
            'Parent', handles.fSecondMenu,...
            'HandleVisibility','callback', ...
            'Label','Include',...
            'Callback', @fObjectTypeCallback);
                
        function fObjectTypeCallback(~,~)
            % --- RADIO BUTTONS -------------------------------------
            handles.figure2 = figure( ...
                'Tag', 'figure1', ...
                'Units', 'characters', ...
                'Position', [103.8 44.6153846153846 28.4 9.7692307692308], ...
                'Name', 'Include or not', ...
                'MenuBar', 'none', ...
                'NumberTitle', 'off', ...
                'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
                'Resize', 'on');
            
            handles.uipanel3 = uibuttongroup( ...
                'Parent', handles.figure2, ...
                'Tag', 'uipanel3', ...
                'Units', 'normalized', ...
                'Position', [0.0141823161189358 0.0334928229665072 0.965277777777778 0.961722488038277], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'Title', 'Include finger');
            
            set(handles.uipanel3,'SelectionChangeFcn',@uibuttongroup1_SelectionChangeFcn);
                        
            function uibuttongroup1_SelectionChangeFcn(~,eventdata)
                switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
                    case 'Exclude'
                        handles.object = 1;
                        strcylinder.yes = 'no';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                        strout = clearoutput(strout);
                    case 'Include'
                        handles.object = 2;
                        strcylinder.yes = 'yes';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();       
                        strout = clearoutput(strout);
                end
            end
                        
            handles.radiobutton1 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'Exclude', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.877777777777778 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Exclude');
                        
            handles.radiobutton2 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'Include', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.666666666666667 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Include');          
            
            switch handles.object
                case 1
                    set(handles.uipanel3,'Selectedobject',handles.radiobutton1);
                case 2
                    set(handles.uipanel3,'Selectedobject',handles.radiobutton2);       
            end
        end        
        
        handles.fParametersMenu      =   uimenu(...       % File menu
            'Parent', handles.fSecondMenu,...
            'HandleVisibility','callback', ...
            'Label','Parameters',...
            'Callback', @fSecondParametersMenuCallback);
        
        function fSecondParametersMenuCallback(~, ~)
            strcylinder0 = strcylinder;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Finger tip setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String',{'parameters of cylinder setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Finger present(yes/no)',strcylinder.yes;...
                'Cylinder radius,m',strcylinder.R;...
                'Cylinder height,m',strcylinder.H;...
                'Center position x,m',strcylinder.X;...
                'Center position y,m',strcylinder.Y;...
                'Center position z,m',strcylinder.Z;...
                'Rotate about x-axis, deg',strcylinder.ax;...
                'Rotate about y-axis, deg',strcylinder.ay;...
                'Rotate about z-axis, deg',strcylinder.az;...          
                'Approximate # of triangles - top face',strcylinder.Tr;...
                'Triangle size at boundaries vs. center size',strcylinder.par};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'You have chosen the cylinder',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback});
            function edit_callback(hObject, ~)
                data = get(hObject,'Data');
                if (data{13}<0.05*data{14})&(data{21}>100);
                    errordlg('The cylinder will have too many faces. Consider changing!', 'Unexpected geometry');
                    return;
                end
                strcylinder.yes       = data{12};
                strcylinder.R         = data{13};
                strcylinder.H         = data{14};
                strcylinder.X         = data{15};
                strcylinder.Y         = data{16};
                strcylinder.Z         = data{17};
                strcylinder.ax        = data{18};
                strcylinder.ay        = data{19};
                strcylinder.az        = data{20};
                strcylinder.Tr        = data{21};
                strcylinder.par       = data{22};
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strcylinder = strcylinder0;
                delete(local);
            end
            
        end        
        
        %----- The Menu of dieectric brick --------------        
        handles.fDielectricMenu   =   uimenu(...       % File menu
            'Parent', handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Dielectric');            
                
        handles.fParametersMenu      =   uimenu(...       % File menu
            'Parent', handles.fDielectricMenu,...
            'HandleVisibility','callback', ...
            'Label','Parameters',...
            'Callback', @fPreParametersCallback);
        
        function fPreParametersCallback(~,~)
            fConductingBrickCallback;
        end        
        %% Brick Case  
        function fConductingBrickCallback(~, ~)
            strbrick0 = strbrick;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Dielectric setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'Dielectric brick setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...              
                'Brick height,m',                               strbrick.H;...           
                'Center position z,m',                          strbrick.Z;...
                'Relative dielectric constant of the brick',    strbrick.epsr;...
                'Approximate # of triangles - top face',        strbrick.Tr;...
                'Triangle size at boundaries vs. center size',  strbrick.par};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'You have chosen the brick case',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback});
            function edit_callback(hObject, ~)
                data = get(hObject,'Data');                           
                strbrick.H              = data{6};         
                strbrick.Z              = data{7};          
                strbrick.epsr           = data{8};
                strbrick.Tr             = data{9};
                strbrick.par            = data{10};                
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strbrick = strbrick0;
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
                'Plot E-field in a plane',strop.yes;...
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
                'CellEditCallback', {@edit_callback});
            function edit_callback(hObject, ~)
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
                if handles.simulate
                    outputgraphics;
                else
                    cleaning;
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
            R_temp = R;
            gauss_temp = gauss;
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
                '# of intergration points',gauss};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'Setup Modeling parameters',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback});
            
            function edit_callback(hObject, ~)
                data = get(hObject,'Data');
                R                       = data{3};
                gauss                   = data{4};
                
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);
                delete(local);
            end
             
            function B2_Callback(~, ~, ~)
                R = R_temp;
                gauss = gauss_temp;
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
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', 'Charge Qij of each pad, pC (or Cij, pF)', 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            handles.outputtable = uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', strout.SelfCapacitance*1e12, 'Visible', 'on', 'BackgroundColor', [1 1 1],...                                
                'ColumnWidth',{50});
            
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
                        view(53, 23);
                    case 'Front'
                        set(0, 'CurrentFigure', handles.figure1);
                        view(90,0);
                    case 'Rear'
                        set(0, 'CurrentFigure', handles.figure1);
                        view(270,0);
                    case 'Top'
                        set(0, 'CurrentFigure', handles.figure1);
                        view(0,90);
                    case 'Bottom'
                        set(0, 'CurrentFigure', handles.figure1);
                        view(180,-90);
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
%             HelpPath = 'index_s34.html';
%             web(HelpPath);
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
            'String', 'This module simulates a capacitive 2D touchpad with a lattice of individual electrodes ("mutual capacitance" method). The default scheme: All electrodes are assigned voltage V1=1V. Ground plane and the human finger phantom are assigned voltage V2=0V. Induced charge Q on every electrode is reported in pC which gives “mutual capacitance” in the form C=Q(V2-V1).');
         

        % --- PUSHBUTTONS -------------------------------------
        handles.pushbutton2 = uicontrol( ...
            'Parent', handles.figure1, ...
            'Tag', 'pushbutton2', ...
            'Style', 'pushbutton', ...
            'Units', 'normalized', ...
            'Position', [0.55 0.0271317829457364 0.105413105413105 0.0542635658914729], ...
            'FontSize', 7.6666666666666, ...
            'FontUnits', 'pixels', ...
            'String', 'Reload',...
            'callback',@GeometryCallback);
        
        function GeometryCallback(~,~)
            handles.simulate = 0;
            geometry();
            cleaning;
        end
        
        handles.pushbutton3 = uicontrol( ...
            'Parent', handles.figure1, ...
            'Tag', 'pushbutton3', ...
            'Style', 'pushbutton', ...
            'Units', 'normalized', ...
            'Position', [0.7 0.0271317829457364 0.106837606837607 0.0542635658914729], ...
            'FontSize', 7.6666666666666, ...
            'FontUnits', 'pixels', ...
            'String', 'Exit',...
            'Callback',@ExitCallback);
        
        function ExitCallback(~,~)
            
            button = questdlg('Save project before closing?');
            
            switch button
                case 'Yes',
                   save('E32project');
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
            'Position', [0.4 0.0271317829457364 0.105413105413105 0.0542635658914729], ...
            'FontSize', 7.6666666666666, ...
            'FontUnits', 'pixels', ...
            'String', 'Simulate',...
            'callback',@SimulateCallback);
        
        function SimulateCallback(~,~)
            simulate();
        end
        
        handles.pushbutton4 = uicontrol( ...
            'Parent', handles.figure1, ...
            'Tag', 'pushbutton3', ...
            'Style', 'pushbutton', ...
            'Units', 'normalized', ...
            'Position', [0.2 0.0271317829457364 0.15 0.0542635658914729], ...
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
                'MenuBar', 'none', ...
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

function [P, t] = brick(L, W, H, Tr, par)
    %   Create a uniform or non-uniform mesh (ratio 1:5) for a 
    %   base brick with the length L, width W, and height H 
    %   Tri is approximate number of triangles for top/bottom plates
    %   SNM Summer 2012

    %   Generate the plate mesh first
    [P0, t0] = plate(L, W, Tr, par, 0, 0, 1);
    P0(:, 3) = [];
    
    %   Apply extrusion
    [P, t] = quadric(H, P0, t0); 
end

function [P, t] = circle(R, Tr, par, nx, ny, nz);
    %   Create a uniform or nonuniform (ratio 1:5) mesh for a base circle 
    %   of radius R, normal vector nx, ny, nz, and with approximately M 
    %   triangles
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson
    %   Adopted and simplified: SNM Summer 2012 

    h    = waitbar(0.5, 'Please wait - computing triangular mesh'); 
    
    if Tr<20
        %   Create a simple structured triangular mesh
        M = Tr;
        x = R*[cos(2*pi*[0:M]/(M+1)) 0];
        y = R*[sin(2*pi*[0:M]/(M+1)) 0];
        P(:, 1) = x';
        P(:, 2) = y';
        P(:, 3) = 0;
        t(:, 1) = [1:M+1  ]';
        t(:, 2) = [2:M+1 1]';
        t(:, 3) = (M+2);
        close(h);
        return;
    end
    
    h0    = 2.5*R/sqrt(Tr);                             %   Approximate desired edge length
    xmin  = -R; xmax = +R; ymin = -R; ymax = +R;        %   Bounding box

    dptol   = 0.001;                                    %   Termination criterion
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
        if count>2^12 break; end;
    end
    
    close(h);
    
       %   Rotate mesh
    P(:, 3) = 0;
    theta   = acos(nz)*sign(ny + eps);           %   rotation about the x-axis following the right-hand rule
    phi     = asin(nx/sqrt(nx^2+ny^2+eps^2));    %   rotation about the z-axis (the same way) 
    theta   = theta/pi*180;
    phi     = phi/pi*180; 
    P       = rotatex(P, theta);
    P       = rotatez(P, phi);
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
        xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); grid on; axis equal; view(-30, 30);
    else
        disp('Electric potential is nearly constant in this plane - contour plot should not be used');       
        xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); grid on; axis equal; view(-30, 30);
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
function [P, t] = cylinder(R, H, Tr, par);
    %   Create a uniform or non-uniform (ratio 1:5) mesh for a base 
    %   cylinder with the radius R, height h, and orientation along 
    %   the z-axis. Tri is approximate number of triangles
    %   SNM Summer 2012

     %   Generate the circle mesh first
    [P0, t0] = circle(R, Tr, par, 0, 0, 1);
    P0(:, 3) = [];
    
    %   Apply extrusion
    [P, t] = quadric(H, P0, t0); 
 
end

function d = drectangle(P, L, W, x, y)
    %   Compute signed distance function for a rectangle with length L and
    %   width W
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson

    x1 = -L/2+x;  x2 = +L/2+x;
    y1 = -W/2+y;  y2 = +W/2+y;
    d1 = y1 - P(:, 2);
    d2 =-y2 + P(:, 2);
    d3 = x1 - P(:, 1);
    d4 =-x2 + P(:, 1);

    d5 = sqrt(d1.^2 + d3.^2);
    d6 = sqrt(d1.^2 + d4.^2);
    d7 = sqrt(d2.^2 + d3.^2);
    d8 = sqrt(d2.^2 + d4.^2);

    d = -min(min(min(-d1, -d2), -d3), -d4);

    ix = d1>0 & d3>0;
    d(ix) = d5(ix);
    ix = d1>0 & d4>0;
    d(ix) = d6(ix);
    ix = d2>0 & d3>0;
    d(ix) = d7(ix);
    ix = d2>0 & d4>0;
    d(ix) = d8(ix);
    end

function d = dcircle(P, R)
    %   Compute signed distance function for a circle with radius R 
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson (PhD thesis, pp. 38-39)

    d = sqrt(P(:, 1).^2 + P(:, 2).^2) - R;
end

function d = dsphere(P, R)
    %   Compute signed distance function for a sphere with radius R 
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson (PhD thesis, pp. 38-39)

    d = sqrt(P(:, 1).^2 + P(:, 2).^2 + P(:, 3).^2) - R;
end

function E = efield(Points, c, P, t, centers, areas, normals, sizes, R, eps, msg)
%   Electric field from free/polarization charges
%   Vectorized for an arbitrary number of observation points
%   SNM Summer 2012
    if msg ==1
        h    = waitbar(0, 'Please wait - computing E-fields in a plane');
    else
        h    = waitbar(0, 'Please wait - computing E-field at the obs. point');
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

function [ ] = graphics_S32(io, strop, stroc, strei, P, t, c, Area, strsc, strbrick);
    %%   Program I/O graphics
    %    variable io is zero for the input graphics and one for the output graphics 

    %%  Prepare common variables
    tM = t((t(:, 4)>0)&(t(:, 4)<1000), :);
    tF = t(t(:, 4)==1000, :);
    tD = t(t(:, 4)==0,  :);
    
    cM = c((t(:, 4)>0)&(t(:, 4)<1000), :);
    cF = c(t(:, 4)==1000, :);
    cD = c(t(:, 4)==0,  :);
    
    Ptemp = P'; ttemp = t'; ind = size(t, 1);  IND = 10000;
    
    ctempM = cM'; ttempM = tM'; 
    XM = reshape(Ptemp(1, ttempM(1:3, :)),[3, size(ttempM, 2)]);
    YM = reshape(Ptemp(2, ttempM(1:3, :)),[3, size(ttempM, 2)]);
    ZM = reshape(Ptemp(3, ttempM(1:3, :)),[3, size(ttempM, 2)]);  
    
    ctempD = cD'; ttempD = tD'; 
    XD = reshape(Ptemp(1, ttempD(1:3, :)),[3, size(ttempD, 2)]);
    YD = reshape(Ptemp(2, ttempD(1:3, :)),[3, size(ttempD, 2)]);
    ZD = reshape(Ptemp(3, ttempD(1:3, :)),[3, size(ttempD, 2)]);  
    
    ctempF = cF'; ttempF = tF'; 
    XF = reshape(Ptemp(1, ttempF(1:3, :)),[3, size(ttempF, 2)]);
    YF = reshape(Ptemp(2, ttempF(1:3, :)),[3, size(ttempF, 2)]);
    ZF = reshape(Ptemp(3, ttempF(1:3, :)),[3, size(ttempF, 2)]);  
    
    
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
        colormap jet;        
        patch('Vertices', P, 'Faces', tD(:, 1:3), 'FaceColor',  'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');         
        patch('Vertices', P, 'Faces', tM(:, 1:3), 'FaceColor',  'flat', 'FaceVertexCData', tM(:, 4), 'EdgeColor', 'k'); 
        if ~isempty(tF) patch('Vertices', P, 'Faces', tF(:, 1:3), 'FaceColor',  [1 0.75 0.65], 'FaceAlpha', 1.0, 'EdgeColor', 'k'); end;            
        box(strbrick.L, strbrick.W, strbrick.H, strbrick.X, strbrick.Y, strbrick.Z, 'y', 0, 1, 'k')
        xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); 
        title('Problem geometry: columns - from left to right; rows - from top to bottom');
    end
    
    %%   Output graphics (separate window)
    if io == 1
        
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
        CF = cv(ttempF(1:3, :));   
     
        if (~strcmp(strop.yes, 'yes'))&(~strcmp(strop.potential, 'yes'))
            ec = 'none';
            patch(XD, YD, ZD, CD, 'FaceAlpha', 0.1, 'EdgeColor', ec, 'FaceLighting', 'none');
            patch(XM, YM, ZM, CM, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');
            patch(XF, YF, ZF, CF, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');
            xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); colorbar;
            title('Solution: Surface charge density on all objects in C/m^2')
        end        
        if strcmp(strop.yes, 'yes')|strcmp(strop.potential, 'yes')
            ec = 'none';
            patch(XD, YD, ZD, CD, 'FaceAlpha', 0.1, 'EdgeColor', ec, 'FaceLighting', 'none');
            patch(XM, YM, ZM, CM, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');
            patch(XF, YF, ZF, CF, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');
            xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); colorbar;
            title('Solution: Surface charge density in C/m^2 and the field in the observation plane'); 
        end
        if strcmp(strop.yes, 'yes')
            %%   Draw electric field (final) in the observation plane
            strop.E = strop.E + 1e-12*max(max(strop.E))*rand(size(strop.E)); %   randomize a bit
            normalize = sqrt(sum(strop.E.*strop.E, 2));
            normalize = max(normalize);
            normalize = normalize/strop.arrow;
            fv = cone(strop.Points(:, 1), strop.Points(:, 2), strop.Points(:, 3), strop.E(:, 1), strop.E(:, 2), strop.E(:, 3), 36, 0.25, 1, normalize/stroc.size);            
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
    view(0, 90);
end

function h = hcircle(P, R, par)
    %   Compute element size function for a circle with radius R 
    %   0<par<1
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson (PhD thesis, pp. 38-39)
    h = 1/par - (1/par-1)*sqrt(sum(P.^2, 2))/R;     %   Relative edge length
end

function h = hrectangle(P, L, W, par)
    %   Compute element size function for a rectangle
    %   0<par<1
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson

    d = drectangle(P, L, W, 0, 0);
    h = 1/par - (1/par-1)*(L/2 + d).*(W/2 + d)/(L*W/4);    %   Relative edge length
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

function [P, t] = plate(L, W, Tr, par, nx, ny, nz)
    %   Create a uniform or nonuniform (ratio 1:5) mesh  for a base plate with 
    %   length L and width W, normal vector nx, ny, nz, and with approximately 
    %   Tr triangles
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson
    %   Adopted and simplified: SNM Summer 2012 
    warning off;
    
    h    = waitbar(0.5, 'Please wait - computing triangular mesh');
    
    h0    = 1.40*sqrt(L*W)/sqrt(Tr);                    %   Approximate desired edge length
    xmin  = -L/2; xmax = +L/2; ymin = -W/2; ymax = +W/2;%   Bounding box

    dptol   = 0.001;                                    %   Termination criterion
    ttol    = 0.1;                                      %   Criterion for retriangulation  
    Fscale  = 1.2;                                      %   Inner "internal pressure" (original)
    if Tr<=200 & par<0.3
        Fscale = 1.10;
    end
    deltat  = 0.2;                                      %   Ratio between force and movement
    geps    = 0.001*h0;                                 %   Relative boundary tolerance (original)
    if par ==1
        geps    = 0.1*h0;                               %   Relative boundary tolerance (relaxed here)
    end  
    deps    = sqrt(eps)*h0;                             %   Accuracy of gradient calculation
    
    %   Create initial distribution in bounding box 
    %   Equilateral triangles - faster convergence
    [x, y] = meshgrid(xmin:h0:xmax, ymin:h0*sqrt(3)/2:ymax);
    x(2:2:end, :) = x(2:2:end, :) + h0/2;                 %   Shift even rows
    P = [x(:), y(:)];                                     %   List of node coordinates

    %   Remove vertices outside the region
    dist = drectangle(P, L, W, 0, 0);                         %  A column of signed distances
    P    = P(dist<geps, :);                             %  Keep only d<0 (geps) points  
    N    = size(P, 1);                                  %  Number of vertices N
    Pold = inf;                                         %  For first iteration

    count   = 0;    
    while 1      
        %   Retriangulation by the Delaunay algorithm
        if max(sqrt(sum((P-Pold).^2,2))/h0)>ttol        %   Any large movement?
            Pold    = P;                                %   Save current positions
            dt      = DelaunayTri(P);                   %   2D Delaunay triangulation
            t       = dt.Triangulation;                 %   List of triangles
            ic      = incenters(dt);                    %   Centroids of triangles
            dist    = drectangle(ic, L, W, 0, 0);             %   A column of signed distances
            t       = t(dist<-geps, :);                 %   Keep interior trianges
            bars    = edges(TriRep(t, P));              %   All sorted mesh edges [M, 2]
            M       = size(bars, 1);                    %   Length of the array of edges        
        end
        count = count + 1;     

        %   Move mesh points based on edge lengths L and forces F
        barvec  = P(bars(:, 1),:)-P(bars(:, 2),:);              %   List of edge vectors [M, 2]
        L1      = sqrt(sum(barvec.^2,2));                       %   Edge lengths [M, 1]       
        barsc   = 0.5*(P(bars(:, 1), :) + P(bars(:, 2), :));    %   Edge centers [M, 2]
        hbars   = hrectangle(barsc, L, W, par);                 %   Element sizes for edge centers [M, 1]
        hbars   = hbars/sqrt(sum(hbars.^2)/M);                  %   Element sizes normalized by its average [M, 1] 
        L0      = Fscale*sqrt(sum(L1.^2)/M)*hbars;              %   Desired (non-uniform) edge lengths [M, 1]        
        F       = max(L0-L1, 0);                                %   Edge forces [M, 1]
        Fvec    = repmat(F./L1, 1, 2).*barvec;                  %   Edge forces (x,y components) [M, 2]

        Ftot = full(sparse(bars(:, [1,1,2,2]), ones(size(F))*[1,2,1,2], [Fvec, -Fvec], N, 2));    
        P = P + deltat*Ftot;                            %   Update node positions

        %   Bring outside vertices back to the boundary
        dist  = drectangle(P, L, W, 0, 0);                    %   Find distances
        ind   = dist>0;                                 %   Find vertices outside (d>0)
        dgradx = (drectangle([P(ind, 1)+deps, P(ind, 2)], L, W, 0, 0)-dist(ind))/deps;   %   Numerical
        dgrady = (drectangle([P(ind, 1), P(ind, 2)+deps], L, W, 0, 0)-dist(ind))/deps;   %   gradient
        dgrad2 = dgradx.^2 + dgrady.^2;
        P(ind, :) = P(ind, :) - ...
            [dist(ind).*dgradx./dgrad2, dist(ind).*dgrady./dgrad2];     %   Projection

        %   Termination criterion: All interior nodes move less than dptol (scaled)
        factor = max(sqrt(sum(deltat*Ftot(dist<-geps,:).^2, 2))/h0);
    
        if factor<dptol, break; end
        if count>2^13 break; end;
    end
    close(h);    
    
    %   Rotate mesh
    P(:, 3) = 0;
    theta   = acos(nz)*sign(ny + eps);           %   rotation about the x-axis following the right-hand rule
    phi     = asin(nx/sqrt(nx^2+ny^2+eps^2));    %   rotation about the z-axis (the same way) 
    theta   = theta/pi*180;
    phi     = phi/pi*180; 
    P       = rotatex(P, theta);
    P       = rotatez(P, phi);
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

function [P, t] = quadric(H, P0, t0)
    %   Create a cylinder (quadric) of height H from any planar mesh using
    %   extrusion
    %   SNM Summer 2012

    %   Find boundary vertices
    bars0     = freeBoundary(TriRep(t0, P0));
    barvec    = P0(bars0(:, 1), :)-P0(bars0(:, 2),: );    %   List of edge vectors [M, 2]
    Lavg      = mean(sqrt(sum(barvec.^2, 2)));            %   Average edge length
    Pb = P0(bars0', :); Pb(2:2:end, :) = [];              %   Boundary vertices
    pb = size(Pb, 1);   barsb = [1:pb; [2:pb 1]]';        %   Boundary connectivity
   
    %  Determine z-divisions   
    NN      = ceil(H/Lavg);                               %   Number of intervals along the height
    hsize   = H/NN;                                       %   Size of the z-divisions
    
    %   Replicate the cap mesh and its boundary; establish the connectivity
    P = P0; p0 = size(P0, 1); P(:, 3) = - H/2;  %   Bottom cap - vertices                                                                              
    t = t0;                                     %   Bottom cap - faces
    
    if NN==1              
        %   Connectivity: bottom to top
        t1(:, 1:2) = bars0;                     %   Lower nodes        
        t1(:, 3)   = bars0(:, 1) + p0;          %   Upper node
        t2(:, 1:2) = bars0       + p0;          %   Upper nodes
        t2(:, 3)   = bars0(:, 2);               %   Lower node
        t = [t' t1' t2']';
    else
        %   Connectivity: bottom to layer
        P(end+1:end+pb, 1:2) = Pb;
        P(end+1-pb:end, 3) = -H/2 + hsize;
        t1(:, 1:2) = bars0;                     %   Lower nodes        
        t1(:, 3)   = barsb(:, 1) + p0;          %   Upper node
        t2(:, 1:2) = barsb       + p0;          %   Upper nodes
        t2(:, 3)   = bars0(:, 2);               %   Lower node
        t = [t' t1' t2']';
        %   Connectivity: layer to layer
        for m =2:NN-1                                 
            P(end+1:end+pb, 1:2) = Pb;
            P(end+1-pb:end, 3) = -H/2 + hsize*m;
            t1(:, 1:2) = barsb       + p0 + pb*(m-2);   %   Lower nodes        
            t1(:, 3)   = barsb(:, 1) + p0 + pb*(m-1);   %   Upper node
            t2(:, 1:2) = barsb       + p0 + pb*(m-1);   %   Upper nodes
            t2(:, 3)   = barsb(:, 2) + p0 + pb*(m-2);   %   Lower node
            t = [t' t1' t2']';
        end
        %   Connectivity: layer to top
        t1(:, 1:2) = barsb       + p0 + pb*(NN-2);   %   Lower nodes        
        t1(:, 3)   = bars0(:, 1) + p0 + (NN-1)*pb;   %   Upper node
        t2(:, 1:2) = bars0       + p0 + (NN-1)*pb;   %   Upper nodes
        t2(:, 3)   = barsb(:, 2) + p0 + pb*(NN-2);   %   Lower node
        t = [t' t1' t2']';
    end
    %   Add the top cap
    P0(:, 3) = H/2; P = [P' P0']';
    t = [t' t0'+(NN-1)*pb+p0]';
end

function [Pnew] = rotatex(P, xangle)
    %%   Rotate the mesh about the x-axis
    %   The local coordinate system is used with the origin at the center of gravity
    %   SNM Summer 2012
    anglex = xangle/180*pi;
    %   Rotation about the x-axis
    LCY  = mean(P(:, 2));
    LCZ  = mean(P(:, 3));
    Pnew(:, 1) = +P(:, 1);
    Pnew(:, 2) = +(P(:, 2) - LCY)*cos(anglex) - (P(:, 3) - LCZ)*sin(anglex);
    Pnew(:, 3) = +(P(:, 2) - LCY)*sin(anglex) + (P(:, 3) - LCZ)*cos(anglex);
    Pnew(:, 2) = Pnew(:, 2) + LCY;
    Pnew(:, 3) = Pnew(:, 3) + LCZ;
end

function [Pnew] = rotatey(P, yangle)
    %%   Rotate the mesh about the y-axis
    %   The local coordinate system is used with the origin at the center of gravity
    %   SNM Summer 2012
    angley = yangle/180*pi;
    %   Rotation about the y-axis
    LCX  = mean(P(:, 1));
    LCZ  = mean(P(:, 3));    
    Pnew(:, 1) = +(P(:, 1) - LCX)*cos(angley) + (P(:, 3) - LCZ)*sin(angley);
    Pnew(:, 2) = +P(:, 2);
    Pnew(:, 3) = -(P(:, 1) - LCX)*sin(angley) + (P(:, 3) - LCZ)*cos(angley);
    Pnew(:, 1) = Pnew(:, 1) + LCX;
    Pnew(:, 3) = Pnew(:, 3) + LCZ;
end

function [Pnew] = rotatez(P, zangle)
    %%   Rotate the mesh about the z-axis
    %   The local coordinate system is used with the origin at the center of gravity
    %   SNM Summer 2012
    anglez = zangle/180*pi;
    %   Rotation about the z-axis
    LCX  = mean(P(:, 1));
    LCY  = mean(P(:, 2));       
    Pnew(:, 1)  = +(P(:, 1) - LCX)*cos(anglez) - (P(:, 2) - LCY)*sin(anglez);
    Pnew(:, 2)  = +(P(:, 1) - LCX)*sin(anglez) + (P(:, 2) - LCY)*cos(anglez);
    Pnew(:, 3)  = +P(:, 3); 
    Pnew(:, 1) = Pnew(:, 1) + LCX;
    Pnew(:, 2) = Pnew(:, 2) + LCY;
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


function [P, t] = sphere(R, Tr);
    %   Create a uniform mesh for a base sphere of radius R with approximately M triangles
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson
    %   Adopted and simplified: SNM Summer 2012 
   
    h    = waitbar(0.5, 'Please wait - computing triangular mesh');  
    
    h0    = 5.5*R/sqrt(Tr);                             %   Approximate desired edge length
    xmin  = -R; xmax = +R; ymin = -R; ymax = +R;        %   3D Bounding box
    zmin  = -R; zmax = +R; 

    dptol   = 0.005;                                    %   Termination criterion
    ttol    = 0.1;                                      %   Criterion for retriangulation
    Fscale  = 1.2;                                      %   Inner "internal pressure"
    deltat  = 0.2;                                      %   Ratio between force and movement
    geps    = 0.001*h0;                                 %   Relative boundary tolerance
    deps    = sqrt(eps)*h0;                             %   Accuracy of gradient calculation

    %   Create initial uniform distribution in bounding box (MATLAB isosurface from grid)
    [x, y, z] = ndgrid(xmin:h0:xmax, ymin:h0:ymax, zmin:h0:zmax);
    pv        = isosurface(x, y, z, reshape(dsphere([x(:), y(:), z(:)], R), size(x)), 0);
    P         = pv.vertices;
    N    = size(P, 1);                                  %  Number of vertices N
    Pold = inf;                                         %  For first iteration

    count = 0;
    while 1
        %   Retriangulation by the Delaunay algorithm
        if max(sqrt(sum((P-Pold).^2,2))/h0)>ttol        %   Any large movement?
            Pold    = P;                                %   Save current positions
            DT  = DelaunayTri(P);                       %   3D Delaunay triangulation
            t   = freeBoundary(DT);                     %   3D Triangular mesh
            bars    = edges(TriRep(t, P));              %   All sorted mesh edges [M, 2]
            M       = size(bars, 1);                    %   Length of the array of edges
            count = count + 1;  
        end

        %   Move mesh points based on edge lengths L1 and forces F
        barvec    = P(bars(:, 1),:)-P(bars(:, 2),:);            %   List of edge vectors [M, 2]
        L1         = sqrt(sum(barvec.^2,2));                    %   Edge lengths [M, 1]
        L0        = Fscale*sqrt(sum(L1.^2)/M)*ones(M, 1);       %   Desired (average) edge lengths [M, 1]
        F         = max(L0-L1, 0);                              %   Edge forces [M, 1]
        Fvec      = repmat(F./L1, 1, 3).*barvec;                %   Edge forces (x,y,z components) [M, 3]

        Ftot = full(sparse(bars(:,[1,1,1,2,2,2]), ones(size(F))*[1,2,3,1,2,3], [Fvec,-Fvec], N, 3));
        P = P + deltat*Ftot;                            %   Update node positions

        %   Bring all vertices back to the boundary
        dist   = dsphere(P, R);                         %   Find distances
        dgradx = (dsphere([P(:, 1)+deps, P(:, 2), P(:, 3)], R)-dist)/deps;   %   Numerical
        dgrady = (dsphere([P(:, 1), P(:, 2)+deps, P(:, 3)], R)-dist)/deps;   %   gradient
        dgradz = (dsphere([P(:, 1), P(:, 2), P(:, 3)+deps], R)-dist)/deps;   %   gradient
        dgrad2 = dgradx.^2 + dgrady.^2 + dgradz.^2;
        P = P - [dist.*dgradx./dgrad2, dist.*dgrady./dgrad2, dist.*dgradz./dgrad2];     
                                                                               %   Projection
        %   Termination criterion: All interior nodes move less than dptol (scaled)
        factor = max(sqrt(sum((P-Pold).^2, 2))/h0);
        if factor<dptol, break; end
        if count>2^10 break; end;
    end
    if exist('h') 
        close(h);
    end
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

function [P, t] = lattice(input)
    %   SNM Fall 2012
    NumberOfPads  = input.NumberOfPads;
    PadSize       = input.PadSize;
    PadSpacing    = input.PadSpacing;
    PadGap        = input.PadGap;   
    Triangles     = input.Triangles;    %   Triangles per one row rectangle    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    if floor(NumberOfPads/2)==NumberOfPads/2
        temp                = [-(floor(NumberOfPads/2)-1):1:floor(NumberOfPads/2)]-0.5;
    else
        temp                = [-floor(NumberOfPads/2):1:floor(NumberOfPads/2)];
    end  
    strc.R = [];
    strc.x = [];
    strc.y = [];
    strc.l = [];
    strc.no = [];
    for m = 1:size(temp, 2) %   loop over rows
        strc.R   = [strc.R PadSize*ones(1, size(temp, 2))/2];
        strc.x   = [strc.x PadSpacing*temp];                    %   Rows of x (x increases from left to right)
        strc.y   = [strc.y PadSpacing*temp(size(temp, 2)-m+1)*ones(size(temp))];%   Columns of y (y increases from top to bottom)
        strc.l   = [strc.l 0.5*ones(size(temp))];                   %   Rel. edge length         - fourth column
        strc.no  = [strc.no [1:length(temp)]+(m-1)*length(temp)];       %   Domain # (0=diel.)       - fifth column       
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Outer shape
    strbrick.L      = NumberOfPads*PadSpacing +PadSpacing/2;     %   brick length, m
    strbrick.W      = strbrick.L;             %   brick width, m
    strbrick.Tr     = Triangles;              %   approximate # of triangles for top face
    
    str0.type    = 'r';
    str0.Tr      = strbrick.Tr;
    str0.L       = strbrick.L;
    str0.W       = strbrick.W;
    str0.x       = 0;
    str0.y       = 0;
    str0.N       = ceil(sqrt(str0.Tr*str0.L/(2.2*str0.W)));
    str0.l       = str0.L/(str0.N-1);       %   Edge length of the base mesh   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %   Create mesh 
    [P, ttop, strc] = combmesh(str0, strc, input);
    %   Assign domain numbers
    t = combdomain(P, ttop, strc);
    t = cutdomain(P, t, strc, input);

    Ptotal = P;
    ttotal = t; 

    %   Move whole mesh
    if floor(NumberOfPads/2)==NumberOfPads/2
        Ptotal(:, 1) = Ptotal(:, 1) - PadSpacing/2; 
    else
        Ptotal(:, 2) = Ptotal(:, 2) - PadSpacing/2; 
    end

    P = Ptotal;
    t = ttotal; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [t] = combdomain(P, t, strc)
    %   Assign domain number for every subdomain

        EdgeCenter1 = 1/2*(P(t(:, 1), :) + P(t(:, 2), :));
        EdgeCenter2 = 1/2*(P(t(:, 1), :) + P(t(:, 3), :));
        EdgeCenter3 = 1/2*(P(t(:, 2), :) + P(t(:, 3), :));    
        TriCenter = 1/3*(P(t(:, 1), :) + P(t(:, 2), :) + P(t(:, 3), :));    
        t(:, 4) = 0;
        tol = 1e-6;
        for m = 1:size(strc.R, 2)
            temp1        = (EdgeCenter1(:, 1) - strc.x(m)).^2 + ...
                           (EdgeCenter1(:, 2) - strc.y(m)).^2 ;
            temp2        = (EdgeCenter2(:, 1) - strc.x(m)).^2 + ...
                           (EdgeCenter2(:, 2) - strc.y(m)).^2 ;
            temp3        = (EdgeCenter3(:, 1) - strc.x(m)).^2 + ...
                           (EdgeCenter3(:, 2) - strc.y(m)).^2 ;
            temp = (temp1<strc.R(m)^2)&(temp2<strc.R(m)^2)&(temp3<strc.R(m)^2);
            t(temp, 4) = strc.no(m); 
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [t] = cutdomain(P, t, strc, input)
        EdgeCenter1 = 1/2*(P(t(:, 1), :) + P(t(:, 2), :));
        EdgeCenter2 = 1/2*(P(t(:, 1), :) + P(t(:, 3), :));
        EdgeCenter3 = 1/2*(P(t(:, 2), :) + P(t(:, 3), :));    
        TriCenter = 1/3*(P(t(:, 1), :) + P(t(:, 2), :) + P(t(:, 3), :));    
        
        tol = 1e-6;
        TEMP = [];
        for m = 1:size(strc.R, 2)
            temp1        = (EdgeCenter1(:, 1) - strc.x(m)).^2 + ...
                           (EdgeCenter1(:, 2) - strc.y(m)).^2 ;
            temp2        = (EdgeCenter2(:, 1) - strc.x(m)).^2 + ...
                           (EdgeCenter2(:, 2) - strc.y(m)).^2 ;
            temp3        = (EdgeCenter3(:, 1) - strc.x(m)).^2 + ...
                           (EdgeCenter3(:, 2) - strc.y(m)).^2 ;
            RR = strc.R(m)+input.PadGap;
            R = strc.R(m);
            temp1 = ((temp1<RR^2)&(temp2<RR^2))|((temp2<RR^2)&(temp3<RR^2))|((temp1<RR^2)&(temp3<RR^2));
            temp2 = ((temp1>R^2)&(temp2>R^2))|((temp2>R^2)&(temp3>R^2))|((temp1>R^2)&(temp3>R^2));
            TEMP = [TEMP; find(temp1&temp2)];            
        end        
        t(TEMP, :) = [];
    end

    function [P, t, strc] = combmesh(str0, strc, input)
        %   Create a mesh for a planar base with subdomains   
        %   str0    - structure of a main boundary
        %   str1    - structure of internal boundaries
        %   SNM Summer 2012

        %%   Outer (base) boundary - outputs P, P0, edges0        
        [P, P0, edges0] = rectangle0(str0);
        free   = size(P, 1);  %   # of non-boundary vertices

        %%  All internal boundaries - outputs P1, edges1, updates str1, and removes P0, edges0 if necessary    
        [Pr, edgesr, strc, P0, edges0] = combcircle(str0, strc, P0, edges0);
        strc.R = strc.R + input.PadGap;
        [PrO, edgesrO, strc, P0, edges0] = combcircle(str0, strc, P0, edges0);
        strc.R = strc.R - input.PadGap;
        
        %%  Combine meshes (the boundary nodes are at the end here)   
        P      = [P; P0; Pr; PrO];                                %   all vertices
        edges  = [edges0+free; edgesr+free+size(P0, 1); edgesrO+free+size(P0, 1)+size(Pr, 1)];     %   all boundary edges

        %   Create and view the initial structured mesh (with possibly
        %   intersecting boundaries)
        warning off; tol = 1e-5; 
        dt      = DelaunayTri(P, edges);
        t       = dt.Triangulation;
        P       = dt.X;                             %   Adding vertices for intersecting boundaries
        edges   = dt.Constraints;                   %   Updating boundary constraints   
        dt      = DelaunayTri(P, edges);            %   Redoing triangulation
        t       = dt.Triangulation;
        q       = simpqual(P, t);
        t(q<tol, :) = [];                           %   Removing triangles of zero quality
        %   Perform Laplacian smoothing
        iter = 20;
        quality = zeros(1, iter);    
        for m = 1:iter
            P = laplace0(P, free, dt);
            quality(m) = min(simpqual(P, t));
            dt  = DelaunayTri(P, edges);
            t   =  dt.Triangulation;
            q   = simpqual(P, t);
            t(q<tol, :) = [];                        %   Removing triangles of zero quality
            if m>1&(quality(m)/quality(m-1)<=1+tol)
                break;
            end
        end
    end

    function [P1, edges1, str1, P0, edges0] = combrectangle(str0, str1, P0, edges0)
    %   This function:
    %   Accumulates boundary nodes and boundary connectivity for all internal
    %   rectangular boundaries
    %   Deletes wrong boundaries and updates the structure str1 for inner boundaries 
    %   if necessary
    %   Deletes the external boundary and replaces it by an internal one if
    %   they coincide
    %   SNM Summer 2012

        edges1 = [];
        P1     = []; 
        tol    = 1e-6;
        M = size(str1.L, 1);
        for m = 1:M
            L = str1.L(m);
            W = str1.W(m);
            x0 = str1.x(m);
            y0 = str1.y(m);
            l = str0.l*str1.l(m);        
            N = ceil(L/l) + 1;      
            %   Add boundary nodes (clockwise)
            sizex   = L/(N-1); 
            sizey   = sizex*sqrt(3)/2; 
            x = linspace(-L/2, L/2, max(ceil(L/sizex), 2)) + x0;
            y = linspace(-W/2, W/2, max(ceil(W/sizey), 2)) + y0;
            temp1 = [(-L/2+x0)*ones(1, length(y));   y                ];
            temp2 = [x(2:end-1);    (+W/2+y0)*ones(1, length(x)-2)    ];
            temp3 = [(+L/2+x0)*ones(1, length(y));   y(end:-1:1)      ];
            temp4 = [x(end-1:-1:2);     (-W/2+y0)*ones(1, length(x)-2)]; 
            Pm    = [temp1 temp2 temp3 temp4]';
            if m<M
                Pm    = rotatez(Pm, 45);
            end
            Im    = size(Pm, 1);
            %   Add boundary edges (connectivity)
            edgesm   = [[1:Im]; [2:Im 1]]';

            if strcmp(str0.type, 'c')                           %   External bounday is a circle
                temp    = str0.R^2 - dot(Pm, Pm, 2)<-tol*str0.R^2;
                if sum(temp)                                    %   A part of the internal boundary is outside the substrate    
                    str1.no(m) = 0;                             %   Domain number is that of substrate   
                    str1.c(m)  = 'y';                           %   Color is that of substrate 
                    h = warndlg('An internal boundary is outside the substrate - boundary deleted', 'Wrong geometry', 'modal');  
                    uiwait(h);
                    Pm = []; edgesm = [];
                end
            end

           if strcmp(str0.type, 'r')                            %   External bounday is a rectangle
                d = drectangle(Pm, str0.L, str0.W, 0, 0);       %   Signed distance for the rectangle
                temp    = abs(d)<+tol*max(str0.L, str0.W);     
                if sum(temp)==size(Pm, 1)                       %    All internal boundary nodes are on the base boundary
                    P0 = []; edges0 = [];
                end   
                temp    = d>+tol*max(str0.L, str0.W);
                if sum(temp)                                    %   A part of the internal boundary is outside the substrate    
                    str1.no(m) = 0;                             %   Domain number is that of substrate   
                    str1.c(m)  = 'y';                           %   Color is that of substrate 
                    h = warndlg('An internal boundary is outside the substrate - boundary deleted', 'Wrong geometry', 'modal');  
                    uiwait(h);
                    Pm = []; edgesm = [];
                end
           end

           P1       = [P1; Pm];
           edges1   = [edges1; edgesm+size(edges1, 1)];
        end
    end

    function [P1, edges1, str1, P0, edges0] = combcircle(str0, str1, P0, edges0)
    %   This function:
    %   Accumulates boundary nodes and boundary connectivity for all internal
    %   circular boundaries; Deletes wrong boundaries and updates the structure str1 for inner boundaries 
    %   if necessary; Deletes the external boundary and replaces it by an internal one if they coincide
    
        edges1 = [];
        P1     = []; 
        tol    = 1e-6;
        for m = 1:size(str1.R, 2)
            R = str1.R(m);
            x = str1.x(m);
            y = str1.y(m);
            l = str0.l*str1.l(m);
            N = round(2*pi*R/l);
            % Add boundary nodes
            N = max(N, 4);
            phi     = [0:2*pi/N:2*pi*(N-1)/N];
            bx      = R*cos(phi) + x; 
            by      = R*sin(phi) + y;
            Im      = size(bx, 2);
            Pm      = [bx; by]'; 
            edgesm  = [[1:Im]; [2:Im 1]]';

            if strcmp(str0.type, 'c')                           %   External bounday is a circle
                temp    = str0.R^2 - dot(Pm, Pm, 2)<+tol*str0.R^2;
                if sum(temp)==size(Pm, 1)                       %    All internal boundary nodes are on the base boundary
                    P0 = []; edges0 = [];
                end   
                temp    = str0.R^2 - dot(Pm, Pm, 2)<-tol*R^2;
                if sum(temp)                                    %   A part of the internal boundary is outside the substrate    
                    str1.no(m) = 0;                             %   Domain number is that of substrate   
                    str1.c(m)  = 'y';                           %   Color is that of substrate 
                    h = warndlg('An internal boundary is outside the substrate - boundary deleted', 'Wrong geometry', 'modal');  
                    uiwait(h);
                    Pm = []; edgesm = [];
                end
            end

           if strcmp(str0.type, 'r')                            %   External bounday is a rectangle
                d = drectangle(Pm, str0.L, str0.W, 0, 0);       %   Signed distance for the rectangle
                temp    = d>+tol*max(str0.L, str0.W);
                if sum(temp)                                    %   A part of the internal boundary is outside the substrate    
                    str1.no(m) = 0;                             %   Domain number is that of substrate   
                    str1.c(m)  = 'y';                           %   Color is that of substrate 
                    h = warndlg('An internal boundary is outside the substrate - boundary deleted', 'Wrong geometry', 'modal');  
                    uiwait(h);
                    Pm = []; edgesm = [];
                end
           end
           P1       = [P1; Pm];
           edges1   = [edges1; edgesm+size(edges1, 1)];
        end
    end


    function d = drectangle(P, L, W, x, y)
        %   Compute signed distance function for a rectangle with length L and
        %   width W
        %   Original code: DISTMESH 2004-2012 Per-Olof Persson

        x1 = -L/2+x;  x2 = +L/2+x;
        y1 = -W/2+y;  y2 = +W/2+y;
        d1 = y1 - P(:, 2);
        d2 =-y2 + P(:, 2);
        d3 = x1 - P(:, 1);
        d4 =-x2 + P(:, 1);

        d5 = sqrt(d1.^2 + d3.^2);
        d6 = sqrt(d1.^2 + d4.^2);
        d7 = sqrt(d2.^2 + d3.^2);
        d8 = sqrt(d2.^2 + d4.^2);

        d = -min(min(min(-d1, -d2), -d3), -d4);

        ix = d1>0 & d3>0;
        d(ix) = d5(ix);
        ix = d1>0 & d4>0;
        d(ix) = d6(ix);
        ix = d2>0 & d3>0;
        d(ix) = d7(ix);
        ix = d2>0 & d4>0;
        d(ix) = d8(ix);
    end

    function [Pnew] = laplace0(P, free, dt);
    %   Iterative grid smoothing based on existing Delaunay connectivity 
    %   P - array of vertices
    %   free - only first FREE vertices will be moved
    %   Centroid Voronoi Tessellation (CVT) smoothing
    %   SNM Summer 2012

        si      = vertexAttachments(dt);    %   Triangles attached to every vertex (cell array)
        ic      = incenters(dt);            %   Centroids of attached triangles

        Pnew = P;        
        for m = 1:free
            icenters = ic(si{m}', :);
            Pnew(m, :) = sum(icenters, 1)/size(icenters, 1);
        end       
    end

    function [P, P0, edges0] = rectangle0(str0)
        %   Outputs inner nodes P, boundary nodes P0, and boundary connectivity
        %   edges0  for a rectangle
        %   SNM Summer 2012

        %   Create a lattice for equilateral triangles within the enclosing rectangle
        L = str0.L; W = str0.W; N = str0.N;
        sizex   = L/(N-1); 
        sizey   = sizex*sqrt(3)/2; 
        x       = -L/2+sizex/2:sizex:L/2-sizex/2; 
        y       = -W/2+sizey/2:sizey:W/2-sizey/2;
        [X, Y]  = meshgrid(x, y);
        X(1:2:end, :) = X(1:2:end, :) + sizex/4;     %   Shift odd rows
        X(2:2:end, :) = X(2:2:end, :) - sizex/4;     %   Shift even rows
        P = [X(:), Y(:)];                            %   List of inner vertices 

        %   Add boundary nodes (clockwise)
        x = linspace(-L/2, L/2, max(ceil(L/sizex), 2));    
        y = linspace(-W/2, W/2, max(ceil(W/sizey), 2));
        temp1 = [-L/2*ones(1, length(y));   y                ];
        temp2 = [x(2:end-1);    +W/2*ones(1, length(x)-2)    ];
        temp3 = [+L/2*ones(1, length(y));   y(end:-1:1)      ];
        temp4 = [x(end-1:-1:2);     -W/2*ones(1, length(x)-2)]; 
        P0    = [temp1 temp2 temp3 temp4]';
        I0    = size(P0, 1);
        %   Add boundary edges (connectivity)
        edges0   = [[1:I0]; [2:I0 1]]';
    end

    function [Pnew] = rotatez(P, zangle)
        %%   Rotate the mesh about the z-axis
        %   The local coordinate system is used with the origin at the center of gravity
        %   SNM Summer 2012
        anglez = zangle/180*pi;
        %   Rotation about the z-axis
        LCX  = mean(P(:, 1));
        LCY  = mean(P(:, 2));       
        Pnew(:, 1)  = +(P(:, 1) - LCX)*cos(anglez) - (P(:, 2) - LCY)*sin(anglez);
        Pnew(:, 2)  = +(P(:, 1) - LCX)*sin(anglez) + (P(:, 2) - LCY)*cos(anglez);   
        Pnew(:, 1) = Pnew(:, 1) + LCX;
        Pnew(:, 2) = Pnew(:, 2) + LCY;
    end

    function [Pnew] = rotatez0(P, zangle)
        %%   Rotate the mesh about the z-axis
        %   The local coordinate system is used with the origin at the center of gravity
        %   SNM Summer 2012
        anglez = zangle/180*pi;
        %   Rotation about the z-axis            
        Pnew(:, 1)  = +P(:, 1)*cos(anglez) - P(:, 2)*sin(anglez);
        Pnew(:, 2)  = +P(:, 1)*sin(anglez) + P(:, 2)*cos(anglez);   
        Pnew(:, 1) = Pnew(:, 1);
        Pnew(:, 2) = Pnew(:, 2);
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
end

