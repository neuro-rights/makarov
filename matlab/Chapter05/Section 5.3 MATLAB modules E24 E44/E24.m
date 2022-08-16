function fig_hdl = E24
%   SYNTAX
%   E24
%   DESCRIPTION
%   This module is an accurate 3D electrostatic MoM solution for
%   multi-conductor structures. Three objects may be placed inside a
%   capacitor: a conducting brick, a sphere, or a cylinder. The capacitor
%   itself may be formed by objects of arbitrary shape. The module also
%   computes charge distributions on conductors, electric field
%   distribution in a plane, and electric field at a point of interest
%   (marked by a cross)
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

%   Parameters for conducting objects (uimenu)
objecttype{1}           = 'circle';     %   'none' or 'plate', 'sphere', 'brick', or 'cylinder' is allowed
objecttype{2}           = 'circle';     %   'none' or 'plate', 'sphere', 'brick', or 'cylinder' is allowed
objecttype{3}           = 'sphere';     %   'none' or 'plate', 'sphere', 'brick', or 'cylinder' is allowed
objectpotential{1}      = 0.0;          %    potential of the first conducting object, V
objectpotential{2}      = 1.0;          %    potential of the second conducting object, V
objectpotential{3}      = 0.5;          %    potential of the third conducting object, V
objectcharge           = 'yes';         %    enforces (yes/no) total charge zero for the third conductor

%   Data for the plate (if 'plate' is selected as a first object) (uitable)
strge.L(1)              = 0.2;         %   plate length, m
strge.W(1)              = 0.2;         %   plate width, m
strge.X(1)              = 0;           %   center position x, m
strge.Y(1)              = 0;           %   center position y, m
strge.Z(1)              = -0.03;       %   center position z, m
strge.ax(1)             = 0;           %   rotate about the x-axis, deg
strge.ay(1)             = 0;           %   rotate about the y-axis, deg
strge.az(1)             = 0;           %   rotate about the z-axis, deg
strge.par(1)            = 1;           %   uniform (0) or non-uniform (1) grid
strge.Tr(1)             = 400;         %   approximate number of triangular patches

%   Data for the plate (if 'plate' is selected as a second object) (uitable)
strge.L(2)              = 0.2;         %   plate length, m
strge.W(2)              = 0.2;         %   plate width, m
strge.X(2)              = 0;           %   center position x, m
strge.Y(2)              = 0;           %   center position y, m
strge.Z(2)              = +0.03;       %   center position z, m
strge.ax(2)             = 0;           %   rotate about the x-axis, deg
strge.ay(2)             = 0;           %   rotate about the y-axis, deg
strge.az(2)             = 0;           %   rotate about the z-axis, deg
strge.par(2)            = 1;           %   uniform (0) or non-uniform (1) grid
strge.Tr(2)             = 400;         %   approximate number of triangular patches

%   Data for the plate (if 'plate' is selected as a third object) (uitable)
strge.L(3)              = 0.05;         %   plate length, m
strge.W(3)              = 0.05;         %   plate width, m
strge.X(3)              = -0.06;        %   center position x, m
strge.Y(3)              = 0;            %   center position y, m
strge.Z(3)              = 0;            %   center position z, m
strge.ax(3)             = 45;           %   rotate about the x-axis, deg
strge.ay(3)             = 0;            %   rotate about the y-axis, deg
strge.az(3)             = 0;            %   rotate about the z-axis, deg
strge.par(3)            = 1;            %   uniform (0) or non-uniform (1) grid
strge.Tr(3)             = 100;          %   approximate number of triangular patches

%   Data for the sphere (if 'sphere' is selected as a first object) (uitable)
strsphere.a(1)          = 0.02;         %   sphere radius, m
strsphere.X(1)          = 0;            %   center position x, m
strsphere.Y(1)          = 0.0;          %   center position y, m
strsphere.Z(1)          = -0.03;        %   center position z, m
strsphere.Tr(1)         = 300;          %   approximate number of triangular patches

%   Data for the sphere (if 'sphere' is selected as a second object) (uitable)
strsphere.a(2)          = 0.02;        %   sphere radius, m
strsphere.X(2)          = 0.0;         %   center position x, m
strsphere.Y(2)          = 0;           %   center position y, m
strsphere.Z(2)          = +0.03;       %   center position z, m
strsphere.Tr(2)         = 300;         %   approximate number of triangular patches

%   Data for the sphere (if 'sphere' is selected as a third object) (uitable)
strsphere.a(3)          = 0.02;        %   sphere radius, m
strsphere.X(3)          = -0.06;       %   center position x, m
strsphere.Y(3)          = 0;           %   center position y, m
strsphere.Z(3)          = 0;           %   center position z, m
strsphere.Tr(3)  = 500;                %   approximate number of triangular patches

%   Data for the brick (if 'brick' is selected as a first object) (uitable)
strbrick.L(1)          = 0.1;         %    brick length, m
strbrick.W(1)          = 0.1;         %    brick width, m
strbrick.H(1)          = 0.005;       %    brick height, m
strbrick.X(1)          = 0;           %    center position x, m
strbrick.Y(1)          = 0;           %    center position y, m
strbrick.Z(1)          = -0.03;       %    center position z, m
strbrick.ax(1)         = 0;           %    rotate about the x-axis, deg
strbrick.ay(1)         = 0;           %    rotate about the y-axis, deg
strbrick.az(1)         = 0;           %    rotate about the z-axis, deg
strbrick.par(1)  = 1;                 %    uniform (0) or non-uniform (1) grid
strbrick.Tr(1)  = 200;                %   approximate number of triangular patches

%   Data for the brick (if 'brick' is selected as a second object) (uitable)
strbrick.L(2)          = 0.2;         %    brick length, m
strbrick.W(2)          = 0.2;         %    brick width, m
strbrick.H(2)          = 0.01;        %    brick height, m
strbrick.X(2)          = 0;           %    center position x, m
strbrick.Y(2)          = 0;           %    center position y, m
strbrick.Z(2)          = +0.03;       %    center position z, m
strbrick.ax(2)         = 0;           %    rotate about the x-axis, deg
strbrick.ay(2)         = 0;           %    rotate about the y-axis, deg
strbrick.az(2)         = 0;           %    rotate about the z-axis, deg
strbrick.par(2)  = 1;                 %    uniform (0) or non-uniform (1) grid
strbrick.Tr(2)  = 200;                %   approximate number of triangular patches

%   Data for the brick (if 'brick' is selected as a third object) (uitable)
strbrick.L(3)          = 0.05;        %    brick length, m
strbrick.W(3)          = 0.05;        %    brick width, m
strbrick.H(3)          = 0.01;        %    brick height, m
strbrick.X(3)          = -0.06;       %    center position x, m
strbrick.Y(3)          = 0;           %    center position y, m
strbrick.Z(3)          = 0;           %    center position z, m
strbrick.ax(3)         = 0;           %    rotate about the x-axis, deg
strbrick.ay(3)         = 0;           %    rotate about the y-axis, deg
strbrick.az(3)         = 0;           %    rotate about the z-axis, deg
strbrick.par(3)  = 1;                 %    uniform (0) or non-uniform (1) grid
strbrick.Tr(3)  = 100;

%   Data for the cylinder (if 'cylinder' is selected as a first object) (uitable)
strcylinder.R(1)          = 0.02;        %     cylinder radius, m
strcylinder.H(1)          = 0.2;         %     cylinder height, m
strcylinder.X(1)          = 0;           %     center position x, m
strcylinder.Y(1)          = 0;           %     center position y, m
strcylinder.Z(1)          = -0.06;       %     center position z, m
strcylinder.ax(1)         = 0;           %     rotate about the x-axis, deg
strcylinder.ay(1)         = 90;          %     rotate about the y-axis, deg
strcylinder.az(1)         = 0;           %     rotate about the z-axis, deg
strcylinder.Tr(1)         = 50;          %     approximate number of triangular patches
strcylinder.par(1)        = 1;

%   Data for the cylinder (if 'cylinder' is selected as a second object) (uitable)
strcylinder.R(2)          = 0.02;        %     cylinder radius, m
strcylinder.H(2)          = 0.2;         %     cylinder height, m
strcylinder.X(2)          = 0;           %     center position x, m
strcylinder.Y(2)          = 0;           %     center position y, m
strcylinder.Z(2)          = 0.06;        %     center position z, m
strcylinder.ax(2)         = 0;           %     rotate about the x-axis, deg
strcylinder.ay(2)         = 90;          %     rotate about the y-axis, deg
strcylinder.az(2)         = 0;           %     rotate about the z-axis, deg
strcylinder.Tr(2)         = 50;          %     approximate number of triangular patches
strcylinder.par(2)        = 1;

%   Data for the cylinder (if 'cylinder' is selected as a third object) (uitable)
strcylinder.R(3)          = 0.02;        %     cylinder radius, m
strcylinder.H(3)          = 0.01;        %     cylinder height, m
strcylinder.X(3)          = -0.06;       %     center position x, m
strcylinder.Y(3)          = 0;           %     center position y, m
strcylinder.Z(3)          = 0;           %     center position z, m
strcylinder.ax(3)         = 0;           %     rotate about the x-axis, deg
strcylinder.ay(3)         = 0;           %     rotate about the y-axis, deg
strcylinder.az(3)         = 0;           %     rotate about the z-axis, deg
strcylinder.Tr(3)         = 50;          %   approximate number of triangular patches
strcylinder.par(3)        = 1;

% This is the circle part

% Data for the circle (if 'circle' is selected as a first object) (uitable)
strcircle.R(1)              = 0.10;        %   circle radius, m
strcircle.X(1)              = 0;           %   center position x, m
strcircle.Y(1)              = 0;           %   center position y, m
strcircle.Z(1)              = -0.03;       %   center position z, m
strcircle.ax(1)             = 0;           %   rotate about the x-axis, deg
strcircle.ay(1)             = 0;           %   rotate about the y-axis, deg
strcircle.az(1)             = 0;           %   rotate about the z-axis, deg
strcircle.par(1)            = 1;           %   uniform (0) or non-uniform (1) grid
strcircle.Tr(1)             = 200;         %   approximate number of triangular patches

% Data for the circle (if 'circle' is selected as a second object) (uitable)
strcircle.R(2)              = 0.10;        %   circle radius, m
strcircle.X(2)              = 0;           %   center position x, m
strcircle.Y(2)              = 0;           %   center position y, m
strcircle.Z(2)              = 0.03;        %   center position z, m
strcircle.ax(2)             = 0;           %   rotate about the x-axis, deg
strcircle.ay(2)             = 0;           %   rotate about the y-axis, deg
strcircle.az(2)             = 0;           %   rotate about the z-axis, deg
strcircle.par(2)            = 1;           %   uniform (0) or non-uniform (2) grid
strcircle.Tr(2)             = 200;         %   approximate number of triangular patches

% Data for the circle (if 'circle' is selected as a third object) (uitable)
strcircle.R(3)              = 0.025;       %   circle radius, m
strcircle.X(3)              = -0.06;       %   center position x, m
strcircle.Y(3)              = 0;           %   center position y, m
strcircle.Z(3)              = 0;           %   center position z, m
strcircle.ax(3)             = 45;          %   rotate about the x-axis, deg
strcircle.ay(3)             = 0;           %   rotate about the y-axis, deg
strcircle.az(3)             = 0;           %   rotate about the z-axis, deg
strcircle.par(3)            = 1;           %   uniform (0) or non-uniform (1) grid
strcircle.Tr(3)             = 100;         %   approximate number of triangular patches

%%  Output graphics parameters
%   Surface charge and general visualization (uitable)
%   Parameters to scale charge distribution for better visual inspection
strsc.positive = 1;         %   positive charge densities higher than this number times the average
%   positive charge density are assigned the same value
strsc.negative = 1;         %   negative charge densities smaller than this number times the average
%   negative charge density are assigned the same value

%   Visualization of the E-field in the observation plane (uitable)
strop.yes = 'no';           %   include (yes) or not (no) plot of the E-field
strop.potential = 'yes';    %   include (yes) or not (no) plot of the electric potential
strop.planetype   = 'yz';   %   xy, xz, or yz plane
strop.planex = -0.06;       %   plane center: x in m
strop.planey = 0.0;         %   plane center: y in m
strop.planez = 0.0;         %   plane center: z in m
strop.planesizex   = 0.25;  %   plane length in m
strop.planesizey   = 0.25;  %   plane width in m
strop.divisionsx   = 16;    %   divisions with respect to length
strop.divisionsy   = 16;    %   divisions with respect to width
strop.arrow = 1.5;          %   relative arrow size versus default size

%   Parameters for the observation point(s) (used to obtain exact values of the field) (uitable)
stroc.yes         = 'no';     %   'yes' - present; 'no' - absent
stroc.x           = 0.0;       %   x position in m
stroc.y           = 0.0;       %   y position in m
stroc.z           = 0.0;       %   z position in m
stroc.size        = 1;         %   relative marker size versus default size

%% global parameters of output results

strout.Cnum             = 0;                     %    Numerical capacitance in pF
strout.ChargeLower      = 0;                     %    Total charge of the lower plate, C
strout.ChargeUpper      = 0;                     %    Total charge of the upper plate, C
strout.ChargeTotal      = 0;                     %    Total charge of the entire structure, C
strout.ObjectPotential  = 0;                     %    Average object potential, V  
strout.E(1,1)           = 0;
strout.E(1,2)           = 0;
strout.E(1,3)           = 0;                     %    E-field at the observation point (if any)
strout.phi              = 0;                     %    Electric potential at the observation point, V (if any)
strout.AvgPotential1    = 0;                     %    Average potential of object #1
strout.PatchesTotal     = 0;                     %    Total number of triangular patches in the mesh
strout.quality          = 0;                     %    Minimum triangle quality
strout.time1            = 0;                     %    CPU time in sec for filling the MoM matrix
strout.time2            = 0;                     %    CPU time in sec for solving the system of MoM eqs.
strout.ChargeTotal      = 0;
strout.PatchesCap       = 0;
strout.ChargeObject = 0;

%%   Numerical parameters (uitable)
R   = 5;        %    dimensionless radius of an enclosing sphere for precise integration
%    R=0 - only self integrals are calculated precisely
%    R=10- integrals for all neighbor triangles whose
%    center-to-center distances from the observation triangle
%    are less than ten times the effective triangle size are calculated
%    precisely
gauss = 7;      %    Number of integration points in the Gaussian quadrature
%    Numbers 1, 4, 7, 13, 25 are permitted

%%  End of GUI window - input parameters
% E24_BUILD
%-------------------------------------------------------------------------------
% File name   : E24_build.m
% Generated on:
% Description :
%-------------------------------------------------------------------------------

%   Input Parameters--------
%%  End of GUI window - output parameters
% Initialize handles structure
handles = struct();
handles.geometry = 0;   %   The real handles
handles.simulate = 0;   %   Only an indicator that the simulations are complete
handles.object = 2;
handles.object2 = 2;
handles.object3 = 3;
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
            E = efield(strop.Points, c, P, t, Center, Area, normals, Size, R, eps0, msg);
            strop.E = E;
        end
        %%  Find the potential in a plane
        if strcmp(strop.potential, 'yes')
            msg = 1;            
            strop.Potential = potential(Points, c, P, t, Center, Area, normals, Size, R, eps0, msg)...
                -strout.AvgPotential1+objectpotential{1};
        end
        %%  Find the E-field/potential at an observation point(s)
        strout.E          = 'none selected';
        strout.phi        = 'none selected';
        if strcmp(stroc.yes, 'yes')
            msg = 0;
            points = [stroc.x stroc.y stroc.z];
            strout.E = efield(points, c, P, t, Center, Area, normals, Size, R, eps0, msg);            
            strout.phi = potential(points, c, P, t, Center, Area, normals, Size, R, eps0, msg)...
                -strout.AvgPotential1+objectpotential{1};
        end      
        %%   Output graphics
        set(0, 'CurrentFigure', handles.figure1);
        if handles.geometry ~= 0; delete(handles.geometry); end
        io = 1;
        handles.geometry = subplot(1,1,1,'Parent',handles.uipanel2);
        graphics_S24(io, strop, stroc, P, t, c, Area, strsc);
    end

    function cleaning
        %   cleaning figure without re-plotting
        strout = clearoutput(strout); io = 0;
        if handles.geometry ~= 0; delete(handles.geometry); handles.geometry = 0; end
        handles.geometry = subplot(1,1,1,'Parent',handles.uipanel2);
        graphics_S24(io, strop, stroc, P, t, c, Area, strsc);
    end

    function strout = clearoutput(strout)
        %   Clear results anyways
        strout.Cnum             = 0;                     %    Numerical capacitance in pF
        strout.ChargeLower      = 0;                     %    Total charge of the lower plate, C
        strout.ChargeUpper      = 0;                     %    Total charge of the upper plate, C
        strout.ChargeTotal      = 0;                     %    Total charge of the entire structure, C
        strout.ObjectPotential  = 0;                     %    Average object potential, V  
        strout.E(1,1)           = 0;
        strout.E(1,2)           = 0;
        strout.E(1,3)           = 0;                     %    E-field at the observation point (if any)
        strout.phi              = 0;                     %    Electric potential at the observation point, V (if any)
        strout.AvgPotential1    = 0;                     %    Average potential of object #1          
        strout.time1            = 0;                     %    CPU time in sec for filling the MoM matrix
        strout.time2            = 0;                     %    CPU time in sec for solving the system of MoM eqs.
        strout.ChargeTotal      = 0;    
        strout.ChargeObject = 0;
    end

    function geometry()        
        for m = 1:size(objecttype, 2)
            if strcmp(objecttype{m}, 'plate')
                [PD{m}, tD{m}] = plate(strge.L(m), strge.W(m), strge.Tr(m), strge.par(m), 0, 0, 1);
                [PD{m}] = rotatex(PD{m}, strge.ax(m));
                [PD{m}] = rotatey(PD{m}, strge.ay(m));
                [PD{m}] = rotatez(PD{m}, strge.az(m));
                [centersD{m}, normalsD{m}] = normcenters(PD{m}, tD{m});
                PD{m}(:, 1) = PD{m}(:, 1) + strge.X(m);
                PD{m}(:, 2) = PD{m}(:, 2) + strge.Y(m);
                PD{m}(:, 3) = PD{m}(:, 3) + strge.Z(m);
                centersD{m}(:, 1) = centersD{m}(:, 1) + strge.X(m);
                centersD{m}(:, 2) = centersD{m}(:, 2) + strge.Y(m);
                centersD{m}(:, 3) = centersD{m}(:, 3) + strge.Z(m);
                tD{m}(:, 4) = m;
            end
            if strcmp(objecttype{m}, 'sphere')
                [PD{m}, tD{m}] = sphere(1, strsphere.Tr(m));
                PD{m} = PD{m}*strsphere.a(m);
                [centersD{m}, normalsD{m}] = normcenters(PD{m}, tD{m});
                PD{m}(:, 1) = PD{m}(:, 1) + strsphere.X(m);
                PD{m}(:, 2) = PD{m}(:, 2) + strsphere.Y(m);
                PD{m}(:, 3) = PD{m}(:, 3) + strsphere.Z(m);
                tD{m}(:, 4) = m;
            end
            if strcmp(objecttype{m}, 'brick')
                [PD{m}, tD{m}] = brick(strbrick.L(m), strbrick.W(m), strbrick.H(m), strbrick.Tr(m), strbrick.par(m));
                [PD{m}] = rotatex(PD{m}, strbrick.ax(m));
                [PD{m}] = rotatey(PD{m}, strbrick.ay(m));
                [PD{m}] = rotatez(PD{m}, strbrick.az(m));
                [centersD{m}, normalsD{m}] = normcenters(PD{m}, tD{m});
                PD{m}(:, 1) = PD{m}(:, 1) + strbrick.X(m);
                PD{m}(:, 2) = PD{m}(:, 2) + strbrick.Y(m);
                PD{m}(:, 3) = PD{m}(:, 3) + strbrick.Z(m);
                centersD{m}(:, 1) = centersD{m}(:, 1) + strbrick.X(m);
                centersD{m}(:, 2) = centersD{m}(:, 2) + strbrick.Y(m);
                centersD{m}(:, 3) = centersD{m}(:, 3) + strbrick.Z(m);
                tD{m}(:, 4) = m;
            end
            if strcmp(objecttype{m}, 'cylinder')
                [PD{m}, tD{m}]  = cylinder(strcylinder.R(m), strcylinder.H(m), strcylinder.Tr(m), strcylinder.par(m));
                [PD{m}] = rotatex(PD{m}, strcylinder.ax(m));
                [PD{m}] = rotatey(PD{m}, strcylinder.ay(m));
                [PD{m}] = rotatez(PD{m}, strcylinder.az(m));
                [centersD{m}, normalsD{m}] = normcenters(PD{m}, tD{m});
                PD{m}(:, 1) = PD{m}(:, 1) + strcylinder.X(m);
                PD{m}(:, 2) = PD{m}(:, 2) + strcylinder.Y(m);
                PD{m}(:, 3) = PD{m}(:, 3) + strcylinder.Z(m);
                centersD{m}(:, 1) = centersD{m}(:, 1) + strcylinder.X(m);
                centersD{m}(:, 2) = centersD{m}(:, 2) + strcylinder.Y(m);
                centersD{m}(:, 3) = centersD{m}(:, 3) + strcylinder.Z(m);
                tD{m}(:, 4) = m;
            end
            if strcmp(objecttype{m}, 'circle')
                [PD{m}, tD{m}] = circle(strcircle.R(m), strcircle.Tr(m), strcircle.par(m), 0, 0, 1);
                [PD{m}] = rotatex(PD{m}, strcircle.ax(m));
                [PD{m}] = rotatey(PD{m}, strcircle.ay(m));
                [PD{m}] = rotatez(PD{m}, strcircle.az(m));
                [centersD{m}, normalsD{m}] = normcenters(PD{m}, tD{m});
                PD{m}(:, 1) = PD{m}(:, 1) + strcircle.X(m);
                PD{m}(:, 2) = PD{m}(:, 2) + strcircle.Y(m);
                PD{m}(:, 3) = PD{m}(:, 3) + strcircle.Z(m);
                centersD{m}(:, 1) = centersD{m}(:, 1) + strcircle.X(m);
                centersD{m}(:, 2) = centersD{m}(:, 2) + strcircle.Y(m);
                centersD{m}(:, 3) = centersD{m}(:, 3) + strcircle.Z(m);
                tD{m}(:, 4) = m;
            end
            
        end
        
        %%  Combining all meshes together (while keeping the fourth index)
        P = [];
        t = [];
        normals = [];
        for m = 1:size(objecttype, 2)
            tD{m}(:, 1:3)   = tD{m}(:, 1:3) + size(P, 1);
            P               = [P' PD{m}']';
            t               = [t' tD{m}']';
            normals         = [normals' normalsD{m}']';     %   normals
        end
        NM = size(t, 1);
        strout.PatchesTotal = NM;
        strout.PatchesCap   = sum(size(tD{1}, 1) + size(tD{2}, 1));
        c = zeros(size(t, 1));
        
        %%  Input graphics
        strout.quality  = min(simpqual(P, t));
        set(0, 'CurrentFigure', handles.figure1);
        if handles.geometry ~= 0; delete(handles.geometry); handles.geometry = 0 ;end
        handles.geometry = subplot(1,1,1,'Parent',handles.uipanel2);
        io = 0;
        graphics_S24(io, strop, stroc, P, t, c, Area, strsc);
    end

    function simulate()
        h    = waitbar(0, 'Please wait - filling the MoM matrix');
        time1 = cputime;
        %%  Parameter initialization for the combined mesh
        NM     = size(t, 1);
        Center = zeros(NM, 3);   %   face center
        Area   = zeros(NM, 1);   %   face area
        Size   = zeros(NM, 1);   %   face size defined as distance from center to furthest vertex
        
        %%   Find base parameters for all faces (metal)
        for m = 1:NM
            Vertexes        = P(t(m, 1:3)', :)';
            r1              = Vertexes(:, 1);
            r2              = Vertexes(:, 2);
            r3              = Vertexes(:, 3);
            tempv           = cross(r2-r1, r3-r1);
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
        ObsPoint                    = zeros(IndexS*NM, 3);
        for p =1:IndexS
            ObsPoint(p:IndexS:end, :) = coeffS(1, p)*P(t(:, 1), :) +  coeffS(2, p)*P(t(:, 2), :) +  coeffS(3, p)*P(t(:, 3), :);
            ObsIndex = repmat([1:IndexS], 1, size(t, 1))';
            WeightsS = repmat(weightsS,  size(t, 1), 1);
        end
        
        %%  Filling Z
        %   Prepare distance matrix DIST for Z
        %   First row - distances from the center of face #1 to all other face centers
        %   Second row - distances from the center of face #2 to all other face centers, etc.
        DIST = zeros(NM, NM);
        for m = 1:NM
            temp = Center'- repmat(Center(m, :)', 1, NM);
            DIST(m, :) = sqrt(dot(temp, temp));
        end
        %   MoM matrix Z initialization (center point integration)
        Z = 1./DIST;
        for n = 1:NM
            Z(:, n) = Z(:, n)*Area(n);
        end
        %   Loop over columns of impedance matrix Z!
        %   n is the number of the inner triangle (every column has the only the n-th inner triangle)
        for n =1:NM
            normdummy    = normals(n, :);
            index   = find(DIST(1:NM, n)'./(Size(1:NM)'*Size(n))<=R+1e-16);   % index is local
            r1      = P(t(n, 1), :);    %   row
            r2      = P(t(n, 2), :);    %   row
            r3      = P(t(n, 3), :);    %   row
            
            m1 = repmat(index'*IndexS-IndexS, [1 IndexS])';
            m3 = m1(:) + ObsIndex(1:IndexS*length(index));
            [I, IRho] = potint(r1, r2, r3, normdummy, ObsPoint(m3, :)); %   I was calculated with the area Area(n)
            Z(index, n) = sum(WeightsS(1:length(index), :).*reshape(I, IndexS, length(index))', 2);
            
            waitbar(n/NM);
        end
        %   Full matrix Z
        Z = +Z/(4*pi*eps0);
        close(h);
        
        %%   The full voltage vector (+/- 0.5V charge)
        V               =  zeros(NM, 1);
        V(t(:, 4)==1)   = objectpotential{1};
        V(t(:, 4)==2)   = objectpotential{2};
        V(t(:, 4)==3)   = objectpotential{3};
        
        %% Charge conservation law (critical here)
        %   Subtract the last row from the others and replace the
        %   last equation by the conservation law
        if strcmp(objectcharge, 'yes')
            %   enforce the charge conservation law for the capacitor
            for m = 1:strout.PatchesCap-1
                V(m)    = V(m)    - V(strout.PatchesCap);
                Z(m, :) = Z(m,: ) - Z(strout.PatchesCap, :);
            end
            % Last row
            Z(strout.PatchesCap, 1:strout.PatchesCap) = ones(1, strout.PatchesCap).*...
                Area(1:strout.PatchesCap)';
            Z(strout.PatchesCap, strout.PatchesCap+1:end) = 0;
            V(strout.PatchesCap) = 0;
            %   enforce the charge conservation law (resulting in zero charge) for the conducting object
            for m = strout.PatchesCap+1:strout.PatchesTotal-1
                V(m)    = V(m)    - V(strout.PatchesTotal);
                Z(m, :) = Z(m,: ) - Z(strout.PatchesTotal, :);
            end
            % Last row
            Z(strout.PatchesTotal, strout.PatchesCap+1:strout.PatchesTotal) = ones(1, strout.PatchesTotal-strout.PatchesCap).*...
                Area(strout.PatchesCap+1:strout.PatchesTotal)';
            Z(strout.PatchesTotal, 1:strout.PatchesCap) = 0;
            V(strout.PatchesTotal) = 0;
        else
            %   enforce charge consevation law for the entire structure
            for m = 1:NM-1
                V(m)    = V(m) - V(end);
                Z(m, :) = Z(m, :) - Z(NM, :);
            end
            % Last metal row
            Z(NM, :) = ones(1, NM).*Area';
            V(end) = 0;        
        end
        
        time1 = cputime - time1;
        %%   Solve MoM equations
        if (NM>1000)
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
        if (NM>1000)
            close(h);
        end
        %%  Charge calculation              
        Charge          = c'.*Area';                  %   charges for every face
        ChargeLower     = sum(Charge(t(:, 4)==1));    %   total charge of the lower plate, C
        ChargeUpper     = sum(Charge(t(:, 4)==2));    %   total charge of the upper plate, C
        ChargeObject    = sum(Charge(t(:, 4)==3));    %   total charge of the body, C
        ChargeTotal     = ChargeLower + ChargeUpper + ChargeObject;
        Cnum            = sum(Charge(t(:, 4)==2))*1e12/(objectpotential{2}-objectpotential{1});
        Ctheory         = 'Undefined';
        if strcmp(objecttype{1},'plate') && strcmp(objecttype{2},'plate')
            temp        = max([abs(strge.X(2)-strge.X(1)) abs(strge.Y(2)-strge.Y(1)) abs(strge.Z(2)-strge.Z(1))]);
            Ctheory     = min(strge.L(1:2))*min(strge.W(1:2))*eps0/temp*1e12;
        else
            if strcmp(objecttype{1},'circle') && strcmp(objecttype{2},'circle')
                temp        = max([abs(strcircle.X(2)-strcircle.X(1)) abs(strcircle.Y(2)-strcircle.Y(1)) abs(strcircle.Z(2)-strcircle.Z(1))]);
                Ctheory     = pi*(min(strcircle.R(1),strcircle.R(2))^2)*eps0/temp*1e12;
            end
        end
        
        %%  Fields and output graphics
        %---------------------------------------------------------------------
        msg = 2;
        Points = Center(t(:, 4)==1, :);
        temp = potential(Points, c, P, t, Center, Area, normals, Size, R, eps0, msg);
        strout.AvgPotential1 = sum(temp.*Area(t(:, 4)==1))/sum(Area(t(:, 4)==1));
        %---------------------------------------------------------------------
        %   Finding average potential of the object
        msg = 3;
        Points = Center(t(:, 4)==3, :);
        temp = potential(Points, c, P, t, Center, Area, normals, Size, R, eps0, msg);
        ObjectPotential = sum(temp.*Area(t(:, 4)==3))/sum(Area(t(:, 4)==3));
        ObjectPotential =  ObjectPotential  - strout.AvgPotential1 + objectpotential{1};        
        %---------------------------------------------------------------------
        handles.simulate = 1;
        outputgraphics;
        
        %%  GUI window - output parameters
        strout.Cnum     = Cnum;                     %    Numerical capacitance in pF
        strout.ChargeLower     = ChargeLower;       %    Total charge of the lower plate, C
        strout.ChargeUpper     = ChargeUpper;       %    Total charge of the upper plate, C
        strout.ChargeObject    = ChargeObject;      %    Total charge of the body, C
        strout.ChargeTotal     = ChargeTotal;       %    Total charge of the entire structure, C
        strout.ObjectPotential = ObjectPotential;   %    Average object potential, V
        strout.E;                                   %    E-field at the observation point (if any)
        strout.PatchesCap;                          %    Total number of triangular patches in the capacitor mesh
        strout.PatchesTotal;                        %    Total number of triangular patches in the body mesh
        strout.quality  = min(simpqual(P, t));
        strout.time1    = time1;                    %    CPU time in sec for filling the MoM matrix
        strout.time2    = time2;                    %    CPU time in sec for solving the system of MoM eqs.
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
            'Position', [103.8 23.0769230769231 145.6 38.4615384615385], ...
            'Name', 'E24-A capacitor with a metal or other conducting object inside: the three-conductor problem', ...
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
        
        %table_position
        Position_local = get(handles.figure1,'position');
        Position_local(1) = 0.5*Position_local(1) ;
        Position_local(2) = 0.87*Position_local(2);
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
            E24;
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
                strsc               = M.strsc;
                strge               = M.strge;
                strop               = M.strop;
                strout              = M.strout;
                R                   = M.R;
                gauss               = M.gauss;
                objectpotential     = M.objectpotential;
                objecttype          = M.objecttype;
                strsphere           = M.strsphere;
                strbrick            = M.strbrick;
                strcylinder         = M.strcylinder;
                strcircle           = M.strcircle;
                stroc               = M.stroc;
                switch objecttype{1}
                    case 'plate'
                        handles.object = 1;
                    case 'circle'
                        handles.object = 2;
                    case 'sphere'
                        handles.object = 3;
                    case 'brick'
                        handles.object = 4;
                    case 'cylinder'
                        handles.object = 5;
                end               
                switch objecttype{2}
                    case 'plate'
                        handles.object2 = 1;
                    case 'circle'
                        handles.object2 = 2;
                    case 'sphere'
                        handles.object2 = 3;
                    case 'brick'
                        handles.object2 = 4;
                    case 'cylinder'
                        handles.object2 = 5;
                end
                switch objecttype{3}
                    case 'plate'
                        handles.object3 = 1;
                    case 'circle'
                        handles.object3 = 2;
                    case 'sphere'
                        handles.object3 = 3;
                    case 'brick'
                        handles.object3 = 4;
                    case 'cylinder'
                        handles.object3 = 5;
                end
                
                geometry();
            end
        end
        
        uimenu(...
            'Parent', handles.fFileMenu,...
            'HandleVisibility','callback', ...
            'Label','Save',...
            'Callback', {@Save_Callback});
        
        function Save_Callback(~, ~)
            save('E24project');
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
        
        
        
        % --- First Plate Menu ---------------------------
        handles.fFirstMenu   =   uimenu(...       % File menu
            'Parent', handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','1st Conductor');
        
        handles.fObjectTypeMenu      =   uimenu(...       % File menu
            'Parent', handles.fFirstMenu,...
            'HandleVisibility','callback', ...
            'Label','Object Type',...
            'Callback', @fObjectTypeCallback);
        
        
        function fObjectTypeCallback(~,~)
            % --- RADIO BUTTONS -------------------------------------
            handles.figure2 = figure( ...
                'Tag', 'figure1', ...
                'Units', 'characters', ...
                'Position', [103.8 38.6153846153846 28.4 15.7692307692308], ...
                'Name', 'Object Type', ...
                'MenuBar', 'none', ...
                'NumberTitle', 'off', ...
                'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
                'Resize', 'on');
            
            handles.uipanel3 = uibuttongroup( ...
                'Parent', handles.figure2, ...
                'Tag', 'uipanel1', ...
                'Units', 'normalized', ...
                'Position', [0.0141823161189358 0.0334928229665072 0.965277777777778 0.961722488038277], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'Title', 'Object Type');
            
            set(handles.uipanel3,'SelectionChangeFcn',@uibuttongroup1_SelectionChangeFcn);
            
            
            function uibuttongroup1_SelectionChangeFcn(~,eventdata)
                switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
                    case 'plate'
                        handles.object = 1;
                        objecttype{1} = 'plate';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                        strout = clearoutput(strout);
                    case 'circle'
                        handles.object = 2;
                        objecttype{1} = 'circle';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                        strout = clearoutput(strout);
                    case 'sphere'
                        handles.object = 3;
                        objecttype{1} = 'sphere';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                        strout = clearoutput(strout);
                    case 'brick'
                        handles.object = 4;
                        objecttype{1} = 'brick';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                        strout = clearoutput(strout);
                    case 'cylinder'
                        handles.object = 5;
                        objecttype{1} = 'cylinder';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                        strout = clearoutput(strout);
                    otherwise
                        errordlg('Fatal error;Contact your Administrator');
                end
            end
            
            
            handles.radiobutton7 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'plate', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.877777777777778 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Plate');
            
            
            handles.radiobutton8 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'circle', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.666666666666667 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Circle');
            
          
            handles.radiobutton9 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'sphere', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.455555555555556 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Sphere');
            
            
            handles.radiobutton10 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'brick', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.244444444444444 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Brick');
            
            
            
            handles.radiobutton11 = uicontrol( ...
                'Parent', handles.uipanel3, ...
                'Tag', 'cylinder', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.0333333333333333 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Cylinder');
            
            switch handles.object
                case 1
                    set(handles.uipanel3,'Selectedobject',handles.radiobutton7);
                case 2
                    set(handles.uipanel3,'Selectedobject',handles.radiobutton8);
                case 3
                    set(handles.uipanel3,'Selectedobject',handles.radiobutton9);
                case 4
                    set(handles.uipanel3,'Selectedobject',handles.radiobutton10);
                case 5
                    set(handles.uipanel3,'Selectedobject',handles.radiobutton11);
                otherwise
                    errordlg('Fatal error;Contact your Administrator');
            end
            
        end
        
        
        handles.fPotentialMenu      =   uimenu(...       % File menu
            'Parent', handles.fFirstMenu,...
            'HandleVisibility','callback', ...
            'Label','Potential',...
            'Callback', @fFirstPotentialCallback);
        
        function fFirstPotentialCallback(~, ~)
            % Callback function run when the Potential menu item is selected
            prompt = {'Value in V(volts)'};
            dlg_title = '1st Conductor V';
            num_lines = 1;
            def = {num2str( objectpotential{1})};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            user_entry = str2double(answer);
            if isnan(user_entry)
                errordlg('You must enter a numeric value','Bad Input','modal');
                return;
            end
            if ~isempty(user_entry)
                objectpotential{1} = user_entry;
                strge.phi  = objectpotential{1};
                cleaning;
            end
        end
        
        
        
        handles.fParametersMenu      =   uimenu(...       % File menu
            'Parent', handles.fFirstMenu,...
            'HandleVisibility','callback', ...
            'Label','Parameters',...
            'Callback', @fPreParametersCallback);
        
        function fPreParametersCallback(~,~)
            
            switch handles.object
                case 1
                    fConductingPlateCallback;
                case 2
                    fConductingCircleCallback;
                case 3
                    fConductingSphereCallback;
                case 4
                    fConductingBrickCallback;
                case 5
                    fConductingCylinderCallback;
                otherwise
                    errordlg('Fatal error;Contact your Administrator');
            end
            
            
        end
        
        
        function fConductingPlateCallback(~, ~)
            strge0 = strge;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'parameters of the plate setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Plate length,m',strge.L(1);...
                'Plate Width,m',strge.W(1);...
                'Center position x,m',strge.X(1);...
                'Center position y,m',strge.Y(1);...
                'Center position z,m',strge.Z(1);...
                'Rotate about x-axis, deg',strge.ax(1);...
                'Rotate about y-axis, deg',strge.ay(1);...
                'Rotate about z-axis, deg',strge.az(1);...
                'Approximate number of triangles',strge.Tr(1);...
                'Triangle size at boundaries vs. center size', strge.par(1)};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'You have chosen the plate case',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback});
            function edit_callback(hObject, ~)
                data = get(hObject,'Data');
                strge.L(1)         = data{11};
                strge.W(1)         = data{12};
                strge.X(1)         = data{13};
                strge.Y(1)         = data{14};
                strge.Z(1)         = data{15};
                strge.ax(1)        = data{16};
                strge.ay(1)        = data{17};
                strge.az(1)        = data{18};
                strge.Tr(1)        = data{19};
                strge.par(1)       = data{20};                
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strge = strge0;
                delete(local);
            end
            
        end
        
        %% Circle case skipped
        function fConductingCircleCallback(~, ~)
            strcircle0 = strcircle;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String',{'parameters of circle setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Circle radius,m',strcircle.R(1);...
                'Center position x,m',strcircle.X(1);...
                'Center position y,m',strcircle.Y(1);...
                'Center position z,m',strcircle.Z(1);...
                'Rotate about x-axis, deg',strcircle.ax(1);...
                'Rotate about y-axis, deg',strcircle.ay(1);...
                'Rotate about z-axis, deg',strcircle.az(1);...
                'Approximate numbers of triangles',strcircle.Tr(1);...
                'Triangle size at boundaries vs. center size', strcircle.par(1)};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'You have chosen the circle',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback});
            function edit_callback(hObject, ~)
                data = get(hObject,'Data');
                strcircle.R(1)         = data{10};
                strcircle.X(1)         = data{11};
                strcircle.Y(1)         = data{12};
                strcircle.Z(1)         = data{13};
                strcircle.ax(1)        = data{14};
                strcircle.ay(1)        = data{15};
                strcircle.az(1)        = data{16};
                strcircle.Tr(1)         = data{17};
                strcircle.par(1) = data{18};
                                
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strcircle = strcircle0;
                delete(local);
            end
            
        end
        
        %% Sphere case
        
        function fConductingSphereCallback(~, ~)
            strsphere0 = strsphere;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units',Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'parameters of sphere setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Sphere radius, m',strsphere.a(1);...
                'Center position x,m',strsphere.X(1);...
                'Center position y,m',strsphere.Y(1);...
                'Center position z,m',strsphere.Z(1);...
                'Approx.# of triangular patches',strsphere.Tr(1)};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'You have chosen the sphere case ',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback});
            function edit_callback(hObject, ~)
                data = get(hObject,'Data');
                strsphere.a(1)                 = data{6};
                strsphere.X(1)                 = data{7};
                strsphere.Y(1)                 = data{8};
                strsphere.Z(1)                 = data{9};
                strsphere.Tr(1)         = data{10};
                
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strsphere = strsphere0;
                delete(local);
            end
            
        end
        %% Brick Case
        
        function fConductingBrickCallback(~, ~)
            strbrick0 = strbrick;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'Parameters of brick setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Brick length,m',strbrick.L(1);...
                'Brick width,m',strbrick.W(1);...
                'Brick height,m',strbrick.H(1);...
                'Center position x,m',strbrick.X(1);...
                'Center position y,m',strbrick.Y(1);...
                'Center position z,m',strbrick.Z(1);...
                'Rotate about x-axis, deg',strbrick.ax(1);...
                'Rotate about y-axis, deg',strbrick.ay(1);...
                'Rotate about z-axis, deg',strbrick.az(1);...
                'Approximate # of triangles - top face',strbrick.Tr(1);...
                'Triangle size at boundaries vs. center size', strbrick.par(1)};
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
                strbrick.L(1)         = data{12};
                strbrick.W(1)         = data{13};
                strbrick.H(1)        = data {14};
                strbrick.X(1)         = data{15};
                strbrick.Y(1)         = data{16};
                strbrick.Z(1)         = data{17};
                strbrick.ax(1)        = data{18};
                strbrick.ay(1)        = data{19};
                strbrick.az(1)        = data{20};
                strbrick.Tr(1)        = data{21};
                strbrick.par(1)        = data{22};
                
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
        
        %% Cylinder case
        
        function fConductingCylinderCallback(~, ~)
            strcylinder0 = strcylinder;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'parameters of cylinder setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Cylinder radius,m',strcylinder.R(1);...
                'Cylinder height,m',strcylinder.H(1);...
                'Center position x,m',strcylinder.X(1);...
                'Center position y,m',strcylinder.Y(1);...
                'Center position z,m',strcylinder.Z(1);...
                'Rotate about x-axis, deg',strcylinder.ax(1);...
                'Rotate about y-axis, deg',strcylinder.ay(1);...
                'Rotate about z-axis, deg',strcylinder.az(1);...
                'Approximate # of triangles - top face',strcylinder.Tr(1);...
                'Triangle size at boundaries vs. center size',strcylinder.par(1)};
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
                strcylinder.R(1)         = data{11};
                strcylinder.H(1)         = data{12};
                strcylinder.X(1)         = data{13};
                strcylinder.Y(1)         = data{14};
                strcylinder.Z(1)         = data{15};
                strcylinder.ax(1)        = data{16};
                strcylinder.ay(1)        = data{17};
                strcylinder.az(1)        = data{18};
                strcylinder.Tr(1)        = data{19};
                strcylinder.par(1)       = data{20};             
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
        
        % --- Second Conductor Menu ----------------------------
        
        handles.fSecondMenu    =   uimenu(...       % File menu
            'Parent',handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','2nd Conductor');
        
        handles.fSecondObjectTypeMenu      =   uimenu(...       % File menu
            'Parent', handles.fSecondMenu,...
            'HandleVisibility','callback', ...
            'Label','Object Type',...
            'Callback', @fSecondObjectTypeCallback);
        
        
        function fSecondObjectTypeCallback(~,~)
            
            % --- RADIO BUTTONS -------------------------------------
            handles.figure3 = figure( ...
                'Tag', 'figure2', ...
                'Units', 'characters', ...
                'Position', [103.8 38.6153846153846 28.4 15.7692307692308], ...
                'Name', 'Object Type', ...
                'MenuBar', 'none', ...
                'NumberTitle', 'off', ...
                'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
                'Resize', 'on');
            
            handles.uipanel4 = uibuttongroup( ...
                'Parent', handles.figure3, ...
                'Tag', 'uipanel1', ...
                'Units', 'normalized', ...
                'Position', [0.0141823161189358 0.0334928229665072 0.965277777777778 0.961722488038277], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'Title', 'Object Type');
            
            set(handles.uipanel4,'SelectionChangeFcn',@uibuttongroup2_SelectionChangeFcn);
            
            
            function uibuttongroup2_SelectionChangeFcn(~,eventdata)
                switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
                    case 'plate'
                        handles.object2 = 1;
                        objecttype{2} = 'plate';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                        strout = clearoutput(strout);
                    case 'circle'
                        handles.object2 = 2;
                        objecttype{2} = 'circle';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                        strout = clearoutput(strout);
                    case 'sphere'
                        handles.object2 = 3;
                        objecttype{2} = 'sphere';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                        strout = clearoutput(strout);
                    case 'brick'
                        handles.object2 = 4;
                        objecttype{2} = 'brick';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                        strout = clearoutput(strout);
                    case 'cylinder'
                        handles.object2 = 5;
                        objecttype{2} = 'cylinder';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                        strout = clearoutput(strout);
                    otherwise
                        errordlg('Fatal error;Contact your Administrator');
                end
            end
            
            
            handles.radiobutton7 = uicontrol( ...
                'Parent', handles.uipanel4, ...
                'Tag', 'plate', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.877777777777778 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Plate');
            
            
            handles.radiobutton8 = uicontrol( ...
                'Parent', handles.uipanel4, ...
                'Tag', 'circle', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.666666666666667 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Circle');
            
            
            
            handles.radiobutton9 = uicontrol( ...
                'Parent', handles.uipanel4, ...
                'Tag', 'sphere', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.455555555555556 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Sphere');
            
            
            
            handles.radiobutton10 = uicontrol( ...
                'Parent', handles.uipanel4, ...
                'Tag', 'brick', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.244444444444444 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Brick');
            
            
            
            handles.radiobutton11 = uicontrol( ...
                'Parent', handles.uipanel4, ...
                'Tag', 'cylinder', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.0333333333333333 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Cylinder');
            
            switch handles.object2
                case 1
                    set(handles.uipanel4,'Selectedobject',handles.radiobutton7);
                case 2
                    set(handles.uipanel4,'Selectedobject',handles.radiobutton8);
                case 3
                    set(handles.uipanel4,'Selectedobject',handles.radiobutton9);
                case 4
                    set(handles.uipanel4,'Selectedobject',handles.radiobutton10);
                case 5
                    set(handles.uipanel4,'Selectedobject',handles.radiobutton11);
                otherwise
                    errordlg('Fatal error;Contact your Administrator');
            end
            
        end
        
        
        handles.fPotentialMenu      =   uimenu(...       % File menu
            'Parent', handles.fSecondMenu,...
            'HandleVisibility','callback', ...
            'Label','Potential',...
            'Callback', @fSecondPotentialCallback);
        
        function fSecondPotentialCallback(~, ~)
            % Callback function run when the Potential menu item is selected
            prompt = {'Value in V(volts)'};
            dlg_title = 'Second Plate V';
            num_lines = 1;
            def = {num2str( objectpotential{2})};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            user_entry = str2double(answer);
            if isnan(user_entry)
                errordlg('You must enter a numeric value','Bad Input','modal');
                return;
            end
            if ~isempty(user_entry)
                objectpotential{2} = user_entry;
                cleaning;
            end
        end
        
        
        
        handles.fSecondParametersMenu      =   uimenu(...       % File menu
            'Parent', handles.fSecondMenu,...
            'HandleVisibility','callback', ...
            'Label','Parameters',...
            'Callback', @fSecondPreParametersCallback);
        
        function fSecondPreParametersCallback(~,~)
            
            switch handles.object2
                case 1
                    fSecondConductingPlateCallback;
                case 2
                    fSecondConductingCircleCallback;
                case 3
                    fSecondConductingSphereCallback;
                case 4
                    fSecondConductingBrickCallback;
                case 5
                    fSecondConductingCylinderCallback;
                otherwise
                    errordlg('Fatal error;Contact your Administrator');
            end
            
            
        end
        
        
        function fSecondConductingPlateCallback(~, ~)
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units',Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'parameters of the plate setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Plate length,m',strge.L(2);...
                'Plate width,m',strge.W(2);...
                'Center position x,m',strge.X(2);...
                'Center position y,m',strge.Y(2);...
                'Center position z,m',strge.Z(2);...
                'Rotate about x-axis, deg',strge.ax(2);...
                'Rotate about y-axis, deg',strge.ay(2);...
                'Rotate about z-axis, deg',strge.az(2);...
                'Approximate number of triangles',strge.Tr(2);...
                'Triangle size at boundaries vs. center size', strge.par(2)};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'You have chosen the plate case',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback});
            function edit_callback(hObject, ~)
                data = get(hObject,'Data');
                strge.L(2)         = data{11};
                strge.W(2)         = data{12};
                strge.X(2)         = data{13};
                strge.Y(2)         = data{14};
                strge.Z(2)         = data{15};
                strge.ax(2)        = data{16};
                strge.ay(2)        = data{17};
                strge.az(2)        = data{18};
                strge.Tr(2)        = data{19};
                strge.par(2) =       data{20};
                
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);
                delete(local);
            end
            
        end
        
        function fSecondConductingCircleCallback(~, ~)
            strcircle0 = strcircle;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'parameters of circle setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Circle radius,m',strcircle.R(2);...
                'Center position x,m',strcircle.X(2);...
                'Center position y,m',strcircle.Y(2);...
                'Center position z,m',strcircle.Z(2);...
                'Rotate about x-axis, deg',strcircle.ax(2);...
                'Rotate about y-axis, deg',strcircle.ay(2);...
                'Rotate about z-axis, deg',strcircle.az(2);...
                'Approximate numbers of triangles',strcircle.Tr(2);...
                'Triangle size at boundaries vs. center size', strcircle.par(2)};
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
                strcircle.R(2)         = data{10};
                strcircle.X(2)         = data{11};
                strcircle.Y(2)         = data{12};
                strcircle.Z(2)         = data{13};
                strcircle.ax(2)        = data{14};
                strcircle.ay(2)        = data{15};
                strcircle.az(2)        = data{16};
                strcircle.Tr(2)         = data{17};
                strcircle.par(2) = data{18};
                
                
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strcircle = strcircle0;
                delete(local);
            end
            
        end
        %% Sphere case
        
        function fSecondConductingSphereCallback(~, ~)
            strsphere0 = strsphere;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'parameters of sphere setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Sphere radius, m',strsphere.a(2);...
                'Center position x,m',strsphere.X(2);...
                'Center position y,m',strsphere.Y(2);...
                'Center position z,m',strsphere.Z(2);...
                'Approx.# of triangular patches',strsphere.Tr(2)};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'You have chosen the sphere case ',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback});
            function edit_callback(hObject, ~)
                data = get(hObject,'Data');
                strsphere.a(2)                 = data{6};
                strsphere.X(2)                 = data{7};
                strsphere.Y(2)                 = data{8};
                strsphere.Z(2)                 = data{9};
                strsphere.Tr(2)         = data{10};
                
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strsphere = strsphere0;
                delete(local);
            end
            
        end
        %% Brick Case
        
        function fSecondConductingBrickCallback(~, ~)
            strbrick0 = strbrick;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String',{'Parameters of brick setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Brick length,m',strbrick.L(2);...
                'Brick width,m',strbrick.W(2);...
                'Brick height,m',strbrick.H(2);...
                'Center position x,m',strbrick.X(2);...
                'Center position y,m',strbrick.Y(2);...
                'Center position z,m',strbrick.Z(2);...
                'Rotate about x-axis, deg',strbrick.ax(2);...
                'Rotate about y-axis, deg',strbrick.ay(2);...
                'Rotate about z-axis, deg',strbrick.az(2);...
                'Approximate # of triangles - top face',strbrick.Tr(2);...
                'Triangle size at boundaries vs. center size', strbrick.par(2)};
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
                strbrick.L(2)         = data{12};
                strbrick.W(2)         = data{13};
                strbrick.H(2)         = data{14};
                strbrick.X(2)         = data{15};
                strbrick.Y(2)         = data{16};
                strbrick.Z(2)         = data{17};
                strbrick.ax(2)        = data{18};
                strbrick.ay(2)        = data{19};
                strbrick.az(2)        = data{20};
                strbrick.Tr(2)        = data{21};
                strbrick.par(2) = data{22};
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);
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
        
        %% Cylinder case
        
        function fSecondConductingCylinderCallback(~, ~)
            strcylinder0 = strcylinder;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'parameters of cylinder setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Cylinder radius,m',strcylinder.R(2);...
                'Cylinder height,m',strcylinder.H(2);...
                'Center position x,m',strcylinder.X(2);...
                'Center position y,m',strcylinder.Y(2);...
                'Center position z,m',strcylinder.Z(2);...
                'Rotate about x-axis, deg',strcylinder.ax(2);...
                'Rotate about y-axis, deg',strcylinder.ay(2);...
                'Rotate about z-axis, deg',strcylinder.az(2);...
                'Approximate # of triangles - top face',strcylinder.Tr(2);...
                'Triangle size at boundaries vs. center size',strcylinder.par(2)};
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
                strcylinder.R(2)         = data{11};
                strcylinder.H(2)         = data{12};
                strcylinder.X(2)         = data{13};
                strcylinder.Y(2)         = data{14};
                strcylinder.Z(2)         = data{15};
                strcylinder.ax(2)        = data{16};
                strcylinder.ay(2)        = data{17};
                strcylinder.az(2)        = data{18};
                strcylinder.Tr(2)        = data{19};
                strcylinder.par(2)       = data{20};
                
                
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
        
        
        %% Conducting body
        
        % --- Conducting body Menu ----------------------------
        
        handles.fSecondMenu    =   uimenu(...       % File menu
            'Parent',handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Conducting Body');
        
        handles.fSecondObjectTypeMenu      =   uimenu(...       % File menu
            'Parent', handles.fSecondMenu,...
            'HandleVisibility','callback', ...
            'Label','Object Type',...
            'Callback', @fThirdObjectTypeCallback);
        
        
        function fThirdObjectTypeCallback(~,~)
            
            % --- RADIO BUTTONS -------------------------------------
            handles.figure3 = figure( ...
                'Tag', 'figure2', ...
                'Units', 'characters', ...
                'Position', [103.8 38.6153846153846 28.4 15.7692307692308], ...
                'Name', 'Object Type', ...
                'MenuBar', 'none', ...
                'NumberTitle', 'off', ...
                'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
                'Resize', 'on');
            
            handles.uipanel4 = uibuttongroup( ...
                'Parent', handles.figure3, ...
                'Tag', 'uipanel1', ...
                'Units', 'normalized', ...
                'Position', [0.0141823161189358 0.0334928229665072 0.965277777777778 0.961722488038277], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'Title', 'Object Type');
            
            set(handles.uipanel4,'SelectionChangeFcn',@uibuttongroup2_SelectionChangeFcn);
            
            
            function uibuttongroup2_SelectionChangeFcn(~,eventdata)
                switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
                    case 'plate'
                        handles.object3 = 1;
                        objecttype{3} = 'plate';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                        strout = clearoutput(strout);
                    case 'circle'
                        handles.object3 = 2;
                        objecttype{3} = 'circle';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                        strout = clearoutput(strout);
                    case 'sphere'
                        handles.object3 = 3;
                        objecttype{3} = 'sphere';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                        strout = clearoutput(strout);
                    case 'brick'
                        handles.object3 = 4;
                        objecttype{3} = 'brick';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                        strout = clearoutput(strout);
                    case 'cylinder'
                        handles.object3 = 5;
                        objecttype{3} = 'cylinder';
                        set(0, 'CurrentFigure', handles.figure1);
                        geometry();
                        strout = clearoutput(strout);
                    otherwise
                        errordlg('Fatal error;Contact your Administrator');
                end
            end
            
            
            handles.radiobutton7 = uicontrol( ...
                'Parent', handles.uipanel4, ...
                'Tag', 'plate', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.877777777777778 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Plate');
            
            
            handles.radiobutton8 = uicontrol( ...
                'Parent', handles.uipanel4, ...
                'Tag', 'circle', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.666666666666667 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Circle');
            
            
            
            handles.radiobutton9 = uicontrol( ...
                'Parent', handles.uipanel4, ...
                'Tag', 'sphere', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.455555555555556 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Sphere');
            
            
            
            handles.radiobutton10 = uicontrol( ...
                'Parent', handles.uipanel4, ...
                'Tag', 'brick', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.244444444444444 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Brick');
            
            
            
            handles.radiobutton11 = uicontrol( ...
                'Parent', handles.uipanel4, ...
                'Tag', 'cylinder', ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Position', [0.165413533834586 0.0333333333333333 0.654135338345865 0.127777777777778], ...
                'FontSize', 10.6666666666667, ...
                'FontUnits', 'pixels', ...
                'String', 'Cylinder');
            
            switch handles.object3
                case 1
                    set(handles.uipanel4,'Selectedobject',handles.radiobutton7);
                case 2
                    set(handles.uipanel4,'Selectedobject',handles.radiobutton8);
                case 3
                    set(handles.uipanel4,'Selectedobject',handles.radiobutton9);
                case 4
                    set(handles.uipanel4,'Selectedobject',handles.radiobutton10);
                case 5
                    set(handles.uipanel4,'Selectedobject',handles.radiobutton11);
                otherwise
                    errordlg('Fatal error;Contact your Administrator');
            end
            
        end
        
        
        handles.fPotentialMenu      =   uimenu(...       % File menu
            'Parent', handles.fSecondMenu,...
            'HandleVisibility','callback', ...
            'Label','Potential',...
            'Callback', @fThirdPotentialCallback);
        
        function fThirdPotentialCallback(~, ~)
            % Callback function run when the Potential menu item is selected
            prompt = {'Value in V(volts)'};
            dlg_title = 'Conducting body V';
            num_lines = 1;
            def = {num2str( objectpotential{3})};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            user_entry = str2double(answer);
            if isnan(user_entry)
                errordlg('You must enter a numeric value','Bad Input','modal');
                return;
            end
            if ~isempty(user_entry)
                objectpotential{3} = user_entry;
                cleaning;
            end
        end
        
        handles.fSecondParametersMenu      =   uimenu(...       % File menu
            'Parent', handles.fSecondMenu,...
            'HandleVisibility','callback', ...
            'Label','Parameters',...
            'Callback', @fThirdPreParametersCallback);
        
        function fThirdPreParametersCallback(~,~)
            
            switch handles.object3
                case 1
                    fThirdConductingPlateCallback;
                case 2
                    fThirdConductingCircleCallback;
                case 3
                    fThirdConductingSphereCallback;
                case 4
                    fThirdConductingBrickCallback;
                case 5
                    fThirdConductingCylinderCallback;
                otherwise
                    errordlg('Fatal error;Contact your Administrator');
            end
            
            
        end
        
        
        function fThirdConductingPlateCallback(~, ~)
            strge0 = strge;
            objectcharge0 = objectcharge;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'parameters of the plate setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Plate length,m',strge.L(3);...
                'Plate width,m',strge.W(3);...
                'Center position x,m',strge.X(3);...
                'Center position y,m',strge.Y(3);...
                'Center position z,m',strge.Z(3);...
                'Rotate about x-axis, deg',strge.ax(3);...
                'Rotate about y-axis, deg',strge.ay(3);...
                'Rotate about z-axis, deg',strge.az(3);...
                'Approximate number of triangles',strge.Tr(3);...
                'Triangle size at boundaries vs. center size', strge.par(3);...
                'Zero the charge of the body(yes/no)',objectcharge};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'You have chosen the plate case',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback});
            function edit_callback(hObject, ~)
                data = get(hObject,'Data');
                strge.L(3)         = data{12};
                strge.W(3)         = data{13};
                strge.X(3)         = data{14};
                strge.Y(3)         = data{15};
                strge.Z(3)         = data{16};
                strge.ax(3)        = data{17};
                strge.ay(3)        = data{18};
                strge.az(3)        = data{19};
                strge.Tr(3)        = data{20};
                strge.par(3)       = data{21};
                objectcharge       = data{22};                
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strge = strge0;
                objectcharge = objectcharge0;
                delete(local);
            end
            
        end
        
        %% Circle case skipped
        function fThirdConductingCircleCallback(~, ~)
            strcircle0 = strcircle;
            objectcharge0 = objectcharge;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'parameters of circle setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Circle radius,m',strcircle.R(3);...
                'Center position x,m',strcircle.X(3);...
                'Center position y,m',strcircle.Y(3);...
                'Center position z,m',strcircle.Z(3);...
                'Rotate about x-axis, deg',strcircle.ax(3);...
                'Rotate about y-axis, deg',strcircle.ay(3);...
                'Rotate about z-axis, deg',strcircle.az(3);...
                'Approximate numbers of triangles',strcircle.Tr(3);...
                'Triangle size at boundaries vs. center size', strcircle.par(3);...
                'Zero the charge of the body(yes/no)',objectcharge};
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
                strcircle.R(3)         = data{11};
                strcircle.X(3)         = data{12};
                strcircle.Y(3)         = data{13};
                strcircle.Z(3)         = data{14};
                strcircle.ax(3)        = data{15};
                strcircle.ay(3)        = data{16};
                strcircle.az(3)        = data{17};
                strcircle.Tr(3)         = data{18};
                strcircle.par(3)        = data{19};
                objectcharge           = data{20};
                
                
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strcircle = strcircle0;
                objectcharge = objectcharge0;
                delete(local);
            end
            
        end
        %% Sphere case
        
        function fThirdConductingSphereCallback(~, ~)
            strsphere0 = strsphere;
            objectcharge_temp = objectcharge;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'parameters of sphere setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Sphere radius, m',strsphere.a(3);...
                'Center Position x,m',strsphere.X(3);...
                'Center Position y,m',strsphere.Y(3);...
                'Center Position z,m',strsphere.Z(3);...
                'Approx.# of triangular patches',strsphere.Tr(3);...
                'Zero the charge of the body(yes/no)',objectcharge};
            ColumnName =   {'Name',   'Value'};
            columnformat    = {'char'};
            columneditable  =  [true];
            uitable('Units', 'normalized', 'Position', Position_table,...
                'Data', data, 'Visible', 'on', 'BackgroundColor', [1 1 1],...
                'RowName', [], 'ColumnName', ColumnName, 'ColumnFormat', columnformat, ...
                'ColumnEditable', columneditable, 'ToolTipString', 'You have chosen the sphere case ',...
                'ColumnWidth',{200,100},...
                'CellEditCallback', {@edit_callback});
            function edit_callback(hObject, ~)
                data = get(hObject,'Data');
                strsphere.a(3)                 = data{7};
                strsphere.X(3)                 = data{8};
                strsphere.Y(3)                 = data{9};
                strsphere.Z(3)                 = data{10};
                strsphere.Tr(3)                = data{11};
                objectcharge                   = data{12};
                
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strsphere = strsphere0;
                objectcharge = objectcharge_temp;
                delete(local);
            end
            
        end
        %% Brick Case
        
        function fThirdConductingBrickCallback(~, ~)
            strbrick0 = strbrick;
            objectcharge_temp = objectcharge;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'Parameters of brick setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Brick length,m',strbrick.L(3);...
                'Brick width,m',strbrick.W(3);...
                'Brick height,m',strbrick.H(3);...
                'Center position x,m',strbrick.X(3);...
                'Center position y,m',strbrick.Y(3);...
                'Center position z,m',strbrick.Z(3);...
                'Rotate about x-axis, deg',strbrick.ax(3);...
                'Rotate about y-axis, deg',strbrick.ay(3);...
                'Rotate about z-axis, deg',strbrick.az(3);...
                'Approximate # of triangles - top face',strbrick.Tr(3);...
                'Triangle size at boundaries vs. center size', strbrick.par(3);...
                'Zero the charge of the body(yes/no)',objectcharge};
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
                strbrick.L(3)         = data{13};
                strbrick.W(3)         = data{14};
                strbrick.H(3)         = data{15};
                strbrick.X(3)         = data{16};
                strbrick.Y(3)         = data{17};
                strbrick.Z(3)         = data{18};
                strbrick.ax(3)        = data{19};
                strbrick.ay(3)        = data{20};
                strbrick.az(3)        = data{21};
                strbrick.Tr(3)        = data{22};
                strbrick.par(3)       = data{23};
                objectcharge            = data{24};
                
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strbrick = strbrick0;
                objectcharge = objectcharge_temp;
                delete(local);
            end
            
        end
        
        %% Cylinder case
        
        function fThirdConductingCylinderCallback(~, ~)
            strcylinder0 = strcylinder;
            objectcharge_temp = objectcharge;
            % Callback function run when_______________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            local = figure('Units', Position_units, 'Position', Position_local);
            set(local, 'Name', 'Conductor setting', 'NumberTitle', 'off', 'MenuBar', 'none');
            b1 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Update','Position', Position_b1, 'Parent', local, 'Callback', @B1_Callback);
            b2 = uicontrol('Style', 'Pushbutton', 'Units', 'normalized', 'String', 'Cancel', 'Position', Position_b2, 'Parent', local, 'Callback', @B2_Callback);
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', {'parameters of cylinder setup';'Press ENTER after changing each parameter'}, 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = {...
                'Cylinder radius,m',strcylinder.R(3);...
                'Cylinder height,m',strcylinder.H(3);...
                'Center position x,m',strcylinder.X(3);...
                'Center position y,m',strcylinder.Y(3);...
                'Center position z,m',strcylinder.Z(3);...
                'Rotate about x-axis, deg',strcylinder.ax(3);...
                'Rotate about y-axis, deg',strcylinder.ay(3);...
                'Rotate about z-axis, deg',strcylinder.az(3);...
                'Approximate # of triangles - top face',strcylinder.Tr(3);...
                'Triangle size at boundaries vs. center size',strcylinder.par(3);...
                'Zero the charge of the body(yes/no)',objectcharge};
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
                strcylinder.R(3)         = data{12};
                strcylinder.H(3)         = data{13};
                strcylinder.X(3)         = data{14};
                strcylinder.Y(3)         = data{15};
                strcylinder.Z(3)         = data{16};
                strcylinder.ax(3)        = data{17};
                strcylinder.ay(3)        = data{18};
                strcylinder.az(3)        = data{19};
                strcylinder.Tr(3)        = data{20};
                strcylinder.par(3)       =data{21};
                objectcharge             = data{22};
                
            end
            
            function B1_Callback(~, ~, ~)
                set(0, 'CurrentFigure', handles.figure1);
                geometry();
                strout = clearoutput(strout);
                delete(local);
            end
            function B2_Callback(~, ~, ~)
                strcylinder = strcylinder0;
                objectcharge = objectcharge_temp;
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
        
        % --- Modeling Setup ------
        
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
            b3 = uicontrol('Style', 'Text', 'Units', 'normalized', 'String', 'Output data and results', 'Position', Position_b3, 'Parent', local, 'BackgroundColor', [0.8, 0.8, 0.8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ColumnName =   {'Name',   'Value'};
            columneditable  =  [true];
            data = {...
                'Numerical capacitance,pF',strout.Cnum;...           
                'Total charge of the first conductor, C',strout.ChargeLower;...
                'Total charge of the second conductor, C',strout.ChargeUpper;...
                'Total charge of the body, C',strout.ChargeObject;...
                'Total charge of the entire structure, C',strout.ChargeTotal;...
                'E field x component at obs. point, V/m',strout.E(1);...
                'E field y component at obs. point, V/m',strout.E(2);...
                'E field z component at obs. point, V/m',strout.E(3);...
                'Electric potential at obs. point, V',strout.phi;...
                'Average electric potential of the third object, V',strout.ObjectPotential;...
                'Total number of triangular patches in the capacitor mesh',strout.PatchesCap;...
                'Total number of triangular patches in all conductors',strout.PatchesTotal;...
                'Minimum triangle quality in the mesh',strout.quality;...
                'CPU time in sec for filling the MoM matrix',strout.time1;...
                'CPU time in sec for solving the system of MoM eqs.',strout.time2};
               
            if (strcmp(stroc.yes,'no'))== 1; data{20} = 'N.A.'; data{21} = 'N.A.'; data{22} = 'N.A.'; data{23} = 'N.A.';  end;
            
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
                        view(-47, 14);
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
                'String', 'Current');
            
            
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
        
        %--- help menu-----
        handles.fHelpMenu    =   uimenu(...       % File menu
            'Parent',handles.figure1,...
            'HandleVisibility','callback', ...
            'Label','Help',...
            'callback',@fHelpSetupCallback);
        
        function fHelpSetupCallback(~, ~)
            HelpPath = 'index_s24.html';
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
            'string','This module is an accurate 3D electrostatic MoM solution for multi-conductor structures. Three objects may be placed inside a capacitor: a conducting brick, a sphere, or a cylinder. The capacitor itself may be formed by objects of arbitrary shape. The module also computes charge distributions on conductors, electric field distribution in a plane, and electric field at a point of interest (marked by a cross).');
        
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
                    Save_Callback();
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

function d = dcircle(P, R)
    %   Compute signed distance function for a circle with radius R 
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson (PhD thesis, pp. 38-39)

    d = sqrt(P(:, 1).^2 + P(:, 2).^2) - R;
end

function d = drectangle(P, L, W)
    %   Compute signed distance function for a rectangle with length L and
    %   width W
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson

    x1 = -L/2;  x2 = +L/2;
    y1 = -W/2;  y2 = +W/2;
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

function [ ] = graphics_S22(io, strop, stroc, strei, P, t, c);
    %%   Program I/O graphics
    %    variable io is zero for the input graphics and one for the output graphics 
  
    %%  Prepare common variables
    Ptemp = P'; ttemp = t'; ctemp = c'; ind = size(t, 1);  IND = 10000;
    X = reshape(Ptemp(1, ttemp(1:3, :)),[3, size(ttemp, 2)]);
    Y = reshape(Ptemp(2, ttemp(1:3, :)),[3, size(ttemp, 2)]);
    Z = reshape(Ptemp(3, ttemp(1:3, :)),[3, size(ttemp, 2)]);       

    %%   Draw the observation plane
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
    %%   Draw electric field (incident or final) in the observation plane
    if strcmp(strop.yes, 'yes')
        strop.E = strop.E + 1e-12*max(max(strop.E))*rand(size(strop.E)); %   randomize a bit
        normalize = strei.normalize/strop.arrow;
        fv = cone(strop.Points(:, 1), strop.Points(:, 2), strop.Points(:, 3), strop.E(:, 1), strop.E(:, 2), strop.E(:, 3), 36, 0.25, 1, normalize);
        patch(fv, 'FaceColor', 'y', 'EdgeColor', 'none'); 
    end
    %%  Draw the observation point(s)
%     stroc.yes         = 'yes';     %   'yes' - present; 'no' - absent
%     stroc.x           = 0.0;       %   x position in m
%     stroc.y           = 0.0;       %   y position in m
%     stroc.z           = 0.0;       %   z position in m
%     stroc.size = 1;                %   relative marker size versus default size 
    if strcmp(stroc.yes, 'yes')
        scaling = 0.075*max(strop.planesizex, strop.planesizey)*stroc.size;
        fv      = kreuz(stroc.x, stroc.y, stroc.z, 4, 0.1*scaling, 1*scaling);
        patch(fv, 'FaceColor', 'g', 'EdgeColor', 'k'); 
    end
    hold on
    %%   Input graphics (within the main GUI window)    
    if io == 0
        colormap jet;
        if ind<IND ec='y'; else ec='none'; end
        h = fill3(X, Y, Z, [0 0 1], 'FaceAlpha', 0.5, 'EdgeColor', ec); xlabel('x, m'); ylabel('y, m');
        zlabel('z, m'); title('Problem geometry: conducting object under study and the incident E-field');
        colorbar; colorbar('YTickLabel', {'0', '0', '0', '0'});
    else
        %   Interpolate charge density for vertexes
        colormap jet;        
        cv = zeros(1, length(Ptemp));
        for m = 1:size(Ptemp, 2)
            [p, q] = find(ttemp(1:3, :)-m==0);
            if (length(q))
                cv(m) = sum(ctemp(q), 2)/length(q);
            end    
        end          
        if strcmp(strop.yes, 'yes')
            C = cv(ttemp(1:3, :));   
            if ind<IND ec='w'; else ec='none'; end
            h = fill3(X, Y, Z, C, 'FaceAlpha', 0.5, 'EdgeColor', ec); xlabel('x, m'); ylabel('y, m');
            zlabel('z, m'); title('Solution: free charge density in C/m^2 and the resulting E-field')
            colorbar;
        else
            C = cv(ttemp(1:3, :));  
            if ind<IND ec='w'; else ec='none'; end
            h = fill3(X, Y, Z, C, 'FaceAlpha', 1.0, 'EdgeColor', ec);
            zlabel('z, m'); title('Solution: free charge density in C/m^2')
            colorbar;            
        end
        if strcmp(strop.potential, 'yes')
            contourca(strop, strop.Potential);
        end
    end
    
    %%   General settings 
    axis('equal'); axis('tight'); grid on; 
    if strcmp(strop.yes, 'yes')
        camproj('orthographic');
        light('Position',[1 1 1]); light('Position',[-3 -1 3]); material shiny; 
    end
    view(165, 15);
end

function [ ] = graphics_S23(io, strop, stroc, P, t, c, Area, strsc);
    %%   Program I/O graphics
    %    variable io is zero for the input graphics and one for the output graphics 

    %%  Prepare common variables
    Ptemp = P'; ttemp = t'; ind = size(t, 1);  IND = 10000;
    X = reshape(Ptemp(1, ttemp(1:3, :)),[3, size(ttemp, 2)]);
    Y = reshape(Ptemp(2, ttemp(1:3, :)),[3, size(ttemp, 2)]);
    Z = reshape(Ptemp(3, ttemp(1:3, :)),[3, size(ttemp, 2)]);      
    
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
%     stroc.yes         = 'yes';     %   'yes' - present; 'no' - absent
%     stroc.x           = 0.0;       %   x position in m
%     stroc.y           = 0.0;       %   y position in m
%     stroc.z           = 0.0;       %   z position in m
%     stroc.size = 1;                %   relative marker size versus default size 
    if strcmp(stroc.yes, 'yes')
        scaling = 0.05*max(strop.planesizex, strop.planesizey)*stroc.size;
        fv      = kreuz(stroc.x, stroc.y, stroc.z, 4, 0.1*scaling, 1*scaling);
        patch(fv, 'FaceColor', 'g', 'EdgeColor', 'k'); 
    end
    hold on
    %%   Input graphics (within the main GUI window)    
    if io == 0
        colormap jet;
        if ind<IND ec='y'; else ec='none'; end
        patch('Vertices', P, 'Faces', t(:, 1:3), 'FaceColor',  'b', 'FaceAlpha', 1.00, 'EdgeColor', ec)     
        xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); 
        title('Problem geometry: conducting objects');
    end
    %%   Output graphics (separate window)
    if io == 1
        %%  Charge scaling (scaled variable is "charge")  
        cplus       = max(c);
        cminus      = min(c);        
        charge      = c;
        charge(charge>=+cplus*strsc.positive)  = cplus*strsc.positive;
        charge(charge<=+cminus*strsc.negative) = cminus*strsc.negative;
        cplus       = max(charge);
        cminus      = min(charge);
        if cplus-cminus<=0; charge = c; end;
        
        %%   Interpolate charge density for vertexes - global   
        ctemp = charge';
        cv = zeros(1, length(Ptemp));
        for m = 1:size(Ptemp, 2)
            [p, q] = find(ttemp(1:3, :)-m==0);
            if (length(q))
                cv(m) = sum(ctemp(q), 2)/length(q);
            end    
        end                 
        C = cv(ttemp(1:3, :)); 
        if (~strcmp(strop.yes, 'yes'))&(~strcmp(strop.potential, 'yes'))
            if ind<IND ec='k'; else ec='none'; end
            patch(X, Y, Z, C, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');          
            xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); colorbar;                
            title('Solution: Surface charge density on all objects in C/m^2');
        end        
        if strcmp(strop.yes, 'yes')|strcmp(strop.potential, 'yes')
            if ind<IND ec='k'; else ec='none'; end
            patch(X, Y, Z, C, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');
            xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); colorbar;
            title('Solution: Surface charge density in C/m^2 and field data in the observation plane'); 
        end
        if strcmp(strop.yes, 'yes')
            %%   Draw electric field (final) in the observation plane
            strop.E = strop.E + 1e-12*max(max(strop.E))*rand(size(strop.E)); %   randomize a bit
            normalize = sqrt(sum(strop.E.*strop.E, 2));
            normalize = max(normalize);
            normalize = normalize/strop.arrow;
            fv = cone(strop.Points(:, 1), strop.Points(:, 2), strop.Points(:, 3), strop.E(:, 1), strop.E(:, 2), strop.E(:, 3), 36, 0.25, 1, normalize);
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
    view(-47, 20);
end

function [ ] = graphics_S24(io, strop, stroc, P, t, c, Area, strsc);
    %%   Program I/O graphics
    %    variable io is zero for the input graphics and one for the output graphics 

    %%  Prepare common variables
    tM = t((t(:, 4)==1)|(t(:, 4)==2), :);
    tD = t(t(:, 4)==3,  :);
    cM = c((t(:, 4)==1)|(t(:, 4)==2), :);
    cD = c(t(:, 4)==3,  :);
    
    Ptemp = P'; ttemp = t'; ind = size(t, 1);  IND = 10000;
    
    ctempM = cM'; ttempM = tM'; 
    XM = reshape(Ptemp(1, ttempM(1:3, :)),[3, size(ttempM, 2)]);
    YM = reshape(Ptemp(2, ttempM(1:3, :)),[3, size(ttempM, 2)]);
    ZM = reshape(Ptemp(3, ttempM(1:3, :)),[3, size(ttempM, 2)]);  
    
    ctempD = cD'; ttempD = tD'; 
    XD = reshape(Ptemp(1, ttempD(1:3, :)),[3, size(ttempD, 2)]);
    YD = reshape(Ptemp(2, ttempD(1:3, :)),[3, size(ttempD, 2)]);
    ZD = reshape(Ptemp(3, ttempD(1:3, :)),[3, size(ttempD, 2)]);  
    
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
%     stroc.yes         = 'yes';     %   'yes' - present; 'no' - absent
%     stroc.x           = 0.0;       %   x position in m
%     stroc.y           = 0.0;       %   y position in m
%     stroc.z           = 0.0;       %   z position in m
%     stroc.size = 1;                %   relative marker size versus default size 
    if strcmp(stroc.yes, 'yes')
        scaling = 0.05*max(strop.planesizex, strop.planesizey)*stroc.size;
        fv      = kreuz(stroc.x, stroc.y, stroc.z, 4, 0.1*scaling, 1*scaling);
        patch(fv, 'FaceColor', 'g', 'EdgeColor', 'k'); 
    end
    hold on
    %%   Input graphics (within the main GUI window)    
    if io == 0
        colormap jet;
        if ind<IND ec='b'; else ec='none'; end
        patch('Vertices', P, 'Faces', tD(:, 1:3), 'FaceColor',  'y', 'FaceAlpha', 1.00, 'EdgeColor', ec);
        if ind<IND ec='y'; else ec='none'; end
        patch('Vertices', P, 'Faces', tM(:, 1:3), 'FaceColor',  'b', 'FaceAlpha', 1.00, 'EdgeColor', ec); 
        xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); 
        title('Problem geometry: two conductors and the conducting object');
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
        
        if (~strcmp(strop.yes, 'yes'))&(~strcmp(strop.potential, 'yes'))
            if ind<IND ec='k'; else ec='none'; end
            patch(XD, YD, ZD, CD, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');
            patch(XM, YM, ZM, CM, 'FaceAlpha', 1.0, 'EdgeColor', ec,  'FaceLighting', 'none');
            xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); colorbar;
            title('Solution: Surface charge density on all objects in C/m^2')
        end        
        if strcmp(strop.yes, 'yes')|strcmp(strop.potential, 'yes')
            if ind<IND ec='k'; else ec='none'; end
            patch(XD, YD, ZD, CD, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');
            patch(XM, YM, ZM, CM, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');
            xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); colorbar;
            title('Solution: Surface charge density in C/m^2 and the field in the observation plane'); 
        end
        if strcmp(strop.yes, 'yes')
            %%   Draw electric field (final) in the observation plane
            strop.E = strop.E + 1e-12*max(max(strop.E))*rand(size(strop.E)); %   randomize a bit
            normalize = sqrt(sum(strop.E.*strop.E, 2));
            normalize = max(normalize);
            normalize = normalize/strop.arrow;
            fv = cone(strop.Points(:, 1), strop.Points(:, 2), strop.Points(:, 3), strop.E(:, 1), strop.E(:, 2), strop.E(:, 3), 36, 0.25, 1, normalize);
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
    view(-47, 14);
end

function [ ] = graphics_S31(io, strop, stroc, strei, P, t, c);
    %%   Program I/O graphics
    %    variable io is zero for the input graphics and one for the output graphics 

    %%  Prepare common variables
    Ptemp = P'; ttemp = t'; ctemp = c'; ind = size(t, 1);  IND = 10000;
    X = reshape(Ptemp(1, ttemp(1:3, :)),[3, size(ttemp, 2)]);
    Y = reshape(Ptemp(2, ttemp(1:3, :)),[3, size(ttemp, 2)]);
    Z = reshape(Ptemp(3, ttemp(1:3, :)),[3, size(ttemp, 2)]);       

    %%   Draw the observation plane
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
    %%   Draw electric field (incident or final) in the observation plane
    if strcmp(strop.yes, 'yes')
        strop.E = strop.E + 1e-12*max(max(strop.E))*rand(size(strop.E)); %   randomize a bit
        normalize = strei.normalize/strop.arrow;
        fv = cone(strop.Points(:, 1), strop.Points(:, 2), strop.Points(:, 3), strop.E(:, 1), strop.E(:, 2), strop.E(:, 3), 36, 0.25, 1, normalize);
        patch(fv, 'FaceColor', 'y', 'EdgeColor', 'none'); 
    end
    %%  Draw the observation point(s)
%     stroc.yes         = 'yes';     %   'yes' - present; 'no' - absent
%     stroc.x           = 0.0;       %   x position in m
%     stroc.y           = 0.0;       %   y position in m
%     stroc.z           = 0.0;       %   z position in m
%     stroc.size = 1;                %   relative marker size versus default size 
    if strcmp(stroc.yes, 'yes')
        scaling = 0.075*max(strop.planesizex, strop.planesizey)*stroc.size;
        fv      = kreuz(stroc.x, stroc.y, stroc.z, 4, 0.1*scaling, 1*scaling);
        patch(fv, 'FaceColor', 'g', 'EdgeColor', 'k'); 
    end
    hold on
    %%   Input graphics (within the main GUI window)    
    if io == 0
        colormap jet;
        if ind<IND ec='w'; else ec='none'; end
        h = fill3(X, Y, Z, [0 0 1], 'FaceAlpha', 0.5, 'FaceColor', [0 0.7 0], 'EdgeColor', ec); 
        xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); title('Problem geometry: dielectric object under study and the incident E-field');
        colorbar; colorbar('YTickLabel', {'0', '0', '0', '0'});
    else
        %   Interpolate charge density for vertexes
        colormap jet;        
        cv = zeros(1, length(Ptemp));
        for m = 1:size(Ptemp, 2)
            [p, q] = find(ttemp(1:3, :)-m==0);
            if (length(q))
                cv(m) = sum(ctemp(q), 2)/length(q);
            end    
        end          
        if strcmp(strop.yes, 'yes')
            C = cv(ttemp(1:3, :));    
            if ind<IND ec='w'; else ec='none'; end
            h = fill3(X, Y, Z, C, 'FaceAlpha', 0.5, 'EdgeColor', ec);
            xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); title('Solution: polarization charge density in C/m^2 and the resulting E-field')
            colorbar;
        else
            C = cv(ttemp(1:3, :));    
            h = fill3(X, Y, Z, C, 'FaceAlpha', 1.0, 'EdgeColor', 'w');
            xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); title('Solution: polarization charge density in C/m^2')
            colorbar;            
        end
        if strcmp(strop.potential, 'yes')
            contourca(strop, strop.Potential);
        end
    end
    
    %%   General settings 
    axis('equal'); axis('tight'); grid on; 
    if strcmp(strop.yes, 'yes')
        camproj('orthographic');
        light('Position',[1 3 2]); light('Position',[-3 -1 3]); material shiny; 
    end
    view(165, 15);
end

function [ ] = graphics_S32(io, strop, stroc, strei, P, t, c, Area, strsc);
    %%   Program I/O graphics
    %    variable io is zero for the input graphics and one for the output graphics 

    %%  Prepare common variables
    tM = t(t(:, 4)<3, :);
    tD = t(t(:, 4)>2,  :);
    cM = c(t(:, 4)<3, :);
    cD = c(t(:, 4)>2,  :);
    
    Ptemp = P'; ttemp = t'; ind = size(t, 1);  IND = 10000;
    
    ctempM = cM'; ttempM = tM'; 
    XM = reshape(Ptemp(1, ttempM(1:3, :)),[3, size(ttempM, 2)]);
    YM = reshape(Ptemp(2, ttempM(1:3, :)),[3, size(ttempM, 2)]);
    ZM = reshape(Ptemp(3, ttempM(1:3, :)),[3, size(ttempM, 2)]);  
    
    ctempD = cD'; ttempD = tD'; 
    XD = reshape(Ptemp(1, ttempD(1:3, :)),[3, size(ttempD, 2)]);
    YD = reshape(Ptemp(2, ttempD(1:3, :)),[3, size(ttempD, 2)]);
    ZD = reshape(Ptemp(3, ttempD(1:3, :)),[3, size(ttempD, 2)]);  
    
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
        if ind<IND ec='b'; else ec='none'; end
        patch('Vertices', P, 'Faces', tD(:, 1:3), 'FaceColor',  'g', 'FaceAlpha', 1.00, 'EdgeColor', 'b')
        if ind<IND ec='y'; else ec='none'; end
        patch('Vertices', P, 'Faces', tM(:, 1:3), 'FaceColor',  'b', 'FaceAlpha', 1.00, 'EdgeColor', 'y'); 
        xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); 
        title('Problem geometry: capacitor plates and the dielectric object');
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
     
        if (~strcmp(strop.yes, 'yes'))&(~strcmp(strop.potential, 'yes'))
            if ind<IND ec='k'; else ec='none'; end
            patch(XD, YD, ZD, CD, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');
            patch(XM, YM, ZM, CM, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');
            xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); colorbar;
            title('Solution: Surface charge density on all objects in C/m^2')
        end        
        if strcmp(strop.yes, 'yes')|strcmp(strop.potential, 'yes')
            if ind<IND ec='k'; else ec='none'; end
            patch(XD, YD, ZD, CD, 'FaceAlpha', 0.1, 'EdgeColor', ec, 'FaceLighting', 'none');
            patch(XM, YM, ZM, CM, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');
            xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); colorbar;
            title('Solution: Surface charge density in C/m^2 and the field in the observation plane'); 
        end
        if strcmp(strop.yes, 'yes')
            %%   Draw electric field (final) in the observation plane
            strop.E = strop.E + 1e-12*max(max(strop.E))*rand(size(strop.E)); %   randomize a bit
            normalize = sqrt(sum(strop.E.*strop.E, 2));
            normalize = max(normalize);
            normalize = normalize/strop.arrow;
            fv = cone(strop.Points(:, 1), strop.Points(:, 2), strop.Points(:, 3), strop.E(:, 1), strop.E(:, 2), strop.E(:, 3), 36, 0.25, 1, normalize);
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
    view(53, 23);
end

function [ ] = graphics_S40(io, strop, stroc, P, t, c, Area, strsc);
    %%   Program I/O graphics
    %    variable io is zero for the input graphics and one for the output graphics 
  
    %%  Prepare common variables
    Ptemp = P'; ttemp = t'; ctemp = c'; ind = size(t, 1);  IND = 3000;
    X = reshape(Ptemp(1, ttemp(1:3, :)),[3, size(ttemp, 2)]);
    Y = reshape(Ptemp(2, ttemp(1:3, :)),[3, size(ttemp, 2)]);
    Z = reshape(Ptemp(3, ttemp(1:3, :)),[3, size(ttemp, 2)]);       
    
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
    %     stroc.yes         = 'yes';     %   'yes' - present; 'no' - absent
    %     stroc.x           = 0.0;       %   x position in m
    %     stroc.y           = 0.0;       %   y position in m
    %     stroc.z           = 0.0;       %   z position in m
    %     stroc.size = 1;                %   relative marker size versus default size 
    if strcmp(stroc.yes, 'yes')
        scaling = 0.05*max(strop.planesizex, strop.planesizey)*stroc.size;
        fv      = kreuz(stroc.x, stroc.y, stroc.z, 4, 0.1*scaling, 1*scaling);
        patch(fv, 'FaceColor', 'g', 'EdgeColor', 'k'); 
    end
    hold on
    %%   Input graphics (within the main GUI window)    
    if io == 0 
        if ind<IND ec='y'; else ec='none'; end
        patch(X, Y, Z, [1 0.75 0.65], 'FaceAlpha', 1.00, 'EdgeColor', ec)     
        axis('equal'); axis('tight'); camproj('orthographic');
        light('Position',[1 3 2]); light('Position',[-3 -1 3]); material shiny; grid on;
        xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); 
        title('Problem geometry: conducting body');
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
        cv = zeros(1, length(Ptemp));
        for m = 1:size(Ptemp, 2)
            [p, q] = find(ttemp(1:3, :)-m==0);
            if (length(q))
                cv(m) = sum(ctemp(q), 2)/length(q);
            end    
        end                 
        C = cv(ttemp(1:3, :)); 
        if (~strcmp(strop.yes, 'yes'))&(~strcmp(strop.potential, 'yes'))
            if ind<IND ec='k'; else ec='none'; end
            patch(X, Y, Z, C, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');          
            xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); colorbar;
            title('Solution: Surface charge density on the body surface in C/m^2')
        end        
        if strcmp(strop.yes, 'yes')|strcmp(strop.potential, 'yes')
            if ind<IND ec='k'; else ec='none'; end
            patch(X, Y, Z, C, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');
            xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); colorbar;
            title('Solution: Surface charge density in C/m^2 and field data in the observation plane'); 
        end
        if strcmp(strop.yes, 'yes')
            %%   Draw electric field (final) in the observation plane
            strop.E = strop.E + 1e-12*max(max(strop.E))*rand(size(strop.E)); %   randomize a bit
            normalize = sqrt(sum(strop.E.*strop.E, 2));
            normalize = max(normalize);
            normalize = normalize/strop.arrow;
            fv = cone(strop.Points(:, 1), strop.Points(:, 2), strop.Points(:, 3), strop.E(:, 1), strop.E(:, 2), strop.E(:, 3), 36, 0.25, 1, normalize);
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
    view(151, 20);
end

function [ ] = graphics_S42(io, strop, stroc, P, t, c, Area, strsc, strei);
    %%   Program I/O graphics
    %    variable io is zero for the input graphics and one for the output graphics 
  
    %%  Prepare common variables
    Ptemp = P'; ttemp = t'; ctemp = c'; ind = size(t, 1);  IND = 3000;
    X = reshape(Ptemp(1, ttemp(1:3, :)),[3, size(ttemp, 2)]);
    Y = reshape(Ptemp(2, ttemp(1:3, :)),[3, size(ttemp, 2)]);
    Z = reshape(Ptemp(3, ttemp(1:3, :)),[3, size(ttemp, 2)]);       
    
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
    
      %%   Draw electric field (incident or final) in the observation plane
    if strcmp(strop.yes, 'yes')
        strop.E = strop.E + 1e-12*max(max(strop.E))*rand(size(strop.E)); %   randomize a bit
        normalize = strei.normalize/strop.arrow;
        fv = cone(strop.Points(:, 1), strop.Points(:, 2), strop.Points(:, 3), strop.E(:, 1), strop.E(:, 2), strop.E(:, 3), 36, 0.25, 1, normalize);
        patch(fv, 'FaceColor', 'y', 'EdgeColor', 'none'); 
    end
    
    %%  Draw the observation point(s)
    %     stroc.yes         = 'yes';     %   'yes' - present; 'no' - absent
    %     stroc.x           = 0.0;       %   x position in m
    %     stroc.y           = 0.0;       %   y position in m
    %     stroc.z           = 0.0;       %   z position in m
    %     stroc.size = 1;                %   relative marker size versus default size 
    if strcmp(stroc.yes, 'yes')
        scaling = 0.05*max(strop.planesizex, strop.planesizey)*stroc.size;
        fv      = kreuz(stroc.x, stroc.y, stroc.z, 4, 0.1*scaling, 1*scaling);
        patch(fv, 'FaceColor', 'g', 'EdgeColor', 'k'); 
    end
    hold on
    %%   Input graphics (within the main GUI window)    
    if io == 0 
        if ind<IND ec='y'; else ec='none'; end
        patch(X, Y, Z, [1 0.75 0.65], 'FaceAlpha', 1.00, 'EdgeColor', ec)     
        axis('equal'); axis('tight'); camproj('orthographic');
        light('Position',[1 3 2]); light('Position',[-3 -1 3]); material shiny; grid on;
        xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); 
        title('Problem geometry: conducting body');
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
        cv = zeros(1, length(Ptemp));
        for m = 1:size(Ptemp, 2)
            [p, q] = find(ttemp(1:3, :)-m==0);
            if (length(q))
                cv(m) = sum(ctemp(q), 2)/length(q);
            end    
        end                 
        C = cv(ttemp(1:3, :)); 
        if (~strcmp(strop.yes, 'yes'))&(~strcmp(strop.potential, 'yes'))
            if ind<IND ec='k'; else ec='none'; end
            patch(X, Y, Z, C, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');          
            xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); colorbar;
            title('Solution: Surface charge density on the body surface in C/m^2')
        end        
        if strcmp(strop.yes, 'yes')|strcmp(strop.potential, 'yes')
            if ind<IND ec='k'; else ec='none'; end
            patch(X, Y, Z, C, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');
            xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); colorbar;
            title('Solution: Surface charge density in C/m^2 and field data in the observation plane'); 
        end
        if strcmp(strop.yes, 'yes')
            %%   Draw electric field (final) in the observation plane
            strop.E = strop.E + 1e-12*max(max(strop.E))*rand(size(strop.E)); %   randomize a bit
            normalize = sqrt(sum(strop.E.*strop.E, 2));
            normalize = max(normalize);
            normalize = normalize/strop.arrow;
            fv = cone(strop.Points(:, 1), strop.Points(:, 2), strop.Points(:, 3), strop.E(:, 1), strop.E(:, 2), strop.E(:, 3), 36, 0.25, 1, normalize);
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
    view(45, 15);
end

function [ ] = graphics_S44(io, strop, stroc, P, t, c, Area, strsc);
    %%   Program I/O graphics
    %    variable io is zero for the input graphics and one for the output graphics 
     %%  Prepare common variables
    tM = t((t(:, 4)==1)|(t(:, 4)==2), :);
    tD = t(t(:, 4)==3,  :);
    cM = c((t(:, 4)==1)|(t(:, 4)==2), :);
    cD = c(t(:, 4)==3,  :);
   
    Ptemp = P'; ttemp = t'; ind = size(t, 1);  IND = 3000;
    
    ctempM = cM'; ttempM = tM'; 
    XM = reshape(Ptemp(1, ttempM(1:3, :)),[3, size(ttempM, 2)]);
    YM = reshape(Ptemp(2, ttempM(1:3, :)),[3, size(ttempM, 2)]);
    ZM = reshape(Ptemp(3, ttempM(1:3, :)),[3, size(ttempM, 2)]);  
    
    ctempD = cD'; ttempD = tD'; 
    XD = reshape(Ptemp(1, ttempD(1:3, :)),[3, size(ttempD, 2)]);
    YD = reshape(Ptemp(2, ttempD(1:3, :)),[3, size(ttempD, 2)]);
    ZD = reshape(Ptemp(3, ttempD(1:3, :)),[3, size(ttempD, 2)]);    
        
    X = reshape(Ptemp(1, ttemp(1:3, :)),[3, size(ttemp, 2)]);
    Y = reshape(Ptemp(2, ttemp(1:3, :)),[3, size(ttemp, 2)]);
    Z = reshape(Ptemp(3, ttemp(1:3, :)),[3, size(ttemp, 2)]);      
    
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
%     stroc.yes         = 'yes';     %   'yes' - present; 'no' - absent
%     stroc.x           = 0.0;       %   x position in m
%     stroc.y           = 0.0;       %   y position in m
%     stroc.z           = 0.0;       %   z position in m
%     stroc.size = 1;                %   relative marker size versus default size 
    if strcmp(stroc.yes, 'yes')
        scaling = 0.05*max(strop.planesizex, strop.planesizey)*stroc.size;
        fv      = kreuz(stroc.x, stroc.y, stroc.z, 4, 0.1*scaling, 1*scaling);
        patch(fv, 'FaceColor', 'g', 'EdgeColor', 'k'); 
    end
    hold on
    %%   Input graphics (within the main GUI window)    
    if io == 0
        ec='none';                
        patch(XD, YD, ZD, [1 0.75 0.65], 'FaceAlpha', 1.00, 'EdgeColor', ec)     
        ec='y';
        patch('Vertices', P, 'Faces', tM(:, 1:3), 'FaceColor',  'b', 'FaceAlpha', 1.00, 'EdgeColor', ec); 
        axis('equal'); axis('tight'); camproj('orthographic');
        light('Position',[1 3 2]); light('Position',[-3 -1 3]); material shiny; grid on;
        xlabel('x, m'), ylabel('y, m'), zlabel('z, m');  
        title('Problem geometry: two conductors and the conducting body');
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
        
        if (~strcmp(strop.yes, 'yes'))&(~strcmp(strop.potential, 'yes'))
            ec='none';
            patch(XD, YD, ZD, CD, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');
            patch(XM, YM, ZM, CM, 'FaceAlpha', 1.0, 'EdgeColor', ec,  'FaceLighting', 'none');
            xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); colorbar;
            title('Solution: Surface charge density on all objects in C/m^2')
        end        
        if strcmp(strop.yes, 'yes')|strcmp(strop.potential, 'yes')
            ec='none';
            patch(XD, YD, ZD, CD, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');
            patch(XM, YM, ZM, CM, 'FaceAlpha', 1.0, 'EdgeColor', ec, 'FaceLighting', 'none');
            xlabel('x, m'), ylabel('y, m'), zlabel('z, m'); colorbar;
            title('Solution: Surface charge density in C/m^2 and the field in the observation plane'); 
        end
        if strcmp(strop.yes, 'yes')
            %%   Draw electric field (final) in the observation plane
            strop.E = strop.E + 1e-12*max(max(strop.E))*rand(size(strop.E)); %   randomize a bit
            normalize = sqrt(sum(strop.E.*strop.E, 2));
            normalize = max(normalize);
            normalize = normalize/strop.arrow;
            fv = cone(strop.Points(:, 1), strop.Points(:, 2), strop.Points(:, 3), strop.E(:, 1), strop.E(:, 2), strop.E(:, 3), 36, 0.25, 1, normalize);
            h = patch(fv, 'FaceColor', 'y', 'EdgeColor', 'none'); 
            camproj('orthographic');
            light('Position',[1 3 2]); light('Position',[-3 -1 3]); material shiny;   
        end
        if strcmp(strop.potential, 'yes')
            contourca(strop, strop.Potential);
        end
    end
    
    %%   General settings 
    axis('equal'); axis('tight'); grid on 
    view(-132, 19);
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

    d = drectangle(P, L, W);
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
    dist = drectangle(P, L, W);                         %  A column of signed distances
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
            dist    = drectangle(ic, L, W);             %   A column of signed distances
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
        dist  = drectangle(P, L, W);                    %   Find distances
        ind   = dist>0;                                 %   Find vertices outside (d>0)
        dgradx = (drectangle([P(ind, 1)+deps, P(ind, 2)], L, W)-dist(ind))/deps;   %   Numerical
        dgrady = (drectangle([P(ind, 1), P(ind, 2)+deps], L, W)-dist(ind))/deps;   %   gradient
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
    if msg ==0
        h    = waitbar(0, 'Please wait - computing electric potential at the obs. point');
    end
    if msg ==1
        h    = waitbar(0,  'Please wait - computing electric potentials in a plane');
    end
    if msg ==2
        h    = waitbar(0, 'Please wait - computing average electric potential of object #1');
    end
    if msg ==3
        h    = waitbar(0, 'Please wait - computing average electric potential of object #3');
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


