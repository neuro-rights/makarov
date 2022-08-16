%   SYNTAX
%   main
%   DESCRIPTION
%   This script finds the static capacitance/inductance matrices for an
%   arbitrary multiconductor transmission line
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

%******************************************************************
clear all
global STRUCT2D_BASE_DIRECTORY;
STRUCT2D_BASE_DIRECTORY     = fileparts(mfilename('fullpath'));
addpath(STRUCT2D_BASE_DIRECTORY);
addpath(fullfile(STRUCT2D_BASE_DIRECTORY, 'codes'));
%******************************************************************

const.epsilon       = 8.85418782e-012;                  %  ANSOFT HFSS value 
const.mu            = 1.25663706e-006;                  %  ANSOFT HFSS value
const.c             = 1/sqrt(const.epsilon*const.mu);
const.eta           = sqrt(const.mu/const.epsilon); 


%------------------------------------------------------------------
% MoM geometry first
% load geometry data from struct2d
val = 1;
matfilename = ['struct2d_var'  num2str(val)];
matfilename = [matfilename    '.mat'];
load(matfilename);
dimensions;
% and create (initial) edge mesh
mesh        = mesh(P, t, re_diel, re_diel_r, set_formula);
mesh0       = mesh;
%------------------------------------------------------------------
M = 4; % number of conductors (except the ground plane)
for m = 1:M    
    mesh = mesh0;
    m
    %------------------------------------------------------------------
    % Voltage vector(s) next
    % define (initial) voltage vector
    V           = zeros(1, length(mesh.e));
    % give potential of +1.0V on the 1st conductor
    index1       = find(mesh.ei == m); 
    V(index1)    = +1.0;
    %------------------------------------------------------------------

    %------------------------------------------------------------------
    % Solve MoM equations and obtain TL parameters
    N = 8; % number of adaptive mesh refinement steps    
    [output, mesh, V]   = mom2d_md(mesh, V, N);
    for n =1:M
        index1              = find(mesh.ei == n); 
        C(m, n)             = sum(output.charge(index1).*mesh.e2(index1));
                              % Maxwell's capacitance matrix 
    end
    % solve MoM equations without dielectric
    mesh1       = mesh;
    mesh1.e1( mesh1.e1>1 )  = 1;
    mesh1.e2( mesh1.e2>1 )  = 1;
    [output1, mesh1, V]     = mom2d_md(mesh1, V, 1);
    for n =1:M
        index1              = find(mesh.ei == n); 
        C0(m, n)             = sum(output1.charge(index1));
    end
end

% Circuit capacitances
C1 = - C;
for m = 1:M
    C1 (m, m) = + sum(C(m, :));
end

C01 = - C0;
for m = 1:M
    C01 (m, m) = + sum(C0(m, :));
end
C
C1
C0
C01

%------------------------------------------------------------------
% Save output data (solution and mesh)
matfilename = ['solution_var'  num2str(val)];
matfilename = [matfilename    '.mat'];
save(matfilename, 'output', 'mesh'); 

%******************************************************************
rmpath(STRUCT2D_BASE_DIRECTORY);
rmpath(fullfile(STRUCT2D_BASE_DIRECTORY, 'codes'));
clear MAT05_STRUCT2D_BASE_DIRECTORY;
%******************************************************************