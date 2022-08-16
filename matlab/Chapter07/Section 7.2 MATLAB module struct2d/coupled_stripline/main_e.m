%   SYNTAX
%   main_e
%   DESCRIPTION
%   This script compares the analytical solution for the coupled stripline
%   with the corresponding MoM result (even mode)
%
%   Outputs:
%   Z0 - characteristic line impedance
%   C - characteristic capacitance per unit length
%   L - characteristic inductance per unit length
%   eps_eff - effctive dielectric constant 
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

% Analytical (exact) solution 
%   S. B. Cohn, “Shielded coupled strip transmission line,” IRE Trans., 
%   vol. MTT-5, Oct. 1955, pp. 29-37.
w       = 1.00e-3;             % Strip width (m)
d       = 1.00e-3;             % Half separation distance between ground plates(m)
S       = 0.50e-3;             % Separation between conductors (m)
eps_r   = 2.2;                 % Relative dielectric constant   
[Z_e, Z_o]      = imp_stripline_coupled(w, d, S, eps_r);
Cohn.Z_e = Z_e; Cohn.Z_o = Z_o; Cohn
%------------------------------------------------------------------

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
%------------------------------------------------------------------

%------------------------------------------------------------------
% Voltage vector(s) for the EVEN mode
% define (initial) voltage vector
V           = zeros(1, length(mesh.e));
% give potential of +1.0V on the 1st conductor
index1       = find(mesh.ei == 1); 
V(index1)    = +1.0;
% give potential of -1.0V on the 2nd conductor
index2       = find(mesh.ei == 2); 
V(index2)    = +1.0;
%------------------------------------------------------------------

%------------------------------------------------------------------
% Solve MoM equations and obtain TL parameters for the even mode
N = 10; % number of adaptive mesh refinement steps
[output, mesh, V]   = mom2d_md(mesh, V, N);
index1              = find(mesh.ei == 1); 
C                   = sum(output.charge(index1).*mesh.e2(index1));
% solve MoM equations without dielectric
mesh1       = mesh;
mesh1.e1( mesh1.e1>1 )  = 1; mesh1.e2( mesh1.e2>1 )  = 1;
[output1, mesh1, V]     = mom2d_md(mesh1, V, 1);
C0                      = sum(output1.charge(index1));
% find characteristic inductance and impedance
L           = const.mu*const.epsilon/C0;
eps_eff_e   = C/C0;
Z           = sqrt(eps_eff_e)/(const.c*C);
%------------------------------------------------------------------

MoM.Z_e = Z; MoM

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