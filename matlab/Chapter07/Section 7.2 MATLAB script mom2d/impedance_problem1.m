%   SYNTAX
%   impedance_problem1
%   DESCRIPTION
%   This script computes static capacitance and characteristic impedance of
%   a parallel-plate air-filled TL via the MoM method
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

%   Vacuum parameters
const.epsilon       = 8.85418782e-012;                  
const.mu            = 1.25663706e-006;                  
const.c             = 1/sqrt(const.epsilon*const.mu);
const.eta           = sqrt(const.mu/const.epsilon); 
 
d       = 0.01;             % Separation distance (m)
w       = 0.01;             % Plate width (m)
h       = w/20;             % Edge length (variable parameter)
eps_r   = 1;                % Dielectric constant (in the entire space)
 
% X- and Y-coordinates for starting and ending points
% Bottom plate
xs_b  = [-w/2:h:w/2-h]; 
xe_b  = xs_b + h;
ys_b  = -d/2*ones(size(xs_b)); 
ye_b  = ys_b;
% Add top plate
xs_t = [-w/2:h:w/2-h]; 
xe_t = xs_t + h;
ys_t = d/2*ones(size(xs_t)); 
ye_t = ys_t;
% Combine two meshes together
xs = [xs_b  xs_t]; xe = [xe_b  xe_t]; 
ys = [ys_b  ys_t]; ye = [ye_b  ye_t]; 
% Mesh size
N = length(xs)
% Potential for two conductors (V)
V(1:N/2)    = 0.0;    % ground plane
V(N/2+1:N)  = 1.0;    % top plate
% Solve the electrostatic problem (C, Z0)
[charge, a] = mom2d(xs, ys, xe, ye, V);
C0          = sum(charge(N/2+1:N))          % In F/m
Z0          = sqrt(eps_r)/(const.c*C0)      % In Ohm
Error       = abs(sum(charge))/max(abs(charge))
C_exact     = 18.73e-12*sqrt(eps_r)
Z0_exact    = sqrt(eps_r)/(const.c*C_exact) 
