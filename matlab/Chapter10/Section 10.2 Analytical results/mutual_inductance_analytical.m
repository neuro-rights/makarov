clear all;
%   SYNTAX
%   mutual_inductance_analytical
%   DESCRIPTION
%   This script finds mutual inductance, M, given by Eq. (10.42) for two
%   equal air-core coaxial coils
%    
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

mu0     = 4*pi*1e-7;        %   permeability of vacuum (air)
omega   = 2*pi*1e6;         %   angular frequency, rad/s
r       = 0.02;             %   coil radius, m
N       = 7;                %   number of turns
d       = 0.4;              %   separation distance, m
M0      = pi*mu0*r^4*N^2./...
        (2*d.^3)            %   mutual inductance (no mag. core)
