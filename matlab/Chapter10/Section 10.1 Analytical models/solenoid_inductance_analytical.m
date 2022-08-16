clear all
%   SYNTAX
%   solenoid_inductance_analytical
%   DESCRIPTION
%   Analytical solution for inductance of an air-core solenoid, L, given by
%   Eq. (10.10)
%    
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

mu0         = 1.25663706e-006;       %   magnetic permeability of vacuum(~air)
strcoil.length = 0.10;               %   coil length, m
strcoil.radius = 0.02;               %   coil radius, m
strcoil.turns  = 17;                 %   number of turns
strcoil.I       = 1;                 %   coil current, A
A = pi*strcoil.radius^2;
w = strcoil.radius/strcoil.length;

Hcenter = strcoil.turns*strcoil.I/strcoil.length;
L = mu0*A*strcoil.turns^2/strcoil.length*(1-8*w/(3*pi)+w^2/2-w^4/4) % H


