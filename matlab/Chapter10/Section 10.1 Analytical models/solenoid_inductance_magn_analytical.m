clear all
%   SYNTAX
%   solenoid_inductance_magn_analytical
%   DESCRIPTION
%   Analytical solution for inductance of a solenoid, L, with a magnetic
%   core given by Eq. (10.14). Core length is equal to the coil length. 
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

L = 0.25*pi*mu0*strcoil.length*strcoil.turns^2/log(strcoil.length/strcoil.radius-1)


