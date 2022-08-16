clear all
%   SYNTAX
%   loop_inductance_analytical
%   DESCRIPTION
%   Analytical solution for inductance of a single loop, L, given by Eq.
%   (10.11)
%    
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

mu0         = 1.25663706e-006;          %   magnetic permeability of vacuum(~air)
r           = 0.02;                     %   loop radius, m
a           = 0.002;                    %   wire radius, m
Y           = 1/2;
L = mu0*r*(log(8*r/a) -2 +Y/2)          %   H


