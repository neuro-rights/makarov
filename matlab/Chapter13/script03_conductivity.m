clear all;
%   SYNTAX
%   script03_conductivity
%   DESCRIPTION
%   Conductivity in 1/(ohm cm) for Si pn-junction as a function of doping
%   concentration
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

q       = 1.60218e-19;          %   electron charge [C]
mun     = 1450;                 %   electron mobility, cm^2/(Vs)
mup     =  500;                 %   hole mobility, cm^2/(Vs)
ND      =  1e17;
sigma = q*ND*mun
