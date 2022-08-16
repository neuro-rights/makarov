clear all
%   SYNTAX
%   script02_built_in_voltage_depl_width
%   DESCRIPTION
%   Built-in voltage (V) for Si pn-junction with terminal concentrations ND0,
%   NA0 and the deplection region width (cm) for the abrupt junction
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

k       = 1.38066e-23;          %   Boltzmann constant [J/K]
q       = 1.60218e-19;          %   electron charge [C]
T       = 273+100;                  %   temperature [K]
VT      = k*T/q;                %   thermal voltage at 300 K [V]
ni      = 1e10;
ND0             = 1.0e15;     
NA0             = 1.0e15;
Vbi             = VT*log(NA0*ND0/ni^2)

eps     = 1.05e-12;             %   dielectric constant of silicon [F/cm]  
W = sqrt(2*eps/q*(ND0+NA0)/(ND0*NA0)*Vbi)