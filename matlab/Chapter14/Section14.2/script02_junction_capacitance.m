clear all
%   SYNTAX
%   script02_junction_capacitance
%   DESCRIPTION
%   Junction capacitance/stored charge for a biased abrupt Si pn-junction
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

k       = 1.38066e-23;          %   Boltzmann constant [J/K]
q       = 1.60218e-19;          %   electron charge [C]
T       = 273+25;               %   temperature [K]
VT      = k*T/q;                %   thermal voltage at 300 K [V]
ni      = 1e10;
V       = [-5:0.1:0.2];
ND0             = 1.0e18;     
NA0             = 1.0e18;
Vbi             = VT*log(NA0*ND0/ni^2)
eps     = 1.05e-12;             %   dielectric constant of silicon [F/cm]  
A       = 1e-4;
W       = sqrt(2*eps/q*(ND0+NA0)/(ND0*NA0)*(Vbi-V));
C       = eps*A./W*1e12;        %   junction capacitance, pF
Q       = sqrt(2*eps*q*(ND0*NA0)/(ND0+NA0)*(Vbi-V));
plot(V, C); grid on;
xlabel('Bias voltage, V'); ylabel('Capacitance, pF')

