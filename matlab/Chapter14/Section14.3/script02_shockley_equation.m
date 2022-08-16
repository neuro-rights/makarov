clear all
%   SYNTAX
%   script02_shockley_equation
%   DESCRIPTION
%   Shockley equation for a biased Si pn-junction (junction current)
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

k       = 1.38066e-23;          %   Boltzmann constant [J/K]
q       = 1.60218e-19;          %   electron charge [C]
T       = 298;                  %   temperature [K]
VT      = k*T/q;                %   thermal voltage at 300 K [V]
ni      = 1.0e10;               %   intrinsic concentration of holes/electrons 

%   Material parameters for the diffusion region (Si)
NA0     = 1.5e16;
ND0     = 3e16;
mun     = 1450;                 %   electron mobility cm^2/(V*sec)
mup     = 500;                  %   hole mobility cm^2/(V*sec)
taun    = 1e-6;                 %   electron lifetime, sec
taup    = 1e-8;                 %   hole lifetime, sec
Dn      = VT*mun;
Dp      = VT*mup;
Ln      = sqrt(Dn*taun);
Lp      = sqrt(Dp*taup);

%   Shockley equation
Jsn          = q*Dn*ni^2./(Ln*NA0)
Jsp          = q*Dp*ni^2./(Lp*ND0)

Js          = ( q*Dn*ni^2./(Ln*NA0)+ ...
                q*Dp*ni^2./(Lp*ND0) )
A = 1e-4;
Is = Js*A
V = 0.55;       %   Bias voltage, V
I = Is*(exp(V/VT)-1)
