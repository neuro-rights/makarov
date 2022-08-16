clear all
%   SYNTAX
%   pn_junction3
%   DESCRIPTION
%   This script models a pn-junction capacitance of a silicon pn junction 
%   with a variety of doping profiles for bias voltages -5.0V < V < +0.7V
%
%   This script models a pn-junction capacitance of a silicon pn junction 
%   with a variety of doping profiles  for bias voltages -5.0V < V < +0.8V.
%   SNM, ECE Dept., WPI Rev. Oct. 2008

%   Material parameters
k       = 1.38066e-23;          %   Boltzmann constant [J/K]
q       = 1.60218e-19;          %   electron charge [C]
T       = 300;                  %   temperature [K]
VT      = k*T/q;                %   thermal voltage at 300 K [V]
eps     = 1.05e-12;             %   dielectric constant of silicon [F/cm]  
ni      = 9.65e09;              %   intrinsic concentration of holes/electrons [a.u.]
L       = 1.0e-5;               %   characteristic length of exponential acceptor's 
                                %   profile [cm]
%   Geometry
x       = [-10*L:L/100:10*L];       %   computational domain [cm]
dx      = x(2) - x(1);              %   discretization step 

%   Doping profiles
ND0             = 1.0e16;     
ND              = ND0*(1 - exp(-x/L));
ND(find(x<0))   = 0;
NA0             = 1.0e16;     
NA              = NA0*(1 - exp(+x/L));
NA(find(x>0))   = 0;
NE              = ND - NA;                  %   effective doping concentration
Vbi             = VT*log(NA0*ND0/ni^2)      %   built-in voltage
                              
%   Applied bias voltage
V       = 0.6;                     %   bias voltage [V]      
dV      = 0.01;                   %   voltage variation, V
W       = ...
        sqrt((2*eps/q)*(Vbi-V)*(ND0 + NA0)/(ND0*NA0))
                                    %   width of the depletion region 
                                    %   for an equivalent abrupt junction

%   Space-charge neutrality model (SCNM) as an initial guess
[ns, ps, phis, Vbi]  = fun_scnm2(ni, NE, VT, V, x); 
%   Iterative solution to the Poisson's equation
[n, p, phi, rho,  E] = fun_ccflm2(ni, q, eps, NE, VT, V, x);
%Q1 = sum(rho(find(rho>0)))*(x(2) - x(1));
Q1 = q*sum(p)*(x(2) - x(1));
V = V + dV;
%   Space-charge neutrality model (SCNM) as an initial guess
[ns, ps, phis, Vbi]  = fun_scnm2(ni, NE, VT, V, x); 
%   Iterative solution to the Poisson's equation
[n, p, phi, rho,  E] = fun_ccflm2(ni, q, eps, NE, VT, V, x);
%Q2 = sum(rho(find(rho>0)))*(x(2) - x(1));
Q2 = q*sum(p)*(x(2) - x(1));

Cexact = (Q2 - Q1)/dV        %   in F per cm^2 
Cdepl  = eps/W               %   in F per cm^2  
