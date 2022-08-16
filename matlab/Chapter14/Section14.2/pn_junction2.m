clear all
%   SYNTAX
%   pn_junction2
%   DESCRIPTION
%   This script models an electric field region of the silicon pn junction with a variety of doping
%   profiles for bias voltages -5.0V < V < +0.7V
%
%   The script calculates and plots 
%   ND(x), NA(x), NE(x);        [a.u.]      %   doping profiles
%   p(x), n(x), ps(x), ns(x);   [a.u.]      %   carrier concentrations
%   phi(x), phis(x);            [V]         %   electric potentials
%   rho(x)                      [C/cm^3]    %   charge density
%   E(x)                        [V/m]       %   electric field

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
V       = -5                      %   bias voltage [V]  - try -1, 0, 0.2, 0.4, and 0.6                                  
W       = ...
        sqrt((2*eps/q)*(Vbi-V)*(ND0 + NA0)/(ND0*NA0))
                                    %   width of the depletion region 
                                    %   for an equivalent abrupt junction

%   Space-charge neutrality model (SCNM) as an initial guess
[ns, ps, phis, Vbi]  = fun_scnm2(ni, NE, VT, V, x);
%   Iterative solution to the Poisson's equation
[n, p, phi, rho,  E] = fun_ccflm2(ni, q, eps, NE, VT, V, x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Graphics
%   Graphics
h1 = subplot (321); plot(x, ND, 'b', x, NA, 'r', x, NE, 'm'); 
grid on; subplot(h1); 
axis tight;
title('Doping profiles, 1/cm^{3}');
xlabel('distance, cm');

h2 = subplot (322); plot(x, n, 'b', x, p, 'r'); 
hold on; plot(x, ns, '--k', x, ps, '--k'); 
grid on; subplot(h2); scale = max(max(n), max(p)); 
axis tight;
title('n(x) and p(x), n_{s}(x) and p_{s}(x)');
xlabel('distance, cm');

h3 = subplot (323); plot(x, phi, 'k', x, phis, '--k');
grid on; subplot(h3);
axis([min(x) max(x) min(min(phi), min(phis)) max(max(phi), max(phis))]);
title('\phi(x) and \phi_{s}(x), V');
xlabel('distance, cm');

h4 = subplot (324); plot(x, rho, 'b', x, mean(rho)*ones(size(rho)), '--k'); 
grid on; subplot(h4);
axis([min(x) max(x) -max(abs(rho)) +max(abs(rho))])
title('\rho(x), C/cm^{3}');
xlabel('distance, cm');

h5 = subplot (325); plot(x, E, 'b'); 
grid on; subplot(h5);
axis([min(x) max(x) -max(abs(E)) +max(abs(E))])
title('Electric field, V/cm');
xlabel('distance, cm');
