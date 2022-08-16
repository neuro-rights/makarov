clear all
%   SYNTAX
%   [ns, ps, phis, Vbi] = fun_scnm1(ni, NE, VT, x)
%   DESCRIPTION
%   This script models a silicon pn junction with a variety of doping
%   profiles at zero bias voltage - demonstrates solution convergence.
%   Press ENTER to start the iterative process.  
%
%   The script calculates and plots 
%   ND(x), NA(x), NE(x);        [a.u.]      %   doping profiles
%   p(x), n(x), ps(x), ns(x);   [a.u.]      %   carrier concentrations
%   phi(x), phis(x);            [V]         %   electric potentials
%   rho(x)                      [C/cm^3]    %   charge density
%   E(x)                        [V/m]       %   electric field
%   for a zero bias voltage V 

%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

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
NA0             = 1.0e16     
NA              = NA0*(1 - exp(+x/L));
NA(find(x>0))   = 0;
NE              = ND - NA;                  %   effective doping concentration
Vbi             = VT*log(NA0*ND0/ni^2)      %   built-in voltage
W       = ...
        sqrt((2*eps/q)*Vbi*(ND0 + NA0)/(ND0*NA0))
                                    %   width of the depletion region 
                                    %   for an equivalent abrupt junction
%   Space-charge neutrality model (SCNM) as an initial guess
[ns, ps, phis, Vbi]  = fun_scnm1(ni, NE, VT, x); 
%   Iterative solution to the Poisson's equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Space-charge neutrality model: ps(x), ns(x), and phis(x) 
[ns, ps, phis, Vbi] = fun_scnm1(ni, NE, VT, x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Iterative solution to the Poisson's equation for the potential 
dx          = x(2) - x(1);          %   discretization step 
phi_past    = phis;                 %   starting guess
J           = 1e-2                  %   Jacobi iteration factor (lower if necessary)
temp        = 1/(max(x)-min(x));
temp1       = (x-min(x));
for m = 1:100        
    charge      = fun_poisson1(x, phi_past, ni, q, eps, VT, NE);    %   right-hand side
    inner_int   = dx*cumsum(charge);                                %   inner integral   
    outer_int   = dx*cumsum(inner_int);                             %   outer integral                  
    C1          = phis(1); 
    C2          = temp*(-C1-outer_int(end)+phis(end));
    phi_next    = outer_int + C1 + C2*temp1;
    error       = norm(phi_next-phi_past)/norm(phi_next);                 
    phi_past      = J*phi_next + (1-J)*phi_past;    
    if (error <0.01*J) 
        iterations = m
        break; 
    end;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cla;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Error
    phi   = phi_next;                     %   total electric potential
    %   Output parameters
    n   = ni.*exp(+phi/VT);               %   n(x)
    p   = ni.*exp(-phi/VT);               %   p(x)
    rho = q*(p - n + NE);                 %   rho(x)
    E   = - fun_derivative(x, phi);       %   electric field as a negative grad of phi 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Graphics
    h1 = subplot (321); plot(x, ND, 'b', x, NA, 'r', x, NE, 'm'); 
    grid on; subplot(h1); 
    axis tight;
    title('Doping profiles, 1/cm^{3}');
    xlabel('distance, cm');

    h2 = subplot (322); plot(x, n, 'b', x, p, 'r', x, ns, '--k', x, ps, '--k');     
    grid on; subplot(h2); scale = max(max(n), max(p)); 
    axis tight;
    title('n(x) and p(x), n_{s}(x) and p_{s}(x)');
    xlabel('distance, cm');

    h3 = subplot (323); plot(x, phi, 'k', x, phis, '--k');
    grid on; subplot(h3);
    axis([min(x) max(x) -0.5 +0.5])
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
    drawnow; pause(0.5/m);
    if m==1 pause; end
        
end

