function  [n, p, phi, rho, E] = fun_ccflm1(ni, q, eps, NE, VT, x)
%   SYNTAX
%   [n, p, phi, rho, E] = fun_ccflm1(ni, q, eps, NE, VT, x)
%   DESCRIPTION
%   This function implements an iterative solution of Poisson's equation
%   at zero bias voltage with Jacobi iterations
%   The function uses the approach of  B. R. Chawla and H. K. Gummel, "Transition 
%   region capacitance of diffused p-n junctions," IEEE Trans. on Electron Devices, 
%   vol. ED18, no.3, March 1971.
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

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
for m = 1:1e6        
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
end
% error
phi   = phi_next; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Output parameters
n   = ni.*exp(+phi/VT);               %   n(x)
p   = ni.*exp(-phi/VT);               %   p(x)
phi = phi;                            %   total electric potential
rho = q*(p - n + NE);                 %   rho(x)
E   = - fun_derivative(x, phi);       %   electric field as a negative grad of phi 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

