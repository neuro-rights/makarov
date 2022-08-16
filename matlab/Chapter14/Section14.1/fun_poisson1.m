function  rhs = fun_poisson1(x, phi, ni, q, eps, VT, NE)
%   SYNTAX
%   rhs = fun_poisson1(x, phi, ni, q, eps, VT, NE)
%   DESCRIPTION
%   This function computes the right-hand side of the nonlinear Poisson
%   equation, rhs = -(q/eps)*(p - n + NE) at zero bias voltage
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

temp    = exp(+phi/VT);
rhs     = -(q/eps)*ni*(1./temp - temp + NE/ni);

