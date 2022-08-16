function  rhs = fun_poisson2(x, phi, ni, q, eps, VT, V, NE)
%   SYNTAX
%   rhs = fun_poisson2(x, phi, ni, q, eps, VT, V, NE)
%   DESCRIPTION
%   This function computes the right-hand side of the nonlinear Poisson
%   equation, rhs = -(q/eps)*(p - n + NE) for arbitrary bias voltages
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

temp    = exp(+phi/VT);
p       = ni*exp(V/VT)./temp;
n       = ni*temp;
rhs     = -(q/eps)*(p - n + NE);

