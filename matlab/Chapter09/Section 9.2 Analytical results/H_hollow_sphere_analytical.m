clear all
%   SYNTAX
%   H_hollow_sphere_analytical
%   DESCRIPTION
%   Analytical solution (H field) for a hollow magnetic sphere given by Eq.
%   (9.35)
%    
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

R = 0.5;
R1 = 0.45;
mu1 = 5000; 
mu2 = 1;
Factor = (R1/R)^3;
FactorMu= (mu1+2*mu2)*(2*mu1+mu2)/(2*(mu1-mu2)^2)
H = 1 - (1-Factor)/(FactorMu - Factor)