clear all
%   SYNTAX
%   H_yoke_analytical
%   DESCRIPTION
%   Analytical solution (H field) in the gap of a magnetic yoke given by
%   Eq. (10.58) (Section 10.3 of Chapter 10)
%    
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

mu0         = 1.25663706e-006;       %   magnetic permeability of vacuum(~air)
h           = 0.080;        %   half of height, m
t           = 0.030;        %   bar width, m
l           = 0.10;         %   arm length, m  
g           = 0.03;         %   half of gap width
delta       = 0.01;         %   smooth edges  

L = 4*h + 4*t + 2*l -2*g;

N       = 13;
i       = 1;  % in A
mur     = 1000; 
A = t*t;
Aeff = (sqrt(A)+2*g)^2;

B = N*i/(L/(mu0*mur)+2*g/mu0);
B = B*A/Aeff


