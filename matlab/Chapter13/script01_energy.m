clear all
%   SYNTAX
%   script01_energy
%   DESCRIPTION
%   Photon energy in eV depeding on the wavelength
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

lambda  = 1.9e-6;
h       = 6.62606957e-34;
eV      = 1.60217657e-19 
c = 3e8;
f = c/lambda
e = h*f/eV