clear all 
%   SYNTAX
%   Two_current_electrodes
%   DESCRIPTION
%   This script computes voltage between the centers of two current
%   electrodes placed on top of a conducting layer backed by a
%   semi-infinite substrate - programs Eq. (8.48) of Chapter 8
%
%   Author: S.N. Makarov
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

%   Problem parameters
R           = 0.002;           %   Electrode radius, m
d           = 0.010;           %   Distance between electrodes, m
h           = 0.002;           %   Thickness of cond. layer, m
I0          = 0.02;            %   Electrode current, A
sigma1      = 0.1;             %   Conductivity of the layer, S/m
sigma2      = 0.5;             %   Conductivity of substrate, S/m
step1       = 0.05/R;          %   Step size of alpha
alpha       = 0+R:step1:50/R;
sigma       = (sigma1-sigma2)/(sigma1+sigma2);    

%   Analytical solution
SIGMA       = sigma*exp(-2*alpha*h);
J0          = I0/(pi*R^2);
Integral    = step1*sum((1-besselj(0,alpha*d)).*...
        ((besselj(1,alpha*R)).*R./(alpha)).*...
        ((1+sigma*exp(-2*alpha*h))./(1-sigma*exp(-2*alpha*h))));
    
%   Resulting voltage
V = 2*J0/sigma1*Integral