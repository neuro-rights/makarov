clear all
%   SYNTAX 
%   test1
%   DESCRIPTION 
%   This script illustrates the use of Gaussian quadratures on triangles to
%   find a single potential MoM integral, I = Is(1/r) 
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

%   Apply Gaussian quadrature
[coeff, weights, IndexF] = tri(25, 10);
%   Generate a triangle (generate its vertices)
p1 = [0 0 0];
p2 = [0 1 0]; 
p3 = [1 0 0];
ObsPoint  = [0.5 0.5 0];
Area = 0.5;

%   Find single potential MoM integral I = Is(1/r) using the Gaussian
%   quadrature
I = 0;
for m = 1:length(weights)
    IntPoint = coeff(1, m)*p1 +  coeff(2, m)*p2 +  coeff(3, m)*p3;    
    r = ObsPoint-IntPoint;
    I = I + weights(m)*1/sqrt(dot(r, r));
end
I = I*Area