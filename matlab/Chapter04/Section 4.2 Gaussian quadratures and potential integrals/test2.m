clear all
%   SYNTAX 
%   test2
%   DESCRIPTION 
%   This script illustrates the use of barycentric quadratures on triangles
%   to find a single potential MoM integral I = Is(1/r) 
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

%   Apply barycentric quadrature
[coeff, weights, IndexF] = tri(100); % number of subtriangles is 100*100=10000
%   Generate a triangle (generate its vertices)
p1 = [0 0 0];
p2 = [0 1 0]; 
p3 = [1 0 0];
ObsPoint  = [0.5 0.5 0];
Area = 0.5;

%   Find single potential MoM integral I = Is(1/r) using the barycentric
%   quadrature
I = 0;
for p =1:length(weights)
    IntPoint = coeff(1, p)*p1 +  coeff(2, p)*p2 +  coeff(3, p)*p3;    
    I = I + weights(p)*1/sqrt(dot(ObsPoint-IntPoint, ObsPoint-IntPoint));
end
I = I*Area