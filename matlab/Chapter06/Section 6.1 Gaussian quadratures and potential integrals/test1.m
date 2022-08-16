clear all
%   SYNTAX 
%   test1
%   DESCRIPTION 
%   This scrip illustrates the use of Gaussian quadratures on triangles to
%   find a single potential MoM integral I = Is(grad(1/r)) 
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

%   Apply Gaussian quadrature
[coeff, weights, IndexS] = tri(25, 10);
%   Generate a triangle (generate its vertices)
p1 = [0 0 0];
p2 = [0 1 0]; 
p3 = [1 0 0];
ObsPoint  = [3 0.25 0];
Area = 0.5;

%   Find single potential MoM integral Int = Is(grad(1/r)) using the Gaussian
%   quadrature
I = [0 0 0];
for p = 1:length(weights)     
    IntPoint = coeff(1, p)*p1 +  coeff(2, p)*p2 +  coeff(3, p)*p3; 
    R = IntPoint - ObsPoint;
    I = I + weights(p)*R/(sqrt(sum(R.*R)))^3;
end
I = Area*I
%   Find the same integral using analytical integration
I = potint2(p1, p2, p3, [0 0 -1], ObsPoint)

