function [charge, a] = mom2d(xs, ys, xe, ye, phi);
%MOM2D Compute charge distribution for 2D electrostatics by MoM
% for an air-filled transmission line
%
% Syntax:
%   [charge, a] = mom2d(xs, ys, xe, ye, phi);
%
% Description:
% INPUTS:
%   xs      = x-coordinates for edge starting points (row)
%   ys      = y-coordinates for edge starting points (row)
%   xe      = x-coordinates for edge ending points (row)
%   ye      = y-coordinates for edge ending points (row)
%   phi     = potential for given edges (row)
% OUTPUTS:
%   a       = MoM coefficients a(n)
%   charge  = total line charge on each element (C/m)
%   SNM Summer 2014

%   Vacuum parameters
const.epsilon       = 8.85418782e-012;                  
const.mu            = 1.25663706e-006;                  
const.c             = 1/sqrt(const.epsilon*const.mu);
const.eta           = sqrt(const.mu/const.epsilon); 
 
x = 0.5*(xs + xe);                  % Observation/integration points
y = 0.5*(ys + ye);                  % Observation/integration points
h = sqrt((xe-xs).^2 + (ye-ys).^2);  % Edge length 
 
% Impedance matrix filling
N = length(x);                      % Number of edges
Z = zeros(N, N);                    % Impedance matrix   
V = zeros(N, 1);                    % Right-hand side
Distance = zeros(1, N);             % Distance vector
Kernel   = zeros(1, N);             % Integral kernel
eps      = 1e-6*min(h);             % Dummy variable
 
for m = 1:N
    Distance    = sqrt((x-x(m)).^2 + (y-y(m)).^2) + eps;
    Kernel      = log(Distance);
    Z(m,:)      = -1/(2*pi*const.epsilon)*Kernel.*h*h(m);
    V(m)        = phi(m)*h(m);
end
% MoM matrix (divide every row of MoM matrix and V(m) by h(m))
% Non-diagonal elements- single-point approximation
for m = 1:N
    Distance    = sqrt((x-x(m)).^2 + (y-y(m)).^2) + eps;
    Kernel      = log(Distance);
    Z(m,:)      = -1/(2*pi*const.epsilon)*Kernel.*h;
    V(m)        = phi(m);
end
% Diagonal (singular elements)
for m = 1:N
    Distance    = sqrt((x-x(m)).^2 + (y-y(m)).^2) + eps;
    Kernel      = log(Distance);
    Z(m,m)      = -1/(2*pi*const.epsilon)...
                  *h(m)*(log(h(m))-1.5);
end
% Charge conservation law:
% Subtract the last metal row from the others and replace the 
% last metal edge equation by the conservation law
for m = 1:N-1
    Z(m,:) = Z(m,:) - Z(N,:);
    V(m)   = V(m)   - V(N);
end
% Last metal row
Z(N,:) = ones(1, N).*h;
V(N)   = 0;
% Solution of linear equations
a       = (Z\V)';      % MoM coefficients a(n) from matrix equation
charge  = h.*a;        % Total line charge per edge (C/m)