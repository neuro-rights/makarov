function  dy = fun_derivative(x, y)
%   SYNTAX
%   dy = fun_derivative(x, y)
%   DESCRIPTION
%   Function derivative
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

dx = x(2) - x(1);
dy = diff(y)/dx;   
dy = [2*dy(1)-dy(2) dy]; 

