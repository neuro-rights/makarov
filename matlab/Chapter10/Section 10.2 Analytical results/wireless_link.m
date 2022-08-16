clear all
%   SYNTAX
%   wireless_link
%   DESCRIPTION
%   This script  will estimate voltage signal (voltage amplitude) induced
%   in the second coaxial coil (RX) if the periodic current (current
%   amplitude) in the first coaxial coil (TX) is known (Example 10.6 of
%   Chapter 10)
%    
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

mu0 = 4*pi*1e-7;    % permeability of vacuum (air)
omega = 2*pi*1e6;   % angular frequency, rad/s
i1 = 0.1;           % amplitude of exciting current i1, A
r = 1e-2;           % coil radius, m
l = 0.1;            % coil length, m
N = 100;            % number of turns
d = [0.1:0.01:2];   % separation distance, m
M0 = pi*mu0*r^4*N^2./(2*d.^3);      % mutual inductance (no mag. core)
v2 = M0*omega*i1;                   % received voltage, V
semilogy(d, v2*1000); grid on;
xlabel('distance d, m'); ylabel('received voltage, mV')