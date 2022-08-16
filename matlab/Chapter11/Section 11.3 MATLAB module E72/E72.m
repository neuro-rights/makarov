clear all
%   SYNTAX
%   E72
%   DESCRIPTION
%   This module programs an analytical solution given by Eq. (11.82) of
%   Chapter 11 for eddy currents generated in  a conducting sphere by an
%   external uniform AC magnetic field H0 directed along the z-axis.
%   C. P. Bidinosti, E. M. Chapple, and M. E. Hayden, “The sphere in a
%   uniform RF field –revisited,” Concepts in Magnetic Resonance Part B
%   (Magnetic Resonance Engineering), vol. 31 B(3), pp. 191-202, 2007.
%    
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

%%   Define EM parameters
eps0        = 8.85418782e-012;                      %   Dielectric permittivity of vacuum(~air) F/m
mu0         = 1.25663706e-006;                      %   Magnetic permeability of vacuum(~air) H/m
c0          = 1/sqrt(eps0*mu0);                     %   Speed of light in vacuum(~air) m/s        
eta0        = sqrt(mu0/eps0);                       %   Vacuum/air 

%%   Sphere parameters:
H0          = 0.01/mu0;                     %   Incident magnetic field, A/m
a           = 0.5;                          %   Radius, m
sigma       = 1.0;                          %   Medium conductivity, S/m
f           = 1e6;                          %   Frequency, Hz  
skindepth   = sqrt(2/(2*pi*f*mu0*sigma))    %   Skin depth, m
parameter    = a/skindepth

%%   Define domain size and grid size
%   Will consider xz-plane
sizex = 2*a;   %   in m
sizez = 2*a;   %   in m
N = 200;
x = linspace(-sizex/2, +sizex/2, N);
z = linspace(-sizez/2, +sizez/2, N);

%%  Compute the solution
[XZ, ZX]    = meshgrid(x, z);
R           = sqrt(XZ.^2 + ZX.^2);
ARG1        = (1+j)*R/skindepth;
ARG0        = (1+j)*a/skindepth;
J1          = sin(ARG1)./ARG1.^2 - cos(ARG1)./ARG1;
J0          = sin(ARG0)./ARG0;
sinTHETA    = XZ./R; 
Jphi        = - 3/2*(1+j)*H0/skindepth*J1/J0.*sinTHETA;
Jphi(R>a)   = 0;
scaleXZ     = max(max(abs(Jphi)));
Jphi        = abs(Jphi)/scaleXZ;

%% Contour plot
v             = [1e-6:0.2:1];                     %   create contour plot levels
[C, h]        = contourf(XZ, ZX, Jphi, v);       %   create contour plot
colormap(summer);
% text_handle = clabel(C,h);
% set(text_handle,'BackgroundColor',[1 1 .6],'Edgecolor',[.7 .7 .7])
phi = linspace(0, 2*pi, N);
x   = a*cos(phi);
y   = a*sin(phi);
line(x, y, 'LineWidth', 5, 'Color', 'm');

s = 1.2;
axis([-s*a +s*a -s*a s*a]);
axis equal; colorbar;
xlabel('x, m'); ylabel('z, m');
title(strcat('Normalized eddy current density in the xz-plane. Max dens.=', num2str(1e3*scaleXZ/1e4), 'mA/cm^2'));

