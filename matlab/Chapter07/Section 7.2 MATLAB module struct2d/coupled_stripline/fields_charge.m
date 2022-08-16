%   SYNTAX
%   fields_charge
%   DESCRIPTION
%   This script plots the charge distribution in a plane
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

clear all
factor = 32;                                            %  scale factor for charge scaling

const.epsilon       = 8.85418782e-012;                  %  ANSOFT HFSS value 
const.mu            = 1.25663706e-006;                  %  ANSOFT HFSS value
const.c             = 1/sqrt(const.epsilon*const.mu);
const.eta           = sqrt(const.mu/const.epsilon); 

% Load data
val = 1;
matfilename = ['struct2d_var'  num2str(val)];
matfilename = [matfilename    '.mat'];
if exist(matfilename)
    load(matfilename);
else
    error('geometry structure does not exist');
end
dimensions;
matfilename = ['solution_var'  num2str(val)];
matfilename = [matfilename    '.mat'];
if exist(matfilename)
    load(matfilename);
else
    error('solution structure does not exist');
end

% Scale and plot surface charge distribution
a = figure;
colorbar;

ch_density  = output.charge./mesh.eh;           % charge density per unit length
variance    = sqrt(dot(ch_density, ch_density))...
             /length(ch_density);
scale       = factor*variance;
ch_density(ch_density<(-scale)) = -scale;
ch_density(ch_density>(+scale)) = +scale;

charge      =(ch_density + scale)/(2*scale);    % normalized to 1
N      = 100;
C      = colormap(jet(N));                      % colormap
charge = min(round(N*charge) + 1, N);           % normalized to colormap

hold on; grid on; axis equal;
for m = 1:length(mesh.e)
    line( [mesh.P(1, mesh.e(1:2, m))], [mesh.P(2, mesh.e(1:2, m))], ...
          'LineWidth', 5, 'Color', [C(charge(m),1:3)] )
end
axis('equal'); xlabel('x, m'); ylabel('y, m'); zlabel('z, m');
title(''); 
title({'Contour color - total charge density;';...
       'blue - most negative charge density; red - most positive charge density'; ''})
