%   SYNTAX
%   viewer_plain
%   DESCRIPTION
%   This script displays a mesh from a *.mat P-t file
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

FileName = uigetfile('*.mat','Select the tissue mesh file to open');
load(FileName, '-mat');

fv.faces = t(:, 1:3); fv.vertices = P; patch(fv, 'FaceColor', [1 0.75 0.65]); 
axis equal; axis tight
view(130, 32); grid on; xlabel('x, mm'); ylabel('y, mm'); zlabel('z, mm');