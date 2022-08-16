clear all
%   SYNTAX
%   mesh_polygon
%   DESCRIPTION
%   This scipt creates and saves a polygonal planar mesh obtained with 
%   distmesh2d
%
%   To display the mesh use: fv.faces = t; fv.vertices = P;
%   patch(fv, 'FaceColor', 'y'); axis equal; view(160, 60); grid on;
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

%   Example: a yoke (start with the upper left corner)
h           = 0.080;        %   half of height, m
t           = 0.030;        %   bar width, m
l           = 0.10;         %   arm length, m  

edgelength  = t/3;        %   desired edge length

%   polygonal shape
pv=         [-t/2   h+t;...
             +t/2+l h+t;...
             +t/2+l +h;...
             +t/2   +h;...
             +t/2   -h;...
             +t/2+l -h;...
             +t/2+l -h-t;...             
             -t/2   -h-t;...
             -t/2   +h+t];       
         
[p, t]=distmesh2d(@dpoly, @huniform, edgelength, [-1,-1; 1,1], pv, pv); 

P = p; P(:, 3) = 0;

save('polygon2.mat', 'P', 't');