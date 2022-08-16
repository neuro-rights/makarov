clear all
%   SYNTAX
%   mesh_cylinder_polygon2
%   DESCRIPTION
%   This scipt creates and saves a cylindrical arbitrary-base mesh obtained by extrusion
%   including outer normal vectors
%
%   To display the mesh use: fv.faces = t; fv.vertices = P;
%   patch(fv, 'FaceColor', 'y'); axis equal; view(160, 60); grid on;
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.


load polygon2.mat;      %   base
H = 0.030;              %   structure height in m
layers = 3;             %   number of layers along the height

%   Determine boundary edges of the cap mesh  
[NonManifoldAttached, edges, edgesb, edgesnm] = meshedges(P, t);

%   Start with the bottom
P0 = P;         
t0 = t;
P(:, 3) = -H/2;
TriCap  = size(P, 1);

%   Do extrusion but keep inner nodes
for n = 1:layers
    ttemp1 = [edgesb(:, 1)+TriCap*(n-1) edgesb(:, 2)+TriCap*(n-1) edgesb(:, 2)+TriCap*(n-0)];
    ttemp2 = [edgesb(:, 1)+TriCap*(n-1) edgesb(:, 1)+TriCap*(n-0) edgesb(:, 2)+TriCap*(n-0)];
    t = [t; ttemp1; ttemp2];
    Ptemp = P0; Ptemp(:, 3) = -H/2 + n*H/layers;
    P  = [P; Ptemp];
end
t = [t; t0+size(P, 1)-size(P0, 1)];

%   Now remove inner nodes
[P, t] = fixmesh(P, t);

P       = meshrotate(P, [1 0 0], pi/2);
snumber = [];

fv.faces = t(:, 1:3); fv.vertices = P; patch(fv, 'FaceColor', [0.8 0.9 1.0]); 
axis equal; axis tight
view(166, 26); grid on;
xlabel('x, m'); ylabel('y, m'); zlabel('z, m');

normals = meshnormals(P, t);
save('mesh_cylinder_polygon2.mat', 'P', 't', 'normals', 'snumber');
