clear all
%   SYNTAX
%   mesh_hollow_sphere
%   DESCRIPTION
%   This script creates a mesh (P, t, normals, snumber) for a hollow
%   sphere. Variable snumber is the surface number or domain number for
%   multidomain surfaces (triangle indicator). Variable snumber may have
%   values of 1, 2, 3; length(snumber) = size(t, 1). For shells, the inner
%   surface must have the value zero.
%
%   To display the mesh use: fv.faces = t; fv.vertices = P;
%   patch(fv, 'FaceColor', 'y'); axis equal; view(160, 60); grid on;
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

R1      = 0.050;  % Radius of the outer shell, m
R2      = 0.045;  % Radius of the inner shell, m

Tr1     = 400;  % Approximate number of triangles for the outer sphere (default 400)
Tr2     = 300;  % Approximate number of triangles for the inner sphere (defgault 300)

[P1, t1]    = sphere(R1, Tr1);
centers1    = meshtricenter(P1, t1);
normals1    = meshnormals(P1, t1);

[P2, t2]    = sphere(R2, Tr2);
centers2    = meshtricenter(P2, t2);
normals2    = meshnormals(P2, t2);

P           = [P1; P2];
t           = [t1; t2+size(P1, 1)];
normals     = [normals1; -normals2];
snumber     = [1*ones(size(t1, 1), 1); 0*ones(size(t2, 1), 1)]; %   zeros for the cavity

save('hollow_sphere', 'P', 't', 'normals', 'snumber');

patch('vertices', P, 'faces', t, 'EdgeColor', 'k', 'FaceAlpha', 0.2,'FaceColor', [0.8 0.9 1.0]);
axis 'equal';  axis 'tight'; grid on;
xlabel('x, m'); ylabel('y, m'); zlabel('z, m');
Triangles = size(t)
Nodes = size(P)
Quality         = min(simpqual(P, t))
