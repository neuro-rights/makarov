function [h] = viewer(arg1, arg2, arg3);
%   SYNTAX
%   viewer(P,t)	viewer(P,t,S)	viewer(geom)
%   h = viewer(...)
%   DESCRIPTION
%   Similar to original MATLAB function TRISURF, the VIEWER function can be
%   used to view surface meshes within MATLAB.  It accepts mesh input in a
%   variety of ways:
% 
%   1.	VIEWER(P,t) displays a surface mesh defined by the point matrix P
%   and the triangle matrix t, as described above.  The mesh is displayed in
%   MATLAB's standard 3-dimensional plot window, and it can be rotated and
%   scaled with the usual tools.  Each triangle will be assigned a color
%   based on its domain number.
%
%   2.	VIEWER(P,t,S) also displays a surface mesh defined by the point
%   matrix P and the triangle matrix t, but the coloring is specified by the
%   matrix S.  S is an array with one element per column in P, each element
%   specifying a value for its corresponding point.  These values are
%   displayed as the relative brightness of each point in the surface mesh.
%
%   3.	VIEWER(geom) displays a surface mesh defined by data stored in the
%   structure geom.  This structure must contain a point matrix P and a
%   triangle matrix t, in which case the command will function as in case 1
%   above.  If the structure also contains an array S, the command will
%   function as in case 2.
%
%   4.	VIEWER('filename') displays a surface mesh defined by data saved in
%   the file filename.mat.  As in case 3, this file must contain a point
%   matrix P and a triangle matrix t.  If it contains an array called S, this
%   will be used to provide coloring for the surface's vertices.
%
%   h = viewer(...) returns a vector of handles to patch graphics objects,
%   one handle per patch.
%   Inputs:
%   P, array of nodes
%   t, array of triangles
%   
%   Authors: M. Liffiton, S.N. Makarov
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

nargs = nargin; 

if  nargs < 1
    error('Requires at least 1 input.');
end

if nargs == 1
    if (strcmp(class(arg1), 'char')) & (exist(arg1) == 2)
        load(arg1);
        
        if (exist('P') ~= 1) | (exist('t') ~= 1)
            error('.mat file must contain ''P'' point matrix and ''t'' triangle matrix.');
        elseif exist('S') ~= 1
            nargs = 2;
        elseif exist('geom') ~= 1
            nargs = 3;
        else
            nargs = 4;
        end
    elseif strcmp(class(arg1), 'struct')
        if ~isfield(arg1,'P') | ~isfield(arg1,'t')
            error('Mesh structure must contain ''P'' point matrix and ''t'' triangle matrix.');
        elseif ~isfield(arg1,'S')
            P = geom.P;
            t = geom.t;
            nargs = 2;
        else
            P = geom.P;
            t = geom.t;
            S = geom.S;
            nargs = 3;
        end
    else
        error('First argument must be a path to a .mat file or a struct containing mesh data.');
    end
else
    P = arg1;
    t = arg2;
    if nargs == 3
        S = arg3;
    end
end

X = reshape(P(1,t(1:3, :)),[3, size(t,2)]);
Y = reshape(P(2,t(1:3, :)),[3, size(t,2)]);
Z = reshape(P(3,t(1:3, :)),[3, size(t,2)]);

if nargs == 2
    if size(t,1) ==4        
        C = t(4, :);
        colormap(summer);
    else
        C = [0.6, 0.6, 0.6];
    end
else
    C = S(t(1:3, :));
    colormap bone;
end

h = fill3(X, Y, Z, C, 'FaceAlpha', 1.0);
axis('equal');
xlabel('x, m');
ylabel('y, m');
zlabel('z, m');
max_domain  = max(t(4, :));
title({'Domain number (from low to high) is displayed by the domain color -'; ...
       'from green (low) to yellow (high)'; strcat('Total domains: ', num2str(max_domain))});

