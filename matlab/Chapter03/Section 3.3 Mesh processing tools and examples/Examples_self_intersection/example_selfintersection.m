clear all; close all;
warning off; %  to avoid globals
%   SYNTAX 
%   example_self_intersection
%   DESCRIPTION
%   This script corrects mesh non-manifold errors and/or mesh
%   self-intersection errors by removing triangles (as many of them as
%   necessary) from the suspicious areas and then filling resulting holes
%   with the good (non-intersecting) triangles manually. It uses MATLAB
%   script select3D.m by Joe Conti to select triangle vertices and add
%   triangles.
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

[FileName, PathName] = uigetfile('*.mat','Select the mesh file for the self-intersection check');
load(FileName);
[str.MasterMeshQ str.MasterMeshIndex] = min(simpqual(P, t));
str

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Check the manifoldness condition first and collect suspicious triangles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edges             = meshconnee(t);
AttachedTriangles = meshconnet(t, edges, 'nonmanifold');
NonManifoldAttached = [];
for m = 1:size(edges, 1)
    if length(AttachedTriangles{m})~=2        
        NonManifoldAttached = [NonManifoldAttached; AttachedTriangles{m}];
    end
end        
triangles = NonManifoldAttached;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Check self-intersections next and collect suspicious triangles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[si11, si12, si2] = meshedgeintersect(P, t, P, edges);
for m = 1:size(si11, 1)
    triangles = [triangles; si11{m}];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Step 1 Return if the mesh is good
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(triangles)
    display('The mesh is 2 manifold and has no self-intersections')
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Step 2 Remove all intersected/non-manifold triangles from the mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t(triangles, :) = [];    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Step 3 Remove all neighbor triangles - as many as necessary (control by M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[NonManifoldAttached, edges, edgesb, edgesnm] = meshedges(P, t);
M = 1;
for m =1:M
    t(NonManifoldAttached, :) = [];
    [NonManifoldAttached, edges, edgesb, edgesnm] = meshedges(P, t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Step 4 Remove unused nodes, display border triangles and border edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[P, t] = fixmesh(P, t);
[NonManifoldAttached, edges, edgesb, edgesnm] = meshedges(P, t);
nodes = reshape(t(NonManifoldAttached, :), 1, 3*length(NonManifoldAttached));
nodes = unique(nodes);
for m = 1:size(edgesb, 1)
    marker1(m) = line('xdata', P(edgesb(m, :),1) ,'ydata', P(edgesb(m, :),2) ,'zdata', P(edgesb(m, :),3),...
    'color', 'b', 'linewidth', 4);           
end
patch('Faces', t(NonManifoldAttached, :), 'Vertices', P, 'FaceColor', 'c', 'EdgeColor', 'k', 'FaceAlpha', 1.0);
grid on; axis equal; axis tight;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Step 5 Add triangles manually by selecting three vertices and hitting enter
%   Hit outside the figure to stop the process 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selection = [];
while 1   
    % p - clicked point; % v - nearest vertex; vi - index into the nearest
    % vertex; facei - index into the clicked face
    triangle = []; flag = 0;
    for m = 1:3
        pause;
        [p v vi face facei] = select3d;        
        if ~isempty(vi)
            triangle  = [triangle vi]; flag = 1;
            marker(m) = line('xdata',v(1),'ydata',v(2),'zdata',v(3),'marker','o',...
                    'erasemode','xor','markerfacecolor','b');         
        else
            flag = 0; break;
        end       
    end
    if flag ==0 break; end;     
    triangle
    t = [t; triangle];
    delete(marker);
    patch('Faces', triangle, 'Vertices', P, 'FaceColor', 'g', 'EdgeColor', 'k', 'FaceAlpha', 1.0);
    drawnow;
end
%   Show node numbers
text(P(nodes, 1), P(nodes, 2), P(nodes, 3), num2str(nodes'), 'FontSize', 20, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');

NewName =  strcat(FileName(1:end-4), '_mod', '.mat');
save(NewName, 'P', 't');

