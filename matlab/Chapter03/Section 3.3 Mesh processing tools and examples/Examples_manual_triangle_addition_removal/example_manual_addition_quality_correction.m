clear all; close all; warning off;
%   SYNTAX 
%   example_manual_addition_quality_correction
%   DESCRIPTION 
%   This script demonstrates the use of manual triangle addition for
%   correcting the mesh quality. Lowest-quality triangle and its neighbors
%   are removed. The resulting hole is filled out manually. The script uses
%   MATLAB script select3D.m by Joe Conti to select triangle vertices and
%   add triangles. To add triangles manually, select three vertices and hit
%   enter. Hit outside the mesh to stop the process
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

[FileName, PathName] = uigetfile('*.mat','Select the mesh file');
load(FileName);
[str.MasterMeshQ str.MasterMeshIndex] = min(simpqual(P, t));
str

% Find triangle with the lowest quality and remove it from the mesh
triangle = str.MasterMeshIndex;
t(triangle, :) = [];    

%   Remove all neighbor triangles - as many as necessary (control by M)
[NonManifoldAttached, edges, edgesb, edgesnm] = meshedges(P, t);
M = 1;
for m =1:M
    t(NonManifoldAttached, :) = [];
    [NonManifoldAttached, edges, edgesb, edgesnm] = meshedges(P, t);
end

%   Remove unused nodes, display border triangles and border edges
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

%   Add triangles manually by selecting three vertices and hitting enter
%   Hit outside the mesh to stop the process 
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

[str.MasterMeshQ str.MasterMeshIndex] = min(simpqual(P, t)); str

NewName =  strcat(FileName(1:end-4), '_mod', '.mat');
save(NewName, 'P', 't');

