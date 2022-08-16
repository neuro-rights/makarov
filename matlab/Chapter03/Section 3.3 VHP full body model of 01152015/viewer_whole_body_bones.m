clear all; clc; warning off; 
%   Full-body viewer (bones only) for the VHP model of 01/15/2015
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.
%   Find edges (indexes) attached to every triangle (Nx3 array)

count = 0;
scrsz = get(0, 'ScreenSize');
a = figure('Position', [scrsz(4)/10 scrsz(4)/10 scrsz(3)/1.2 scrsz(4)/1.2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
load('VHP_Calcaneous_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 0.5,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Calcaneous_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 0.5,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Average_body.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'none', 'FaceAlpha', 0.15,'FaceColor', [1 0.75 0.65]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Feet1Phalange_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Feet2Phalange_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Feet3Phalange_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Feet4Phalange_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Feet5Phalange_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Feet1Phalange_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Feet2Phalange_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Feet3Phalange_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Feet4Phalange_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Feet5Phalange_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Femur_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'none', 'FaceAlpha', 0.2,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Femur_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'none', 'FaceAlpha', 0.2,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Fibula_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Fibula_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Hands1Phalange_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Hands2Phalange_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Hands3Phalange_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Hands4Phalange_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Hands5Phalange_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Hands1Phalange_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Hands2Phalange_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Hands3Phalange_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Hands4Phalange_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Hands5Phalange_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Humerus_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Humerus_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Jaw_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
load('VHP_Navicular_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1.0,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Navicular_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1.0,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Patella_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Patella_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Pelvic_Girdle_M.mat');   % for pelvic girdle
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 0.8,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Ribs_left11_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Ribs_left12_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Ribs_right10_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Ribs_right11_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Ribs_right12_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Ribs_Cartilage_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Scapula_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Scapula_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Skin.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'none', 'FaceAlpha', 0.15,'FaceColor', [1 0.75 0.65]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Skull.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1.0,'FaceColor', 'r');
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Spine_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
load('VHP_Talus_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Talus_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_lower17_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_lower18_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_lower19_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_lower28_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_lower29_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_lower30_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_lower31_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_upper11_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_upper11_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_upper12_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_upper13_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_upper14_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_upper15_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_upper16_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_upper2_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_upper3_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_upper4_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_upper5_6_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_upper7_8_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Teeth_upper9_10_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Tibia_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 0.5,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Tibia_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 0.5,'FaceColor', [0.87 0.49 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Trabecular_lower_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'none', 'FaceAlpha', 0.35,'FaceColor', 'y');
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Trabecular_lower_right.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'none', 'FaceAlpha', 0.35,'FaceColor', 'y');
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Trabecular_upper_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'none', 'FaceAlpha', 0.35,'FaceColor', 'y');
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Trabecular_upper_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'none', 'FaceAlpha', 0.35,'FaceColor', 'y');
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Ulna_Radius_left_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VHP_Ulna_Radius_right_M.mat');
PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1,'FaceColor', [1 0 0]);
count = count + size(t,1);

disp('Total number of triangles:'); count
axis 'equal';  axis 'tight';
xlabel('x, mm'); ylabel('y, mm'); zlabel('z, mm');
view(0,0);
set(gcf,'Color','White');
axis off;
camroll(90);
cameratoolbar('Show');
camlight('headlight');




