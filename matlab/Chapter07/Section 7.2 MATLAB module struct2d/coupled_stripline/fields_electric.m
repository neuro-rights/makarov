%   SYNTAX
%   fields_electric
%   DESCRIPTION
%   This script plots the the E-field distribution in a plane
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

clear all
factor = 128;                                           %  scale factor for E-field scaling

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

% Calculate electric field
c = zeros(3, length(t));
E = zeros(3, length(t));
for m =1:length(t)
    c(:,m) = (P(:, t(1, m)) + P(:, t(2, m)) + P(:, t(3, m)))/3;
                                             % triangle center
end
Ntotal = length(mesh.e);
for m = 1:Ntotal
    dist_x      = c(1, :) - mesh.ec(1, m);
    dist_y      = c(2, :) - mesh.ec(2, m);    
    Distance    = dist_x.^2 + dist_y.^2 + eps;
    Kernel      = 1./Distance;
    Kernel      = 1/(2*pi*const.epsilon)*Kernel*output.charge(m);
    E(1, :)     = E(1,:) + dist_x.*Kernel;
    E(2, :)     = E(2,:) + dist_y.*Kernel;
end
e   = sqrt(E(1, :).^2 + E(2, :).^2);     % field magnitude

count = 0;                              % field within metal is zero
for i = 1:length(re_diel)    
    if ~isempty(findstr(set_formula, num2str(i)))
        count = count + 1;
        if ~re_diel(i)
            index = find(t(4, :)==count);
            e(index) = 0;
        end
    end
end

% Scale and plot electric field
a = figure;
variance    = sqrt(dot(e, e))/length(e);
scale       = factor*variance;
e(e<(-scale)) = -scale;
e(e>(+scale)) = +scale;

e      =(e + scale)/(2*scale);                  % normalized to 1
N      = 100;
C      = colormap(jet(N));                      % colormap
brighten(0.5);
e      = min(round(N*e) + 1, N);                % normalized to colormap

X = reshape(P(1,t(1:3, :)),[3, size(t, 2)]);
Y = reshape(P(2,t(1:3, :)),[3, size(t, 2)]);
Z = reshape(P(3,t(1:3, :)),[3, size(t, 2)]);
h = fill3(X, Y, Z, e, 'FaceAlpha', 1.0);
axis('equal'); xlabel('x, m'); ylabel('y, m'); zlabel('z, m');
title({'Total E-field - magnitude distribution;';...
       'blue - smallest field magnitude; red - largest field magnitude'; ''});
view(0, 90);