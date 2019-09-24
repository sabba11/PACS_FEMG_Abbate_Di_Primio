% Star graph generator
%
% This script generates a star graph in a .txt (giving edges with initial and
% final points coordinates with boundary conditions afterwards).
% It is a 2 dimensional graph in a 3d points setting put in the box
% [0,L;0,L].
% Then the star will be translate of (x_trans,y_trans).
%
% Input parameters:
%      n      = points of the star (at least 3 to obtain a Y conf)
% out_filename = name of the output file. Put it in the ./data folder.
%     ndim    = dimension of the problem (in this case is 3 even if the
% final star will be only on one plane)
%     L       = dimension of the box in which the star will be contained
% x_trans, y_trans = setting for translation

close all
clear all

addpath('./matlab_functions')

% Number of points of the star
n = 4;
% Setting dimension of the problem
ndim = 2;
% ndim = 3;
% Output file name
out_filename = '../data/star.txt';

% Box lenght
L = 2;
% Setting factor for translation
x_trans = -1;
y_trans = -1;
% Total point = points(tips) of the star+ central
ntot = n+1;
points = zeros(ntot,2);

% COORDINATES
% Central point
points(1,:) = [L/2,L/2];
% Tips of the star
for i = 2:ntot
   points(i,:) = [L/2+L/2*cos((i-2)*2*pi/n),L/2+L/2*sin((i-2)*2*pi/n)];
end

%Augmenting dimension
if ndim == 3
    points = FEMG_augment_dim(points);
end
points(:,2) = points(:,2) - ones(length(points(:,2)),1);
points(:,1) = points(:,1) - ones(length(points(:,2)),1);

% CONNECTION OF THE VERTEXES
con_vectors = zeros(ntot,ntot);
% Central point connected wiht all
con_vectors(1,:) = [2:ntot,0];
% Tips connected with central
con_vectors(2:ntot,1) = 1;

% Building adjacency
[Matrix_adj] = FEMG_build_adjacency(ntot, con_vectors);
FEMG_check_graph(Matrix_adj,n)

% BOUNDARY CONDITIONS
% We give only inflow in one tips and the outflow and int for central
BC = zeros(ntot,5);
BC(2,:) = ['DIR ', '1'];
% Assign outflow and internal
[BC] = FEMG_assign_cond(BC,ntot,Matrix_adj);



% Writing output and plotting the graph
[Edges,Finaltext] = FEMG_build_graphtext(ndim,ntot,points,Matrix_adj,BC);

if ndim == 2
    FEMG_plot2d_graph(Edges,points)
else
    FEMG_plot3d_graph(Edges,points)
end


dlmwrite(out_filename, Finaltext, 'delimiter', '');
