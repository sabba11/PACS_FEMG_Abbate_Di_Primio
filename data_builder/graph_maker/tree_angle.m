% Tree graph generator
%
% This script generates a tree graph in a *.txt (source and target coordinates, then their
% buondary conditions). It is a 3 dimensional  graph set only on a x-y plane.
% You can build it 2-dimensional just by commenting the signaled line.
% This tree is made by changing the bifurcation angle (tightening), but manatining the edge length.
%
% Input parameters:
%      n      = number of the order of bifurcation (at least 2 to obtain
%		        a Y configuration.)
% out_filename = name of the output file. Put it in the ./data folder.
%     ndim    = dimension of the problem (in this case is 2)
%    length   = length of all Edges
%    theta    = initial angle

close all
clear all

addpath('./matlab_functions')

% SETTING PARAMTERS
% Number of order of bifurcation
nbif = 4;
% Setting dimension of the problem (comment what you want)
ndim = 2;
% ndim = 3;
% Output file name
out_filename = '../data/txt_files/tree_angle.txt';
% Edge lenght
length_edge = 1;
% Initial angle of bifurcation with respect to the axis which is x-axis
% Will be divided by the order of biforcation minus 1
theta = pi/4;
% Starting point coordinates
x_start = 0;
y_start = 1;

% COORDINATES & CONNECTION OF THE VERTEXES
% Total number of points
n = 2^nbif;
% Creating points vector
points = zeros(n,2);
con_vectors = zeros(n,n);
% First two points
points(1,:) = [x_start,y_start];
points(2,:) = [x_start+length_edge,y_start];
con_vectors(1,1) = 2;
con_vectors(2,1) = 1;
% Counter of points that are already done
counter = 2;
% Creating all other points with right connections:
for j = 2:nbif
        % After every order of bifurcation nodes are written in the opposite
        % way (bottom-top/top-bottom of the iteration before.(j index for)
    for i = 1:2^(j-2)
        counter = counter+1;
        % Every iteration assign the coordinates of two nodes from one of
        % the j-index iteration before this one
        points(counter,:) = [points(counter-1-3*(i-1),1)+cos(theta/(j-1))*length_edge,
            points(counter-1-3*(i-1),2)+sin(theta/(j-1))*length_edge];
        points(counter+1,:) = [points(counter-1-3*(i-1),1)+cos(theta/(j-1))*length_edge,
            points(counter-1-3*(i-1),2)-sin(theta/(j-1))*length_edge];
        % And then connects them:
        con_vectors(counter-1-3*(i-1),2:3) = [counter,counter+1];
        con_vectors(counter,1) = counter-1-3*(i-1);
        con_vectors(counter+1,1) = counter-1-3*(i-1);
        counter = counter+1;
    end
end
%Augmenting dimension
if ndim == 3
    points = FEMG_augment_dim(points);
end
% Building adjacency
[Matrix_adj] = FEMG_build_adjacency(n, con_vectors);
FEMG_check_graph(Matrix_adj,n)

% BOUNDARY CONDITIONS
% We give only inflow in one tips and the outflow and int for central
BC = zeros(n,5);
BC(1,:) = ['DIR ', '1'];
% Assign outflow and internal
[BC] = FEMG_assign_cond(BC,n,Matrix_adj);

% Writing output and plotting the graph
[Edges,Finaltext] = FEMG_build_graphtext(ndim,n,points,Matrix_adj,BC);
% Plot
if ndim == 2
    FEMG_plot2d_graph(Edges,points)
else
    FEMG_plot3d_graph(Edges,points)
end
dlmwrite(out_filename, Finaltext, 'delimiter', '');
