% Tree graph generator
%
% This script generates a tree graph in a .txt (giving edges with initial and
% final points coordinates with boundary conditions afterwards). It is a 3 dimensional
% graph set only on a x-y plane in the box [0,1;0,1];
%
% Input parameters:
%      n      = number of the order of bifurcation (at least 2 to obtain
%		        a Y configuration.)
% out_filename = name of the output file. Put it in the ./data folder.
%     ndim    = dimension of the problem (in this case is 2)



close all
clear all

addpath('./matlab_functions')

% Number of order of bifurcation
nbif = 4;
% Setting dimension of the problem
ndim = 3;
% Output file name
out_filename = '../data/simple_tree_edge.txt';

% Box lenght: unitary
L = 1;
% Total number of points
n = 2^nbif;

%Angle of bifurcation with respect to the axis which is x-axis
theta = pi/4;

% Length of the edge
% Will be divided by the order of biforcation minus 1
length_edge = L/(1+cos(theta)*(nbif-1));


% COORDINATES & CONNECTION OF THE VERTEXES
points = zeros(n,2);
con_vectors = zeros(n,n);

% First two points

points(1,:) = [L,L/2];
points(2,:) = [L-length_edge,L/2];


con_vectors(1,1) = 2;
con_vectors(2,1) = 1;
% Counter of points that are already done
counter = 2;

for j = 2:nbif
        % after every order of bifurcation nodes are written in the opposite
        % way (bottom-top/top-bottom of the iteration before.(j index for)
    for i = 1:2^(j-2)
        counter = counter+1;
        % every iteration assign the coordinates of two nodes from one of
        % the j-index iteration before this one
        points(counter,:) = [points(counter-1-3*(i-1),1)-cos(theta)*length_edge/(j-1),
            points(counter-1-3*(i-1),2)+sin(theta)*length_edge/(j-1)];
        points(counter+1,:) = [points(counter-1-3*(i-1),1)-cos(theta)*length_edge/(j-1),
            points(counter-1-3*(i-1),2)-sin(theta)*length_edge/(j-1)];

        % And then connects them

        con_vectors(counter-1-3*(i-1),2:3) = [counter,counter+1];
        con_vectors(counter,1) = counter-1-3*(i-1);
        con_vectors(counter+1,1) = counter-1-3*(i-1);

        counter = counter+1;
    end
end

% Checking that the graph makes sense and stay in the box
if points>1 & points<0
    error('Graph out of the box, change edge length or angle')
end
points = FEMG_augment_dim(points);
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

%FEMG_plot2d_graph(Edges,points)
FEMG_plot3d_graph(Edges,points)

dlmwrite(out_filename, Finaltext, 'delimiter', '');
