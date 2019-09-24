% Graphene graph generator
%
% This script generates a simple graphene graph in a .txt (giving edges with
% initial and final points coordinates with boundary conditions afterwards).
% It is a 2 dimensional graph in the box [0,L;0,L] made by two exagon
% linked by an edge
% The setting and the points are saved as like as it is a 3d graph
%
% Input parameters:
% out_filename = name of the output file. Put it in the ./data folder.
%     ndim    = dimension of the problem (in this case is 2)
%       L      = length of the box containing the graphene



close all
clear all

addpath('./matlab_functions')

% Number of vertexes
n = 12;
% Setting dimension of the problem
ndim = 3;
% Output file name
out_filename = '../data/graphene.txt';

% Box lenght
L = 5;
% Vettore coordinate
points = zeros(n,2);

% COORDINATES
% The two exagons

for i= 1:6
    points(i,:) = [L/5+L/5*cos((i-1)*2*pi/6),L/2+L/5*sin((i-1)*2*pi/6)];
    points(i+6,:) = [4*L/5+L/5*cos((i-1)*2*pi/6),L/2+L/5*sin((i-1)*2*pi/6)];
end

points = FEMG_augment_dim(points);

points(:,2) = points(:,2) - ones(length(points(:,2)),1)*points(5,2) ;

% CONNECTION OF THE VERTEXES
con_vectors = zeros(n,n);
% Connection middle points
for i = 2:5
    con_vectors(i,1:2) = [i-1,i+1];
    con_vectors(i+6,1:2) = [i+5,i+7];
end
% Connection first and last point of the two exagon
con_vectors(1,1:2) = [2,6];
con_vectors(6,1:2) = [1,5];
con_vectors(7,1:2) = [8,12];
con_vectors(12,1:2) = [7,11];

% Central line
con_vectors(1,3) = 10;
con_vectors(10,3) = 1;
% Building adjacency
[Matrix_adj] = FEMG_build_adjacency(n, con_vectors);
FEMG_check_graph(Matrix_adj,n)

% BOUNDARY CONDITIONS
% We give only inflow in one tips and the outflow and int for central
BC = zeros(n,5);
% Assign outflow and internal
[BC] = FEMG_assign_cond(BC,n,Matrix_adj);



% Writing output and plotting the graph
[Edges,Finaltext] = FEMG_build_graphtext(ndim,n,points,Matrix_adj,BC);

%FEMG_plot2d_graph(Edges,points)
FEMG_plot3d_graph(Edges,points)

dlmwrite(out_filename, Finaltext, 'delimiter', '');
