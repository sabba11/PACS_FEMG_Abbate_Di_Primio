% Graphene graph generator
%
% This script generates a simple graphene graph in a .txt (giving edges with
% initial and final points coordinates with boundary conditions afterwards).
% It is a 2 dimensional graph in the box [0,L;0,L] made by two exagon
% linked by an edge
% The setting and the points are saved as like as it is a 3d graph;
% You can build it 2-dimensional just by commenting the signaled line.
%
% Input parameters:
% out_filename = name of the output file. Put it in the ./data folder.
%     ndim     = dimension of the problem (in this case is 2)
%       L      = length of the box containing the graphene

close all
clear all
addpath('./matlab_functions')

% SETTING PARAMETERS
% Number of vertexes
n = 11;
% Setting dimension of the problem
ndim = 2;
% ndim = 3;
% Output file name
out_filename = '../data/txt_files/graphene_noline.txt';
% Box lenght: unitary
L = 4;

% COORDINATES
% Coordinates Vector
points = zeros(n,2);
% The two exagons
for i = 1:6
    points(i,:) = [L/4+L/4*cos((i-1)*2*pi/6),L/2+L/4*sin((i-1)*2*pi/6)];
    %second exagon
    if i ~=6
    points(i+6,:) = [3*L/4+L/4*cos((i+3)*2*pi/6),L/2+L/4*sin((i+3)*2*pi/6)];
    end
end
%Augmenting dimension
if ndim == 3
    points = FEMG_augment_dim(points);
end
points(:,2) = points(:,2) - ones(length(points(:,2)),1)*points(5,2) ;

% CONNECTION OF THE VERTEXES
con_vectors = zeros(n,n);
% Connection middle points
for i = 2:5
    con_vectors(i,1:2) = [i-1,i+1];
    %second exagon
    if i~=5
    con_vectors(i+6,1:2) = [i+5,i+7];
    end
end
% Connection with first central point fixed
con_vectors(1,1:4) = [2,6,7,11];
con_vectors(6,1:2) = [1,5];
con_vectors(7,1:2) = [1,8];
con_vectors(11,1:2) = [1,10];
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
% Plot 
if ndim == 2
    FEMG_plot2d_graph(Edges,points)
else
    FEMG_plot3d_graph(Edges,points)
end
dlmwrite(out_filename, Finaltext, 'delimiter', '');
