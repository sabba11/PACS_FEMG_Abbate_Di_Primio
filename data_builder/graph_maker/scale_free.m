% Barabasi-Albert Graph
% 
% This file builds a 2-dimensional scale-free graph with the
% Barabasi-Albert algorithm found in the Contest MATLAB toolkit.
% After building the adjacency  matrix the points will be taken on a
% circle.
% If you want the points to be 3 dimensional with a zero z-coordinate
% untoggle the comments at the signaled line.

close all
clear all

% Adding needed function (FEMG & Contest) and output directory
addpath('./matlab_functions')
out_filename = '../data/txt_files/scale_free_prova.txt';

% Dimension (comment the unwanted one)
ndim = 2;
% ndim = 3;

% Circle information
% Circle radius
L = 1;
% Center of the circle
x = 0;
y = 0;

% Number of vertices required (even)
N = 20;

% Building geometrical points of the graph
points = zeros(N,2);
theta = pi/(N/2);
scanning_x = x-L;
% Extrema on the circle
points(1,:) = [scanning_x, y];
points(N,:) = [x+L, y];
% Building other points on the circle equispace
scanning_index = 2;
for i = 1:(N-2)/2
    scanning_x = x-L*cos(theta*i);
    points(scanning_index,:) = [scanning_x, y + sqrt(L^2 - (scanning_x-x)^2)];
    scanning_index = scanning_index + 1;
    points(scanning_index,:) = [scanning_x, y - sqrt(L^2 - (scanning_x-x)^2)];
    scanning_index = scanning_index + 1;
end

%Augmenting dimension
if ndim == 3
    points = FEMG_augment_dim(points);
end

% Building adjacency matrix and natural Boundary Condition for the graph
[Matrix_adj] = CONTEST_pref(N, 2);
BC = zeros(N,5);
Matrix_adj = full(Matrix_adj);
[BC] = FEMG_assign_cond(BC,N,Matrix_adj);

% Writing output and plotting the graph
% Builing .txt
[Edges,Finaltext] = FEMG_build_graphtext(ndim,N,points,Matrix_adj,BC);
% Plot
if ndim == 2
    FEMG_plot2d_graph(Edges,points)
else
    FEMG_plot3d_graph(Edges,points)
end
%Exporting .txt
dlmwrite(out_filename, Finaltext, 'delimiter', '');