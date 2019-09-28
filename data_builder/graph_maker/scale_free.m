close all
clear all

addpath('./matlab_functions')
out_filename = '../data/scale_free.txt';

ndim = 2;
% ndim = 3;

%Circle radius
L = 4;

%Center of the circle
x = 4;
y = 4;

%Number of vertices required (even)
N = 20;

points = zeros(N,2);

theta = pi/(N/2);
scanning_x = x-L;

points(1,:) = [scanning_x, y];
points(N,:) = [x+L, y];

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

[Matrix_adj] = CONTEST_pref(N, 2);
BC = zeros(N,5);
Matrix_adj = full(Matrix_adj);
[BC] = FEMG_assign_cond(BC,N,Matrix_adj);
% Writing output and plotting the graph

[Edges,Finaltext] = FEMG_build_graphtext(ndim,N,points,Matrix_adj,BC);

if ndim == 2
    FEMG_plot2d_graph(Edges,points)
else
    FEMG_plot3d_graph(Edges,points)
end

dlmwrite(out_filename, Finaltext, 'delimiter', '');
