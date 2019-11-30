% FEMG_plot3d_graph

% Plotting graph in the 3 dimensional case
% It takes as inputs:
% Edges = matrix containg the edges 6*n
% pointsd = matrix 3*n containing the coordiantes of the points

function FEMG_plot3d_graph(Edges,points)
%Plotting edges
for i = 1:size(Edges,1)
    plot3(Edges(i,[1,4]), Edges(i,[2,5]), Edges(i,[3,6]), '-k', 'LineWidth', 2)
    hold on
end
% Plotting points
plot3(points(:,1), points(:,2),points(:,3), 'r*', 'MarkerSize', 10, 'LineWidth', 10)
end