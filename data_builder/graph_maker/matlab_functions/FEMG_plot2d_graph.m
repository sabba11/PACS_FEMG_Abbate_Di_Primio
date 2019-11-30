% FEMG_plot2d_graph

% Plotting graph in the 2 dimensional case
% It takes as inputs:
% Edges = matrix containg the edges 4*n
% pointsd = matrix 2*n containing the coordiantes of the points

function FEMG_plot2d_graph(Edges,points)
%Plotting edges
for i = 1:size(Edges,1)
    plot(Edges(i,[1,3]), Edges(i,[2,4]), '-k', 'LineWidth', 2)
    hold on
end
% Plotting points
plot(points(:,1), points(:,2), 'r*', 'MarkerSize', 10, 'LineWidth', 10)
end