% FEMG_check_graph

% This function simply checks if the graph has no direction in the edge
% (symmetric adjacency matrix) and if all vertexes are connected.
% It takes as inputs:
%       n     = number of vertex of the graph
%  Matrix_adj = adjacency matrix of the graph

function FEMG_check_graph(Matrix_adj,n)
%Checking no direction graph
if issymmetric(Matrix_adj) == 0    
    error('Adjacency matrix not symmetric, edges should not have direction')    
end
%Checking all vertexes are connected
for i = 1:n    
    if sum(Matrix_adj(i,:))<1        
        error('One vertex is not connected')        
    end
end
end