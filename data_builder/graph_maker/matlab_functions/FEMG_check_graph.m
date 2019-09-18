function FEMG_check_graph(Matrix_adj,n)

% This function simply checks if the graph has no direction in the edge
% (symmetric adjacency matrix) and if all vertexes are connected.


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