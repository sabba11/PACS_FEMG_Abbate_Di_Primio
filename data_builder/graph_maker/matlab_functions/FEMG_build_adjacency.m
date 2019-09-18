function [Matrix_adj] = FEMG_build_adjacency(n, con_vectors)

% This function takes in:
%
%     n       = number of ordered vertex of the graph
%
% con_vectors = matrix of row vectors that say for each vertex to which
%               ones is linked (with their numbers)
%
% It builds the adjacency matrix of the graph



%Adjancency matrix
Matrix_adj = zeros(n,n);

%Builder
for i = 1:n
    
    for j = 1:n
        
        if con_vectors(i,j) ~= 0
            
            Matrix_adj(i,con_vectors(i,j)) = 1;
            
        else
            
            break
            
        end
    end
end


end