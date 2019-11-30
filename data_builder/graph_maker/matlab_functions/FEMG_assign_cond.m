% FEMG_assign_cond

% This function automatically assign natural condition to the missing point
% of the graph.
% Internal points will have INT condition (Neumann-Kirkhoff) with a dummy 0
% value.
% External point will have MIX conditon (Neumann) automatically set at 0.

% This function takes in:
%    BC      = matrix n*2 to save BC and their values
%     n      = number of ordered vertex of the graph
% Matrix_adj = adjacency matrix of the graph

function [BC] = FEMG_assign_cond(BC,n,Matrix_adj)
% This function works well only if you run check_graph.m before.
% It takes in input only undirectional graphs with all vertexes connected.
for i = 1:n
    if BC(i,1) == 0
        %If outflow vertex put mixed condition with default zero value
        if sum(Matrix_adj(i,:))==1
            BC(i,:) = ['MIX ','0'];
        %If internal vertex put internal condition with default zero value
        else
            BC(i,:) = ['INT ','0'];
        end
    end
end
end