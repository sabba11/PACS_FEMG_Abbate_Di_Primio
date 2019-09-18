function [BC] = FEMG_assign_cond(BC,n,Matrix_adj)

% This function works well only if you run check_graph.m before (It takes in
% input only undirectional graphs with all vertexes connected


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
