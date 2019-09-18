function [Edges,Finaltext] = FEMG_build_graphtext(ndim,n,points,Matrix_adj,BC)

%  This function takes in:
%
%      ndim   = dimension of the problem
%
%       n     = number of vertex of the graph
%
%     points  = matrix n*ndim with the vector of coordinates of the points
%
%  Matrix_adj = adjacency matrix of the graph
%
%      BC     = matrix n*2 of type of BC and value
%
% This function works only with a graph with no direction!

%Generation of the empty text
Edges = zeros(sum(sum(Matrix_adj)/2), 2*ndim);

%Counter for the row
counter = 1;


for i = 1:n-1
    
    for j = i+1:n
        
        if Matrix_adj(i,j) == 1
            
            % Put coordinates of beginning and ending vertex of the edge
            % and their boundary conditions
            
            Edges(counter,:) = [points(i,:), points(j,:)];
            
            % Costruction in function of number of coordinates,
            % for now only 2 or 3 could be better implemented
            
            
            if ndim == 2
                % Those two cicles are implemented to handle the minus sign
                for k = 1:ndim
                    if points(i,k)>=0
                        ch(k) = ' ';
                    else
                        ch(k) = '-';
                    end
                end
                
                for k = 1:ndim
                    if points(j,k)>=0
                        ch(2+k) = ' ';
                    else
                        ch(2+k) = '-';
                    end
                end
                
                Finaltext(counter, :) = [ch(1), num2str(abs(points(i,1)),'%f'), ' ', ch(2), num2str(abs(points(i,2)),'%f'),...
                    ' ', ch(3), num2str(abs(points(j,1)),'%f'), ' ', ch(4), num2str(abs(points(j,2)),'%f'), ...
                    ' ', BC(i,:), ' ', BC(j,:)];
            end
            
            if ndim == 3
                % Those two cicles are implemented to handle the minus sign
                for k = 1:ndim
                    if points(i,k)>=0
                        ch(k) = ' ';
                    else
                        ch(k) = '-';
                    end
                end
                
                for k = 1:ndim
                    if points(j,k)>=0
                        ch(3+k) = ' ';
                    else
                        ch(3+k) = '-';
                    end
                end

                Finaltext(counter, :) = [ ch(1), num2str(abs(points(i,1)),'%f'), ' ', ch(2), num2str(abs(points(i,2)),'%f'),...
                    ' ', ch(3), num2str(abs(points(i,3)),'%f'), ' ', ch(4), num2str(abs(points(j,1)),'%f'), ...
                    ' ', ch(5), num2str(abs(points(j,2)),'%f'), ' ', ch(6), num2str(abs(points(j,3)),'%f'), ...
                    ' ', BC(i,:), ' ', BC(j,:)];
            end
            
            %Counter of the row
            counter = counter+1;
            
        end
    end
end

end
