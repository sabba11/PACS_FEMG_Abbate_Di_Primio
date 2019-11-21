% FEMG_check_BC.m
% Checks Kirchhoff-Neumann boundary conditions at the nodes of the graph.
%
% INPUT:
% solution     a vector containing the solution of a differential problem.
%
% mf           a GetFEM++ MeshFem object. It represents the mesh on which
%              finite elements are built.
%
% mesh         a GetFEM++ Mesh object. It represents the mesh on which the
%              the eigenvectors are computed. It may possibly be equal
%              to mf.

function FEMG_check_BC(solution, mf, mesh)
   access_idxs = FEMG_assign_vertices(mf, mesh);
   RIDs = gf_mesh_get(mesh, 'regions');
   CVFIDs = gf_mesh_get(mesh, 'region', RIDs(end)); % get convexes in last region
   [DOFs, IDx] = gf_mesh_fem_get(mf, 'basic dof from cvid');
   PTs = gf_mesh_get(mesh, 'pts');
   
   n_branches = size(RIDs,2) - 1; % number of branches
   total_nodes = 0;
   for i = 1 : n_branches
       cvf = gf_mesh_get(mesh, 'region', RIDs(i)); % get convexes in region i
       dof = gf_mesh_fem_get(mf, 'basic dof from cvid', cvf(1,:));
       dof = unique(dof);
       total_nodes = total_nodes + length(dof);
   end
   node_idx = total_nodes + 1; % first real node index
   points = [];
   for i = 1 : size(CVFIDs,2)
       points = [points, DOFs(IDx(CVFIDs(1,i)):IDx(CVFIDs(1,i)+1)-1)'];
   end
   bc_check = cell(2, length(unique(DOFs)) - node_idx + 1);
   bc_check(1, :) = num2cell([node_idx:length(unique(DOFs))]);
   for j = 1 : size(points,2)
       current = points(:, j);
       distance = norm(PTs(:,points(1,j)) - PTs(:, points(2,j)));
       check = (current >= node_idx);
       if (sum(check) ~= 1)
           disp('Convex is not valid to check boundary conditions');
           return;
       end
       id = find(check == 1);
       aux = current(id) - node_idx + 1;
       bc_check{2,aux} = [bc_check{2,aux}, current];
   end
   for k = 1 : size(bc_check,2)
       for l = 1 : size(bc_check(2,k))
           value = 0;
           current2 = bc_check{2,k};
           current3 = current2(:,l);
           id = find(current3 >= node_idx);
           eval = solution(current3);
           if id == 1
               value = value + (eval(1)-eval(2))/distance;
           elseif id == 2
               value = value + (eval(2)-eval(1))/distance;
           else
               disp('Convex is not one-dimensional');
               return;
           end
       end
       bc_check{2,k} = [];
       bc_check{2,k} = value;
   end
   figure()
   hold on
   grid on
   for m = 1 : size(bc_check,2)
       plot(bc_check{1,m}, bc_check{2,m}, 'ko');
   end
   title('Boundary conditions at each node');
end





