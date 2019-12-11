% FEMG_check_BC.m
% Checks Kirchhoff-Neumann and Dirichlet boundary conditions at the 
% original nodes of the graph, plotting results.
%
% OUTPUT:
% bc_check     a cell array containing original node global indexes in the
%              mesh (first row), type of BC (second row, DIR or NEU),
%              computed value (third row) and coordinate array (fourth row).
%              Values in the third row are directly computed from the
%              solution and are to be confronted with target BCs.
% INPUT:
% solution     a vector containing the solution of a differential problem.
%
% mf           a GetFEM++ MeshFem object. It represents the mesh on which
%              finite elements are built.
%
% mesh         a GetFEM++ Mesh object. It represents the mesh on which the
%              the eigenvectors are computed. It may possibly be equal
%              to mf. BC data should be contained in its regions.
%
% n_bc         number of BC types in the problem. If the problem is fully
%              Dirichlet or Neumann-Kirchhoff n_bc equals 1, otherwise n_bc
%              equals 2. No other values can be assigned.
%
% subtitle     subtitle to add to the plot title. Optional argument,
%              efaults to an ampty string.
%
% Remark. The GetFEM++ Mesh object should have n_branches + 3 regions,
%         where n_branches is the number of branches in the original graph.
%         The first n_branches regions (with indexes running from 0 to
%         n_branches - 1) contain convexes lying on a fixed branch, whose
%         endpoints are discretization nodes (not belonging to the original
%         graph). They are the "interal" convexes of each branch.
%         The region with index n_branches contains convexes whose
%         endpoints are constitued by one real point and one discretization
%         node.
%         The region with index n_branches + 1 contains the faces of the
%         convexes (hence nodes) in which Neumann-Kirchoff boundary 
%         condtions are imposed.
%         The region with index n_branches + 2 contains contains the faces
%         of the convexes (hence nodes) in which Dirichlet boundary
%         conditions are imposed.

function bc_check = FEMG_check_BC(solution, mf, mesh, n_bc, subtitle)
   % Check input args
   if (nargin < 5)
       subtitle = '';
   end
   
   if (n_bc ~= 1 && n_bc ~= 2)
       disp('Invalid number of BC types: should be 1 or 2');
       return;
   end
   bc_label = 'MIX'; % default is mixed problem
   
   % Ask for which BC has been chosen (if 1)
   if (n_bc == 1)
       bc_label = input('Insert BC type (DIR for Dirichlet-type or NEU for Neumann-type):\n','s');
       if (strcmp(bc_label, 'DIR') == 0 && strcmp(bc_label, 'NEU') == 0)
            disp('Invalid BC label (should be DIR or NEU)');
            return;
       end
   end
   
   % Read mesh data
   RIDs = gf_mesh_get(mesh, 'regions'); % get mesh regions
   CVFIDs = gf_mesh_get(mesh, 'region', RIDs(end-n_bc)); % get convexes in region 
   [DOFs, IDx] = gf_mesh_fem_get(mf, 'basic dof from cvid'); % get dofs
   PTs = gf_mesh_get(mesh, 'pts'); % get point coords
   if (bc_label == "MIX")
       neu_CVFIDs = gf_mesh_get(mesh, 'region', RIDs(end-1)); % get Neumann BC data
       dir_CVFIDs = gf_mesh_get(mesh, 'region', RIDs(end)); % get Dirichlet BC data
   elseif (bc_label == "DIR")
       dir_CVFIDs = gf_mesh_get(mesh, 'region', RIDs(end)); % get Dirichlet BC data
   else
       neu_CVFIDs = gf_mesh_get(mesh, 'region', RIDs(end)); % get Neumann BC data   
   end
   
   points = [];
   for i = 1 : size(CVFIDs,2)
       points = [points, DOFs(IDx(CVFIDs(1,i)):IDx(CVFIDs(1,i)+1)-1)'];
   end
   
   % Computing number of branches and number of nodes from mesh data.
   n_branches = size(RIDs,2) - n_bc - 1; % number of branches
   total_nodes = 0;
   for i = 1 : n_branches
       cvf = gf_mesh_get(mesh, 'region', RIDs(i)); % get convexes in region i
       dof = gf_mesh_fem_get(mf, 'basic dof from cvid', cvf(1,:));
       dof = unique(dof);
       total_nodes = total_nodes + length(dof);
   end
   
   % Creating BC data structure
   node_idx = total_nodes + 1; % first real node index
   bc_check = cell(4, length(unique(DOFs)) - node_idx + 1);
   bc_check(1, :) = num2cell([node_idx:length(unique(DOFs))]);
   for i = 1 : size(bc_check, 2)
       bc_check{4, i} = PTs(:, bc_check{1, i});
   end
   
   % Checking Neumann-Kirchhoff conditions
   if (bc_label ~= "DIR") 
       neu_dofs = [];
       for i = 1 : size(neu_CVFIDs,2)
           temp = DOFs(IDx(neu_CVFIDs(1,i)):IDx(neu_CVFIDs(1,i)+1)-1);
           face = neu_CVFIDs(2,i);
           face = -face + 3; % for face numbering issues (changes 1 in 2 and viceversa)...
           neu_dofs = [neu_dofs, temp(face)];
       end
       neu_dofs = unique(neu_dofs); 

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
           bc_check{3,aux} = [bc_check{3,aux}, current];
       end

       for k = 1 : size(neu_dofs,2)
           current_dof = neu_dofs(k);
           aux = current_dof - node_idx + 1;
           value = 0;
           for l = 1 : size(bc_check{3,aux},2)
               current2 = bc_check{3,aux};
               current3 = current2(:,l);
               id = find(current3 >= node_idx);
               eval = solution(current3);
               if id == 1
                   value = value + (eval(2)-eval(1))/distance;
               elseif id == 2
                   value = value + (eval(1)-eval(2))/distance;
               else
                   disp('Convex is not one-dimensional');
                   return;
               end
           end
           bc_check{3,aux} = [];
           bc_check{3,aux} = value;
           bc_check{2,aux} = "NEU";
       end
   end
   
   % Checking Dirichlet boundary conditions
   if (bc_label ~= 'NEU')
       dir_dofs = [];
       for i = 1 : size(dir_CVFIDs,2)
           temp = DOFs(IDx(dir_CVFIDs(1,i)):IDx(dir_CVFIDs(1,i)+1)-1);
           face = dir_CVFIDs(2,i);
           face = -face + 3;
           dir_dofs = [dir_dofs, temp(face)];
       end
       dir_dofs = unique(dir_dofs); 

       for j = 1 : size(dir_dofs, 2)
           aux = dir_dofs(j) - node_idx + 1;
           bc_check{2, aux} = "DIR";
           bc_check{3, aux} = solution(dir_dofs(j));
       end
   end
   
   % Plotting results
   figure()
   hold on
   grid on
   dim = length(bc_check{4,end});
   if (dim == 2)
       warning('3D plot may be unreadable for complex graphs, use the returned data structure instead');
       for m = 1 : size(bc_check, 2)
           coords = bc_check{4,m};
           x = coords(1);
           y = coords(2);
           plot3(x, y, 0, 'ro');
           text(x,y,[bc_check{2,m}, ': ', bc_check{3,m}]);
       end
       title({'Boundary conditions at each node'; subtitle});
   else
       disp('Plot is available only for planar graphs (with 2 spatial dimensions)');    
   end
end





