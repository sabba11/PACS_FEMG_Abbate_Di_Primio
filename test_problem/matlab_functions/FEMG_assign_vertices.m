% FEMG_assign_vertices.m
% Assigns real vertices to regions of graph nodes.
%
% OUTPUT:
% access_idxs  a cell array containing vectors of indices, one for each
%              branch of the graph. Indices are put in such a way that
%              eigvect(access_idxs{i}) computes the values of the
%              eigenvector eigvect on the i-th branch (from source to
%              target vertex).
% INPUT:
% mf           a GetFEM++ MeshFem object. It represents the mesh on which
%              finite elements are built.
%
% mesh         a GetFEM++ Mesh object. It represents the mesh on which the
%              the eigenvectors are computed. It may possibly be equal
%              to mf.

function access_idxs = FEMG_assign_vertices(mf, mesh)
    RIDs = gf_mesh_get(mesh, 'regions'); % get valid region indexes
    [DOFs, IDx] = gf_mesh_fem_get(mf, 'basic dof from cvid'); % get convex and dof enumeration

    n_branches = size(RIDs,2) - 3; % number of branches
    branch_idx = size(IDx,2) - 2*n_branches; % auxiliary counter

    access_idxs = {}; % cell array to store idxs to access eigenvector values on a branch

    % loop over all regions except the last one (i.e. on all branches)
    for i = 1 : n_branches
        CVFIDs = gf_mesh_get(mesh, 'region', RIDs(i)); % get convexes in region RIDs(id)
        vec = zeros(1, 2*(size(CVFIDs,2) + 2)); % allocate space for access vector
        vec(1:2) = DOFs(IDx(branch_idx):IDx(branch_idx + 1)-1); % get source node dofs
        for j = 1 : size(CVFIDs,2)
            vec(2*j+1:2*j+2) = DOFs(IDx(CVFIDs(1,j)):IDx(CVFIDs(1,j) + 1)-1); % get mesh nodes dofs
        end
        branch_idx = branch_idx + 1; % update counter
        vec(end-1:end) = DOFs(IDx(branch_idx):IDx(branch_idx + 1)-1); % get target node dofs
        vec([end-1 end]) = vec([end end-1]); % correct ordering
        vec = unique(vec, 'stable'); % remove duplicates
        access_idxs = [access_idxs, vec]; % push vector in cell array
        branch_idx = branch_idx + 1; % update counter
    end
end
