% compute_oscillations.m 
% Computes number of oscillations with respect to eigenvalue index.
% The function assumes eigenvectors are passed sorted by increasing
% corresponding eigenvalues.

function access_idxs = FEMG_compute_oscillations(eigvects, mf, mesh) 
    RIDs = gf_mesh_get(mesh, 'regions'); % get valid region indexes
    [DOFs, IDx] = gf_mesh_fem_get(mf, 'basic dof from cvid'); % get convex and dof enumeration
    
    n_branches = size(RIDs,2) - 1; % number of branches
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
    
    oscillations = zeros(1, size(eigvects, 1));
    % evaluate eigenvector on each branch
    for i = 1 : size(eigvects,1) % loop on each eigenvector
        eigvect = eigvects(:, i); % consider current eigenvector
        for j = 1 : size(access_idxs,2) % for each access vector
            access = access_idxs{j}; % consider current access vector
            eval = eigvect(access); % evaluate eigenvector on branch
            for k = 1 : size(eval)-1
                if eval(k)*eval(k+1) < 0 
                    oscillations(i) = oscillations(i) + 1;
                end
            end
        end
    end
        
    idxs = [1:size(eigvects, 1)]; % x-values for the final plot
    figure()
    plot(idxs, oscillations);
    title('Number of eigenvector oscillations');
    xlabel('Eigenvalue index');
    ylabel('Number of corresponding eigenvector oscillations');
    
    oscillations(end)
end
