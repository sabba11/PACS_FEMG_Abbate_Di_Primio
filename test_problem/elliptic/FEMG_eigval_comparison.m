% FEMG_eigval_comparison.m
% This script compares the eigenvalue trend of the Laplacian matrix, the
% extended Laplacian matrix and the Hamiltonian matrix.

%% 1. Get eigenvalue data.
% This section should be run to keep eigenvalue data across different
% eigenvalue problems. Eigenvalue data should already be in the workspace
% (eigvals variable in FEMG_eigen_postprocessing.m).

%laplacian_eigvals = ;
%extd_laplacian_eigvals = eigvals; % run one of these lines at a time, comment the other
hamiltonian_eigvals = eigvals; % run one of these lines at a time, comment the other

%% 2. Plot eigenvalue data.

dim = 5; % this should be the number of vertices in the original graph
figure()
idxs = [1:dim];
plot(idxs, extd_laplacian_eigvals(1:dim), 'k-o');
hold on
plot(idxs, hamiltonian_eigvals(1:dim), 'r-o');
plot(idxs, [0,1,1,1,5], 'b-o');
title('Eigenvalue behaviour comparison (ascending order)');
legend('Extended Laplacian', 'Hamiltonian', 'Original Laplacian');


