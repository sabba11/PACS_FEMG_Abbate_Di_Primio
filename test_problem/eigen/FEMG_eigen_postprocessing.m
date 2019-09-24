% view_eigenvectors.m
% MATLAB routine to postprocess eigenvalue problem solutions. 
%
% Requirements:
% (a) GetFEM++-MATLAB interface installed.

%% 0. Clearing workspace.
% Avoid clearing variables if saving data for further elaboration.

clear all

%% 1. Setup. Setting paths.

close all
clc

path_to_files = 'star/export/Hamiltonian/401 point-mesh/QR'; % (!) path to eigenvectors/eigenvalues

addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath(path_to_files) % adds path(s) to .U and mesh files

%% 2. Loading eigenvalues and eigenvector files.

eigen_files = dir([path_to_files, '/*.U']);    % search for .U files
N = length(eigen_files)-1;   % total number of eigenvector files 
eigvals = zeros(N, 1); % vector to store eigenvalues
eigvects = zeros(N, N); % matrix to store eigenvectors
mf = gfMeshFem('load',[path_to_files, '/solution.mf']); % import mesh_fem object 
mesh = gfMesh('load',[path_to_files, '/mesh.mh']); % import mesh object

% eigenvalue file
thisfile = eigen_files(1).name;
U = load(thisfile);
U = U';
eigvals = U;
% loop for each eigenvector file 
for i = 2 : N+1
    thisfile = eigen_files(i).name;
    U = load(thisfile);
    U = U';
    eigvects(:, i-1) = U;
end
fclose('all');

%% 3. Computing eigenvector oscillations

compute_oscillations(eigvects, mf, mesh) ;

%% 4. Eigenvalue plots

figure()
plot([1:N], eigvals, 'r-o', 'LineWidth', 0.1);
title('Eigenvalues of the operator sorted by ascending order');
xlabel('Eigenvalue index');
ylabel('Eigenvalue');

%% 5. Eigenvector plots

n_plots = 10;
for i = 1 : n_plots
    figure()
    val = eigvals(i);
    title(['Eigenvector associated to \lambda = ', num2str(val)]);
    gf_plot(mf, eigvects(:, i)', 'zplot', 'on', 'mesh', 'on', 'title', ['Eigenvector associated to \lambda = ', num2str(val)]);
end

%% 6. Saving eigenvalue data for further comparison

% FEMG_eigval_comparison();

%% 7. 
