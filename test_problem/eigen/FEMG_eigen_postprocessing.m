% view_eigenvectors.m
% MATLAB routine to postprocess eigenvalue problem solutions. 
%
% Requirements:
% (a) GetFEM++-MATLAB interface installed.

%% 0. Clearing workspace.
% Avoid clearing variables if saving data for further elaboration.

clear all

%% 1. Setup. Setting paths.
% The line indicated with (!) possibly requires user input.

close all
clc


path_to_files = 'star/export/1460 point-mesh/QZ'; % (!) path to eigenvectors/eigenvalues


addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath(path_to_files) % adds path(s) to .U and mesh files

%% 2. Loading eigenvalues and eigenvector files.
3
eigen_files = dir([path_to_files, '/*.U']);    % search for .U files
N = length(eigen_files)-1;   % total number of eigenvector files 
eigvals = zeros(N, 1); % vector to store eigenvalues
eigvects = zeros(N, N); % matrix to store eigenvectors
mf = gfMeshFem('load',[path_to_files, '/solution.mf']); % import mesh_fem object 
mesh = gfMesh('load',[path_to_files, '/mesh.mh']); % import mesh object

% importing data
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

FEMG_compute_oscillations(eigvects, mf, mesh) ;

%% 4. Eigenvalue plots

% Plots eigenvalues in ascending order.


figure()
plot([1:N], eigvals, 'r-o', 'LineWidth', 0.1);
title('Eigenvalues of the operator sorted by ascending order');
xlabel('Eigenvalue index');
ylabel('Eigenvalue');
grid on

%% 4. Eigenvector plots
% Plots the first n_plots eigenvectors.

n_plots = 10;
for i = 1 : n_plots
    val = eigvals(i);

    figure()
    gf_plot(mf, 10*eigvects(:, i)', 'zplot', 'on', 'mesh', 'on', 'title', ['Eigenvector associated to \lambda = ', num2str(val)], 'disp_options', 'off');
    title(['Eigenvector associated to \lambda = ', num2str(val)]);
    pbaspect([2 1 1]);
    xlabel('x');
    ylabel('y');
    zlabel('Eigenvector');
    grid on

    colorbar
end

%% 5. Computing eigenvector oscillations
% Computes the number of oscillations of all eigenvectors, plotting the 

FEMG_compute_oscillations(eigvects, mf, mesh) ;

%% 6. Saving eigenvalue data for further comparison
% Saves eigenvalues to avoid variable overwriting upon postprocessing of
% different numerical results.

% FEMG_eigval_comparison();

%% 7. Boundary conditions check
% Checks boundary conditions for the i-th eigenvector. Includes a for cycle
% if multiple eigenvectors are to be checked.

idx = 2;
FEMG_check_BC(eigvects(:,idx), mf, mesh);

% for i = 1 : idx
%   figure()
%   
% end

%% 8. Other postprocessing
% Space to define custom postprocessing routines.

% ...