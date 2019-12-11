% MATLAB routine to see differences between mesh uniform on edge vs on
% whole graph.
%
% Requirements:
% (a) GetFEM++-MATLAB interface installed.

%% 0. Clearing workspace.
clear all
% The line indicated with (!) possibly requires user input.
close all
clc

%% 1. Processing uniform mesh with step h

% path_to_files = 'laplacian/tree/export/93 point-mesh/QZ';
% path_to_files = 'laplacian/tree/export/973 point-mesh/QZ';
path_to_files = 'laplacian/scale_free/export/376 point-mesh/QZ';% (!) path to eigenvectors/eigenvalues
addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath(path_to_files) % adds path(s) to .U and mesh files

eigen_files_h = dir([path_to_files, '/*.U']);    % search for .U files
N_h = length(eigen_files_h)-1;   % total number of eigenvector files
eigvals_h = zeros(N_h, 1); % vector to store eigenvalues

% importing data
% eigenvalue file
thisfile = eigen_files_h(1).name;
U = load(thisfile);
U = U';
eigvals_h = U;

fclose('all');

%% 2. Processing mesh with N sub-edges for each edge and almost same number as before

% path_to_files = 'laplacian/tree/export/91 point-mesh/QZ';
% path_to_files = 'laplacian/tree/export/976 point-mesh/QZ';
path_to_files = 'laplacian/scale_free/export/368 point-mesh/QZ';% (!) path to eigenvectors/eigenvalues
addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath(path_to_files) % adds path(s) to .U and mesh files

eigen_files_N = dir([path_to_files, '/*.U']);    % search for .U files
N_N = length(eigen_files_N)-1;   % total number of eigenvector files
eigvals_N = zeros(N_N, 1); % vector to store eigenvalues

% importing data
% eigenvalue file
thisfile = eigen_files_N(1).name;
U = load(thisfile);
U = U';
eigvals_N = U;

fclose('all');

%% 3. Processing mesh now divided uniformly on whole graph with step h equals to minimum length of the step of the preciding mesh

% path_to_files = 'laplacian/tree/export/263 point-mesh/QZ';
% path_to_files = 'laplacian/tree/export/2849 point-mesh/QZ';
path_to_files = 'laplacian/scale_free/export/1629 point-mesh/QZ';% (!) path to eigenvectors/eigenvalues
addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath(path_to_files) % adds path(s) to .U and mesh files

eigen_files_h_new = dir([path_to_files, '/*.U']);    % search for .U files
N_h_new = length(eigen_files_h_new)-1;   % total number of eigenvector files
eigvals_h_new = zeros(N_h_new, 1); % vector to store eigenvalues

% importing data
% eigenvalue file
thisfile = eigen_files_h_new(1).name;
U = load(thisfile);
U = U';
eigvals_h_new = U;

fclose('all');

%% 4. Eigenvalue plots

% Plots eigenvalues in ascending order for each of the mesh and compare
% them

figure()
plot([1:N_h_new], eigvals_h_new, 'b-', 'LineWidth', 1);
title('Eigenvalues of the operator sorted by ascending order');
xlabel('Eigenvalue index');
ylabel('Eigenvalue');
grid on
hold on
plot([1:N_N], eigvals_N, 'g-', 'LineWidth', 1);
plot([1:N_h], eigvals_h, 'r-', 'LineWidth', 1);
legend('h-type mesh', 'N-type mesh: same dofs', 'h-type mesh: minimum step length','Location','northwest')
