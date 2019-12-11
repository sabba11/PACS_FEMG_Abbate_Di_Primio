% FEMG_eigen_postprocessing.m
% MATLAB routine to postprocess eigenvalue problem solutions.
% We strongly recommend to run this code section by section.
% 
% Requirements:
% (a) GetFEM++-MATLAB interface installed.
% (b) test_problem/eigen as working directory.

%% 0. Clearing workspace.
% Avoid clearing variables needed for further elaboration.

clearvars -except eigenvalue_vector legends % add variables to keep here...

%% 1. Setup. Setting paths.
% The lines indicated with (!) possibly requires user input.

close all
clc

path_to_files = 'laplacian/graphene/export/64 point-mesh/QZ'; % (!) path to eigenvectors/eigenvalues
path_to_data = '../../data_builder/data'; % (!) path to .txt data to build Laplace matrices

addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath('../matlab_functions') % adds custom MATLAB functions
addpath(path_to_files) % adds path(s) to .U and mesh files
addpath(path_to_data) % adds path to .txt data

%% 2. Loading eigenvalues and eigenvector files.
% Loads data in workspace.

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

%% 3. Eigenvalue plots
% Plots eigenvalues in ascending order.

figure()
plot([1:N], eigvals, 'r-o', 'LineWidth', 0.1);
title('Eigenvalues of the operator sorted by ascending order');
xlabel('Eigenvalue index');
ylabel('Eigenvalue');
grid on

%% 4. Eigenvector plots
% Plots the first n_plots eigenvectors.

n_plots = 5;
for i = 1 : n_plots
    val = eigvals(i);
    figure()
    gf_plot(mf, eigvects(:, i)', 'zplot', 'on', 'mesh', 'on', 'title', ['Eigenvector associated to \lambda = ', num2str(val)], 'disp_options', 'off');
    title(['Eigenvector associated to \lambda = ', num2str(val)]);
    %pbaspect([2 1 1]);
    xlabel('x');
    ylabel('y');
    zlabel('Eigenvector');
    grid on
    colorbar
end

%% 5. Computing eigenvector oscillations
% Computes the number of oscillations of all eigenvectors, plotting the
% trend.

FEMG_compute_oscillations(eigvects, mf, mesh);

%% 6. Boundary conditions check
% Checks boundary conditions for the i-th eigenvector. Includes a for cycle
% if multiple eigenvectors are to be checked.

n_bc = 1; % number of types of BC (1 if full Dirichlet, 1 if full Neumann, 2 if mixed)

idx = 2;
val = eigvals(idx);
subtitle = ['(eigenvector associated to \lambda = ', num2str(val), ')'];
bc_check = FEMG_check_BC(eigvects(:,idx), mf, mesh, n_bc, subtitle);

% for i = 1 : idx
%   figure()
%   val = eigvals(idx);
%   subtitle = ['Eigenvector associated to \lambda = ', num2str(val)];
%   bc_check = FEMG_check_BC(eigvects(:,idx), mf, mesh, n_bc, subtitle);
% end

%% 7. Eigenvalue comparison.
% Compares the eigenvalue trend of the Laplacian matrix, and that of the
% extended Laplacian matrix, across possibly different numerical
% experiments.
% This section is intended to be used as follows:
%   1. Every time the code in this section is run, eigenvalues in workspace
%      (they are assumed to be called eigvals) are saved for further
%      comparison in the vector eigenvalue_vector. A legend for plot
%      clarity is requested, to be assigned to the variable legend_ below.
%   2. The first element of eigenvalue_vector consists of the combinatorial
%      Laplace matrix eigenvalues, computed with the MATLAB eig() routine.
%   3. Set the variable plot_ to true once results are to be plotted.
% Consider comparing no more than 5 eigenvalue trends at a time. Otherwise,
% colors should be added in the color_vec structure (RGB notation [R, G, B] 
% can be used).

plot_ = false;
legend_ = 'test legend'; % (!) insert legend here

if exist('eigenvalue_vector','var') == 0
    laplace = FEMG_build_laplace_matrix(".txt"); % (!) specify which text file to read
    original_eig = eig(laplace);
    eigenvalue_vector = {original_eig};
    legends = 'Original eigenvalues';
end
eigenvalue_vector = [eigenvalue_vector, eigvals]; % do not clear this variable!
if (legend_ == '')
    legend_ = input('An empty legend is not allowed. Insert legend:\n');
end
legends = [legends, legend_]; % do not clear this variable!
if (plot_ == true)
    % plotting:
    color_vec = {'k', 'b', 'r', 'g', 'c', 'm'}; % (!) add colors if necessary
    figure()
    idxs = [1:size(original_eig)];
    hold on
    grid on
    for i = 1 : size(eigenvalue_vector)
        plot(idxs, eigenvalue_vector{i}(idxs), 'color', color_vec{i}, 'marker', 'o');
    end
    title('Eigenvalue behaviour comparison (ascending order)');
    legend(legends);
end

%% 8. Other postprocessing
% Space to define custom postprocessing routines.
