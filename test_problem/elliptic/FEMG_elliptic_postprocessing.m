% FEMG_elliptic_postprocessing.m
% MATLAB routine to postprocess results of an elliptic differential
% problem.
% We strongly recommend to run this script section by section.
%
% Requirements:
% (a) GetFEM++-MATLAB interface installed.
% (b) test_problem/elliptic as working directory.

%% 0. Clearing workspace.
% Avoid clearing variables needed for further elaboration.

clearvars -except % add variables here...

%% 1. Setup. Setting paths.
% The lines indicated with (!) may possibly require user input.

close all
clc

path_to_files = 'laplacian/star/export/21 point-mesh/QMR'; % (!) path to solution files

addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath('../matlab_functions') % adds custom MATLAB functions
addpath(path_to_files) % adds path(s) to solution and mesh files

%% 2. Loading solution data.

solution_file = dir([path_to_files, '/*.U']);    % search for .U files (there should be only one)
mf = gfMeshFem('load',[path_to_files, '/solution.mf']); % import mesh_fem object 
mesh = gfMesh('load',[path_to_files, '/mesh.mh']); % import mesh object

thisfile = solution_file(1).name;
U = load(thisfile);
U = U';
solution = U;
fclose('all');

%% 3. Solution plot.
% Plots the solution over the graph.

figure()
gf_plot(mf, solution, 'zplot', 'on', 'mesh', 'on', 'title', 'Approximated solution of the differential problem', 'disp_options', 'off');
title('Approximated solution of the differential problem');
xlabel('x');
ylabel('y');
zlabel('Solution');
grid on
colorbar

%% 4. Boundary conditions check.
% Checks Neumann-Kirchhoff and Dirichlet boundary conditions.

n_bc = 1; % number of types of BC (1 if full Dirichlet, 1 if full Neumann, 2 if mixed)
bc_check = FEMG_check_BC(solution, mf, mesh, n_bc);

%% 5. Other postprocessing.
% Space to define custom postprocessing routines.
