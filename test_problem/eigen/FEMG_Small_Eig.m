
%% 0. Clearing workspace.
clear all
% The line indicated with (!) possibly requires user input.
close all
clc

%% 1. Processing uniform mesh with step h

path_to_files = 'laplacian/tree/export/angle/QZ';% (!) path to eigenvectors/eigenvalues
addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath(path_to_files) % adds path(s) to .U and mesh files

eigen_files_1 = dir([path_to_files, '/*.U']);    % search for .U files
N_1 = length(eigen_files_1)-1;   % total number of eigenvector files 
eigvals_1 = zeros(N_1, 1); % vector to store eigenvalues

% importing data
% eigenvalue file
thisfile = eigen_files_1(1).name;
U = load(thisfile);
U = U';
eigvals_1 = U;

fclose('all');

path_to_files = 'laplacian/tree/export/angle_edge/QZ';% (!) path to eigenvectors/eigenvalues
addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath(path_to_files) % adds path(s) to .U and mesh files

eigen_files_2 = dir([path_to_files, '/*.U']);    % search for .U files
N_2 = length(eigen_files_2)-1;   % total number of eigenvector files 
eigvals_2 = zeros(N_2, 1); % vector to store eigenvalues

% importing data
% eigenvalue file
thisfile = eigen_files_2(1).name;
U = load(thisfile);
U = U';
eigvals_2 = U;

fclose('all');

path_to_files = 'laplacian/tree/export/angle_edge_conv/QZ';% (!) path to eigenvectors/eigenvalues
addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath(path_to_files) % adds path(s) to .U and mesh files

eigen_files_3 = dir([path_to_files, '/*.U']);    % search for .U files
N_3 = length(eigen_files_3)-1;   % total number of eigenvector files 
eigvals_3 = zeros(N_3, 1); % vector to store eigenvalues

% importing data
% eigenvalue file
thisfile = eigen_files_3(1).name;
U = load(thisfile);
U = U';
eigvals_3 = U;

fclose('all');

path_to_files = 'laplacian/tree/export/pi4/QZ';% (!) path to eigenvectors/eigenvalues
addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath(path_to_files) % adds path(s) to .U and mesh files

eigen_files_4 = dir([path_to_files, '/*.U']);    % search for .U files
N_4 = length(eigen_files_4)-1;   % total number of eigenvector files 
eigvals_4 = zeros(N_4, 1); % vector to store eigenvalues

% importing data
% eigenvalue file
thisfile = eigen_files_4(1).name;
U = load(thisfile);
U = U';
eigvals_4 = U;

fclose('all');

path_to_files = 'laplacian/tree/export/pi6/QZ';% (!) path to eigenvectors/eigenvalues
addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath(path_to_files) % adds path(s) to .U and mesh files

eigen_files_5 = dir([path_to_files, '/*.U']);    % search for .U files
N_5 = length(eigen_files_5)-1;   % total number of eigenvector files
eigvals_5 = zeros(N_5, 1); % vector to store eigenvalues

% importing data
% eigenvalue file
thisfile = eigen_files_5(1).name;
U = load(thisfile);
U = U';
eigvals_5 = U;

fclose('all');

path_to_files = 'laplacian/tree/export/pi4_conv/QZ';% (!) path to eigenvectors/eigenvalues
addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath(path_to_files) % adds path(s) to .U and mesh files

eigen_files_6 = dir([path_to_files, '/*.U']);    % search for .U files
N_6 = length(eigen_files_6)-1;   % total number of eigenvector files
eigvals_6 = zeros(N_1, 1); % vector to store eigenvalues

% importing data
% eigenvalue file
thisfile = eigen_files_6(1).name;
U = load(thisfile);
U = U';
eigvals_6 = U;

fclose('all');

path_to_files = 'laplacian/tree/export/pi6_conv/QZ';% (!) path to eigenvectors/eigenvalues
addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath(path_to_files) % adds path(s) to .U and mesh files

eigen_files_7 = dir([path_to_files, '/*.U']);    % search for .U files
N_7 = length(eigen_files_7)-1;   % total number of eigenvector files
eigvals_7 = zeros(N_7, 1); % vector to store eigenvalues

% importing data
% eigenvalue file
thisfile = eigen_files_7(1).name;
U = load(thisfile);
U = U';
eigvals_7 = U;

fclose('all');

figure()
plot([1:15], eigvals_1(1:15), '-o', 'LineWidth', 0.1);
title('First eigenvalues multiplicty comparison');
xlabel('Eigenvalue index');
ylabel('Eigenvalue');
grid on
hold on
plot([1:15], eigvals_2(1:15), '-o', 'LineWidth', 0.1);
plot([1:15], eigvals_3(1:15), '-o', 'LineWidth', 0.1);
plot([1:15], eigvals_4(1:15), '-o', 'LineWidth', 0.1);
plot([1:15], eigvals_5(1:15), '-o', 'LineWidth', 0.1);
plot([1:15], eigvals_6(1:15), '-o', 'LineWidth', 0.1);
plot([1:15], eigvals_7(1:15), '-o', 'LineWidth', 0.1);
legend('const edge',...
    'changing both angle and edge', ...
    'changing both angle and edge,graph length converges',...
    'theta = pi/4', ...
    'theta = pi/6', ...
    'theta = pi/4, graph length converges',...
    'theta = pi/6, graph length converges',...
    'Location', 'northwest')
hold off

figure()
plot([15:30], eigvals_1(15:30), '-o', 'LineWidth', 0.1);
title('First eigenvalues multiplicty comparison-2');
xlabel('Eigenvalue index');
ylabel('Eigenvalue');
grid on
hold on
plot([15:30], eigvals_2(15:30), '-o', 'LineWidth', 0.1);
plot([15:30], eigvals_3(15:30), '-o', 'LineWidth', 0.1);
plot([15:30], eigvals_4(15:30), '-o', 'LineWidth', 0.1);
plot([15:30], eigvals_5(15:30), '-o', 'LineWidth', 0.1);
plot([15:30], eigvals_6(15:30), '-o', 'LineWidth', 0.1);
plot([15:30], eigvals_7(15:30), '-o', 'LineWidth', 0.1);
legend('const edge',...
    'changing both angle and edge', ...
    'changing both angle and edge,graph length converges',...
    'theta = pi/4', ...
    'theta = pi/6', ...
    'theta = pi/4, graph length converges',...
    'theta = pi/6, graph length converges',...
    'Location', 'northwest')
hold off

figure()
plot([30:45], eigvals_1(30:45), '-o', 'LineWidth', 0.1);
title('First eigenvalues multiplicty comparison-3');
xlabel('Eigenvalue index');
ylabel('Eigenvalue');
grid on
hold on
plot([30:45], eigvals_2(30:45), '-o', 'LineWidth', 0.1);
plot([30:45], eigvals_3(30:45), '-o', 'LineWidth', 0.1);
plot([30:45], eigvals_4(30:45), '-o', 'LineWidth', 0.1);
plot([30:45], eigvals_5(30:45), '-o', 'LineWidth', 0.1);
plot([30:45], eigvals_6(30:45), '-o', 'LineWidth', 0.1);
plot([30:45], eigvals_7(30:45), '-o', 'LineWidth', 0.1);
legend('const edge',...
    'changing both angle and edge', ...
    'changing both angle and edge,graph length converges',...
    'theta = pi/4', ...
    'theta = pi/6', ...
    'theta = pi/4, graph length converges',...
    'theta = pi/6, graph length converges',...
    'Location', 'northwest')
hold off


