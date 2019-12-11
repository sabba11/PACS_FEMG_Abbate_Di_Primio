% MATLAB routine to study the behaviour of eigenvalues towards infinity
% First thing we want to exclude the spourious one cause the discretization
% does ruin the bigger one.
% After that we claim that the  trend is quadratic and we try to study the
% coefficient.
%
% Requirements:
% (a) GetFEM++-MATLAB interface installed.

%% 0. Clearing workspace.
clear all
% The line indicated with (!) possibly requires user input.
close all
clc

%% 1. Processing uniform mesh with step h

path_to_files = 'laplacian/tree/export/295 point-mesh/QZ';% (!) path to eigenvectors/eigenvalues
% path_to_files = 'laplacian/tree/export/452 point-mesh/QZ';
% path_to_files = 'laplacian/tree/export/241 point-mesh/QZ';
% path_to_files = 'laplacian/tree/export/366 point-mesh/QZ';
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

%% 2. Processing uniform mesh with step h/2

path_to_files = 'laplacian/tree/export/604 point-mesh/QZ';% (!) path to eigenvectors/eigenvalues
% path_to_files = 'laplacian/tree/export/906 point-mesh/QZ';
% path_to_files = 'laplacian/tree/export/494 point-mesh/QZ';
% path_to_files = 'laplacian/tree/export/742 point-mesh/QZ';
addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath(path_to_files) % adds path(s) to .U and mesh files

eigen_files_h2 = dir([path_to_files, '/*.U']);    % search for .U files
N_h2 = length(eigen_files_h2)-1;   % total number of eigenvector files 
eigvals_h2 = zeros(N_h2, 1); % vector to store eigenvalues

% importing data
% eigenvalue file
thisfile = eigen_files_h2(1).name;
U = load(thisfile);
U = U';
eigvals_h2 = U;

fclose('all');

%% 3. Processing uniform mesh with step h/4

path_to_files = 'laplacian/tree/export/1215 point-mesh/QZ';% (!) path to eigenvectors/eigenvalues
% path_to_files = 'laplacian/tree/export/1824 point-mesh/QZ';
% path_to_files = 'laplacian/tree/export/997 point-mesh/QZ';
% path_to_files = 'laplacian/tree/export/1490 point-mesh/QZ';
addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath(path_to_files) % adds path(s) to .U and mesh files

eigen_files_h4 = dir([path_to_files, '/*.U']);    % search for .U files
N_h4 = length(eigen_files_h4)-1;   % total number of eigenvector files 
eigvals_h4 = zeros(N_h4, 1); % vector to store eigenvalues

% importing data
% eigenvalue file
thisfile = eigen_files_h4(1).name;
U = load(thisfile);
U = U';
eigvals_h4 = U;

fclose('all');

%% 4. Processing uniform mesh with step h/8

path_to_files = 'laplacian/tree/export/2437 point-mesh/QZ';% (!) path to eigenvectors/eigenvalues
% path_to_files = 'laplacian/tree/export/3661 point-mesh/QZ';
% path_to_files = 'laplacian/tree/export/1993 point-mesh/QZ';
% path_to_files = 'laplacian/tree/export/2991 point-mesh/QZ';
addpath('/usr/local/getfem_toolbox') % adds gf_* functions
addpath(path_to_files) % adds path(s) to .U and mesh files

eigen_files_h8 = dir([path_to_files, '/*.U']);    % search for .U files
N_h8 = length(eigen_files_h8)-1;   % total number of eigenvector files 
eigvals_h8 = zeros(N_h8, 1); % vector to store eigenvalues

% importing data
% eigenvalue file
thisfile = eigen_files_h8(1).name;
U = load(thisfile);
U = U';
eigvals_h8 = U;

fclose('all');

%% 5. Eigenvalue plots

% Plots eigenvalues in ascending order for each of the mesh and compare
% them

figure()
plot([1:N_h], eigvals_h, 'r-', 'LineWidth', 0.5);
title('Eigenvalues of the operator sorted in ascending order');
xlabel('Eigenvalue index');
ylabel('Eigenvalue');
grid on
hold on
plot([1:N_h2], eigvals_h2, 'b-', 'LineWidth', 0.5);
plot([1:N_h4], eigvals_h4, 'g-', 'LineWidth', 0.5);
plot([1:N_h8], eigvals_h8, 'y-', 'LineWidth', 0.5);
legend('step h', 'step h/2', 'step h/4','step h/8','Location','northwest')
hold off


%% 6. First approach at studying the trend

%Loglog graph to observe the quadratic behaviour
figure()
loglog([1:N_h], abs(eigvals_h), 'r-', 'LineWidth', 0.5);
title('Logarithmic plot of the eigenvalues');
xlabel('Eigenvalue index');
ylabel('Eigenvalue');
grid on
hold on
loglog([1:N_h2], abs(eigvals_h2), 'b-', 'LineWidth', 0.5);
loglog([1:N_h4], abs(eigvals_h4), 'g-', 'LineWidth', 0.5);
loglog([1:N_h8], abs(eigvals_h8), 'y-', 'LineWidth', 0.5);
loglog([1:N_h8], [1:N_h8].^2,'k-', 'LineWidth', 0.5);
legend('step h', 'step h/2', 'step h/4','step h/8','N^2','Location','northwest')

p = polyfit([1:N_h8], abs(eigvals_h8(1:N_h8)), 2);
quadratic_vec = polyval(p,[1:N_h8]);
figure()
plot([1:N_h8], abs(eigvals_h8(1:N_h8)), 'k-', 'LineWidth', 0.5);
title('Plot of eigenvalues vs their quadratic estimate');
grid on
hold on
plot([1:N_h8], quadratic_vec, 'r-', 'LineWidth', 0.5);
legend('Eigenvalues','Quadratic function','Location','northwest')

%% 7. excluding last eigenvalues

for i = 2:N_h
    if eigvals_h(i)*1.10<eigvals_h2(i) || eigvals_h(i)*0.90>eigvals_h2(i)
        index(1) = i;
        errors(1) = i/N_h;
        break
    end
end
for i = 2:N_h2
    if eigvals_h2(i)*1.10<eigvals_h4(i) || eigvals_h2(i)*0.90>eigvals_h4(i)
        index(2) = i;
        errors(2) = i/N_h2;
        break
    end
end
for i = 2:N_h4
    if eigvals_h4(i)*1.10<eigvals_h8(i) || eigvals_h4(i)*0.90>eigvals_h8(i)
        index(3) = i;
        errors(3) = i/N_h4;
        break
    end
end


figure()
plot([1:floor(min(errors)*N_h)], eigvals_h(1:floor(min(errors)*N_h)), 'r-', 'LineWidth', 0.5);
title('Trimmed (first part) of the eigenvalues');
xlabel('Eigenvalue index');
ylabel('Eigenvalue');
grid on
hold on
plot([1:floor(min(errors)*N_h2)], eigvals_h2(1:floor(min(errors)*N_h2)), 'b-', 'LineWidth', 0.5);
plot([1:floor(min(errors)*N_h4)], eigvals_h4(1:floor(min(errors)*N_h4)), 'g-', 'LineWidth', 0.5);
plot([1:floor(min(errors)*N_h8)], eigvals_h8(1:floor(min(errors)*N_h8)), 'y-', 'LineWidth', 0.5);
legend('step h', 'step h/2', 'step h/4','step h/8','Location','northwest')
hold off

p_new = polyfit([1:floor(min(errors)*N_h8)], abs(eigvals_h8(1:floor(min(errors)*N_h8))), 2);
quadratic_vec_new = polyval(p,[1:floor(min(errors)*N_h8)]);
figure()
plot([1:floor(min(errors)*N_h8)], abs(eigvals_h8(1:floor(min(errors)*N_h8))), 'k-', 'LineWidth', 0.1);
title('Quadratic estimate for the trimmed part');
grid on
hold on
plot([1:floor(min(errors)*N_h8)], quadratic_vec_new, 'r-', 'LineWidth', 0.5);
legend('Eigenvalues','Quadratic function','Location','northwest')

p_last = polyfit([floor(min(errors)*N_h8/3):floor(min(errors)*N_h8)], abs(eigvals_h8(floor(min(errors)*N_h8/3):floor(min(errors)*N_h8))), 2);
quadratic_vec_last = polyval(p,[floor(min(errors)*N_h8/3):floor(min(errors)*N_h8)]);
figure()
plot([floor(min(errors)*N_h8/3):floor(min(errors)*N_h8)], abs(eigvals_h8(floor(min(errors)*N_h8/3):floor(min(errors)*N_h8))), 'k-', 'LineWidth', 0.1);
title('Quadratic estimate for last eigenvalues of the trimmed part');
grid on
hold on
plot([floor(min(errors)*N_h8/3):floor(min(errors)*N_h8)], quadratic_vec_last, 'r-', 'LineWidth', 0.5);
legend('Eigenvalues','Quadratic function','Location','northwest')

