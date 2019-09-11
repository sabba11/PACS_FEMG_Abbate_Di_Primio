addpath('/usr/local/getfem_toolbox')
clear all
close all
clc
eigen_files = dir('*.U') ;    % list of .U files
N = length(eigen_files) ;   % total number of files 
% loop for each file 
for i = 1 : 10
    thisfile = eigen_files(i).name ;
    U = load(thisfile);
    U = U';
    mf=gfMeshFem('load','solution.mf');
    figure()
    gf_plot(mf,U,'mesh','on', 'zplot', 'on');
end