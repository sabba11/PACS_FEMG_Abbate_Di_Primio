addpath('/usr/local/getfem_toolbox')
clear all
close all
clc
eigen_files = dir('*.U') ;    % list of .U files
N = length(eigen_files) ;   % total number of files 
% loop for each file 
 mf=gfMeshFem('load','solution.mf');
for i = 1:10
    thisfile = eigen_files(i).name ;
    U = load(thisfile);
    U = 10*U'
    
   
    figure()
    gf_plot(mf,U,'mesh','on', 'zplot', 'on');
end