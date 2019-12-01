% This file is used to produce a radii.pts file such that all edges found
% in the .txt have random radii with a maximum indication.
close all
clear all
% File infos
out_filename = 'star_rradii.txt';
in_filename = 'star.txt';
% Decide the maximum radius
max_radius = 0.1;
% Writing routine
fin = fopen(in_filename);
M = textscan(fin,'%s','delimiter','\n');
n = length(M{1});
 Filetext = "BEGIN_LIST";
for i=2:n+1
    Filetext(i,1) = num2str(rand*max_radius,'%f');
end
Filetext(n+2,1) = "END_LIST";
writematrix( Filetext,out_filename );
