% This file is used to produce a radii.pts file such that all edges found
% in the .txt are ttreated as they have costant radii.
close all
clear all
% File infos
out_filename = 'txt_files/cradii_star.txt';
in_filename = 'txt_files/star.txt';
% out_filename = 'txt_files/cradii_graphene.txt';
% in_filename = 'txt_files/graphene.txt';
% out_filename = 'txt_files/cradii_tree.txt';
% in_filename = 'txt_files/tree_angle.txt';
% out_filename = 'txt_files/cradii_scale_free.txt';
% in_filename = 'txt_files/scale_free_20.txt';
% Decide the radius
radius = 0.1;
% Writing routine
fin = fopen(in_filename);
M = textscan(fin,'%s','delimiter','\n');
n = length(M{1});
 Filetext = "BEGIN_LIST";
for i=2:n+1
    Filetext(i,1) = num2str(radius,'%f');
end
Filetext(n+2,1) = "END_LIST";
writematrix( Filetext,out_filename );
