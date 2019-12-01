% This file is used to produce a radii.pts file such that all edges found
% in the .txt are ttreated as they have costant radii.
close all
clear all
% File infos
% out_filename = 'cradii_star.txt';
% in_filename = 'star.txt';
% out_filename = 'cradii_graphene.txt';
% in_filename = 'graphene.txt';
% out_filename = 'cradii_tree.txt';
% in_filename = 'tree_edge.txt';
out_filename = 'cradii_scale free.txt';
in_filename = 'scale_free.txt';
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
