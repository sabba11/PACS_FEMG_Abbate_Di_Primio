% FEMG_augment_dim

% This function takes as an  imput the matrix containing the column of
% our 2d points and it adds the third coordinate as a zero.
% It takes as inputs:
% points2d = matrix 2*n containing the coordiantes of the points


function [points3d] = FEMG_augment_dim(points2d)
[r,c] = size(points2d);
points3d = [points2d, zeros(r,1)];
end
