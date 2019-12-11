% FEMG_build_laplace_matrix.m
% Computes the Laplace matrix of a graph from a .txt file.
%
% OUTPUT:
% laplace   Laplace matrix of the graph described in a .txt file. The .txt
%           file is assumed to contain a list of arcs in the form
%                     x1 y1 z1 x2 y2 z2 bc1 val1 bc2 val2
%           where
%           (x1,y1,z1) are the coordinates of the source vertex,
%           (x2,y2,z2) are the coordinates of the target vertex,
%           (bc1, val1) are label and value of BC at the source vertex,
%           (bc2, val2) are label and value of BC at the target vertex.
%           z-coordinates may also be not given.
% INPUT:
% filename  name of the .txt file containing edge data.

function laplace = FEMG_build_laplace_matrix(filename)
	points = []; % counter for points
	arcs = 0; % counter for arcs
    connectivity_data = [1,2]; % storing connectivity info
	fid = fopen(filename);
    tline = fgetl(fid); % read first line
    n_coords = 0;
    sizes = [0, 8, 10]; % auxiliary vector to check .txt is valid
    first = true;
    while(ischar(tline))
        line = regexprep(tline,' +',' '); % removing extra whitespace
        line = strtrim(line); % removing leading whitespace
        line = split(line); % separating line data
        if (first == true) % set correct number of coordinates
            if (size(line,1) == 8)
                n_coords = 2;
            elseif (size(line,1) == 10)
                n_coords = 3;
            else
                disp("Points should have 2 or 3 coordinates.")
                fclose(fid);
                return;
            end
            source = [];
            target = [];
            for i = 1 : n_coords
                source = [source, str2num(line{i})];
                target = [target, str2num(line{i+n_coords})];
            end
            points = [points; source];
            points = [points; target];
            first = false;
        else
            if (size(line,1) ~= sizes(n_coords)) % check points have right number of coordinates
                disp(["Points should have the same number of coordinates (]", n_coords, ")."]);
                fclose(fid);
                return;
            else
                source = [];
                target = [];
                for i = 1 : n_coords
                    source = [source, str2num(line{i})];
                    target = [target, str2num(line{i+n_coords})];
                end
                [sres, sidx] = ismember(source, points, 'rows');
                [tres, tidx] = ismember(target, points, 'rows');
                if (sres == 0) % check if source vertex has already been inserted
                    points = [points; source];
                    connectivity_data = [connectivity_data, size(points, 1)];
                else
                    connectivity_data = [connectivity_data, sidx];
                end
                if (tres == 0) % check if target vertex has already been inserted
                    points = [points; target];
                    connectivity_data = [connectivity_data, size(points, 1)];
                else
                    connectivity_data = [connectivity_data, tidx];
                end
            end
        end
        arcs = arcs + 1;
        tline = fgetl(fid);
    end
    fclose(fid);
	adj = zeros(size(points,1), arcs);
    column_counter = 1;
    % fill adjacency matrix
	for i = 1 : 2 : length(connectivity_data)-1
        adj(connectivity_data(i), column_counter) = -1;
        adj(connectivity_data(i+1), column_counter) = -1;
        column_counter = column_counter + 1;
    end
    % compute laplace matrix
	laplace = zeros(size(points,1), size(points,1));
	laplace = adj * adj';
end
