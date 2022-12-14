%--------------------------------------------------------------------------
% Compute nodal displacement Neumann conditions given user input format
%--------------------------------------------------------------------------
function [fix, fix_vals] = BCtoNodal(X, BC)
	dim = 3;
	n_nodes = size(X,1);
    hits   = zeros(n_nodes, dim);
    fix_vals   = zeros(n_nodes, dim);
	for plane_idx = 1:size(BC, 1)
        axis    = BC(plane_idx,1);
        coord   = BC(plane_idx,2);
        t_axis  = BC(plane_idx,3);
        t_value = BC(plane_idx,4);
        hit_idx = X(:,axis) == coord;
        hits(hit_idx, t_axis) = 1;
        fix_vals(hit_idx, t_axis) = t_value;
	end
    hits = logical(hits);
	fix = hits;
end