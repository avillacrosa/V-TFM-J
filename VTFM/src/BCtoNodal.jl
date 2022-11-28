function BCtoNodal(X, uBC)
    n_nodes = size(X,1);
    dim = 3
    hits = zeros(n_nodes, dim);
    fix_vals = zeros(n_nodes, dim)
    for planeidx in 1:4
        axis    = trunc(Int, uBC[planeidx, 1])
        coord   = trunc(Int, uBC[planeidx, 2])
        t_axis  = trunc(Int, uBC[planeidx, 3])
        t_value = uBC[planeidx, 4]
        hit_idx = X[:, axis] .== coord;
        hits[hit_idx, t_axis] .= 1
        fix_vals[hit_idx,t_axis] .= t_value
    end
    hits = iszero.(hits)
    fix = hits
    return fix, fix_vals
end