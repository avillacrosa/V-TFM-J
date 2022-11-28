function assembleK(X, n)
    n_elem = size(n,1)
    dim = 3
    Kg1 = zeros(Int,8*3, n_elem);
    Kg2 = zeros(Int,8*3, n_elem);
    for e in 1:n_elem
        ne = n[e,:];
        Kg1_e = zeros(8*3, 1);
        Kg2_e = zeros(8*3, 1);
        for a in 1:8
            for b in 1:8
                sl_a_e = (dim*(a-1)+1):dim*a;
                sl_b_e = (dim*(b-1)+1):dim*b;

                sl_a = (dim*(ne[a]-1)+1):dim*ne[a];
                sl_b = (dim*(ne[b]-1)+1):dim*ne[b];

                Kg1_e[sl_a_e] = sl_a;
                Kg2_e[sl_b_e] = sl_b;
            end
        end
        Kg1[:, e] = Kg1_e;
        Kg2[:, e] = Kg2_e;
    end
    return Kg1, Kg2
end