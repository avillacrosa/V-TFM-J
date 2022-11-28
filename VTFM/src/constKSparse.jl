using SparseArrays
function constKSparse(X, x, n, l, mu, Kg1, Kg2)
    gaussPoints  =  [-0.5774  -0.5774  -0.5774
    -0.5774  -0.5774   0.5774
    -0.5774   0.5774  -0.5774
    -0.5774   0.5774   0.5774
     0.5774  -0.5774  -0.5774
     0.5774  -0.5774   0.5774
     0.5774   0.5774  -0.5774
     0.5774   0.5774   0.5774]
    gaussWeights = ones(Int32, 8,1)
    dim = 3
    n_elem = size(n,1)
    n_node = size(X,1)
    c = 1;
    ll = 8*3;
    K_id1 = zeros(Int64,ll^2*n_elem);
    K_id2 = zeros(Int64,ll^2*n_elem);
    K = zeros(Float64, ll^2*n_elem);
    for e in 1:n_elem
        xe   = x[n[e,:],:];
		Xe   = X[n[e,:],:];
        Ke  = zeros(8*3, 8*3);
        for gp in 1:8
            z = gaussPoints[gp,:]
            D = ctensor(l, mu)
            dNdx, J = getdNdx(Xe,z)
            B = getB(dNdx)
            for a = 1:8
                for b = 1:8
                    sl_a_e = (dim*(a-1)+1):(dim*a);
                    sl_b_e = (dim*(b-1)+1):(dim*b);
                    zz = J*gaussWeights[gp,:]
                    zz2 = (B[:,:,a]'*D)*B[:,:,b].*zz
                    Ke[sl_a_e, sl_b_e] = Ke[sl_a_e, sl_b_e] + zz2;
                end
            end
        end
        for aa = 1:8*3
            for bb = 1:8*3
                K_id1[c] = Kg1[aa,e];
                K_id2[c] = Kg2[bb,e];
                K[c] = Ke[aa,bb];
                c = c+1;
            end
        end
    end
    K = sparse(K_id1, K_id2, K);
    return K
end