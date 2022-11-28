function constK(X, x, n, l, mu)
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
    K = zeros(n_node*dim, n_node*dim);
    for e in 1:n_elem
        xe   = x[n[e,:],:];
		Xe   = X[n[e,:],:];
        for gp in 1:8
            z = gaussPoints[gp,:]
            D = ctensor(l, mu)
            dNdx, J = getdNdx(Xe,z)
            B = getB(dNdx)
            for a = 1:8
                for b = 1:8
                    sl_a_e = (dim*(n[e,a]-1)+1):(dim*n[e,a]);
                    sl_b_e = (dim*(n[e,b]-1)+1):(dim*n[e,b]);
                    zz = J*gaussWeights[gp,:]
                    zz2 = (B[:,:,a]'*D)*B[:,:,b].*zz
                    K[sl_a_e, sl_b_e] = K[sl_a_e, sl_b_e] + zz2;
                end
            end
        end
    end
    return K
end