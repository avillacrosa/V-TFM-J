function internalF(X, x, n, l, mu)
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
    T = zeros(n_node*dim, 1);

	for e in 1:n_elem
		xe   = x[n[e,:],:];
		Xe   = X[n[e,:],:];
        for gp in 1:8
	        z = gaussPoints[gp,:];
	        s_e = stress(xe, Xe, z, l, mu);
            dNdx, J = getdNdx(Xe, z);
            B = getB(dNdx);
            for a = 1:8
                sl_a_e = (dim*(n[e,a]-1)+1):(dim*n[e,a]);
                zz = J*gaussWeights[gp,:]
                T[sl_a_e] = T[sl_a_e] + B[:,:,a]'*s_e.*zz;
            end
        end
	end
    return T
end