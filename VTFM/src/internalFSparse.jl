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
    c = 1;
    ll = 8*3;
    T_id1 = zeros(Int32,ll^2*n_elem);
    T = zeros(Float32, ll^2*n_elem);

	for e in 1:n_elem
		xe   = x[n[e,:],:];
		Xe   = X[n[e,:],:];
        Te  = zeros(8*3, 8*3);
        for gp in 1:8
	        z = gaussPoints[gp,:];
	        s_e = stress(xe, Xe, z, l, mu);
            dNdx, J = getdNdx(Xe, z);
            B = getB(dNdx);
            for a = 1:8
                sl_a_e = (dim*(a-1)+1):(dim*a);
                zz = J*gaussWeights[gp,:]
                Te[sl_a_e] = Te[sl_a_e] + B[:,:,a]'*s_e.*zz;
            end
        end
        for aa = 1:8*3
            T_id1[c] = Kg1[aa,e];
            K[c] = Ke[aa,bb];
            c = c+1;
        end
	end
    return T
end