%--------------------------------------------------------------------------
% Computation of the internal in both 2D and 3D
%--------------------------------------------------------------------------
function T = internalF(X, x, n, l, mu)
    n_elem = size(n,1);
    n_nodes = size(X,1);
	dim = 3;
    gaussPoints  =  [-0.5774  -0.5774  -0.5774
                     -0.5774  -0.5774   0.5774
                     -0.5774   0.5774  -0.5774
                     -0.5774   0.5774   0.5774
                      0.5774  -0.5774  -0.5774
                      0.5774  -0.5774   0.5774
                      0.5774   0.5774  -0.5774
                      0.5774   0.5774   0.5774];
    gaussWeights = ones(8,1);
	T = zeros(n_nodes*dim, 1);
	for e = 1:n_elem
		xe   = x(n(e,:),:,:);
		Xe   = X(n(e,:),:);
        for gp = 1:size(gaussPoints,1)
	        z = gaussPoints(gp,:);
	        se = stress(xe, Xe, z, l , mu);
            [dNdx, J] = getdNdx(Xe, z);
            B = getB(dNdx);
            for a = 1:8
                sl_a_e = (dim*(n(e,a)-1)+1):(dim*n(e,a));
                T(sl_a_e) = T(sl_a_e) + B(:,:,a)'*se*J*gaussWeights(gp,:);
            end
        end
	end
end