function K = constKSparse(X,x,n,l,mu,Kg1,Kg2)
	n_nodes_elem = 8;
	dim = 3;
	n_elem = size(n,1);
    ll = n_nodes_elem*dim;
    K_id1 = zeros(ll^2*n_elem,1);
    K_id2 = zeros(ll^2*n_elem,1);
    K = zeros(ll^2*n_elem,1);
    c = 1;
    gaussPoints  =  [-0.5774  -0.5774  -0.5774
                 -0.5774  -0.5774   0.5774
                 -0.5774   0.5774  -0.5774
                 -0.5774   0.5774   0.5774
                  0.5774  -0.5774  -0.5774
                  0.5774  -0.5774   0.5774
                  0.5774   0.5774  -0.5774
                  0.5774   0.5774   0.5774];
    gaussWeights = ones(8,1);
    for e = 1:n_elem
        ne = n(e,:);
        xe = x(ne,:,:);
        Xe = X(ne,:,:);
    	Ke  = zeros(n_nodes_elem*dim, n_nodes_elem*dim);
    	for gp = 1:size(gaussPoints,1)
        	z = gaussPoints(gp,:);
			D = ctensor(l,mu);
	
        	[dNdx, J] = getdNdx(Xe, z);
        	B = getB(dNdx);
        	for a = 1:n_nodes_elem
            	for b = 1:n_nodes_elem
                	sl_a_e = (dim*(a-1)+1):dim*a;
                	sl_b_e = (dim*(b-1)+1):dim*b;
	
                	Ke(sl_a_e, sl_b_e) = Ke(sl_a_e, sl_b_e)+ ...
                       	(B(:,:,a)'*D*B(:,:,b))*gaussWeights(gp,:)*J;
            	end
        	end
    	end
		
		for aa = 1:size(Ke,1)
        	for bb = 1:size(Ke,2)
            	K_id1(c) = Kg1(aa,e);
            	K_id2(c) = Kg2(bb,e);
            	K(c) = Ke(aa,bb);
            	c = c+1;
        	end
    	end
	end
    K = sparse(K_id1, K_id2, K);
end