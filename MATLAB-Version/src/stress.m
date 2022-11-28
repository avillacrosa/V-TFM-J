function sigma = stress(x,X,z,l,mu)
    Fd = deformF(x,X,z);
    D = [l+2*mu		l		l		0	0	0	
    l	  l+2*mu	l		0	0	0
    l		l	  l+2*mu	0	0	0
    0		0		0		mu	0	0	
    0		0		0		0	mu	0
    0		0		0		0	0	mu];
    eps = (Fd'+Fd)/2-eye(size(Fd));
    eps = [eps(1,1), eps(2,2), eps(3,3), 2*eps(1,2), 2*eps(1,3), 2*eps(2,3)];
    sigma = D*eps';
end