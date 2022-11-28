function c = ctensor(l,mu)
    c = [l+2*mu		l		l		0	0	0	
    l	  l+2*mu	l		0	0	0
    l		l	  l+2*mu	0	0	0
    0		0		0	   mu	0	0	
    0		0		0		0   mu	0
    0		0		0		0	0	mu];
end