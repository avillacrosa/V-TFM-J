function runFEM(ns, ds, E, nu, uBC)
    l  = E*nu/((1+nu)*(1-2*nu));
    mu = E/(2*(1+nu));
    
    x, n = meshgen(ns, ds)
    
    X = copy(x)
    Kg1, Kg2 = assembleK(X,n)
    dof, fix_vals = BCtoNodal(X, uBC)
    dofv = dof'[:];
    n_node = size(X,1)

    x = x + fix_vals
    xv = x'[:]
    T = internalF(X, x, n, l, mu);
    # Newton Algorithm
    it = 1
    tol = norm(T[dofv]);
    while tol > 1e-10
        # K = constK(X, x, n, l, mu)
        K = constKSparse(X, x, n, l, mu, Kg1, Kg2)
        Kdf = K[dofv, dofv]
        Tdf = -T[dofv]
        du = Kdf\Tdf
        xv[dofv] = xv[dofv] + du
        x = reshape(xv, 3, n_node)'
        T = internalF(X, x, n, l, mu)
        T = sparse(T)
        tol = norm(T[dofv])
        tolX = norm(du)
        println("ITER = $it, tolR = $tol, tolX = $tolX\n")
        it = it + 1
    end
    return x
end