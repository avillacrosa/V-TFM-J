function meshgen(ns, ds)
    nx = ns[1];
    ny = ns[2];
    nz = ns[3];
    
    dx = ds[1];
    dy = ds[2];
    dz = ds[3];


    dim = 3;
    nelem = (ns[1]-1)*(ns[2]-1)*(ns[3]-1);
    
    x  = zeros(Float32, nx*ny*nz, 3);
    n  = zeros(Int32, nelem, 2^dim);
    for nzi in 1:nz
        for nyi in 1:ny
            for nxi in 1:nx
                xc = (nxi-1)*dx;
                yc = (nyi-1)*dy;
				zc = (nzi-1)*dz;
                n_idx = (nzi-1)*(ny*nx) + nx*(nyi-1) + nxi;
                x[n_idx, :] = [xc, yc, zc];
                bl = n_idx;
                cond = nxi != nx && nyi != ny && nzi != nz;
                if cond
                    e_idx = n_idx + (1-nyi) + (nx+ny-1)*(1-nzi);
                    cs  = [bl, bl + 1, bl + nx + 1, bl + nx];
                    cs = vcat(cs, cs .+ nx*ny);
                    n[e_idx, :] = cs;
                end
            end
        end
    end    
    return x, n
end