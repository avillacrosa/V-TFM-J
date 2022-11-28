function getdNdx(x, z)
    z1 = z[1];
    z2 = z[2];
    z3 = z[3];

    N = [(1-z1)*(1-z2)*(1-z3)
         (1+z1)*(1-z2)*(1-z3)
         (1+z1)*(1+z2)*(1-z3)
         (1-z1)*(1+z2)*(1-z3)
         (1-z1)*(1-z2)*(1+z3)
         (1+z1)*(1-z2)*(1+z3)
         (1+z1)*(1+z2)*(1+z3)
         (1-z1)*(1+z2)*(1+z3)
        ]/8;

    dNdz = [ -(z2-1)*(z3-1) -(z1-1)*(z3-1) -(z1-1)*(z2-1)
              (z2-1)*(z3-1)  (z1+1)*(z3-1)  (z1+1)*(z2-1)
             -(z2+1)*(z3-1) -(z1+1)*(z3-1) -(z1+1)*(z2+1)
              (z2+1)*(z3-1)  (z1-1)*(z3-1)  (z1-1)*(z2+1)
              (z2-1)*(z3+1)  (z1-1)*(z3+1)  (z1-1)*(z2-1)
             -(z2-1)*(z3+1) -(z1+1)*(z3+1) -(z1+1)*(z2-1)
              (z2+1)*(z3+1)  (z1+1)*(z3+1)  (z1+1)*(z2+1)
             -(z2+1)*(z3+1) -(z1-1)*(z3+1) -(z1-1)*(z2+1)
            ]/8;
    dxdz = dNdz'*x;
    dNdx = (dxdz\dNdz')';
    J    = det(dxdz);
    return dNdx, J
end