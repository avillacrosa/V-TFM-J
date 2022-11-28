function getB(dNdx)
    n = size(dNdx,1);
    d = size(dNdx,2);
    B = zeros(6, 3, 8);
    
    for h in 1:8
        B[:,:,h] = [dNdx[h,1]      0          0 
                        0       dNdx[h,2]     0
                        0          0        dNdx[h,3]
                        dNdx[h,2]   dNdx[h,1]    0
                        dNdx[h,3]     0        dNdx[h,1]
                        0        dNdx[h,3]  dNdx[h,2]];
    end
    return B
end