function Q0 = getStructureMatrix( P, Np, r, gl, Ra, F );

circumf = 2.*pi.*r;
area = pi.*r.^2;

A_tmp = 1./circumf;
B_tmp = area./Ra;


Q0 = P *( diag(A_tmp) * F * diag(B_tmp) * F - diag(gl) );





