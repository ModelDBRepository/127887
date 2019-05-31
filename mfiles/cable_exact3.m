function Vexact = cable_exact3(xa,W,r,L,Ra,V0,gl,x,I0,Cm,T)

NumTerms = 2000; 


[B0, Bn] = cable_exact_coeffs(xa,W,NumTerms,r,L,Ra);

circumf = 2.*pi.*r; area = pi.*r.^2; 
D = area./(Ra.*circumf); % assuming uniform properties
II0 = I0.*2./(W.*2.*pi.*r);

Vexact = cable_passive_testing_exactsoln_func2(B0,Bn,V0,gl,D,L,x,II0,Cm,T,NumTerms);
