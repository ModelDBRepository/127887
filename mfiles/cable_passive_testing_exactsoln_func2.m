function V = cable_passive_testing_exactsoln_func2(B0,Bn,V0,g,D,L,x,I0,C,t,NumTerms)

nn = [1:NumTerms];
gn = g + D.*(nn.*pi./L).^2;
egt = exp(-g.*t);
egnt = exp(-gn.*t);

f1 = V0.*egt + I0./(2.*C).*(B0./g).*(1-egt);

f2n = I0./(2.*C).*cos(x*nn.*pi./L)*diag((Bn./gn).*(1-egnt));

V = f1 + sum(f2n,2);







