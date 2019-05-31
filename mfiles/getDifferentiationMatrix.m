function [F,G,x] = getDifferentiationMatrix( space_discretization_method, Np, L)

x0 = 0;
xL = L;

% differentiation matrix and grid

if space_discretization_method==0,
%    fprintf('Using Cheby collocation\n');
    [Fch,Gch,xch] = CHEB(xL,x0,Np);
    x = xch; % includes boundaries
    F = Fch; G = Gch; clear xch Fch Gch
else,
    orderFD = space_discretization_method; % 2,4,6 order FD
%    fprintf('Using FD order %d\n',orderFD);
    dx = (xL-x0)./(Np-1);
    F = FD(Np,orderFD,1,dx);
    G = zeros(Np); 
    G(2:Np,2:Np) = inv(F(2:Np,2:Np));
    x = (x0 + [0:Np-1].*dx)';
end
