function [F, G, x_ch] = CHEB(Th,L,N)
% Chebyshev grid and operator
% L,Th: lower,upper bounds
% N grid size, Vs,r: reset,reversal potential
% n,b: gamma distr. parameters

% differentiation matrix and grid
[F,y]=cheb0(N-1);  % y goes from 1 to -1 and has size N
y = flipud(y);
F = -(2./(Th-L)).*F;

% integration matrix
G = zeros(N); 
G(2:N,2:N) = inv(F(2:N,2:N));

% actual grid from cheby grid
x_ch = (Th-L)./2 .*(y+1) + L;



