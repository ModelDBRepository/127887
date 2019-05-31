function [B0, Bn] = cable_exact_coeffs(xa,W,NumTerms,r,L,Ra);

% current input coeffs (needed for the exact soln)

xb = xa + W;
nn = [1:NumTerms]; % size of sum
pp = 2.*pi;
C = cos(pp.*xa./W);
S = sin(pp.*xa./W);
k1n = pi.*(2./W + nn./L);
k2n = pi.*(2./W - nn./L);
B0 = W./L;
s1n = ( sin(k1n.*xb) - sin(k1n.*xa) )./ k1n;
s2n = ( sin(k2n.*xb) - sin(k2n.*xa) )./ k2n; % div by 0
c1n = ( cos(k1n.*xb) - cos(k1n.*xa) )./ k1n;
c2n = ( cos(k2n.*xb) - cos(k2n.*xa) )./ k2n;
Bn = (2./(nn.*pi)).*(sin(nn.*pi.*xb./L) - sin(nn.*pi.*xa./L)) ...
    -(C./L).*(s1n + s2n) ...
    +(S./L).*(c1n + c2n);
