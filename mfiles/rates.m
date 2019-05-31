function [am,bm,ah,bh,an,bn] = rates(v)

q10 = 1;
am = .1 .* vtrap(-(v+40),10);
bm =  4 .* exp(-(v+65)./18);
ah = .07 .* exp(-(v+65)./20);
bh = 1 ./ (exp(-(v+35)./10) + 1);
an = .01.*vtrap(-(v+55),10);
bn = .125.*exp(-(v+65)./80);
