function RR = getRHSInput( I0, r, W, xa, Np, x, dummyarg )

I0 = I0.*2./(W.*2.*pi.*r); % assumed const radius here

xb = xa + W;
pp = 2.*pi;

% current input
I02C = (I0./2);
na = max(find(x<xa));
nb = max(find(x<xb));

J = zeros(Np-1,1);

if na==nb,
    J(na) = I02C.*( (xb-xa) - W./(pp).*sin(pp./W.*(xb-xa)) );
else,
    J(na) = I02C.*( (x(na+1)-xa) - W./(pp).*sin(pp./W.*(x(na+1)-xa)) );
    J(nb) = I02C.*( (xb-x(nb)) - W./(pp).*(sin(pp./W.*(xb-xa)) - sin(pp./W.*(x(nb)-xa))) );
    NN = [1:Np];
    indx = find( NN>na & NN<nb );
    J(indx) = I02C.*( (x(indx+1)-x(indx)) ...
        - W./(pp).*(sin(pp./W.*(x(indx+1)-xa)) - sin(pp./W.*(x(indx)-xa))) );
end

RR = reduceByOne(J);

