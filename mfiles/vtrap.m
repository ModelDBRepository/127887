function vtrap = vtrap(x,y)

if abs(x./y) < 1e-6
    vtrap = y.*(1 - x./y./2);
else
    vtrap = x./(exp(x./y) - 1);
end
