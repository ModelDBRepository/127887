function RR = calcTreeRHS(numsections,Np,RS)

k = 1;
label = ['sect' num2str(k)];
RR = getfield( RS, label );

for k=2:numsections,
    label = ['sect' num2str(k)];
    R = getfield( RS, label );
    RR = [RR; R];
end




