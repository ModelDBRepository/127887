function PP = calcTreeOper(numsections,Np,PS)

N = Np - 2;
numcols = sum(Np,2);
numrows = numcols - numsections.*2;
PP = zeros(numrows,numcols);


% col start and end marker for BC matrix
colStart(1)=1;
for j=2:numsections, 
    colStart(j) = colStart(j-1) + Np(j-1); 
end
colEnd = colStart + Np - 1;

% row start and end marker for BC matrix
rowStart(1)=1;
for j=2:numsections, 
    rowStart(j) = rowStart(j-1) + N(j-1); 
end
rowEnd = rowStart + N - 1;

% stack them diagonally
for k=1:numsections,
    label = ['sect' num2str(k)];
    P = getfield( PS, label );
    rowindx=rowStart(k):rowEnd(k);
    colindx=colStart(k):colEnd(k);
    PP(rowindx,colindx) = P;
end
PP = sparse(PP);


