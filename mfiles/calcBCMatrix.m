function [BC,ColIndx] = calcBCMatrix(numsections,Np,r0,rL,Ra,FS,BC_branchends,BC_nodecurrents,BC_nodeequality)

% col start and end marker for BC matrix
colStart(1)=1;
for j=2:numsections, colStart(j)=colStart(j-1)+Np(j-1); end
colEnd = colStart + Np - 1;

ColIndx = [];
for k=1:numsections, ColIndx = [ColIndx,colStart(k),colEnd(k)]; end

% section cross section areas at left and right
A0 = r0.^2.*pi;  D0 = A0./Ra;
AL = rL.^2.*pi;  DL = AL./Ra;

% size of BC matrix
numrows = size(BC_branchends,1) + size(BC_nodecurrents,1) + size(BC_nodeequality,1);
numcols = sum(Np);
BC = zeros(numrows,numcols);

countrow = 1;

% branch sealed ends
M = BC_branchends;
for row = 1:size(M,1),
    secNum = M(row,1); % which section
    secEnd = M(row,2); % which end is sealed       
    if secEnd==0, secEnd=1; else, secEnd=Np(secNum); end    
    label = ['sect' num2str(secNum)];
    F = getfield(FS,label);
    BC(countrow,colStart(secNum):colEnd(secNum)) = F(secEnd,:);
    countrow = countrow + 1;
end


% current sum at nodes
M = BC_nodecurrents;
for row = 1:size(M,1),
    numbranches = M(row,1);
    for j = 2:2:numbranches.*2,
        secNum = M(row,j); % which section
        secEnd = M(row,j+1); % which end is the node
        if secEnd==0, 
            secEnd=1; DD = D0(secNum);
        else,
            secEnd=Np(secNum); DD = -DL(secNum);
        end
        label = ['sect' num2str(secNum)];
        F = getfield(FS,label);
        BC(countrow,colStart(secNum):colEnd(secNum)) = DD.*F(secEnd,:);
    end
    countrow = countrow + 1;
end


% mem pot equality at nodes
M = BC_nodeequality;
for row = 1:size(M,1),
    secNum1 = M(row,1); secEnd1 = M(row,2);
    secNum2 = M(row,3); secEnd2 = M(row,4);
    if secEnd1==0, n1 = colStart(secNum1); else n1 = colEnd(secNum1); end
    if secEnd2==0, n2 = colStart(secNum2); else n2 = colEnd(secNum2); end
    BC(countrow,n1) = 1;
    BC(countrow,n2) = -1;
    countrow = countrow + 1;
end









