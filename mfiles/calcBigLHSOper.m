function PPbarbar = calcBigLHSOper(PPbar,Np)
%
% Create the LHS matrix for full system (including channel vars m,h,n)
% Input PPbar is the LHS oper for membrane pot
% Np is the array of section sizes
%

n = size(PPbar,2);

N = Np - 2;

Size = sum( N + 3.*Np );

PPbarbar = speye( Size );

PPbarbar(1:n,1:n) = PPbar;

