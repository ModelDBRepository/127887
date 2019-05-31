function K = getK(BC,ColIndx)

BC = rearrangeCols(BC,ColIndx);

Np = size(BC,2); % size of full system 
N = Np - size(ColIndx,2); % size of reduced system

B_hat = BC(:,1:N);
C = BC(:,N+1:end);
K = -inv(C)*B_hat;
