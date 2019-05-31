function PPbar = reduceByBC(PP,K,ColIndx)

Np = size(PP,2);
N = Np - size(ColIndx,2);

PP = rearrangeCols(PP,ColIndx);

PP_hat = PP(:,1:N);
T = PP(:,N+1:end);

PPbar = PP_hat + T*K;
