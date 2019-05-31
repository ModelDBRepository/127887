function BC = rearrangeCols(BC,ColIndx)
%
% ColIndx contains col numbers which are to be solved for
% Those cols are to become the last cols of BC
%

% find cols which are not in ColIndx
numcols = size(BC,2);
NN = [1:numcols];
indx = setdiff(NN,ColIndx); % the reduced system indices

% assign the non-ColIndx cols
BCtmp = BC;
BC = BCtmp(:,indx);

% assign the ColIndx cols at the end
BC = [BC, BCtmp(:,ColIndx)];

