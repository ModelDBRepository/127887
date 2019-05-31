function PP = reduceByOne( PP );

if 1,
    % Merge the last two equations
    % Gives best results with point input located away from the end
    PP(end-1,:,:) = .5.* PP(end-1,:,:) + .5.* PP(end,:,:) ; PP(end,:,:) = []; return
else
    % Merge each pair of equations
    [nrows,ncols,nn] = size(PP);
    PPP = zeros(nrows-1,ncols,nn);
    for ii=1:nrows-1,
        PPP(ii,:,:) = .5.* PP(ii,:,:) + .5.* PP(ii+1,:,:);
    end
    PP = PPP;
    return;
end
