function [ ainv ] = rankKGeneralInv( q, r, k, e )
rinv = r(1:k,1:k)\eye(k);
ainv = cat(1,rinv*q(:,1:k)',zeros(size(r,1)-k,size(q,1)));
ainv = ainv(e,:);
end

