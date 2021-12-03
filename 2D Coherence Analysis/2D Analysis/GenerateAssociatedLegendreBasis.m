function [B] = GenerateAssociatedLegendreBasis(X,Y,N,m)
y1 = linspace(-1,1,X);
y2 = linspace(-1,1,Y);
% B = zeros([X Y ((N+1)*N/2)]);
B = zeros([X Y (N*N)]);
c = 0;
for k=(1:N);
    temp1 = legendre(k+m-1,y1,'norm')/sqrt(X);
    for q=(1:(N));
        c = c+1;
        temp2 = legendre(q+m-1,y2,'norm')/sqrt(Y);
        B(:,:,c) = temp1(m+1,:)'*temp2(m+1,:);
    end
end
end