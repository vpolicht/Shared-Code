function [ C ] = padpost2( A,padsize )
dim=size(A);
B=zeros(padsize,dim(2));
C=cat(1,A,B);


end

