function [ r, c ] = FindInMatrix( M, x )
%Inputs:
%M is the matrix to be searched
%x is the target to be found
%
%Outputs:
%r the row the target is found on
%c the column the target is found on

y = dsearchn(M(:),x);
r = mod(y,size(M,1));
c = (y-r)/size(M,1)+1;
end

