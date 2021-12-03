function [ x, k ] = qrlinreg( A, b, tol)
%This function performs a qr decomposition based linear regression for the
%problem: Ax = b, solve for x.  Matlab already implements this directly in
%mldivide, but it doesn't return q, r, or the rank of the solution it
%found.
%inputs:
%   A, a struct containing the fields Q, R, and Z corresponding to the
%   output of a qr(.,0) command.
%   b, a vector.
%outputs: 
%   x, the solution
%   k, rank estimate.
%
if nargin<3
    tol = eps*size(A.Q,1);
end
k = sum(diag(abs(A.R))>tol);
t = A.R(1:k,1:k)\(A.Q(:,1:k)'*b);
x = cat(1,t,zeros(size(A.Q,2)-k,1));
x = x(A.Z);

end

