function Lf = LFuncs(x,N,k)

% Lf = LFuncs(x,N,k)
%
% Scaled Laguerre functions of order 0 to N, evaluated at x.
% (The default scaling is k==1.)
% Lf is of size (Nx by N+1), where Nx is the length of x.

  
  if nargin < 3,  k = 1;  end

  xx = k*x;   Lf = diag(exp(-xx/2)) * LPolys(xx,N);
  