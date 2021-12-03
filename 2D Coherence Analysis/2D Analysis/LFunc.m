function Lf = LFunc(x,N,k)

% Lf = LFunc(x,N,k)
%
% Scaled Laguerre function of order N, evaluated at x.
% (The default scaling is k==1.)
% Lf is of size (Nx by N+1), where Nx is the length of x.

  
  if nargin < 3,  k = 1;  end

  xx = k*x;   Lf = exp(-xx/2) .* LPoly(xx,N);
  