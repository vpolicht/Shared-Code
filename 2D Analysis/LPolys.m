function L = LPolys(x,N)

% L = LPolys(x,N)
%
% Unscaled Laguerre polynomials of order 0 to N, evaluated at x.
% L is of size (Nx by N+1), where Nx is the length of x.

  x = x(:);  Nx = length(x);  L = zeros(Nx,N+1);

  f0 = ones(size(x));  L(:,1) = f0;   if N==0,  return, end;
  f1 = 1 - x;          L(:,2) = f1;   if N==1,  return, end;
  
  for n = 2:N
      f2 = ((2*n - 1 - x).*f1 - (n-1)*f0)/n;   L(:,n+1) = f2;
      f0 = f1;  f1 = f2;
  end
