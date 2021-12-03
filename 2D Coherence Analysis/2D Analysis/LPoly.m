function L = LPoly(x,N)

% L = LPoly(x,N)
%
% Unscaled Laguerre polynomial of order N, evaluated at x.


  f0 = ones(size(x));  if N==0,  L = f0;  return, end;
  f1 = 1 - x;          if N==1,  L = f1;  return, end;
  
  for n = 2:N
      f2 = ((2*n - 1 - x).*f1 - (n-1)*f0)/n;
      f0 = f1;  f1 = f2;
  end

  L = f2;
  