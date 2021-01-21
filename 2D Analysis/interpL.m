function C = interpL(x,xx,k)

% C = interpL(x,xx,k)

  x = k*x;  xx = k*xx;   %  revert to unscaled coords

  N = length(x) - 1;   dLf = dLfun(xx,N);   Lf = LFunc(x,N+1);
  
  [X,XX] = meshgrid(x,xx);   [LF,DLF] = meshgrid(Lf,dLf);
  
  Equal = XX==X;   dX = XX - X + Equal;   [rows,cols] = find(Equal);
  
  C = -DLF./(LF.*dX);   C(rows,:) = 0;   C = C + Equal;
  