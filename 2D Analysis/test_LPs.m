function emax = test_LPs(N,xmax,flag)

  [D,x,kappa] = LPs(N,xmax);
  
  xx    = xmax*(0:0.005:1)';
  C     = interpL(x,xx,kappa);
  I     = eye(N+1);
  rhs   = zeros(N+1,1);
  
  switch flag
      case 1
          exact = exp(-xx);    R = D^2 - I;          
      case 2
          exact = exp(-xx) .* (cos(2*xx) + sin(2*xx)/2);
          R     = D^2 + 2*D + 5*I; 
          R(2,:) = D(1,:);   rhs(2) = 0;
      case 3
          exact = exp(-xx) .* (cos(xx) + sin(xx));
          R     = D^4 + 4*I; 
          R(2,:) = D(1,:);   rhs(2) = 0;
  end
  
  R(1,:) = 0;   R(1,1) = 1;  rhs(1) = 1;
  
  f     = R\rhs;
  ff    = C*f;
  err   = ff - exact;
  emax  = max(abs(err(xx < 10)));
  
  subplot(2,1,1)
  plot(x,f,'.',xx,ff,'-')
  axis([0 xmax -0.2 1.05]),  xlabel x,  ylabel f,  grid on
  title(['Maximum error = ', num2str(emax) ])
  
  subplot(2,1,2)
  plot(xx,err)
  