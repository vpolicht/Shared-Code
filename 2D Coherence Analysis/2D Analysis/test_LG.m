function test_LG(N,rho)

  kappa = N/rho;
  xmax  = 8;                  % For plotting
  x     = xmax*(0:0.005:1)';  %   purposes only
  L     = LFuncs(x,N,kappa);
  D     = LG(N,kappa);
  I     = eye(N+1);
  rhs   = zeros(N+1,1);
  
  for i = 1:3 
      
      switch i
          case 1,  exact = exp(-x);
                   R     = D^2 - I;  
          case 2,  exact = exp(-x) .* (cos(2*x) + sin(2*x)/2);
                   R     = D^2 + 2*D + 5*I;
          case 3,  exact = exp(-x) .* (cos(x) + sin(x));
                   R     = D^4 + 4*I;
      end
      
      R(end,:) = 1;   rhs(end) = 1;
      
      if i > 1,  R(end-1,:) = sum(D);  rhs(end-1) = 0;   end
      
      a     = R\rhs;
      f     = L*a;
      err   = f - exact;
      emax  = max(abs(err));
      
      subplot(3,1,i)
      plot(x,f),  axis([0 xmax -0.6 1.2]),  grid on
      title(['Function # ', int2str(i), ...
             '   (max error = ', num2str(emax), ')' ])
  end
  
  xlabel x
  