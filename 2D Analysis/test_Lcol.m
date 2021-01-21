function test_Lcol(N,xmax)

[D,k,x] = LG2(N,xmax);
  xx    = xmax*(0:0.005:1)';
  Lx    = LFuncs(x,N,k);
  Lxx   = LFuncs(xx,N,k);
  I     = eye(N+1);
  rhs   = zeros(N+1,1);

  for i = 1:3 
      
      switch i
          case 1,   exact = exp(-xx);
                    R     = Lx*(D^2 - I);
          case 2,   exact = exp(-xx) .* (cos(2*xx) + sin(2*xx)/2);
                    R     = Lx*(D^2 + 2*D + 5*I);
          case 3,   exact = exp(-xx) .* (cos(xx) + sin(xx));
                    R     = Lx*(D^4 + 4*I);
      end
      
      R(end,:) = 1;   rhs(end) = 1;
      
      if i > 1,  R(end-1,:) = sum(D);  rhs(end-1) = 0;   end
      
      a     = R\rhs;
      f     = Lx*a;
      ff    = Lxx*a;
      err   = ff - exact;
      emax  = max(abs(err(xx < 10)));
      
      subplot(3,1,i)
      plot(x,f,'.',xx,ff,'-'),  axis([0 xmax -0.6 1.2]),  grid on
      title(['Function # ', int2str(i), ...
             '   (max error = ', num2str(emax), ')' ])
  end
  
  xlabel x
