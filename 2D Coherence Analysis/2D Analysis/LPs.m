function [D,x,k] = LPs(N,xmax)

% [D,x,k] = LPs(N,xmax)
%
% Order-N pseudospectral scheme based on scaled Laguerre functions
% Lf(n,x),  0 <= n <= N.  The functions are scaled by a factor "k" so that
% the collocation points lie in the domain [0,xmax].
%
%   x  -  collocation points  (column vector of length N+1)
%   D  -  differentiation matrix  (square, size N+1)
%
%
% EXAMPLE:
%
% Solve the following ODE:
%
%      f'' - 4*f  =  0,    f(0) = kappa,   f -> 0   as  x -> infty.
%
% This, of course, has the solution  f(x) = kappa*exp(-2*x).
% Suppose that we don't know this, but that we do suspect that it decays
% to zero as exp(-2*x).
% In that case, for best results we should probably set xmax to 10.0 (or even
% larger), implying at least a million-fold decay of the function over the
% computational interval.
%
% Solution:
%
%    N = 20;  xmax = 10;  kappa = 1.0;
%
%    [D,x,k] = LPs(N,xmax);   D2 = D^2;   I = eye(N+1);
%
%    Nx = 40;  xx = xmax*((0:Nx)/Nx);   C = interp_L(x,xx,k);
%
%    A = D2 - 4*I;   b = zeros(N+1,1);           %  differential eqn
%    A(1,:) = 0;   A(1,1) = 1;   b(1) = kappa;   %  boundary condition
%
%    f = A\b;   ff = C*f;   plot(x,f,'*', xx,ff,'-')



% Collocation points of unscaled functions are found by solving
%
%         Lfun(t,N+1) == Lfun(t,N).
%
  t_try = [ (1/N)  (2/N):0.1:3  4:(4*N) ];    % sample points for t > 0
  df    = dLfun(t_try,N);
  ind   = find( df(1:end-1).*df(2:end) < 0 );
  t_est = (t_try(ind) + t_try(ind+1))/2;
  
  if (length(t_est) ~= N),  error('Failed to find all points'), end
  
% unscaled points
  t = zeros(N+1,1);
  warning off
      for n = 1:N,  t(n+1) = fzero('dLfun', t_est(n),[],N);  end
  warning on
  
  tmax = t(end);
  k    = tmax/xmax;
  x    = t/k;
  I    = eye(N+1);
  X    = repmat(x,1,N+1);
  dX   = X - X' + I;
  Lf   = LFunc(x,N+1,k);
  LF   = repmat(Lf,1,N+1);
  
% Off-diagonal entries
%   (Diagonal entries are given false values - in fact, set to unity.)
  D = LF./(LF' .* dX);
  
% Set diagonal entries (for x > 0) to zero
%    (the diagonal entry for x==0 has a known non-zero value)
  D = D - diag(diag(D));   D(1,1) = -k*(N+1)/2;
  
  
  
  
  
  