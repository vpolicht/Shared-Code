function [D,k,x] = LG2(N,xmax)

% [D,k,x] = LG2(N,xmax)
%
% Spectral differentiation scheme based on scaled Laguerre functions
% Lf(n,x),  0 <= n <= N.  The functions are scaled by a factor "k" so that
% the coresponding collocation points "x" lie in the domain [0,xmax].
%
% D is a square matrix of size N+1.
% Its (n+1)-th column expresses the derivative of Lf(n,x) as a linear
% combination of  { Lf(k,x),  0 <= k <= n }.
%
% Let f(x) be a function defined on [0,infty), and let "c" be its Laguerre-
% function coefficients, expressed as a column vector of length N+1.
% Then the scalar
%                d_n  =  D(n+1,:) * c
%
% is the spectral coefficient of Lf(n,x) for the function f'(x).
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
%    N = 20;  xmax = 10;  kappa = 1.0;  Nx = 40;
%
%    [D,k] = LG2(N,xmax);   D2 = D^2;    I = eye(N+1);
%
%    xx  = xmax*((0:Nx)'/Nx);  %  interpolating points
%    c2f = LFuncs(xx,N,k);     %  convert spectral coefficients "c"
%                              %    to function values at "xx"
%    A = D2 - 4*I;             %  basic operator - its (n+1)-th row is the
%                              %    coefficient of Lf(n,x) in the residual
%    b = zeros(N+1,1);         %  want residual to be zero-valued ...
%                              %  ... but must sacrifice one eqn to
%                              %      implement BC, ie f(0) = 1 ...
%    A(end,:) = 0;             %  ... so sacrifice the condition of
%                              %      zero Lf(N) in residual ...
%                              %  Since  Lf(n,0) == 1  for all n,
%    A(end,:) = 1;             %    and we require  f(0) = kappa,
%    b(end,1) = kappa;         %    it follows that the spectral 
%                              %    coefficients of f must sum to kappa.
%    c  = A\b;                 %  spectral coefficients of f ...
%    ff = c2f * c;             %  ... and the solution!
%    plot(xx,ff)


  [D_PS,x,k] = LPs(N,xmax);   %  Get scale factor and collocation points
                               %   from pseudospectral code  (ignore D_PS,
                               %   the PS differentiation matrix)                   
  One = ones(N+1);
  Eye = eye(N+1);
  D   = k * (Eye/2 - triu(One));
  
  
  
  