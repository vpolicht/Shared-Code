function D = LG(N,k)

% D = LG(N,k)
%
% Differentiation matrix for Laguerre Galerkin scheme with spectral order N
% and scale factor "k".
%
% This function is deprecated.  A better version is L_Galerkin, which
% chooses "k" internally to match a user-defined interval [0,xmax].

  One = ones(N+1);
  Eye = eye(N+1);
  D   = k * (Eye/2 - triu(One));
  