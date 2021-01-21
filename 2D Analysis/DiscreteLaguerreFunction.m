function [Ft,Ff,f] = DiscreteLaguerreFunction(t,n,a)
% This function implements a discrete laguerre function, as described in 
% "Discrete Orthonormal Sequences" by Paul Broome, 1965.
% http://dl.acm.org/citation.cfm?id=321265
%
% Laguerre Function is defined to be the product of the Laguerre Polynomial
% with the Exponential exp(-ax), where 'a' is means of scaling the
% function.  Scaling is critical for applications, so it is passed as an
% argument.
% 
% Inputs:
% t = the discrete time points you would like the function evaluated at.
% Must be even
% n = the order of the Laguerre Function.  This function will return a
%      matrix of all Laguerre function up to and including the Nth term.
% a = a scalar which scales the time axis. |A|<1 is required.
%
% Outputs: 
% Ft = A matrix of size [T x (N+1)], containing the time domain Laguerre
% functions at the time points requested in T
% Fw = A matrix of size [T x (N+1)], containing the frequency domain
% Laguerre functions at the frequency points 'w', corresponding to T.
% w = a row vector the same length as T, containing associated frequencies
% to the time axis T.

% Formulae implemented:
% Ft(z,a,n) = sqrt(1 - a^2)*((a - z.^(-1) ).^(n))./((1 - a*z.^(-1)).^(n+1));
% Ff(f,a,n) = dT*sqrt(1 - a^2)*((a - exp(-1i*2*pi*f*dT)).^(n))./((1 - a*exp(-1i*2*pi*f*dT)).^(n+1))

assert(mod(length(t),2) == 0);
if size(t,1)>1
    t = t';
end
if size(t,1)>1; error('input t must be a vector'); end;

dT = mean(diff(t));
f = (-length(t)/2:1:((length(t)/2) - 1)).*(1/range(t));

Ft = zeros(length(t),n+1);
Ff = zeros(length(t),n+1);
Ft(:,1) = (sqrt(1 - a^2)*1./((1 - a*t.^(-1))))';
Ff(:,1) = (dT*sqrt(1 - a^2)*1./((1 - a*exp(-1i*2*pi*f*dT)).^(1)))';
for k=2:(n+1);
    Ft(:,k) = Ft(:,k-1).*(((a - t.^(-1) ))./((1 - a*t.^(-1))))';
    Ff(:,k) = Ff(:,k-1).*(((a - exp(-1i*2*pi*f*dT)))./((1 - a*exp(-1i*2*pi*f*dT))))';
end


end