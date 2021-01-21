function [ Phi ] = multi_exp_phi( alpha, beta, t )
f = @(a,tau) exp(-tau/a);
n = length(alpha);
d = length(beta);
m = length(t);
Phi = zeros(m,n+d+1);
for k=1:length(alpha)
    Phi(:,k) = f(alpha(k),t);
end
for k=1:length(beta);
    Phi(:,n+k) = f(beta(k),t);
end
Phi(:,n+d+1) = ones(m,1);
end

