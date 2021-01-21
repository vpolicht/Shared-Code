function [ Phi ] = multi_exp_with_oscillations_phi( alpha, beta, Nosc, t )
f = @(a,b,tau) exp(-tau/a).*sin(2*pi*t/b);
g = @(a,b,tau) exp(-tau/a).*cos(2*pi*t/b);
h = @(a,tau) exp(-tau/a);
n = length(alpha);
n2 = length(beta);
m = length(t);
Phi = zeros(m,n+1);
for k=1:2:(Nosc*2)
    Phi(:,k) = f(alpha(k),alpha(k+1),t);
    Phi(:,k+1) = g(alpha(k),alpha(k+1),t);
end
for k=(Nosc*2+1):1:n
    Phi(:,k) = h(alpha(k),t);
end
c = 0;
for k=(n+1):(n+n2);
    c = c+1;
    Phi(:,k) = h(beta(c),t);
end
Phi(:,n+n2+1) = ones(m,1);
end

