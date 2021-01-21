function [ Phi, dPhi, Ind] = multi_exp_with_oscillations_phi_and_dphi( alpha, beta, Nosc, t )
f = @(a,b,tau) exp(-tau/a).*sin(2*pi*t/b);
df_a = @(a,b,tau) (-1/a).*exp(-tau/a).*sin(2*pi*t/b);
df_b = @(a,b,tau) (2*pi/b).*exp(-tau/a).*cos(2*pi*t/b);
g = @(a,b,tau) exp(-tau/a).*cos(2*pi*t/b);
dg_a = @(a,b,tau) (-1/a).*exp(-tau/a).*cos(2*pi*t/b);
dg_b = @(a,b,tau) -(2*pi/b).*exp(-tau/a).*sin(2*pi*t/b);
h = @(a,tau) exp(-tau/a);
dh_a = @(a,tau) (-1/a)*exp(-tau/a);
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

NdPhi = 4*Nosc+(n-2*Nosc);
dPhi = zeros(m,NdPhi);
Ind = zeros(2,NdPhi);
for k=1:2:(2*Nosc)
    Ind(1,2*(k-1)+1) = k;
    dPhi(:,2*(k-1)+1) = df_a(alpha(k),alpha(k+1),t);
    Ind(1,2*(k-1)+2) = k;
    dPhi(:,2*(k-1)+2) = df_b(alpha(k),alpha(k+1),t);
    Ind(1,2*(k-1)+3) = k+1;
    dPhi(:,2*(k-1)+3) = dg_a(alpha(k),alpha(k+1),t);
    Ind(1,2*(k-1)+4) = k+1;
    dPhi(:,2*(k-1)+4) = dg_b(alpha(k),alpha(k+1),t);
    Ind(2,2*(k-1)+1) = k;
    Ind(2,2*(k-1)+2) = k+1;
    Ind(2,2*(k-1)+3) = k;
    Ind(2,2*(k-1)+4) = k+1;
end
if NdPhi>4*Nosc;
    for k=1:1:(n-2*Nosc);
        Ind(1,4*Nosc+k) = 2*Nosc+k;
        Ind(2,4*Nosc+k) = 2*Nosc+k;
        dPhi(:,4*Nosc+k) = dh_a(alpha(2*Nosc+k),t);
    end
end
end

