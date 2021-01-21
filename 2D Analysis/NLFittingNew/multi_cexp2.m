function [ Phi, dPhi, Ind, dPhi2, Ind2] = multi_cexp2( beta, Nosc, const_rates, const_bool, t )
betaosc = beta(1:(2*Nosc));
betar = beta(2*Nosc+1:end);
f = @(a,b,t) exp(-t*a).*sin(2*pi*t*b);
g = @(a,b,t) exp(-t*a).*cos(2*pi*t*b);
h = @(a,t) exp(-t*a);
n = length(beta);
m = length(betar);
z = length(const_rates);
l = length(t);
const_num = const_bool>0;
if const_num;
    Phi = zeros(l,n+z+1);
else
    Phi = zeros(l,n+z+const_num);
end
for k=1:2:(2*Nosc)
    Phi(:,k) = g(betaosc(k),betaosc(k+1),t);
    Phi(:,k+1) = f(betaosc(k),betaosc(k+1),t);
end
c = 0;
for k=(Nosc*2+1):n
    c = c+1;
    Phi(:,k) = h(betar(c),t);
end
c = 0;
for k=(n+1):(n+z);
    c = c+1;
    Phi(:,k) = h(const_rates(c),t);
end
if const_num;
    Phi(:,n+z+1) = ones(l,1);
end

    if nargout>1
        dfda = @(a,b,t) -t.*exp(-t*a).*sin(2*pi*t*b);
        dgda = @(a,b,t) -t.*exp(-t*a).*cos(2*pi*t*b);
        dfdb = @(a,b,t) 2*pi*t.*exp(-t*a).*cos(2*pi*t*b);
        dgdb = @(a,b,t) -2*pi*t.*exp(-t*a).*sin(2*pi*t*b);
        dhda = @(a,t) -t.*exp(-t*a);
        d2fda2 = @(a,b,t) (t.^2).*exp(-t*a).*sin(2*pi*t*b);
        d2fdadb = @(a,b,t) -2*pi*(t.^2).*exp(-t*a).*cos(2*pi*t*b);
        d2gda2 = @(a,b,t) (t.^2).*exp(-t*a).*cos(2*pi*t*b);
        d2gdadb = @(a,b,t) 2*pi*(t.^2).*exp(-t*a).*sin(2*pi*t*b);
        d2fdb2 = @(a,b,t) -4*pi^2*(t.^2).*exp(-t*a).*sin(2*pi*t*b);
        d2gdb2 = @(a,b,t) -4*pi^2*(t.^2).*exp(-t*a).*cos(2*pi*t*b);
        d2hda2 = @(a,t) (t.^2).*exp(-t*a);
        NdPhi = 4*Nosc + m; %2 functions per oscillation, and 2 derivatives per function + 1 derivative per rate.

        dPhi = zeros(l,NdPhi);
        Ind = zeros(2,NdPhi);
        for k=1:2:(2*Nosc)
            Ind(1,2*(k-1)+1) = k; %derivative of kth function wrt kth value of beta
            dPhi(:,2*(k-1)+1) = dgda(beta(k),beta(k+1),t);
            Ind(1,2*(k-1)+2) = k; %derivative of kth function wrt to k+1 value of beta
            dPhi(:,2*(k-1)+2) = dgdb(beta(k),beta(k+1),t);
            Ind(1,2*(k-1)+3) = k+1; %derivative of k+1th function wrt kth value of beta
            dPhi(:,2*(k-1)+3) = dfda(beta(k),beta(k+1),t);
            Ind(1,2*(k-1)+4) = k+1; %derivative of k+1th function wrt k+1th value of beta
            dPhi(:,2*(k-1)+4) = dfdb(beta(k),beta(k+1),t);
            Ind(2,2*(k-1)+1) = k; %derv wrt k
            Ind(2,2*(k-1)+2) = k+1; %derv wrt k+1
            Ind(2,2*(k-1)+3) = k; %derv wrt k
            Ind(2,2*(k-1)+4) = k+1; %derv wrt k+1
        end
        if NdPhi>4*Nosc;
            for k=1:m;
                Ind(1,4*Nosc+k) = 2*Nosc+k; %derv of 2*Nosc + kth function
                Ind(2,4*Nosc+k) = 2*Nosc+k; %wrt 2*Nosc + kth value of beta
                dPhi(:,4*Nosc+k) = dhda(beta(2*Nosc+k),t);
            end
        end

        if nargout > 3
            NdPhi2 = 8*Nosc + m;%2 functions per oscillation, 4 derivatives per function + 1 derivative per rate
            dPhi2 = zeros(l,NdPhi2);
            Ind2 = zeros(3,NdPhi2);
            for k=1:2:(2*Nosc)
                Ind2(1,4*(k-1)+1) = k; %d2f_(k)/(d(k)d(k))
                Ind2(2,4*(k-1)+1) = k; 
                Ind2(3,4*(k-1)+1) = k; 
                dPhi2(:,4*(k-1)+1) = d2gda2(beta(k),beta(k+1),t);
                
                Ind2(1,4*(k-1)+2) = k; %d2f_(k)/(d(k+1)d(k+1))
                Ind2(2,4*(k-1)+2) = k+1; 
                Ind2(3,4*(k-1)+2) = k+1; 
                dPhi2(:,4*(k-1)+2) = d2gdb2(beta(k),beta(k+1),t);
                
                Ind2(1,4*(k-1)+3) = k+1; %d2f_(k+1)/(d(k)d(k))
                Ind2(2,4*(k-1)+3) = k; 
                Ind2(3,4*(k-1)+3) = k; 
                dPhi2(:,4*(k-1)+3) = d2fda2(beta(k),beta(k+1),t);
                
                Ind2(1,4*(k-1)+4) = k+1; %d2f_(k+1)/(d(k+1)d(k+1))
                Ind2(2,4*(k-1)+4) = k+1; 
                Ind2(3,4*(k-1)+4) = k+1; 
                dPhi2(:,4*(k-1)+4) = d2fdb2(beta(k),beta(k+1),t);
                
                Ind2(1,4*(k-1)+5) = k; %d2f_(k)/(d(k)d(k+1))
                Ind2(2,4*(k-1)+5) = k; 
                Ind2(3,4*(k-1)+5) = k+1; 
                dPhi2(:,4*(k-1)+5) = d2gdadb(beta(k),beta(k+1),t);
                
                Ind2(1,4*(k-1)+6) = k+1; %d2f_(k+1)/(d(k)d(k+1))
                Ind2(2,4*(k-1)+6) = k; 
                Ind2(3,4*(k-1)+6) = k+1; 
                dPhi2(:,4*(k-1)+6) = d2fdadb(beta(k),beta(k+1),t);
                
                Ind2(1,4*(k-1)+7) = k; %d2f_(k)/(d(k)d(k+1))
                Ind2(2,4*(k-1)+7) = k+1; 
                Ind2(3,4*(k-1)+7) = k; 
                dPhi2(:,4*(k-1)+7) = d2gdadb(beta(k),beta(k+1),t);
                
                Ind2(1,4*(k-1)+8) = k+1; %d2f_(k+1)/(d(k)d(k+1))
                Ind2(2,4*(k-1)+8) = k+1; 
                Ind2(3,4*(k-1)+8) = k; 
                dPhi2(:,4*(k-1)+8) = d2fdadb(beta(k),beta(k+1),t);
            end
            if NdPhi2>8*Nosc;
                for k=1:m;
                    Ind2(1,8*Nosc+k) = 4*Nosc+k; %d2f_(k)/(d(k)d(k))
                    Ind2(2,8*Nosc+k) = 4*Nosc+k;
                    Ind2(3,8*Nosc+k) = 4*Nosc+k;
                    dPhi2(:,8*Nosc+k) = d2hda2(beta(2*Nosc+k),t);
                end
            end
        end

    end
end

