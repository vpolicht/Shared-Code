function [obj, grad, J, P] = varpro_gradient_test4(y, beta, ada)
LIM.q = length(beta);
LIM.L = size(y,1);
LIM.M = size(y,2);
[phi0, dphi0, Ind] = feval(ada, beta);
LIM.n1 = size(phi0,2);
indtrix = zeros(LIM.q,size(dphi0,2));
for k=1:LIM.q; indtrix(k,:) = (Ind(2,:) == k); end;
indtrix = sparse(indtrix);

if nargout>3
    [obj, grad, J, ~, P] = obj_fun(beta);
elseif nargout==3
    [obj, grad, J] = obj_fun(beta);
elseif nargout==2
    [obj, grad] = obj_fun(beta);
else
    [obj] = obj_fun(beta);
end

function [o, gradient, J, c, P] = obj_fun(beta)
if nargout>1
    gradient = zeros(1,LIM.q);
    if nargout>2
        J = zeros(LIM.L*LIM.M,LIM.q);
    end
end
P = struct();
P.resid = zeros(LIM.L,LIM.M);
[P.P, P.dP, P.dPI] = feval(ada, beta);
c = zeros(LIM.n1,LIM.M);
tol = eps*norm(P.P,'fro');%eps*size(Phi,1);
[P.Q, P.R, e] = qr(P.P,0);
[~,P.Z] = sort(e);
P.k = sum(diag(abs(P.R))>tol);
P.Pm = rankKGeneralInv(P.Q,P.R,P.k,P.Z);
c(:,1) = P.Pm*y(:,1);
P.resid(:,1) = y(:,1)-(P.P*c(:,1));
if nargout>1;
    P.Sm = (bsxfun(@times,indtrix,c(P.dPI(1,:),1)'))';
    P.T = P.dP*P.Sm;
    gradient = gradient -2*(P.resid(:,1)'*P.T) + 2*(((P.resid(:,1)'*P.P)*P.Pm)*P.T);
    if nargout>2
        P.Sm = (bsxfun(@times,indtrix,c(P.dPI(1,:),1)'))';
        P.T1 = P.dP*P.Sm;
        P.T2 = diag((P.dP*indtrix')'*P.resid(:,1));
        J(((1-1)*LIM.L+1):(1*LIM.L),:) = (P.T1 - (P.P*(P.Pm*P.T1)))+P.Pm(1:LIM.q,:)'*P.T2;
    end
end
if nargout>1
    for i=2:LIM.M;
        c(:,i) = P.Pm*y(:,i);
        P.resid(:,i) = y(:,i)-(P.P*c(:,i));
        P.Sm = (bsxfun(@times,indtrix,c(P.dPI(1,:),i)'))';
        P.T  = P.dP*P.Sm;
        gradient = gradient -2*(P.resid(:,i)'*P.T) + 2*(((P.resid(:,i)'*P.P)*P.Pm)*P.T);
    end
    if nargout>2
        for i=2:LIM.M;
            P.Sm = (bsxfun(@times,indtrix,c(P.dPI(1,:),i)'))';
            P.T1 = P.dP*P.Sm;
            P.T2 = diag((P.dP*indtrix')'*P.resid(:,i));
            J(((i-1)*LIM.L+1):(i*LIM.L),:) = (P.T1 - (P.P*(P.Pm*P.T1)))+P.Pm(1:LIM.q,:)'*P.T2;
        end
    end
else
    for i=2:LIM.M;
        c(:,i) = P.Pm*y(:,i);
        P.resid(:,i) = y(:,i)-(P.P*c(:,i));
    end
end
o = sum(P.resid(:).^2);
end

end

