function [Fit, Stats, ExitCode] = ...
          varpro_multiway_gradient3(y, beta, ada, lb, ub, global_method, wall_clock_time, local_tol, local_method)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialization: Check input, set default parameters and options.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[L,M] = size(y);      % L = number of time points, M = number of observations

[q,ell] = size(beta);   % q = number of nonlinear parameters

if (ell > 1)
   error('beta must be a column vector containing initial guesses for nonlinear parameters')
end

if (nargin < 6)  
    error('a lower bound must be specified');
else
    [q1,ell] = size(lb);
    if (q1 > 0) && (ell > 0)
       if (q1 ~= q) || (ell > 1)
          error('lb must be empty or a column vector of the same length as beta')
       end
    end
end

if (nargin < 7)
    error('an upper bound must be specified');
else
    [q1,ell] = size(ub);
    if (q1 > 0) && (ell > 0)
       if (q1 ~= q) || (ell > 1)
          error('ub must be empty or a column vector of the same length as beta')
       end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make the first call to ada and do some error checking / initialization
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try 
    [Phi, dPhi, Ind] = feval(ada, beta);
catch %#ok<CTCH>
%     fprintf(1,'No derivative information supplied.\n');
    [Phi] = feval(ada, beta);
    dPhi = [];
    Ind = [];
end
[d1,d2] = size(dPhi);
[I1,I2] = size(Ind);

[m1,n1] = size(Phi);     % n1 = number of basis functions Phi.

if (L ~= m1)
    error('In user function ada: number of rows in Phi must be L');
end

if (m1 ~= d1) && (d1 > 0)
   error('In user function ada: Phi and dPhi must have the same number of rows.')
end

if (I2 > 0) && (I1 ~= 2)
   error('In user function ada: Ind must have two rows.')
end

if (d2 > 0) && (d2 ~= I2)
   error('In user function ada: dPhi and Ind must have the same number of columns.')
end

LIM.q = q; %num non-lin params
LIM.L = L; %the number of time points
LIM.M = M; %the number of pixels
LIM.n1 = n1; %the number of linear parameters per pixel (including offsets and fixed rates).

if d2>0;
    indtrix = zeros(LIM.q,d2);
    for k=1:q; indtrix(k,:) = (Ind(2,:) == k); end;
    indtrix = sparse(indtrix);
end

if (q > 0) % The problem is nonlinear.
    if isempty(global_method)
        opt.algorithm = local_method;
        opt.ftol_rel = local_tol;
    else
        opt.algorithm = global_method;
        opt.maxtime = wall_clock_time;
        opt.local_optimizer.algorithm = local_method;
        opt.local_optimizer.ftol_rel = local_tol; %stopping criteria for local optimizer
    end
   opt.min_objective = @(x) obj_fun(x);
   opt.lower_bounds = lb;
   opt.upper_bounds = ub;
   [Fit.beta, ~, retcode] = nlopt_optimize(opt,beta);
   if nargout>1
       [RSS,~,J,Fit.alpha,~] = obj_fun(Fit.beta);
       %Fit.indtrix = indtrix;
       ExitCode = retcode;
       DoF = LIM.L*LIM.M-LIM.q-LIM.n1*LIM.M;
       hess = 2*(J'*J);
       Stats = FitStatistics(hess,RSS,DoF,Fit.beta,Fit.alpha,numel(y));
   end
   
else       % The problem is linear.
    error('The problem is linear.  This shit be way too fancy for you; go away and use the backslash operator.');
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
if nargout==1
    [P.P] = feval(ada, beta);
else
    [P.P, P.dP, P.dPI] = feval(ada, beta);
end
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
        %P.T2 = diag((P.dP*indtrix')'*P.resid(:,1));
        %J(((1-1)*LIM.L+1):(1*LIM.L),:) = (P.T1 - (P.P*(P.Pm*P.T1)))+P.Pm(1:LIM.q,:)'*P.T2;
        A = zeros(size(P.P,2),LIM.q);
        for which_face=1:LIM.q;
            which_cols = (P.dPI(2,:)==which_face);
            which_rows = P.dPI(1,which_cols);
            A(which_rows,which_face) = (P.dP(:,which_cols))'*P.resid(:,1);
        end
        J(((1-1)*LIM.L+1):(1*LIM.L),:) = (P.T1 - (P.P*(P.Pm*P.T1)))+P.Pm'*A;
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
            %P.T2 = diag((P.dP*indtrix')'*P.resid(:,i));
            %J(((i-1)*LIM.L+1):(i*LIM.L),:) = (P.T1 - (P.P*(P.Pm*P.T1)))+P.Pm(1:LIM.q,:)'*P.T2;
            A = zeros(size(P.P,2),LIM.q);
            for which_face=1:LIM.q;
                which_cols = (P.dPI(2,:)==which_face);
                which_rows = P.dPI(1,which_cols);
                A(which_rows,which_face) = (P.dP(:,which_cols))'*P.resid(:,i);
            end
            J(((1-1)*LIM.L+1):(1*LIM.L),:) = (P.T1 - (P.P*(P.Pm*P.T1)))+P.Pm'*A;
        end
    end
else
    for i=2:LIM.M;
        c(:,i) = P.Pm*y(:,i);
        P.resid(:,i) = y(:,i)-(P.P*c(:,i));
    end
end
o = norm(P.resid);
end


end


