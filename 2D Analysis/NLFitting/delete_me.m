function [beta, c, resid, wresid_norm, y_est, Regression] = ...
          varpro_multiway_gradient(y, beta, ada, lb, ub, global_evals, local_tol, method)

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
% Make the first call to ada and do some error checking.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Phi, dPhi, Ind] = feval(ada, beta);
[d1,d2] = size(dPhi);
[I1,I2] = size(Ind);

[m1,n1] = size(Phi);     % n1 = number of basis functions Phi.

if (n1 < q)
   error('In user function ada: The number of columns in Phi must be >= #non_lin_params')
end

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

if (q > 0) % The problem is nonlinear.
   opt.algorithm = NLOPT_GN_MLSL_LDS; %multi level single leakage with quasi-random sampling (low discrepancy sequence)
   %opt.min_objective = @(x) sum(f_lsq(x).^2); %the objective function to be used.
   opt.min_objective = @(x) obj_fun(x);
   opt.maxeval = global_evals; %stopping criteria for global optimizer
   opt.lower_bounds = lb;
   opt.upper_bounds = ub;
   opt.local_optimizer.algorithm = method;
   opt.local_optimizer.ftol_rel = local_tol; %stopping criteria for local optimizer
   [beta, wresid_norm2, retcode] = nlopt_optimize(opt,beta);
   [resid, y_est, c] = f_lsq(beta);
   wresid_norm = sqrt(wresid_norm2);
   Regression.report.exitflag = retcode;

else       % The problem is linear.
    error('The problem is linear.  This shit be way too fancy for you; go away.');
end
fprintf(1,'VARPRO is finished.\n');

function [objective, obj_gradient] = obj_fun(beta_trial)
    [residual, obj_gradient] = f_lsq(beta_trial);
    objective = sum(residual.^2);
end
    

function [residual, gradient, y_est, c] = f_lsq(beta_trial)

[Phi_trial, dPhidB, dPhiInd] = feval(ada, beta_trial);
[c, residual, y_est] = calcResidual(Phi_trial);
if nargout>1
    sdPhiM = zeros(L*M,size(dPhidB,2));
    J = zeros(L*M,q);
    for i=0:(M-1)
        for j=1:q;
            ind_filt = (dPhiInd(1,:)==j);
            sdPhiM((L*i+1):(L*(i+1)),ind_filt) = c(j,i+1)*dPhidB(:,ind_filt);
        end
    end
    for j=1:q
        ind_filt = (dPhiInd(2,:)==j);
        J(:,j) = sum(sdPhiM(:,ind_filt),2);
    end
    gradient = 2*residual'*J;
end
end

function [c, resid, y_est] = calcResidual(Phi)
warning('off','MATLAB:rankDeficientMatrix');
c = zeros(n1,M);
resid = 0*y;
y_est = 0*y;
T = Phi*Phi';
try rankPhi = rank(T);
catch err %#ok<NASGU>
    rankPhi = [];
end
if rankPhi == n1;
    f = @(P,y) (P'*P)\(P'*y); %this is way faster, but is only good if P'P has full rank.
else
    f = @(P,y) P\y; %else we use this slower, but less rank sensitive solver
end
for kkk=1:M;
    c(:,kkk) = f(Phi,y(:,kkk));
    y_est(:,kkk) = Phi*c(:,kkk);
    resid(:,kkk) = y(:,kkk) - y_est(:,kkk);
end
resid = reshape(resid,[L*M 1]);

end

end

