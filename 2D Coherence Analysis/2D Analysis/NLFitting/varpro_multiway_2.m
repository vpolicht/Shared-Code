function [alpha, c, wresid, wresid_norm, y_est, Regression] = ...
          varpro_multiway_2(y, w, alpha, n, ada, lb, ub, global_evals, local_tol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialization: Check input, set default parameters and options.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,z] = size(y);      % m = number of observations
[m1,z1] = size(w);

if (m1 ~= m) || (z ~= z1)
   error('y and w must be of the same size')
end

[q,ell] = size(alpha);   % q = number of nonlinear parameters

if (ell > 1)
   error('alpha must be a column vector containing initial guesses for nonlinear parameters')
end

if (nargin < 6)  
    error('a lower bound must be specified');
else
    [q1,ell] = size(lb);
    if (q1 > 0) && (ell > 0)
       if (q1 ~= q) || (ell > 1)
          error('lb must be empty or a column vector of the same length as alpha')
       end
    end
end

if (nargin < 7)
    error('an upper bound must be specified');
else
    [q1,ell] = size(ub);
    if (q1 > 0) && (ell > 0)
       if (q1 ~= q) || (ell > 1)
          error('ub must be empty or a column vector of the same length as alpha')
       end
    end
end

W = cell(1);
for k=1:z;
    W{k} = spdiags(w(:,k),0,m,m); % Create an m x m diagonal matrix from the vector w
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make the first call to ada and do some error checking.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Phi, dPhi, Ind] = feval(ada, alpha);
[d1,d2] = size(dPhi);
[I1,I2] = size(Ind);

[m1,n1] = size(Phi);     % n1 = number of basis functions Phi.

if (n1 < n)
   error('In user function ada: The number of columns in Phi must be >= n')
end

if (m ~= m1)
    error('In user function ada: number of rows in Phi must be m');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solve the least squares problem using lsqnonlin or, if there
% are no nonlinear parameters, using the SVD procedure in formJacobian.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (q > 0) % The problem is nonlinear.
   opt.algorithm = NLOPT_GN_MLSL_LDS; %multi level single leakage with quasi-random sampling (low discrepancy sequence)
   opt.min_objective = @(x) sum(f_lsq(x).^2); %the objective function to be used.
   opt.maxeval = global_evals; %stopping criteria for global optimizer
   opt.lower_bounds = lb;
   opt.upper_bounds = ub;
   opt.local_optimizer.algorithm = NLOPT_LN_BOBYQA;  %the local optimizer algorithm
   opt.local_optimizer.ftol_rel = local_tol; %stopping criteria for local optimizer
   [alpha, wresid_norm2, retcode] = nlopt_optimize(opt,alpha);
   [wresid, ~, y_est, c] = f_lsq(alpha);
   wresid_norm = sqrt(wresid_norm2);
   Regression.report.exitflag = retcode;

else       % The problem is linear.
    error('The problem is linear.  This shit be way too fancy for you; go away.');
end
fprintf(1,'VARPRO is finished.\n');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The computation is now completed.  
%
% varpro uses the following two functions, f_lsq and formJacobian.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------- Beginning of f_lsq --------------------------

function [wr_trial, Phi_trial, y_est, c] = f_lsq(alpha_trial)

[Phi_trial] = feval(ada, alpha_trial);

[c, wr_trial, y_est] = calcResidual(Phi_trial);

end %--------------------------- End of f_lsq ---------------------------

%----------------------- Beginning of formJacobian ----------------------

function [c, wresid, y_est] = calcResidual(Phi)
warning('off','MATLAB:rankDeficientMatrix');
c = zeros(n1,z);
wresid = 0*y;
y_est = 0*y;
T = Phi*Phi';
try rankPhi = rank(T);
catch err %#ok<NASGU>
    rankPhi = [];
end
if rankPhi == n1;
    f = @(P,y) P'*P\(P'*y);
else
    f = @(P,y) P\y;
end
for kkk=1:z;
    c(:,kkk) = f(Phi,y(:,kkk));
    y_est(:,kkk) = Phi*c(:,kkk);
    wresid(:,kkk) = y(:,kkk) - y_est(:,kkk);
end
wresid = reshape(wresid,[m*z 1]);

end %-------------------------- End of calcResidual ----------------------

end %------------------------------ End of varpro ------------------------

