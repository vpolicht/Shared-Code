function [alpha, x, wresid, wresid_norm, y_est, Regression] = ...
          varpro_multiway_3(y, w, alpha, n, ada, lb, ub)
%VARPRO Solve a separable nonlinear least squares problem.
% [alpha, c, wresid, wresid_norm, y_est, Regression] =
%             VARPRO(y, w, alpha, n, ada, lb, ub, options)
%
% Given a set of m observations y(1),...,y(m)
% this program computes a weighted least squares fit using the model
%
%    eta(alpha,c,t) = 
%            c_1 * phi_1 (alpha,t) + ...  + c_n * phi_n (alpha,t) 
% (possibly with an extra term  + phi_{n+1} (alpha,t) ).
%
% This program determines optimal values of the q nonlinear parameters
% alpha and the n linear parameters c, given observations y at m
% different values of the "time" t and given evaluation of phi and 
% (optionally) derivatives of phi.
%
% On Input:
%
%   y    m x z   vector containing the m observations x z instances of that observation.
%   w    m x 1   vector of weights used in the least squares
%                fit.  We minimize the norm of the weighted residual
%                vector r, where, for i=1:m,
%
%                r(i) = w(i) * (y(i) - eta(alpha, c, t(i,:))).
%
%                Therefore, w(i) should be set to 1 divided by
%                the standard deviation in the measurement y(i).  
%                If this number is unknown, set w(i) = 1.
%   alpha q x 1  initial estimates of the parameters alpha.
%                If alpha = [], Varpro assumes that the problem
%                is linear and returns estimates of the c parameters.
%   n            number of linear parameters c
%   ada          a function handle, described below.
%   lb    q x 1  lower bounds on the parameters alpha. 
%   (Optional)   (Omit this argument or use [] if there are
%                no lower bounds.)
%   ub    q x 1  upper bounds on the parameters alpha. 
%   (Optional)   (Omit this argument or use [] if there are
%                no upper bounds.)
%   options      The Matlab optimization parameter structure,
%   (Optional)   set by "optimset", to control convergence
%                tolerances, maximum number of function evaluations,
%                information displayed in the command window, etc. 
%                To use default options, omit this parameter.
%                To determine the default options, type
%                    options = optimset('lsqnonlin')
%                After doing this, the defaults can be modified;
%                to modify the display option, for example, type
%                    options = optimset('lsqnonlin');
%                    optimset(options,'Display','iter');
%
% On Output:
%
%  alpha  q x 1  estimates of the nonlinear parameters.
%  c      n x 1  estimates of the linear parameters.
%  wresid m x 1  weighted residual vector, with i-th component
%                w(i) * (y(i) - eta(alpha, c, t(i,:))).
%  wresid_norm   norm of wresid.
%  y_est  m x 1  the model estimates = eta(alpha, c, t(i,:)))
%  Regression    a structure containing diagnostics about the model fit.
%                **************************************************
%                *                C a u t i o n:                  *
%                *   The theory that makes statistical            *
%                *   diagnostics useful is derived for            *
%                *   linear regression, with no upper- or         *
%                *   lower-bounds on variables.                   *
%                *   The relevance of these quantities to our     *
%                *   nonlinear model is determined by how well    *
%                *   the linearized model (Taylor series model)   *
%                *         eta(alpha_true, c_true)                *
%                *            +  Phi * (c  - c_true)              *
%                *            + dPhi * (alpha - alpha_true)       *
%                *   fits the data in the neighborhood of the     *
%                *   true values for alpha and c, where Phi       *
%                *   and dPhi contain the partial derivatives     *
%                *   of the model with respect to the c and       *
%                *   alpha parameters, respectively, and are      *
%                *   defined in ada.                              *
%                **************************************************
%  Regression.report:  
%                This structure includes information on the solution
%                process, including the number of iterations, 
%                termination criterion, and exitflag from lsqnonlin.
%                (Type 'help lsqnonlin' to see the exit conditions.)
%                Regression.report.rank is the computed rank of the 
%                matrix for the linear subproblem.  If this equals
%                n, then the linear coefficients are well-determined.
%                If it is less than n, then although the model might
%                fit the data well, other linear coefficients might
%                give just as good a fit.
%  Regression.sigma:        
%                The estimate of the standard deviation is the
%                weighted residual norm divided by the square root
%                of the number of degrees of freedom.
%                This is also called the "regression standard error"
%                or the square-root of the weighted SSR (sum squared
%                residual) divided by the square root of the
%                number of degrees of freedom.
%  Regression.RMS:
%                The "residual mean square" is equal to sigma^2:
%                RMS = wresid_norm^2 / (m-n+q)
%  Regression.coef_determ:
%                The "coefficient of determination" for the fit,
%                also called the square of the multiple correlation
%                coefficient, is sometimes called R^2.
%                It is computed as 1 - wresid_norm^2/CTSS,
%                where the "corrected total sum of squares"
%                CTSS is the norm-squared of W*(y-y_bar),
%                and the entries of y_bar are all equal to
%                (the sum of W_i^2 y_i) divided by (the sum of W_i^2).
%                A value of .95, for example, indicates that 95 per 
%                cent of the CTTS is accounted for by the fit.
%
%  Regression.CovMx: (n+q) x (n+q)
%                This is the estimated variance/covariance matrix for
%                the parameters.  The linear parameters c are ordered
%                first, followed by the nonlinear parameters alpha.
%                This is empty if dPhi is not computed by ada.
%  Regression.CorMx: (n+q) x (n+q)
%                This is the estimated correlation matrix for the 
%                parameters.  The linear parameters c are ordered
%                first, followed by the nonlinear parameters alpha.
%                This is empty if dPhi is not computed by ada.
%  Regression.std_param: (n+q) x 1
%                This vector contains the estimate of the standard 
%                deviation for each parameter.
%                The k-th element is the square root of the k-th main 
%                diagonal element of Regression.CovMx
%                This is empty if dPhi is not computed by ada.
%  Regression.t_ratio:   (n+q) x 1
%                The t-ratio for each parameter is equal to the
%                parameter estimate divided by its standard deviation.
%                (linear parameters c first, followed by alpha)
%                This is empty if dPhi is not computed by ada.
%  Regression.standardized_wresid:
%                The k-th component of the "standardized weighted 
%                residual" is the k-th component of the weighted 
%                residual divided by its standard deviation.
%                This is empty if dPhi is not computed by ada.
%
%---------------------------------------------------------------
% Specification of the function ada, which computes information
% related to Phi:
%
%   function [Phi,dPhi,Ind] = ada(alpha)
%
%     This function computes Phi and, if possible, dPhi.
%
%     On Input: 
%
%        alpha q x 1    contains the current value of the alpha parameters.
%
%        Note:  If more input arguments are needed, use the standard
%               Matlab syntax to accomplish this.  For example, if
%               the input arguments to ada are t, z, and alpha, then
%               before calling varpro, initialize t and z, and in calling 
%               varpro, replace "@ada" by "@(alpha)ada(t,z,alpha)".
%
%     On Output:
%
%        Phi   m x n1   where Phi(i,j) = phi_j(alpha,t_i).
%                       (n1 = n if there is no extra term; 
%                        n1 = n+1 if an extra term is used)
%        dPhi  m x p    where the columns contain partial derivative
%                       information for Phi and p is the number of 
%                       columns in Ind 
%                       (or dPhi = [] if derivatives are not available).
%        Ind   2 x p    Column k of dPhi contains the partial
%                       derivative of Phi_j with respect to alpha_i, 
%                       evaluated at the current value of alpha, 
%                       where j = Ind(1,k) and i = Ind(2,k).
%                       Columns of dPhi that are always zero, independent
%                       of alpha, need not be stored. 
%        Example:  if  phi_1 is a function of alpha_2 and alpha_3, 
%                  and phi_2 is a function of alpha_1 and alpha_2, then 
%                  we can set
%                          Ind = [ 1 1 2 2
%                                  2 3 1 2 ]
%                  In this case, the p=4 columns of dPhi contain
%                          d phi_1 / d alpha_2,
%                          d phi_1 / d alpha_3,
%                          d phi_2 / d alpha_1,
%                          d phi_2 / d alpha_2,
%                  evaluated at each t_i.
%                  There are no restrictions on how the columns of
%                  dPhi are ordered, as long as Ind correctly specifies
%                  the ordering.
%
%        If derivatives dPhi are not available, then set dPhi = Ind = [].
%      
%---------------------------------------------------------------
%
%  Varpro calls lsqnonlin, which solves a constrained least squares
%  problem with upper and lower bounds on the constraints.  What
%  distinguishes varpro from lsqnonlin is that, for efficiency and
%  reliability, varpro causes lsqnonlin to only iterate on the
%  nonlinear parameters.  Given the information in Phi and dPhi, this 
%  requires an intricate but inexpensive computation of partial 
%  derivatives, and this is handled by the varpro function formJacobian.
%
%  lsqnonlin is in Matlab's Optimization Toolbox.  Another solver
%  can be substituted if the toolbox is not available.
%
%  Any parameters that require upper or lower bounds should be put in
%  alpha, not c, even if they appear linearly in the model.
%  
%  The original Fortran implementation of the variable projection 
%  algorithm (ref. 2) was modified in 1977 by John Bolstad 
%  Computer Science Department, Serra House, Stanford University,
%  using ideas of Linda Kaufman (ref. 5) to speed up the 
%  computation of derivatives.  He also allowed weights on
%  the observations, and computed the covariance matrix.
%  Our Matlab version takes advantage of 30 years of improvements
%  in programming languages and minimization algorithms.  
%  In this version, we also allow upper and lower bounds on the 
%  nonlinear parameters.
%
%  It is hoped that this implementation will be of use to Matlab
%  users, but also that its simplicity will make it easier for the
%  algorithm to be implemented in other languages.
%
%  This program is documented in
%  Dianne P. O'Leary and Bert W. Rust,
%  Variable Projection for Nonlinear Least Squares Problems,
%  Computational Optimization and Applications (2012)
%  doi 10.1007/s10589-012-9492-9.
%
%  US National Inst. of Standards and Technology, 2010.
%
%  Main reference:
%
%    0.  Gene H. Golub and V. Pereyra, 'Separable nonlinear least   
%        squares: the variable projection method and its applications,'
%        Inverse Problems 19, R1-R26 (2003).
%         
%  See also these papers, cited in John Bolstad's Fortran code:
%                                                                       
%    1.  Gene H. Golub and V. Pereyra, 'The differentiation of      
%        pseudo-inverses and nonlinear least squares problems whose 
%        variables separate,' SIAM J. Numer. Anal. 10, 413-432      
%         (1973).                                                    
%    2.  ------, same title, Stanford C.S. Report 72-261, Feb. 1972.
%    3.  Michael R. Osborne, 'Some aspects of non-linear least      
%        squares calculations,' in Lootsma, Ed., 'Numerical Methods 
%        for Non-Linear Optimization,' Academic Press, London, 1972.
%    4.  Fred Krogh, 'Efficient implementation of a variable projection
%        algorithm for nonlinear least squares problems,'           
%        Comm. ACM 17:3, 167-169 (March, 1974).                    
%    5.  Linda Kaufman, 'A variable projection method for solving  
%        separable nonlinear least squares problems', B.I.T. 15,   
%        49-57 (1975).                                             
%    6.  C. Lawson and R. Hanson, Solving Least Squares Problems,
%        Prentice-Hall, Englewood Cliffs, N. J., 1974.          
%
%  These books discuss the statistical background:
%
%    7.  David A. Belsley, Edwin Kuh, and Roy E. Welsch, Regression 
%        Diagnostics, John Wiley & Sons, New York, 1980, Chap. 2.
%    8.  G.A.F. Seber and C.J. Wild, Nonlinear Regression,
%        John Wiley & Sons, New York, 1989, Sec. 2.1, 5.1, and 5.2.
%
%  Dianne P. O'Leary, NIST and University of Maryland, February 2011.
%  Bert W. Rust,      NIST                             February 2011.
%  Comments updated 07-2012.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

[Phi] = feval(ada, alpha);

[m1,n1] = size(Phi);     % n1 = number of basis functions Phi.

if (n1 < n) || (n1 > n+1)
   error('In user function ada: The number of columns in Phi must be n or n+1.')
end

if (m ~= m1)
    error('In user function ada: number of rows in Phi must be m');
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
   opt.maxeval = 1e3; %stopping criteria for global optimizer
   opt.lower_bounds = lb;
   opt.upper_bounds = ub;
   opt.local_optimizer.algorithm = NLOPT_LN_BOBYQA;  %the local optimizer algorithm
   opt.local_optimizer.ftol_rel = 1e-6; %stopping criteria for local optimizer
 
%    [alpha, wresid_norm2, wresid, exitflag,output] = ...
%                   lsqnonlin(@f_lsq, alpha, lb, ub, options);
   [alpha, wresid_norm2, retcode] = nlopt_optimize(opt,alpha);
   [wresid, ~, y_est,x] = f_lsq(alpha);
   wresid_norm = sqrt(wresid_norm2);
   Regression.report.exitflag = retcode;

else       % The problem is linear.

   if (~strcmp(options.Display,'off'))
      fprintf(1,'VARPRO problem is linear, since length(alpha)=0.');
   end

   [x, wresid, y_est, Regression.report.rank] =  ...
                                       calcResidual(alpha, Phi, dPhi);
   wresid_norm = norm(wresid);
   wresid_norm2 = wresid_norm * wresid_norm;

end % if q > 0
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Compute some statistical diagnostics for the solution.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate sample variance,  the norm-squared of the residual
%    divided by the number of degrees of freedom.

% sigma2 = wresid_norm2 / (m-n-q);
% 
% % Compute  Regression.sigma:        
% %                square-root of weighted residual norm squared divided 
% %                by number of degrees of freedom.
% 
% Regression.sigma = sqrt(sigma2);
% 
% % Compute  Regression.coef_determ:
% %                The coeficient of determination for the fit,
% %                also called the square of the multiple correlation
% %                coefficient, or R^2.
% %                It is computed as 1 - wresid_norm^2/CTSS,
% %                where the "corrected total sum of squares"
% %                CTSS is the norm-squared of W*(y-y_bar),
% %                and the entries of y_bar are all equal to
% %                (the sum of W_i y_i) divided by (the sum of W_i).
% 
% y_bar = sum(w.*y) / sum(w);
% CTTS = norm(W * (y - y_bar)) ^2;
% Regression.coef_determ = 1 - wresid_norm^2 / CTTS;
% 
% % Compute  Regression.RMS = sigma^2:
% %                the weighted residual norm divided by the number of
% %                degrees of freedom.
% %                RMS = wresid_norm / sqrt(m-n+q)
% 
% Regression.RMS = sigma2;                        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% End of statistical diagnostics computations.
% Print some final information if desired.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  fprintf(1,' ');
  fprintf(1,'VARPRO Results:\n');
  fprintf(1,' Linear Parameters:\n');
  fprintf(1,' %15.7e \n',x);
  
  fprintf(1,' Nonlinear Parameters:\n');
  fprintf(1,' %15.7e \n',alpha);
  fprintf(1,' ');
  
  fprintf(1,' Norm-squared of weighted residual  = %15.7e\n',wresid_norm2);
  fprintf(1,' Norm-squared of data vector        = %15.7e\n',norm(w.*y)^2);
  fprintf(1,' Norm         of weighted residual  = %15.7e\n',wresid_norm);
  fprintf(1,' Norm         of data vector        = %15.7e\n',norm(w.*y));
%   fprintf(1,' Expected error of observations     = %15.7e\n',Regression.sigma);
%   fprintf(1,' Coefficient of determination       = %15.7e\n',Regression.coef_determ);
  fprintf(1,'VARPRO is finished.\n');
  fprintf(1,'-------------------\n');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The computation is now completed.  
%
% varpro uses the following two functions, f_lsq and formJacobian.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------- Beginning of f_lsq --------------------------

function [wr_trial, Phi_trial, y_est, x] = f_lsq(alpha_trial)

[Phi_trial] = feval(ada, alpha_trial);

[x, wr_trial, y_est] = calcResidual(Phi_trial);

end %--------------------------- End of f_lsq ---------------------------

%----------------------- Beginning of formJacobian ----------------------

function [x, wresid, y_est, rankPhi] = calcResidual(Phi)
x = zeros(n1,z);
wresid = 0*y;
y_est = 0*y;
[~,nl] = size(Phi);
[Q,R,P] = qr(Phi);

% Determine rank of A.
% The diagonal entries of R satisfy
%|R(1,1)| >= |R(2,2)| >= |R(3,3)| >= ..
% Find the smallest integer r such that
%|R(r+1,r+1)| < max(size(A))*eps*|R(1,1)|
tol = max(size(Phi))*eps*abs(R(1,1));
Rtrunc = R(1:nl,1:nl);
Rdiag = diag(Rtrunc);
rankPhi = sum(abs(Rdiag)>tol);
Q2 = Q(:,rankPhi+1:end);

for kkk=1:z;
    c = Q'*y(:,kkk);
    y1 = R(1:rankPhi,1:rankPhi) \ c(1:rankPhi);
    y2 = zeros(nl-rankPhi,1);
    xtemp = (P*[y1;y2]);
    x(:,kkk) = xtemp;
    y_est(:,kkk) = Phi*x(:,kkk);
    wresid(:,kkk) = Q2*Q2'*y(:,kkk);
end
wresid = reshape(wresid,[m*z 1]);

end %-------------------------- End of calcResidual ----------------------

end %------------------------------ End of varpro ------------------------

