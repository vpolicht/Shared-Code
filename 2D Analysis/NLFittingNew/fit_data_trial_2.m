clear all; close all; clc;
addpath('C:\Users\fullerf\Downloads\DERIVESTsuite\DERIVESTsuite');
true_rates = [1/100 1/400 1/2000];
true_freqs = [1/50 1/300];
true_beta = [true_rates(1) true_freqs(1) true_rates(2) true_freqs(2) true_rates(3)];
%true_beta = [1/2000];
L = 1000;
t = linspace(100,3000,L);
Nosc = 2;
[true_phi] = multi_cexp( true_beta, Nosc, [], t);
M = 1;
rand_amps = rand(M,length(true_beta)+1);
true_sig = (rand_amps*true_phi')';
SNR = 5;
noise_mat = (rand(length(t),M)/max(max(true_sig)))/SNR;
noisy_sig = true_sig + noise_mat;
beta_lb = [1/200 1/100 1/1000 1/800 1/4000]';
beta_ub = [1/30 1/10 1/100 1/30 1/500]';
%beta_guess = [1/50 1/40 1/100 1/50 1/1000]';
beta_guess = true_beta';
ei = [0 0 0 0 1]';
%ei = [1];
h = 1E-6;
beta_guess_h_p = beta_guess+h*ei;
beta_guess_h_m = beta_guess-h*ei;
q = length(beta_guess);
n1 = q+1;
ada = @(x) multi_cexp(x,Nosc,[],t);
ZZ = 25;
h_lin = h*(-ZZ:1:ZZ);
zero_loc = dsearchn(h_lin',0);
beta_mat = beta_guess(:,ones(1,length(h_lin))) + h_lin(ones(length(ei),1),:).*ei(:,ones(1,length(h_lin)));
objfun = @(x) varpro_gradient_test(noisy_sig,x,ada);
obj_mat = zeros(1,length(h_lin)); for k=1:length(h_lin); obj_mat(k) = varpro_gradient_test(noisy_sig,beta_mat(:,k),ada); end;
[obj_val, obj_grad, obj_hess, J, c, residual, Hnn] = varpro_gradient_test(noisy_sig,beta_guess,ada);
plot(h_lin,obj_mat);
pfit = polyfit(h_lin,obj_mat,2);
panalytic = [obj_hess(boolean(ei),boolean(ei))/2 obj_grad(boolean(ei)) obj_val];
hess_hat = pfit(1)*2;
grad_hat = pfit(2);
hess_hat_robust = hessian(objfun,beta_guess);
hess_hat_diff = (obj_mat(zero_loc+1)-2*obj_mat(zero_loc)+obj_mat(zero_loc-1))/(h^2);
grad_hat_diff = (obj_mat(zero_loc+1)-obj_mat(zero_loc-1))/(2*h);
rel_change = (obj_mat(zero_loc+1)-obj_mat(zero_loc))/obj_mat(zero_loc);
hold all; plot(h_lin,polyval(pfit,h_lin),'o');
plot(h_lin,polyval(panalytic,h_lin),'x');
fprintf('Analytic Gradient: %g\n',obj_grad(boolean(ei)));
fprintf('Estimated Gradient from Fit: %g\n',grad_hat);
fprintf('Estimated Gradient from Diff: %g\n',grad_hat_diff);
fprintf('Analytic Hessian: %g\n',obj_hess(boolean(ei),boolean(ei)));
fprintf('Estimated Hessian from Fit: %g\n',hess_hat);
fprintf('Estimated Hessian from Diff: %g\n',hess_hat_diff);
fprintf('Estimated Hessian from robust estimate: %g\n',hess_hat_robust);
fprintf('Relative Change of Obj: %g\n',rel_change);

% % tic
% % %[beta_opt, lin_terms, resid, resid_norm, y_est, Regression] =  ...
% % %varpro_multiway_2(noisy_sig,0*noisy_sig+1,beta_guess,length(beta_guess),@(x) multi_cexp(x,Nosc,[],t),beta_lb,beta_ub,1e3,1e-4);
% % [beta_opt, lin_terms, resid, resid_norm, y_est, Regression] =  ...
% % varpro_multiway_gradient(noisy_sig,beta_guess,@(x) multi_cexp(x,Nosc,[],t),beta_lb,beta_ub,1e3,1e-4,NLOPT_LD_MMA);
% % toc
% % for k=1:length(beta_guess); fprintf(1,'%g\n',1./beta_opt(k)); end;
% 
% [obj, grad1, hess1] = varpro_gradient_test(noisy_sig,beta_guess,ada);
% [obj_hp, grad2, hess2] = varpro_gradient_test(noisy_sig,beta_guess_h_p,ada);
% [obj_hm, grad3, hess3] = varpro_gradient_test(noisy_sig,beta_guess_h_m,ada);
% 
% grad_hat = (obj_hp-obj_hm)/(2*h); %central difference approximation of first derv
% hess_hat = (obj_hp-2*obj+obj_hm)/(h^2); %central difference approx of 2nd derv
% rel_obj = (obj_hp-obj)/obj;
% 
% fprintf('Analytic Gradient: %g\n',grad1(boolean(ei)));
% fprintf('Estimated Gradient: %g\n',grad_hat);
% fprintf('Analytic Hessian: %g\n',hess1(boolean(ei),boolean(ei)));
% fprintf('Estimated Hessian: %g\n',hess_hat);
% fprintf('Relative Obj Change: %g\n',rel_obj);