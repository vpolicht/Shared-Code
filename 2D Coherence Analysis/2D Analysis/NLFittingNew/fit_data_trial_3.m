clear all; close all; clc;
addpath('C:\Users\fullerf\Downloads\DERIVESTsuite\DERIVESTsuite');
test_scaling = 0;
test_local_minimization = 0;
test_global_minimization = 0;
test_accuracy = 0;
test_hessian = 0;
% matlabpool(4);
true_rates = [1/100 1/400 1/200 1/250 1/2000];
true_freqs = [1/50 1/300 1/70 1/150];
L = 300;
t = linspace(1,3000,L);
NoscTrue = length(true_freqs);
true_beta = zeros(1,length(true_rates)+length(true_freqs));
true_beta(1:2:(2*NoscTrue)) = true_rates(1:NoscTrue);
true_beta(2:2:(2*NoscTrue)) = true_freqs(1:NoscTrue);
true_beta((2*NoscTrue+1):end) = true_rates(NoscTrue+1:end);
Nbeta = length(true_beta);
[true_phi] = multi_cexp2( true_beta, NoscTrue, [], 1, t);
M = 1000;
rand_amps = rand(M,size(true_phi,2))-0.5; %[1 1/3 1 1/3 2 0];
true_sig = (rand_amps*true_phi')';
SNR = 1;
noise_mat = (rand(length(t),M)/max(max(true_sig)))/SNR;
noisy_sig = true_sig + noise_mat;
y = noisy_sig;
Nosc = 4;
ada = @(x) multi_cexp2(x,Nosc,[],1,t);
beta_lb = [1/1000 1/1000 1/1000 1/1000 1/1000 1/1000 1/1000 1/1000 1/4000]';
beta_ub = [1/40   1/40   1/40   1/40   1/40   1/40   1/40   1/40   1/500]';
gen_rand_beta = @(lb,ub) ((ub-lb).*rand(length(ub),1)+lb);
beta_guess = (beta_lb+beta_ub)/2;
q = length(beta_guess);
n1 = q+1;

[phi_guess] = multi_cexp2( beta_guess, Nosc, [], 1,t);
%% Test the accuracy of the gradient against a numerically calculated one
if test_accuracy;
    objfun = @(x) varpro_gradient_test4(y,x,ada);
    %beta_test = gen_rand_beta(beta_lb,beta_ub);
    beta_test = beta_guess;
    g = gradest(objfun,beta_test);
    [~, ga, P] = varpro_gradient_test4(y,beta_test,ada);
    fprintf(1,'%g\n',(g-ga)./g);
    fprintf(1,'----\n');
end


%% Test the scaling of the objective function

if test_scaling;
    tic; for k=1:1000; [o] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
    tic; for k=1:1000; [o] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
    tic; for k=1:1000; [o] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;

    tic; for k=1:1000; [o,g] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
    tic; for k=1:1000; [o,g] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
    tic; for k=1:1000; [o,g] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;

    Nosc = 3;
    ada = @(x) multi_cexp2(x,Nosc,[],1,t);
    beta_lb = [1/1000 1/1000 1/1000 1/1000 1/1000 1/1000 1/4000]';
    beta_ub = [1/40   1/40   1/40   1/40   1/40   1/40   1/500]';

    tic; for k=1:1000; [o] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
    tic; for k=1:1000; [o] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
    tic; for k=1:1000; [o] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;

    tic; for k=1:1000; [o,g] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
    tic; for k=1:1000; [o,g] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
    tic; for k=1:1000; [o,g] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;

    Nosc = 2;
    ada = @(x) multi_cexp2(x,Nosc,[],1,t);
    beta_lb = [1/1000 1/1000 1/1000 1/1000 1/4000]';
    beta_ub = [1/40   1/40   1/40   1/40   1/500]';

    tic; for k=1:1000; [o] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
    tic; for k=1:1000; [o] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
    tic; for k=1:1000; [o] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;

    tic; for k=1:1000; [o,g] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
    tic; for k=1:1000; [o,g] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
    tic; for k=1:1000; [o,g] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;

    Nosc = 1;
    ada = @(x) multi_cexp2(x,Nosc,[],1,t);
    beta_lb = [1/1000 1/1000 1/4000]';
    beta_ub = [1/40   1/40   1/500]';

    tic; for k=1:1000; [o] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
    tic; for k=1:1000; [o] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
    tic; for k=1:1000; [o] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;

    tic; for k=1:1000; [o,g] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
    tic; for k=1:1000; [o,g] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
    tic; for k=1:1000; [o,g] = varpro_gradient_test4(y,gen_rand_beta(beta_lb,beta_ub),ada); end; toc;
end

%% Test local minimization methods
if test_local_minimization;
    Nosc = 3;
    ada = @(x) multi_cexp2(x,Nosc,[],1,t);
    beta_lb = [1/1000 1/1000 1/1000 1/1000 1/1000 1/1000 1/4000]';
    beta_ub = [1/40   1/40   1/40   1/40   1/40   1/40   1/500]';
    local_tol = 1e-6;
    wall_clock_time = 30; %seconds
    global_method = [];
    Ntests = 10;
    for k=1:Ntests;
        beta_test = gen_rand_beta(beta_lb,beta_ub);


        tic
        [Fit1,~] =  ...
        varpro_multiway_gradient3(noisy_sig,beta_test,@(x) multi_cexp2(x,Nosc,[],1,t), ...
                            beta_lb,beta_ub,global_method,wall_clock_time,local_tol,NLOPT_LD_LBFGS);
        toc
%         fprintf(1,'%g\n',1./beta_opt1);
        fprintf(1,'--\n');

        tic
        [Fit2,~] =  ...
        varpro_multiway_gradient3(noisy_sig,beta_test,@(x) multi_cexp2(x,Nosc,[],1,t), ...
                            beta_lb,beta_ub,global_method,wall_clock_time,local_tol,NLOPT_LN_SBPLX);
        toc
%         fprintf(1,'%g\n',1./beta_opt2);
        fprintf(1,'xxxx\n');
    end
end

%% Test Global Minimization Methods
if test_global_minimization;
    Nosc = 3;
    ada = @(x) multi_cexp2(x,Nosc,[],1,t);
    beta_lb = [1/1000 1/1000 1/1000 1/1000 1/1000 1/1000 1/4000]';
    beta_ub = [1/40   1/40   1/40   1/40   1/40   1/40   1/500]';
    local_tol = 1e-6;
    wall_clock_time = 10; %seconds
    global_method = NLOPT_G_MLSL_LDS;
    Ntests = 1;
    for k=1:Ntests;
        beta_test = gen_rand_beta(beta_lb,beta_ub);


        tic
        [Fit1,~] =  ...
        varpro_multiway_gradient3(noisy_sig,beta_test,@(x) multi_cexp2(x,Nosc,[],1,t), ...
                            beta_lb,beta_ub,global_method,wall_clock_time,local_tol,NLOPT_LD_LBFGS);
        toc
%         fprintf(1,'%g\n',1./beta_opt1);
        fprintf(1,'--\n');

        tic
        [Fit2,~] =  ...
        varpro_multiway_gradient3(noisy_sig,beta_test,@(x) multi_cexp2(x,Nosc,[],1,t), ...
                            beta_lb,beta_ub,global_method,wall_clock_time,local_tol,NLOPT_LN_SBPLX);
        toc
%         fprintf(1,'%g\n',1./beta_opt2);
        fprintf(1,'xxxx\n');
    end
end

%% Test Hessian Approximation and Asymtotic Error Near Found Solution

if test_hessian
    Nosc = 3;
    ada = @(x) multi_cexp2(x,Nosc,[],1,t);
    beta_lb = [1/1000 1/1000 1/1000 1/1000 1/1000 1/1000 1/4000]';
    beta_ub = [1/40   1/40   1/40   1/40   1/40   1/40   1/500]';
    local_tol = 1e-6;
    wall_clock_time = 20; %seconds
    global_method = NLOPT_G_MLSL_LDS;
    Ntests = 1;
    objfun = @(x) varpro_gradient_test4(y,x,ada);
    for k=1:Ntests;
        beta_test = gen_rand_beta(beta_lb,beta_ub);
        [Fit,Stats,exitcode] =  ...
        varpro_multiway_gradient3(noisy_sig,beta_test,@(x) multi_cexp2(x,Nosc,[],1,t), ...
                            beta_lb,beta_ub,global_method,wall_clock_time,local_tol,NLOPT_LD_LBFGS);
        ha = 2*(Fit.J'*Fit.J);
        h = hessian(objfun,Fit.beta);
        discrep = (h-ha)./h;
        fprintf(1,'%g\n',diag(discrep));
        fprintf(1,'--\n');
        NumStats = FitStatistics(h,Stats.RSS,Stats.DoF,Fit.beta,Fit.alpha,numel(y));
    end
end
