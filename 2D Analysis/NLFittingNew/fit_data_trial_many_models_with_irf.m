clear all; close all; clc;
%parpool(12)
addpath(fullfile('C:','Users','Franklin','Documents','Work','Matlab','2D Analysis','NLFitting'));
L = 300;
t = linspace(100,3000,L);
fub = 1/(2*mean(diff(t))); % nyquist frequency
flb = 1/(t(end)/(pi)); %require at least 1 half period
gen_rand_beta = @(lb,ub) ((ub-lb).*rand(length(ub),1)+lb);
NoscTrue = 2;
NdynRatesTrue = 2;
StatRatesTrue = 1/10000;
OscBoundsTrue.fub = fub;
OscBoundsTrue.flb = flb;
OscBoundsTrue.rub = fub;
OscBoundsTrue.rlb = flb;
RateBoundsTrue.rub = [1/100 1/1000];
RateBoundsTrue.rlb = [1/1000 1/4000];
IRFBoundsTrue.ub = [10 300];
IRFBoundsTrue.lb = [-10 20];
ModelTrue = CreateSoLModelWithIRF(NoscTrue,NdynRatesTrue,StatRatesTrue,0,OscBoundsTrue,RateBoundsTrue,IRFBoundsTrue,t);
BetaTrue = gen_rand_beta(ModelTrue.beta_lb,ModelTrue.beta_ub);
PhiTrue = ModelTrue.modelfun(BetaTrue);
M = 1000;
rand_amps = rand(M,size(PhiTrue,2))-0.5;
true_sig = (rand_amps*PhiTrue')';
SNR = 5;
noise_mat = (rand(length(t),M)/max(max(true_sig)))/SNR;
noisy_sig = true_sig + noise_mat;
y = noisy_sig;

MaxNOsc = 2;
MinNOsc = 1;
MaxNdynRates = 2;
MinNdynRates = 1;
MaxConstRates = [1/5000 1/10000];
OscBounds = OscBoundsTrue;
RateBounds = RateBoundsTrue;
[ModelCell, Diagnostics] = GenManySoLModelsWithIRF(MaxNOsc, MinNOsc, MaxNdynRates, MinNdynRates, MaxConstRates, OscBounds, RateBounds, IRFBoundsTrue, t);

%constant fitting parameters
local_tol = 1e-5;
wall_clock_time = 5; %seconds
global_method = NLOPT_G_MLSL_LDS;
Ntests = 1;

parfor k=1:length(ModelCell)
    beta_start = gen_rand_beta(ModelCell{k}.beta_lb,ModelCell{k}.beta_ub);
    ada = ModelCell{k}.modelfun;
    [ModelCell{k}.Fit,ModelCell{k}.Stats,ModelCell{k}.exitcode] =  ...
        varpro_multiway_gradient3(y,beta_start,ada, ...
                            ModelCell{k}.beta_lb,ModelCell{k}.beta_ub,global_method,wall_clock_time,local_tol,NLOPT_LD_LBFGS);
end
if ispc
    base = 'Z:';
else
    base = '/mnt/data/';
end
savepath = fullfile(base,'2D','Frank','Matlab','2D Analysis','NLFittingNew','Results');
save([savepath 'fit_trial_results.mat'],'ModelCell');



