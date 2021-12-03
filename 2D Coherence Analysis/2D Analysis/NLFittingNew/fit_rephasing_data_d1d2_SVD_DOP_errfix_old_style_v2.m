%% Preliminaries
close all; clear all; clc
if ispc;
    base = 'Z:';
else
    base = '/mnt/data';
end
%parpool(12)
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','NLFittingNew'));
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','DOPBoxV1_7','DOPBoxV1-7','DOPbox'));
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis'));
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','NLFitting','cm_and_cb_utilities','cm_and_cb_utilities'));
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','freezeColors_v23_cbfreeze','freezeColors'));
c = 299.792458; %nm/fs
w2wn = @(w) 1E7*w/(2*pi*c);
wn2f = @(wn) (c*wn)./(1E7);
wn2w = @(wn) (2*pi*c*wn)./(1E7);
wn2T = @(wn) 1./((c*wn)./(1E7)); %wave number to period
T2wn = @(tau) 1./((c*tau)./(1E7));
w2l = @(w) 2*pi*c./w;
%% Inputs
data_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','PostProcessingData');
results_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','PostProcessingData','NewFits');
figures_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','Figures');
data_file_name_r = 'data_stack_r.mat';
data_file_name_nr = 'data_stack_nr.mat';
%% Load the rephasing data, compress it.
if (~exist('all_data_mat_r','var') || ~exist('T2vec','var'));
    load(fullfile(data_path,data_file_name_r));
end
Trange = 1:185;
T = -T2vec(Trange);
all_data_mat_r = all_data_mat_r(:,:,Trange);
Nbf = 60;
Bx = dop( size(all_data_mat_r,2), Nbf );
By = dop( size(all_data_mat_r,1), Nbf );
[SRX, SRY, SRZ] = size(all_data_mat_r);
Rcompressed = zeros(Nbf,Nbf,SRZ);
for k=1:SRZ; Rcompressed(:,:,k) = By' * squeeze(all_data_mat_r(:,:,k)) * Bx; end;
R = zeros(Nbf*Nbf,SRZ);
for k=1:Nbf
    R(((k-1)*Nbf+1):(k*Nbf),:) = squeeze(Rcompressed(:,k,:));
end
[Ur, Sr, Vr] = svd(R,0);
Yr = (Sr*Vr')';
%% Load the non-rephasing data, compress it.
% if (~exist('all_data_mat_nr','var') || ~exist('T2vec','var'));
%     load(fullfile(data_path,data_file_name_nr));
% end
% all_data_mat_nr = all_data_mat_nr(:,:,Trange);
% Nbf = 60;
% NRcompressed = zeros(Nbf,Nbf,SRZ);
% for k=1:SRZ; NRcompressed(:,:,k) = By' * squeeze(all_data_mat_nr(:,:,k)) * Bx; end;
% NR = zeros(Nbf*Nbf,SRZ);
% for k=1:Nbf
%     NR(((k-1)*Nbf+1):(k*Nbf),:) = squeeze(NRcompressed(:,k,:));
% end
% [Unr, Snr, Vnr] = svd(NR,0);
% Ynr = (Snr*Vnr')';
%% Prepare Nonlinear Fitting
vectorize = @(m) reshape(m,[size(m,1)*size(m,2) size(m,3)]);
matricize = @(m,X,Y,Z) reshape(m,[X Y Z]);
%Y = cat(2,Yr,Ynr);
Y = Yr;

fub = 1/(2*mean(diff(T))); % nyquist frequency
flb = 1/(T(end)/(2*pi*0.75)); %require at least 3/4 full period
rub = 1/100;
rlb = 1/(T(end));
gen_rand_beta = @(lb,ub) ((ub-lb).*rand(length(ub),1)+lb);
OscBounds.fub = fub;
OscBounds.flb = flb;
OscBounds.rub = rub;
OscBounds.rlb = rlb;
% OscBounds.sub = 1.0;
% OscBounds.slb = 0.1;
RateBounds.rub = [1/50]; %#ok<NBRAK>
RateBounds.rlb = [1/1000]; %#ok<NBRAK>
% RateBounds.sub = 1.0;
% RateBounds.slb = 0.1;
MaxNOsc = 0;
MinNOsc = 0;
MaxNdynRates = 1;
MinNdynRates = 1;
MaxConstRates = [1/2000 1/10000];
%All the possible models
[ModelCell, Diagnostics] = GenManySoLModels(MaxNOsc, MinNOsc, MaxNdynRates, MinNdynRates, MaxConstRates, OscBounds, RateBounds, T);

%constant fitting parameters
local_tol = 1e-6;
wall_clock_time = 60; %seconds
global_method = NLOPT_G_MLSL_LDS;
Ntests = 1;
candidate_freqs = [90 130 250 500 730 850 970];
%% Fit Loop
parfor k=1:length(ModelCell)
    beta_start = gen_rand_beta(ModelCell{k}.beta_lb,ModelCell{k}.beta_ub);
    Nosc = ModelCell{k}.Nosc;
%     beta_start(2:2:(Nosc*2)) = 1./T2wn(candidate_freqs(randi([1,length(candidate_freqs)],[1,Nosc]))); %#ok<PFBNS>
    ada = ModelCell{k}.modelfun;
    [ModelCell{k}.Fit,ModelCell{k}.Stats,ModelCell{k}.exitcode] =  ...
        varpro_multiway_gradient3(Y,beta_start,ada, ...
                            ModelCell{k}.beta_lb,ModelCell{k}.beta_ub,global_method,wall_clock_time,local_tol,NLOPT_LN_SBPLX);
end
a = clock;
s = [num2str(a(1)) '-' num2str(a(2)) '-' num2str(a(3)) '-' num2str(a(4)) '-' num2str(a(5))];
%save(fullfile(results_path,['fit_results_' s '.mat']),'ModelCell','-v7.3');
%% Refine Model
aicvec = zeros(1,length(ModelCell)); for k=1:length(ModelCell); aicvec(k) = ModelCell{k}.Stats.AICc; end;
[minval, minind] = min(aicvec);
best_model = ModelCell{minind};
uber_model = best_model;
ada = best_model.modelfun;
[uber_model.Fit,uber_model.Stats,uber_model.exitcode] =  ...
        varpro_multiway_gradient3(Y,best_model.Fit.beta,ada, ...
                            best_model.beta_lb,best_model.beta_ub,global_method,wall_clock_time,1e-10,NLOPT_LD_LBFGS);
save(fullfile(results_path,['fit_results_' s '.mat']),'ModelCell','-v7.3');
uber_model.Stats.relerr

%% print some results to screen
fprintf(1,'------- Results ------- \n');

fprintf(1,'Oscillator Freqs\n');
fprintf(1,'%3.3g\n',T2wn(1./uber_model.Fit.beta(2:2:(2*uber_model.Nosc))));
fprintf(1,'Kinetic Lifetimes\n');
fprintf(1,'%3.3g\n',1./uber_model.Fit.beta((2*uber_model.Nosc+1):end));





