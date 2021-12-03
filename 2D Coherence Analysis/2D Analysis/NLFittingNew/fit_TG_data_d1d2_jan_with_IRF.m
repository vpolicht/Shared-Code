%% Preliminaries
close all; clear all; clc
if ispc;
    base = 'Z:';
else
    base = '/mnt/data';
end
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','NLFittingNew'));
%addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','NLFitting'));
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
% data_path = fullfile(base,'2D','Data','14-02-06','d1d2_tg-batch13','processed_signal');
% results_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','14-02-06','PostProcessingData');
% figures_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','14-02-06','Figures');
% phase_tg_path = fullfile(base,'2D','Data','14-02-06','d1d2_tg_phasing-batch00','processed_signal','Processed-TG.mat');
% phase_pp_path = fullfile(base,'2D','Data','14-02-06','d1d2_pp-batch00','processed_signal','Processed-TG.mat');
% data_file_name = 'Processed-TG.mat';
data_path = fullfile(base,'2D','Data','14-02-02','d1d2_tg-batch08','processed_signal');
results_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','14-02-02','PostProcessingData');
figures_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','14-02-02','Figures');
phase_tg_path = fullfile(base,'2D','Data','14-02-02','d1d2_tg_phasing-batch01','processed_signal','Processed-TG.mat');
phase_pp_path = fullfile(base,'2D','Data','14-02-02','d1d2_pp-batch00','processed_signal','Processed-TG.mat');
data_file_name = 'Processed-TG.mat';
%% Compute the data phase (or load it)
phase_it = 0;
if phase_it
    load(phase_pp_path);
    load(phase_tg_path);
    PP = PPSignal;
    TG = TGSignal;
    W = vWe;
    [Theta, Alpha, Beta] = MinimizationPhasing(TG, PP, W, [1 1340], [0 -200 -500]', [(2*pi-eps) 200 500]', 10);
    clear TGSignal PPSignal vWe T;
    plot(W,real(Alpha*TG.*exp(1i*Theta))); hold all; plot(W,PP);
    save(fullfile(results_path,'PhasingData.mat'),'Theta','Alpha','Beta');
else
    load(fullfile(results_path,'PhasingData.mat'));
end
%% Load the data, truncate/average as desired.
if (~exist('TGSignal','var') || ~exist('T','var'));
    load(fullfile(data_path,data_file_name));
end
T = -T;
minT = 50; %fs
maxT = 830;
minTind = dsearchn(T,minT);
maxTind = dsearchn(T,maxT);
Trange = minTind:maxTind;
T = T(Trange);
Y = (squeeze(mean(TGSignal(:,Trange,:),3)))';
if exist('Theta','var')
    Corr = (Alpha*repmat(exp(1i*Theta),[1 size(Y,1)]))';
    Y = real(Y.*Corr);
else
    Y = abs(Y);
end
% L = 3;
% [u,s,v] = svd(Y);
% Y = u(:,1:L);

% plot_ind = dsearchn(vWe',2.35);
% figure; plot(T,Y(:,plot_ind))

% Prepare Nonlinear Fitting
vectorize = @(m) reshape(m,[size(m,1)*size(m,2) size(m,3)]);
matricize = @(m,X,Y,Z) reshape(m,[X Y Z]);

fub = 1/(2*mean(diff(T))); % nyquist frequency
flb = 1/(T(end)/(pi)); %require at least 1 half period
rub = 1/200;
rlb = 1/(T(end));
gen_rand_beta = @(lb,ub) ((ub-lb).*rand(length(ub),1)+lb);
OscBounds.fub = fub;
OscBounds.flb = flb;
OscBounds.rub = rub;
OscBounds.rlb = rlb;
IRFBounds.ub = [10 200];
IRFBounds.lb = [-10 30];
% RateBounds.rub = [1/200  1/1500  1/10000 1/30000];
% RateBounds.rlb = [1/3000 1/15000 1/50000 1/130000];
RateBounds.rub = [1/50  1/300  1/1500];
RateBounds.rlb = [1/400 1/2000 1/5000];
MaxNOsc = 0;
MinNOsc = 0;
MaxNdynRates = 2;
MinNdynRates = 1;
MaxConstRates = [1/800 1/8000 1/600000];
%All the possible models
[ModelCell, Diagnostics] = GenManySoLModelsWithIRF(MaxNOsc, MinNOsc, MaxNdynRates, MinNdynRates, MaxConstRates, OscBounds, RateBounds, IRFBounds, T);
% [ModelCell, Diagnostics] = GenManySoLModels(MaxNOsc, MinNOsc, MaxNdynRates, MinNdynRates, MaxConstRates, OscBounds, RateBounds, T);

%constant fitting parameters
local_tol = 1e-6;
wall_clock_time = 20; %seconds
global_method = NLOPT_G_MLSL_LDS;
Ntests = 1;
%% Fit Loop
parfor k=1:length(ModelCell)
    beta_start = gen_rand_beta(ModelCell{k}.beta_lb,ModelCell{k}.beta_ub);
    Nosc = ModelCell{k}.Nosc;
    ada = ModelCell{k}.modelfun;
    [ModelCell{k}.Fit,ModelCell{k}.Stats,ModelCell{k}.exitcode] =  ...
        varpro_multiway_gradient3(Y,beta_start,ada, ...
                            ModelCell{k}.beta_lb,ModelCell{k}.beta_ub,global_method,wall_clock_time,local_tol,NLOPT_LD_LBFGS);
end
aicvec = zeros(1,length(ModelCell)); for k=1:length(ModelCell); aicvec(k) = ModelCell{k}.Stats.AICc; end;
[minval, minind] = min(aicvec);
M = ModelCell{minind};
fprintf(1,'---Fit Results---\n');
fprintf(1,'Life Times:\n');
for k=1:(length(M.Fit.beta))
    fprintf(1,'%g\n',1./M.Fit.beta(k));
end
% fprintf(1,'Instrument Response Fit\n');
% fprintf(1,'%g\n',M.Fit.beta(1));
% fprintf(1,'%g\n',M.Fit.beta(2));
a = clock;
s = [num2str(a(1)) '-' num2str(a(2)) '-' num2str(a(3)) '-' num2str(a(4)) '-' num2str(a(5))];
save(fullfile(results_path,['fit_results_' s '.mat']),'M','-v7.3');
fprintf(1,'Results saved to %s',fullfile(results_path,['fit_results_' s '.mat']));




