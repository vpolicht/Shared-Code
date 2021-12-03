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
% addpath('C:\Users\Franklin\Documents\Work\Matlab\2D Analysis\DOPBoxV1_7\DOPBoxV1-7\DOPbox')
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis'));
c = 299.792458; %nm/fs
w2wn = @(w) 1E7*w/(2*pi*c);
wn2f = @(wn) (c*wn)./(1E7);
wn2w = @(wn) (2*pi*c*wn)./(1E7);
wn2T = @(wn) 1./((c*wn)./(1E7)); %wave number to period
T2wn = @(tau) 1./((c*tau)./(1E7));
w2l = @(w) 2*pi*c./w;
%% Inputs
data_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','14-03-11','PostProcessingData');
results_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','14-03-11','PostProcessingData','NewFits');
figures_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','14-03-11','Figures');
data_file_name_r = 'data_stack_r.mat';
data_file_name_nr = 'data_stack_nr.mat';
%% Load the rephasing data, compress it.
if (~exist('all_data_mat_r','var') || ~exist('T2vec','var'));
    load(fullfile(data_path,data_file_name_r));
end
f = @(L,flen) [exp(-linspace(0,flen-1,flen)'.^4/(flen/2)^4); zeros(L-flen,1)];
L1 = 150;
L2 = 150;
f1 = f(size(all_data_mat_r,1),L1);
f2 = f(size(all_data_mat_r,2),L2);
Fmat = f1*f2';
Rcompressed = zeros(150,150,size(all_data_mat_r,3));
for k=1:size(all_data_mat_r,3);
Rcompressed(:,:,k)  = idct(idct(dct(dct(squeeze(all_data_mat_r(:,:,k)))')'.*Fmat,150)',150)';
end
R = zeros(L1*L2,size(all_data_mat_r,3));
for k=1:L2
    R(((k-1)*L1+1):(k*L1),:) = squeeze(Rcompressed(:,k,:));
end
[Ur, Sr, Vr] = svd(R,0);
Yr = (Sr*Vr')';
Trange = 1:140;
T = T2vec(Trange);
%% Load the non-rephasing data, compress it.
if (~exist('all_data_mat_nr','var') || ~exist('T2vec','var'));
    load(fullfile(data_path,data_file_name_nr));
end
L1 = 150;
L2 = 150;
f1 = f(size(all_data_mat_r,1),L1);
f2 = f(size(all_data_mat_r,2),L2);
Fmat = f1*f2';
NRcompressed = zeros(150,150,size(all_data_mat_nr,3));
for k=1:size(all_data_mat_nr,3);
NRcompressed(:,:,k)  = idct(idct(dct(dct(squeeze(all_data_mat_nr(:,:,k)))')'.*Fmat,150)',150)';
end
NR = zeros(L1*L2,size(all_data_mat_nr,3));
for k=1:L2
    NR(((k-1)*L1+1):(k*L1),:) = squeeze(NRcompressed(:,k,:));
end
[Unr, Snr, Vnr] = svd(NR,0);
Ynr = (Snr*Vnr')';
%% Prepare Nonlinear Fitting
vectorize = @(m) reshape(m,[size(m,1)*size(m,2) size(m,3)]);
matricize = @(m,X,Y,Z) reshape(m,[X Y Z]);
Y = cat(2,Yr,Ynr);

fub = 1/(2*mean(diff(T))); % nyquist frequency
flb = 1/(T(end)/(pi)); %require at least 1 half period
rub = 1/85;
rlb = 1/(T(end));
gen_rand_beta = @(lb,ub) ((ub-lb).*rand(length(ub),1)+lb);
OscBounds.fub = fub;
OscBounds.flb = flb;
OscBounds.rub = rub;
OscBounds.rlb = rlb;
RateBounds.rub = [1/50 1/300 1/800];
RateBounds.rlb = [1/500 1/1000 1/2000];
MaxNOsc = 10;
MinNOsc = 1;
MaxNdynRates = 1;
MinNdynRates = 1;
MaxConstRates = [1/20000];
%All the possible models
[ModelCell, Diagnostics] = GenManySoLModels(MaxNOsc, MinNOsc, MaxNdynRates, MinNdynRates, MaxConstRates, OscBounds, RateBounds, T);

%constant fitting parameters
local_tol = 1e-6;
wall_clock_time = 30; %seconds
global_method = NLOPT_G_MLSL_LDS;
Ntests = 1;
% candidate_freqs = [90 130 250 340 500 730 850 950];
%% Fit Loop
parfor k=1:length(ModelCell)
    beta_start = gen_rand_beta(ModelCell{k}.beta_lb,ModelCell{k}.beta_ub);
    Nosc = ModelCell{k}.Nosc;
%     beta_start(2:2:(Nosc*2)) = 1./T2wn(candidate_freqs(randi([1,length(candidate_freqs)],[1,Nosc]))); %#ok<PFBNS>
    ada = ModelCell{k}.modelfun;
    [ModelCell{k}.Fit,ModelCell{k}.Stats,ModelCell{k}.exitcode] =  ...
        varpro_multiway_gradient3(Y,beta_start,ada, ...
                            ModelCell{k}.beta_lb,ModelCell{k}.beta_ub,global_method,wall_clock_time,local_tol,NLOPT_LD_LBFGS);
end
a = clock;
s = [num2str(a(1)) '-' num2str(a(2)) '-' num2str(a(3)) '-' num2str(a(4)) '-' num2str(a(5))];
save(fullfile(results_path,['fit_results_' s '.mat']),'ModelCell','-v7.3');




