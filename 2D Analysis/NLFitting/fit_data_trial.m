close all; clear all; clc
addpath('Z:\2D\Frank\Matlab\2D Analysis\DOPBoxV1_7\DOPBoxV1-7\DOPbox')
%% Load the data, compress it, save it.
if (~exist('all_data_mat','var') || ~exist('T2vec','var'));
    load(['Z:\2D\Frank\LyX_Lab_Notebook\LaserLab\13-09-05\PostProcessingData\' 'data_stack']);
end
Nbf = 60;
Trange = 1:185;
rDm = zeros(Nbf,Nbf,length(Trange));
[~,rDm(:,:,1),Bx,By] = dop2DApprox(all_data_mat(:,:,1),Nbf,Nbf);
for k=2:length(Trange);
    [~,rDm(:,:,k),~,~] = dop2DApprox(all_data_mat(:,:,k),Nbf,Nbf);
end
T = -T2vec(Trange);
save(['Z:\2D\Frank\LyX_Lab_Notebook\LaserLab\13-09-05\PostProcessingData\' 'reduced_stack'],'rDm','Bx','By','T');
clear all_data_mat;
%% Load the reduced data if not already in memory
if (~exist('rDm','var') || ~exist('T','var'));
    load(['Z:\2D\Frank\LyX_Lab_Notebook\LaserLab\13-09-05\PostProcessingData\' 'reduced_stack']);
end
vectorize = @(m) reshape(m,[size(m,1)*size(m,2) size(m,3)]);
matricize = @(m,X,Y,Z) reshape(m,[X Y Z]);
Y = vectorize(rDm)';
tic
non_lin_param_lb = [0]';
non_lin_param_ub = [1000]';
alpha0 = (non_lin_param_lb + non_lin_param_ub)/2;
[alpha, c, wresid, wresid_norm, y_est, Regression] =  ...
    varpro_multiway_2(Y,(0*Y+1),alpha0,1,@(x) multi_exp_phi(x,3000,T),non_lin_param_lb,non_lin_param_ub);
toc
plot(Y); hold all; plot(y_est);
figure;
 Z = Y-y_est;
 imagesc(abs(fft(Z,2^12,1)))
 F = matricize(y_est',size(rDm,1),size(rDm,2),size(rDm,3));
 save(['Z:\2D\Frank\LyX_Lab_Notebook\LaserLab\13-09-05\PostProcessingData\' 'reduced_bkg_stack'],'F','Bx','By','T');