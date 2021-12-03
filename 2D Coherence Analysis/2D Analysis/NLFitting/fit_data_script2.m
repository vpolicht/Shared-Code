%% Preliminaries
close all; clear all; clc
if ispc;
    base = 'Z:';
else
    base = '/mnt/data';
end
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','DOPBoxV1_7','DOPBoxV1-7','DOPbox'))
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis'))
c = 299.792458;
w2wn = @(w) 1E7*w/(2*pi*c);
wn2f = @(wn) (c*wn)./(1E7);
wn2w = @(wn) (2*pi*c*wn)./(1E7);
wn2T = @(wn) 1./((c*wn)./(1E7)); %wave number to period
T2wn = @(tau) 1./((c*tau)./(1E7));
%% Load the data, compress it, save it.
if (~exist('all_data_mat','var') || ~exist('T2vec','var'));
    load(fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','PostProcessingData','data_stack.mat'));
end
Nbf = 60;
Trange = 1:185;
rDm = zeros(Nbf,Nbf,length(Trange));
[~,rDm(:,:,1),Bx,By] = dop2DApprox(all_data_mat(:,:,1),Nbf,Nbf);
for k=2:length(Trange);
    [~,rDm(:,:,k),~,~] = dop2DApprox(all_data_mat(:,:,k),Nbf,Nbf);
end
T = -T2vec(Trange);
% save(fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','PostProcessingData','reduced_stack.mat'),'rDm','Bx','By','T');
clear all_data_mat;
%% Load the reduced data, fit it, save it.
if (~exist('rDm','var') || ~exist('T','var'));
    load(fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','PostProcessingData','reduced_stack.mat'));
end
vectorize = @(m) reshape(m,[size(m,1)*size(m,2) size(m,3)]);
matricize = @(m,X,Y,Z) reshape(m,[X Y Z]);
Y = vectorize(rDm)';
tic
% non_lin_param_lb = [0]';
% non_lin_param_ub = [1000]';
% beta = [2000 10000]';
% alpha0 = (non_lin_param_lb + non_lin_param_ub)/2;
% [alpha, c, wresid, wresid_norm, y_est, Regression] =  ...
%     varpro_multiway_2(Y,(0*Y+1),alpha0,1,@(x) multi_exp_phi(x,3000,T),non_lin_param_lb,non_lin_param_ub,1e3,1e-6);
osc_lb = [100     wn2T(120)   100      wn2T(200)    100      wn2T(299)  100     wn2T(400)   100     wn2T(550)   100     wn2T(800)  100     wn2T(900)   100    wn2T(1000)];
osc_ub = [10000   wn2T(50)    10000    wn2T(100)    10000    wn2T(200)  10000   wn2T(300)   10000   wn2T(400)   10000   wn2T(600)  10000   wn2T(800)   10000  wn2T(900)];
var_rates_lb = [0]; %#ok<NBRAK>
var_rates_ub = [1000]; %#ok<NBRAK>
osc_guesses = [88 127 247 339 495 730 854 974];
osc_guesses = osc_guesses + (rand(1,length(osc_guesses))-0.5)*20;
osc_guesses = T2wn(osc_guesses);
rate_guesses = [2000 2000 2000 2000 2000 2000 2000 2000 213];
non_lin_param_lb = [osc_lb var_rates_lb]';
non_lin_param_ub = [osc_ub var_rates_ub]';
static_rates = [2000 10000]';
assert(length(osc_guesses) == length(osc_lb)/2);
Nosc = length(osc_guesses);
alpha0 = 0*non_lin_param_lb;
alpha0(1:2:(Nosc*2-1)) = rate_guesses(1:Nosc);
alpha0(2:2:(Nosc*2)) = osc_guesses;
alpha0((Nosc*2+1):end) = rate_guesses(Nosc+1:end);
Nnonlin = length(non_lin_param_lb);
[alpha, lin_terms, wresid, wresid_norm, y_est, Regression] =  ...
    varpro_multiway_2(Y,(0*Y+1),alpha0,Nnonlin,@(x) multi_exp_with_oscillations_phi(x,static_rates,Nosc,T),non_lin_param_lb,non_lin_param_ub,1e3,1e-6);
toc
T2wn(alpha(2:2:(2*Nosc)))
fprintf(1,'%f\n',alpha(1:2:(2*Nosc-1)));
 Z = Y-y_est;
 F = matricize(y_est',size(rDm,1),size(rDm,2),size(rDm,3));
 save(fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','PostProcessingData','reduced_bkg_stack_with_oscillations.mat'),'F','Bx','By','T');
 save(fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','PostProcessingData','linear_parameters.mat'),'lin_terms','Bx','By','T');
 save(fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','PostProcessingData','nonlinear_parameters.mat'),'alpha');
 
 %% Uncompress the background fit and save it.
% if (~exist('F','var') || ~exist('T','var') || ~exist('Bx','var') || ~exist('By','var'));
%     load(fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','PostProcessingData','reduced_bkg_stack.mat'));
% end
% init_cube = UncompressData( squeeze(F(:,:,1)), Bx, By);
% bkg_cube = zeros(size(init_cube,1),size(init_cube,2),length(T));
% for k=1:length(T);
%     bkg_cube(:,:,k) = UncompressData( squeeze(F(:,:,k)), Bx, By);
% end
% % save(['Z:\2D\Frank\LyX_Lab_Notebook\LaserLab\13-09-05\PostProcessingData\' 'bkg_cube'],'bkg_cube','T');
% clear bkg_cube;

%% Construct the Compressed 3D spectrum with no background

rZt = rDm - F;
filter_f = @(L,a) [0.5; subsref(tukeywin(2*L,a),struct('type','()','subs',{{(L+2):(2*L),1}}))];
filter_m = @(x,y,L,a) permute(repmat(filter_f(L,a),[1 x y]),[2 3 1]);
rZtf = rZt.*filter_m(size(rZt,1),size(rZt,2),size(rZt,3),0.3);
rZf = fft(rZtf,2^10,3);
rZfaxis = w2wn(MakeFourierOmegaAxis(T,2^10));
frobeniusSpectrum = zeros(size(rZf,3),1);
for k=1:size(rZf,3);
    frobeniusSpectrum(k) = norm(rZf(:,:,k),'fro');
end
figure; plot(rZfaxis,frobeniusSpectrum);
saveas(gcf,fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','Figures',['frobenius_spectrum.png']));

%% Construct the amp maps
% line_diag_eq = @(x,b) x+b;
% for k=1:2:(2*Nosc-1);
%     complex_map = ((UncompressData(squeeze(matricize(lin_terms(k+1,:),60,60,1)),Bx,By)+1i*UncompressData(squeeze(matricize(lin_terms(k,:),60,60,1)),Bx,By)));
%     contourf(x_axis_w_roi,y_axis_w_roi,abs(complex_map),10); %set(gca,'YDir','normal');
%     line(x_axis_w_roi,line_diag_eq(x_axis_w_roi,0),'linewidth',2,'Color','w');
%     line(x_axis_w_roi,line_diag_eq(x_axis_w_roi,-wn2w(T2wn(alpha(k+1)))),'linewidth',2,'Color','w');
%     line(x_axis_w_roi,line_diag_eq(x_axis_w_roi,+wn2w(T2wn(alpha(k+1)))),'linewidth',2,'Color','w');
%     xlim([2.7 3.05]);
%     ylim([2.7 3.05]);
%     saveas(gcf,fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','Figures','Beat Analysis',['beats_at_' num2str(T2wn(alpha(k+1))) '_wn_lifetime_' num2str(alpha(k)) 'fs.png']));
%     save(fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','Figures','Beat Analysis',['complex_map_at_' num2str(T2wn(alpha(k+1))) '.mat']),'complex_map');
% end