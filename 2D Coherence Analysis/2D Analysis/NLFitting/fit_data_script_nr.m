%% Preliminaries
close all; clear all; clc
if ispc;
    base = 'Z:';
else
    base = '/mnt/data';
end
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
figures_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','Figures');
data_file_name = 'data_stack_nr.mat';
%options
save_compressed_data = 1;
save_uncompressed_bkg = 0;
save_fit_results = 1;
save_amp_maps = 0;
%% Load the data, compress it, save it.
if (~exist('all_data_mat_nr','var') || ~exist('T2vec','var'));
    load(fullfile(data_path,'data_stack_nr.mat'));
end
Nbf = 60;
Trange = 1:185;
rDm = zeros(Nbf,Nbf,length(Trange));
[~,rDm(:,:,1),Bx,By] = dop2DApprox(all_data_mat_nr(:,:,1),Nbf,Nbf);
for k=2:length(Trange);
    [~,rDm(:,:,k),~,~] = dop2DApprox(all_data_mat_nr(:,:,k),Nbf,Nbf);
end
T = -T2vec(Trange);
if save_compressed_data;
    save(fullfile(data_path,'reduced_stack_nr.mat'),'rDm','Bx','By','T');
end
% clear all_data_mat_nr;
%% Load the reduced data, fit it, save it.
if (~exist('rDm','var') || ~exist('T','var'));
    load(fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','PostProcessingData','reduced_stack_nr.mat'));
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
% osc_lb = [100     wn2T(120)   100      wn2T(200)    100      wn2T(299)  100     wn2T(400)   100     wn2T(550)   100     wn2T(800)  100     wn2T(900)   100    wn2T(1000)];
% osc_ub = [10000   wn2T(50)    10000    wn2T(100)    10000    wn2T(200)  10000   wn2T(300)   10000   wn2T(400)   10000   wn2T(600)  10000   wn2T(800)   10000  wn2T(900)];
osc_lb = [];
osc_ub = [];
var_rates_lb = [0]; %#ok<NBRAK>
var_rates_ub = [1000]; %#ok<NBRAK>
% osc_guesses = [88 127 247 339 495 730 854 974];
% osc_guesses = osc_guesses + (rand(1,length(osc_guesses))-0.5)*20;
% osc_guesses = T2wn(osc_guesses);
osc_guesses = [];
% rate_guesses = [2000 2000 2000 2000 2000 2000 2000 2000 213];
rate_guesses = [213];%#ok<NBRAK>
non_lin_param_lb = [osc_lb var_rates_lb]';
non_lin_param_ub = [osc_ub var_rates_ub]';
static_rates = [2000 10000]';
assert(length(osc_guesses) == length(osc_lb)/2);
Nosc = length(osc_guesses);
alpha0 = 0*non_lin_param_lb;
% alpha0(1:2:(Nosc*2-1)) = rate_guesses(1:Nosc);
% alpha0(2:2:(Nosc*2)) = osc_guesses;
alpha0((Nosc*2+1):end) = rate_guesses(Nosc+1:end);
N_lin_params = length(non_lin_param_lb)+length(static_rates);
[alpha, lin_terms, wresid, wresid_norm, y_est, Regression] =  ...
    varpro_multiway_2(Y,(0*Y+1),alpha0,N_lin_params,@(x) multi_exp_with_oscillations_phi_and_dphi(x,static_rates,Nosc,T),non_lin_param_lb,non_lin_param_ub,1e4,1e-6);
toc
% T2wn(alpha(2:2:(2*Nosc)))
fprintf(1,'%f\n',alpha(1:2:(2*Nosc-1)));
 Z = Y-y_est;
 F = matricize(y_est',size(rDm,1),size(rDm,2),size(rDm,3));
 R = matricize(wresid',size(rDm,1),size(rDm,2),size(rDm,3));
 save(fullfile(data_path,'reduced_bkg_stack_nr_no_oscillations.mat'),'F','Bx','By','T');
 save(fullfile(data_path,'linear_parameters_nr.mat'),'lin_terms','Bx','By','T');
 save(fullfile(data_path,'nonlinear_parameters_nr.mat'),'alpha');
 
 %% Uncompress the background fit and save it.
if save_uncompressed_bkg;
    if (~exist('F','var') || ~exist('T','var') || ~exist('Bx','var') || ~exist('By','var')); %#ok<UNRCH>
        load(fullfile(data_path,'reduced_bkg_stack_nr.mat'));
    end
    init_cube = UncompressData( squeeze(F(:,:,1)), Bx, By);
    bkg_cube_nr = zeros(size(init_cube,1),size(init_cube,2),length(T));
    for k=1:length(T);
        bkg_cube_nr(:,:,k) = UncompressData( squeeze(F(:,:,k)), Bx, By);
    end
    save(['Z:\2D\Frank\LyX_Lab_Notebook\LaserLab\13-09-05\PostProcessingData\' 'bkg_cube_nr'],'bkg_cube_nr','T');
    clear bkg_cube;
end

%% Construct the Compressed 3D spectrum with no background

if ~exist('F','var');
    load(fullfile(data_path,'reduced_bkg_stack_nr_no_oscillations.mat'));
end
if ~exist('rDm','var')
    load(fullfile(data_path,'reduced_stack_nr.mat'));
end

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
saveas(gcf,fullfile(figures_path,'frobenius_spectrum_nr.png'));
dlmwrite(fullfile(figures_path,'frobenius_spectrum_nr.txt'),[rZfaxis' frobeniusSpectrum],'\t');

%% Construct the amp maps from fit
% if ~exist(fullfile(figures_path,'Beat Analysis'),'dir');
%     system([mkdir ' ' fullfile(figures_path,'Beat Analysis\')]);
% end
% if save_amp_maps;
%     line_diag_eq = @(x,b) x+b; %#ok<UNRCH>
%     for k=1:2:(2*Nosc-1);
%         complex_map = ((UncompressData(squeeze(matricize(lin_terms(k+1,:),60,60,1)),Bx,By)+1i*UncompressData(squeeze(matricize(lin_terms(k,:),60,60,1)),Bx,By)));
%         contourf(x_axis_w_roi,y_axis_w_roi,abs(complex_map),10); %set(gca,'YDir','normal');
%         line(x_axis_w_roi,line_diag_eq(x_axis_w_roi,0),'linewidth',2,'Color','w');
%         line(x_axis_w_roi,line_diag_eq(x_axis_w_roi,-wn2w(T2wn(alpha(k+1)))),'linewidth',2,'Color','w');
%         line(x_axis_w_roi,line_diag_eq(x_axis_w_roi,+wn2w(T2wn(alpha(k+1)))),'linewidth',2,'Color','w');
%         xlim([2.7 3.05]);
%         ylim([2.7 3.05]);
%         saveas(gcf,fullfile(figures_path,'Beat Analysis',['beats_at_' num2str(T2wn(alpha(k+1))) '_wn_lifetime_' num2str(alpha(k)) 'fs.png']));
%         save(fullfile(figures_path,'Beat Analysis',['complex_map_at_' num2str(T2wn(alpha(k+1))) '.mat']),'complex_map');
%     end
% end

%% Construct the fourier amplitude maps at peaks selected from Frobenius Spectrum
if ~exist(fullfile(figures_path,'Beat_Analysis_Fourier'),'dir');
    system(['mkdir' ' ' fullfile(figures_path,'Beat_Analysis_Fourier_NonRephasing')]);
end
if ~exist('y_axis_w_roi','var') || ~exist('x_axis_w_roi','var');
    load(fullfile(data_path,'data_stack_nr.mat'));
    clear all_data_mat;
end
make_exc_mat = @(v) abs(triu(repmat(v',[1 length(v)]) - repmat(v',[1 length(v)])' + diag(v)));
kristin_excitons = make_exc_mat([14661 14719 14778 14806 14832 14854 14874 14905 15085]);
novo_excitons = make_exc_mat([14663 14719 14778 14846 14888 14908 14925 14993 15110]);
renger_excitons = make_exc_mat([14599 14749 14815 14859 14970 15198]);
peaks = [91.24 127.1 250.9 338.9 729.9 853.7 974.3];
line_diag_eq = @(x,b) x+b;
ylim_low = dsearchn(y_axis_w_roi',2.7);
xlim_low = dsearchn(x_axis_w_roi',2.7);
ylim_high = dsearchn(y_axis_w_roi',3.05);
xlim_high = dsearchn(x_axis_w_roi',3.05);
bkg_img = UncompressData(squeeze(mean(rDm,3)),Bx,By);
max_bkg = max(max(bkg_img));
min_bkg = min(min(bkg_img));
D_excitons = diag(kristin_excitons);
D_excitons_r = diag(renger_excitons);
D_excitons_n = diag(novo_excitons);
for k=1:length(peaks);
    fig = figure;
    [r,c] = FindInMatrix(kristin_excitons,peaks(k));
    [r1,c1] = FindInMatrix(renger_excitons,peaks(k));
    [r2,c2] = FindInMatrix(novo_excitons,peaks(k));
    tind = dsearchn(rZfaxis',peaks(k));
    timg = UncompressData(squeeze(rZf(:,:,tind)),Bx,By);
    max_timg = max(max(abs(timg(ylim_low:ylim_high,xlim_low:xlim_high))));
    v = linspace(0.25*max_timg,max_timg,10);
    contourf(w2wn(x_axis_w_roi)*1e-2,w2wn(y_axis_w_roi)*1e-2,abs(timg),v,'linewidth',0.01);
    colormap('jet');
%     colorbar;
%     cbfreeze;
%     freezeColors;
    hold all;
    contour(w2wn(x_axis_w_roi)*1e-2,w2wn(y_axis_w_roi)*1e-2,abs(bkg_img),linspace(min_bkg,0.5*max_bkg,8));
%     colormap('gray');
%     freezeColors;
    line(w2wn(x_axis_w_roi)*1e-2,w2wn(line_diag_eq(x_axis_w_roi,0))*1e-2,'linestyle','--','linewidth',1,'Color','k');
    line(w2wn(x_axis_w_roi)*1e-2,w2wn(line_diag_eq(x_axis_w_roi,-wn2w(peaks(k))))*1e-2,'linestyle','--','linewidth',1,'Color','k');
    line(w2wn(x_axis_w_roi)*1e-2,w2wn(line_diag_eq(x_axis_w_roi,+wn2w(peaks(k))))*1e-2,'linestyle','--','linewidth',1,'Color','k');
    line(w2wn(x_axis_w_roi)*1e-2,w2wn(line_diag_eq(x_axis_w_roi,-2*wn2w(peaks(k))))*1e-2,'linestyle','--','linewidth',1,'Color','k');
    if (peaks(k)<max(max(kristin_excitons - diag(diag(kristin_excitons)))) && abs(kristin_excitons(r,c) - peaks(k))<10);
        line(w2wn(x_axis_w_roi)*1e-2,ones(size(x_axis_w_roi))*D_excitons(r)*1e-2,'linestyle','-.','linewidth',0.75,'Color','r');
        line(w2wn(x_axis_w_roi)*1e-2,ones(size(x_axis_w_roi))*D_excitons(c)*1e-2,'linestyle','-.','linewidth',0.75,'Color','r');
        line(ones(size(y_axis_w_roi))*D_excitons(r)*1e-2,w2wn(y_axis_w_roi)*1e-2,'linestyle','-.','linewidth',0.75,'Color','r');
        line(ones(size(y_axis_w_roi))*D_excitons(c)*1e-2,w2wn(y_axis_w_roi)*1e-2,'linestyle','-.','linewidth',0.75,'Color','r');
    end
%     if peaks(k)<max(max(renger_excitons - diag(diag(renger_excitons))));
%         line(w2wn(x_axis_w_roi)*1e-2,ones(size(x_axis_w_roi))*D_excitons_r(r1)*1e-2,'linestyle','--','linewidth',0.75,'Color','r');
%         line(w2wn(x_axis_w_roi)*1e-2,ones(size(x_axis_w_roi))*D_excitons_r(c1)*1e-2,'linestyle','--','linewidth',0.75,'Color','r');
%         line(ones(size(y_axis_w_roi))*D_excitons_r(r1)*1e-2,w2wn(y_axis_w_roi)*1e-2,'linestyle','--','linewidth',0.75,'Color','r');
%         line(ones(size(y_axis_w_roi))*D_excitons_r(c1)*1e-2,w2wn(y_axis_w_roi)*1e-2,'linestyle','--','linewidth',0.75,'Color','r');
%     end
    if (peaks(k)<max(max(novo_excitons - diag(diag(novo_excitons)))) && abs(novo_excitons(r2,c2) - peaks(k))<10);
        line(w2wn(x_axis_w_roi)*1e-2,ones(size(x_axis_w_roi))*D_excitons_n(r2)*1e-2,'linestyle','-','linewidth',0.75,'Color','r');
        line(w2wn(x_axis_w_roi)*1e-2,ones(size(x_axis_w_roi))*D_excitons_n(c2)*1e-2,'linestyle','-','linewidth',0.75,'Color','r');
        line(ones(size(y_axis_w_roi))*D_excitons_n(r2)*1e-2,w2wn(y_axis_w_roi)*1e-2,'linestyle','-','linewidth',0.75,'Color','r');
        line(ones(size(y_axis_w_roi))*D_excitons_n(c2)*1e-2,w2wn(y_axis_w_roi)*1e-2,'linestyle','-','linewidth',0.75,'Color','r');
    end
    xlim(w2wn([2.7 3.05])*1e-2);
    ylim(w2wn([2.7 3.05])*1e-2);
    xlabel('Excitation Frequency in cm^{-1}x10^{2}','FontSize',14);
    ylabel('Detection Frequency in cm^{-1}x10^{2}','FontSize',14);
    set(gca,'FontSize',14)
    if peaks(k)>100;
        str_precision = 3;
    else
        str_precision = 2;
    end
    title(['\omega_{2}: ' num2str(peaks(k),str_precision) ' cm^{-1}'],'FontSize',14);
%     saveas(gcf,fullfile(figures_path,'Beat_Analysis_Fourier',['beats_at_' num2str(peaks(k)) '_wn.eps']));
    plot2svg(fullfile(figures_path,'Beat_Analysis_Fourier_NonRephasing',['beats_at_' num2str(peaks(k)) '_wn.svg']),fig);
%     setenv('PATH','C:\Program Files (x86)\Inkscape');
    export_string = ['inkscape ' fullfile(figures_path,'Beat_Analysis_Fourier_NonRephasing',['beats_at_' num2str(peaks(k)) '_wn.svg'])  ...
        ' --export-png=' fullfile(figures_path,'Beat_Analysis_Fourier_NonRephasing',['beats_at_' num2str(peaks(k)) '_wn.png']) ...
        ' --export-dpi=600'];
    system(export_string);
    close all;
end

%% Create Pixel Traces at specified points
% if ~exist(fullfile(figures_path,'Pixel_Traces'),'dir');
%     system(['mkdir' ' ' fullfile(figures_path,'Pixel_Traces')]);
% end
% 
% a90 = [678.06,682.314]; %D- 90wn
% b90 = [678.6773,674.3371]; %D+ 90wn
% c90 = [678.7,678.7]; %D 90wn
% 
% a127 = [677.45,683.07]; %D-
% b127 = [678.88,672.973]; %D+
% c127 = [676.8,676.8]; %D
% 
% a250 = [666.19,677.61]; %D-
% b250 = [674.4,663.2]; %D+
% d250 = [660.5,683.2]; %D--
% 
% a339 = [668.37,683.824]; %D-
% b339 = [673.8,658.7]; %D+
% c339 = [676.2,676.2]; %D
% d339 = [656.3,686.8]; %D--
% 
% a730 = [644.23,676.023]; %D-
% b730 = [675.2,643.5]; %D+
% c730 = [676.8,676.8]; %D
% d730 = [619.87,681.5617]; %D--
% 
% a853 = [641.1,678.2]; %D-
% b853 = [672.8,636.2]; %D+
% c853 = [676.6,676.6]; %D
% 
% a974 = [632.6,674.2]; %D-
% b974 = [667.8,627]; %D+
% c974 = [674.8,674.8]; %D
% 
% if (~exist('data','var') || ~exist('residual','var') || ~exist('fit','var'));
% data = zeros(1340,672,size(rDm,3));
% residual = zeros(1340,672,size(rDm,3));
% fit = zeros(1340,672,size(rDm,3));
%     for k=1:size(rDm,3);
%         data(:,:,k) = UncompressData(squeeze(rDm(:,:,k)),Bx,By);
%         residual(:,:,k) = UncompressData(squeeze(rZt(:,:,k)),Bx,By);
%         fit(:,:,k) = UncompressData(squeeze(F(:,:,k)),Bx,By);
%     end
% end
% 
% find_point = @(d,x) squeeze(d(dsearchn(y_axis_w_roi',w2l(x(1,2))),dsearchn(x_axis_w_roi',w2l(x(1,1))),:));
% 
% points = [a90; b90; c90; a127; b127; c127; a250; b250; d250; a339; b339; c339; d339; a730; b730; c730; d730; a853; b853; c853; a974; b974; c974];
% for k=1:size(points,1);
%     dataT = find_point(data,squeeze(points(k,:)));
%     dlmwrite(fullfile(figures_path,'Pixel_Traces',['time_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.txt']),[T dataT],'\t');
%     fitT = find_point(fit,points(k,:));
%     dlmwrite(fullfile(figures_path,'Pixel_Traces',['time_domain_fit_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.txt']),[T fitT],'\t');
%     dataW = abs(fft(find_point(residual,points(k,:)),2^10));
%     dlmwrite(fullfile(figures_path,'Pixel_Traces',['freq_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.txt']),[rZfaxis' dataW],'\t');
%     figure; plot(T,dataT); hold all;
%     plot(T,fitT); hold off;
%     saveas(gcf,fullfile(figures_path,'Pixel_Traces',['time_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.png']));
%     figure; plot(rZfaxis,dataW);
%     xlim([0 1000]);
%     saveas(gcf,fullfile(figures_path,'Pixel_Traces',['freq_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.png']));
%     
% end
% close all;

%% Create Pixel traces on boxes
% filter_f = @(L,a) [0.5; subsref(tukeywin(2*L,a),struct('type','()','subs',{{(L+2):(2*L),1}}))];
% 
% if ~exist(fullfile(figures_path,'Pixel_Traces_On_Box'),'dir');
%     system(['mkdir' ' ' fullfile(figures_path,'Pixel_Traces_On_Box')]);
% end
% 
% if (~exist('data','var') || ~exist('residual','var') || ~exist('fit','var'));
% data = zeros(1340,672,size(rDm,3));
% residual = zeros(1340,672,size(rDm,3));
% fit = zeros(1340,672,size(rDm,3));
%     for k=1:size(rDm,3);
%         data(:,:,k) = UncompressData(squeeze(rDm(:,:,k)),Bx,By);
%         residual(:,:,k) = UncompressData(squeeze(rZt(:,:,k)),Bx,By);
%         fit(:,:,k) = UncompressData(squeeze(F(:,:,k)),Bx,By);
%     end
% end
% 
% lower_box_corners = ...
%     [15829.7018, 14794.7042; ...
%      15598.2049, 14746.1591; ...
%      15513.6195, 14785.4575; ...
%      15095.1442, 14771.5875; ...
%      14997.2032, 14753.4058; ...
%      14899.2622, 14766.9641; ...
%      14801.3210, 14720.7307];
%  
%  upper_diag_box_corners = [lower_box_corners(:,1) lower_box_corners(:,1)];
%  lower_diag_box_corners = [lower_box_corners(:,2) lower_box_corners(:,2)];
%  upper_box_corners = [lower_box_corners(:,2) lower_box_corners(:,1)];
%  
%  find_point = @(d,x) RetrieveAveragedPixelTrace(d,w2wn(x_axis_w_roi),w2wn(y_axis_w_roi),x,15);
%  
%  points = upper_diag_box_corners;
%  for k=1:size(lower_box_corners,1)
%     dataT = find_point(data,squeeze(points(k,:)));
%     dlmwrite(fullfile(figures_path,'Pixel_Traces_On_Box',['time_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.txt']),[T dataT],'\t');
%     fitT = find_point(fit,points(k,:));
%     dlmwrite(fullfile(figures_path,'Pixel_Traces_On_Box',['time_domain_fit_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.txt']),[T fitT],'\t');
%     dataW = abs(fft(find_point(residual,points(k,:)).*filter_f(185,0.3),2^10));
%     dlmwrite(fullfile(figures_path,'Pixel_Traces_On_Box',['freq_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.txt']),[rZfaxis' dataW],'\t');
%     figure; plot(T,dataT); hold all;
%     plot(T,fitT); hold off;
%     saveas(gcf,fullfile(figures_path,'Pixel_Traces_On_Box',['time_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.png']));
%     figure; plot(rZfaxis,dataW);
%     xlim([0 1000]);
%     saveas(gcf,fullfile(figures_path,'Pixel_Traces_On_Box',['freq_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.png']));
%  end
%  
%   points = lower_diag_box_corners;
%  for k=1:size(lower_box_corners,1)
%     dataT = find_point(data,squeeze(points(k,:)));
%     dlmwrite(fullfile(figures_path,'Pixel_Traces_On_Box',['time_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.txt']),[T dataT],'\t');
%     fitT = find_point(fit,points(k,:));
%     dlmwrite(fullfile(figures_path,'Pixel_Traces_On_Box',['time_domain_fit_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.txt']),[T fitT],'\t');
%     dataW = abs(fft(find_point(residual,points(k,:)).*filter_f(185,0.3),2^10));
%     dlmwrite(fullfile(figures_path,'Pixel_Traces_On_Box',['freq_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.txt']),[rZfaxis' dataW],'\t');
%     figure; plot(T,dataT); hold all;
%     plot(T,fitT); hold off;
%     saveas(gcf,fullfile(figures_path,'Pixel_Traces_On_Box',['time_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.png']));
%     figure; plot(rZfaxis,dataW);
%     xlim([0 1000]);
%     saveas(gcf,fullfile(figures_path,'Pixel_Traces_On_Box',['freq_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.png']));
%  end
%  
%   points = lower_box_corners;
%  for k=1:size(lower_box_corners,1)
%     dataT = find_point(data,squeeze(points(k,:)));
%     dlmwrite(fullfile(figures_path,'Pixel_Traces_On_Box',['time_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.txt']),[T dataT],'\t');
%     fitT = find_point(fit,points(k,:));
%     dlmwrite(fullfile(figures_path,'Pixel_Traces_On_Box',['time_domain_fit_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.txt']),[T fitT],'\t');
%     dataW = abs(fft(find_point(residual,points(k,:)).*filter_f(185,0.3),2^10));
%     dlmwrite(fullfile(figures_path,'Pixel_Traces_On_Box',['freq_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.txt']),[rZfaxis' dataW],'\t');
%     figure; plot(T,dataT); hold all;
%     plot(T,fitT); hold off;
%     saveas(gcf,fullfile(figures_path,'Pixel_Traces_On_Box',['time_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.png']));
%     figure; plot(rZfaxis,dataW);
%     xlim([0 1000]);
%     saveas(gcf,fullfile(figures_path,'Pixel_Traces_On_Box',['freq_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.png']));
%  end
%  
%   points = upper_box_corners;
%  for k=1:size(lower_box_corners,1)
%     dataT = find_point(data,squeeze(points(k,:)));
%     dlmwrite(fullfile(figures_path,'Pixel_Traces_On_Box',['time_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.txt']),[T dataT],'\t');
%     fitT = find_point(fit,points(k,:));
%     dlmwrite(fullfile(figures_path,'Pixel_Traces_On_Box',['time_domain_fit_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.txt']),[T fitT],'\t');
%     dataW = abs(fft(find_point(residual,points(k,:)).*filter_f(185,0.3),2^10));
%     dlmwrite(fullfile(figures_path,'Pixel_Traces_On_Box',['freq_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.txt']),[rZfaxis' dataW],'\t');
%     figure; plot(T,dataT); hold all;
%     plot(T,fitT); hold off;
%     saveas(gcf,fullfile(figures_path,'Pixel_Traces_On_Box',['time_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.png']));
%     figure; plot(rZfaxis,dataW);
%     xlim([0 1000]);
%     saveas(gcf,fullfile(figures_path,'Pixel_Traces_On_Box',['freq_domain_at_' num2str(points(k,1)) '_nm_' num2str(points(k,2)) '_nm.png']));
%  end
%  
%  close all;
%  
 %% Make a pretty 2D plot
 
 clear data residual fit;
 if (~exist('all_data_mat_nr','var') || ~exist('T2vec','var'));
    load(fullfile(data_path,'data_stack_nr.mat'));
 end
ind = dsearchn(T,180);
img_data = squeeze(all_data_mat_nr(:,:,ind));
Trange = 1:185;
T = -T2vec(Trange);
min_img = min(min(img_data));
max_img = max(max(img_data));
tol = mean(diff(linspace(min_img,max_img,100)));


img_data_pos = img_data;
img_data_pos(img_data<=0) = NaN;
img_data_neg = img_data;
img_data_neg(img_data>0-tol) = NaN;

v = linspace(min_img,max_img,100);
v0 = dsearchn(v',0);
v((v0):1:(v0)) = [];
dec_fac = 1;
img_data_neg_dec = img_data_neg(1:dec_fac:end,1:dec_fac:end);
x_axis_w_roi_dec = x_axis_w_roi(1:dec_fac:end);
y_axis_w_roi_dec = y_axis_w_roi(1:dec_fac:end);

figure('units','inches','position',[1 1 3 5]);
[~,h1] = contour(w2wn(x_axis_w_roi)*1e-2,w2wn(y_axis_w_roi)*1e-2,img_data_pos,v,'linestyle','-');
colormap('jet');
hold all;
[~,h2] = contourf(w2wn(x_axis_w_roi_dec)*1e-2,w2wn(y_axis_w_roi_dec)*1e-2,img_data_neg_dec,v,'linestyle','-');
set(gca,'YDir','normal');
% [~,h3] = contour(w2wn(x_axis_w_roi_dec)*1e-2,w2wn(y_axis_w_roi_dec)*1e-2,img_data_neg_dec,v,'linestyle','-');
axis square;
line(w2wn(x_axis_w_roi)*1e-2,w2wn(x_axis_w_roi)*1e-2,'linewidth',0.75,'Color','k');
colormap('jet');
xlabel('Excitation Frequency in cm^{-1}x10^{2}','FontSize',9);
ylabel('Detection Frequency in cm^{-1}x10^{2}','FontSize',9);
set(gca,'FontSize',7)
title('2DES of d1d2 at t_{2} = 170 fs','FontSize',9);
set(h1,'LineWidth',.3);
set(h2,'LineWidth',.3);
xlim([142 167]);
% set(h3,'LineWidth',0.5);
 