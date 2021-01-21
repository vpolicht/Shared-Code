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
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','plot2svg'));
c = 299.792458; %nm/fs
w2wn = @(w) 1E7*w/(2*pi*c);
wn2f = @(wn) (c*wn)./(1E7);
wn2w = @(wn) (2*pi*c*wn)./(1E7);
wn2T = @(wn) 1./((c*wn)./(1E7)); %wave number to period
T2wn = @(tau) 1./((c*tau)./(1E7));
w2l = @(w) 2*pi*c./w;
%% Inputs
data_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-07-06','PostProcessingData');
figures_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-07-06','Figures_New');
data_file_name = 'data_stack_nr.mat';
%options
save_compressed_data = 1;
save_uncompressed_bkg = 0;
save_fit_results = 1;
save_amp_maps = 0;
%% Load the data, compress it, save it.
if (~exist('all_data_mat_r','var') || ~exist('T2vec','var'));
    load(fullfile(data_path,'data_stack_nr.mat'));
end
Nbf = 60;
Trange = 1:66;
rDm = zeros(Nbf,Nbf,length(Trange));
[~,rDm(:,:,1),Bx,By] = dop2DApprox(all_data_mat_nr(:,:,1),Nbf,Nbf);
for k=2:length(Trange);
    [~,rDm(:,:,k),~,~] = dop2DApprox(all_data_mat_nr(:,:,k),Nbf,Nbf);
end
T = -T2vec(Trange);
if save_compressed_data;
    save(fullfile(data_path,'reduced_stack_nr.mat'),'rDm','Bx','By','T');
end
% clear all_data_mat_r;
%% Load the reduced data, fit it, save it.
if (~exist('rDm','var') || ~exist('T','var'));
    load(fullfile(data_path,'reduced_stack_nr.mat'));
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
var_rates_lb = [0];
var_rates_ub = [900];
% osc_guesses = [88 127 247 339 495 730 854 974];
% osc_guesses = osc_guesses + (rand(1,length(osc_guesses))-0.5)*20;
% osc_guesses = T2wn(osc_guesses);
osc_guesses = [];
% rate_guesses = [2000 2000 2000 2000 2000 2000 2000 2000 213];
rate_guesses = [600];
non_lin_param_lb = [osc_lb var_rates_lb]';
non_lin_param_ub = [osc_ub var_rates_ub]';
static_rates = [2000]'; %#ok<NBRAK>
assert(length(osc_guesses) == length(osc_lb)/2);
Nosc = length(osc_guesses);
alpha0 = 0*non_lin_param_lb;
% alpha0(1:2:(Nosc*2-1)) = rate_guesses(1:Nosc);
% alpha0(2:2:(Nosc*2)) = osc_guesses;
alpha0((Nosc*2+1):end) = rate_guesses(Nosc+1:end);
N_lin_params = length(non_lin_param_lb)+length(static_rates);
[alpha, lin_terms, wresid, wresid_norm, y_est, Regression] =  ...
    varpro_multiway_2(Y,(0*Y+1),alpha0,N_lin_params,@(x) multi_exp_with_oscillations_phi_and_dphi(x,static_rates,Nosc,T),non_lin_param_lb,non_lin_param_ub,1e3,1e-6);
toc
% T2wn(alpha(2:2:(2*Nosc)))
fprintf(1,'%f\n',alpha(1:2:(2*Nosc-1)));
 Z = Y-y_est;
 F = matricize(y_est',size(rDm,1),size(rDm,2),size(rDm,3));
 R = matricize(wresid',size(rDm,1),size(rDm,2),size(rDm,3));
 save(fullfile(data_path,'reduced_bkg_stack_no_oscillations_non_rephasing.mat'),'F','Bx','By','T');
 save(fullfile(data_path,'linear_parameters_non_rephasing.mat'),'lin_terms','Bx','By','T');
 save(fullfile(data_path,'nonlinear_parameters_non_rephasing.mat'),'alpha');
 
 %% Uncompress the background fit and save it.
if save_uncompressed_bkg;
    if (~exist('F','var') || ~exist('T','var') || ~exist('Bx','var') || ~exist('By','var')); %#ok<UNRCH>
        load(fullfile(data_path,'reduced_bkg_stack_no_oscillations_non_rephasing.mat'));
    end
    init_cube = UncompressData( squeeze(F(:,:,1)), Bx, By);
    bkg_cube = zeros(size(init_cube,1),size(init_cube,2),length(T));
    for k=1:length(T);
        bkg_cube(:,:,k) = UncompressData( squeeze(F(:,:,k)), Bx, By);
    end
    save(['Z:\2D\Frank\LyX_Lab_Notebook\LaserLab\13-09-05\PostProcessingData\' 'bkg_cube_nr'],'bkg_cube_r','T');
    clear bkg_cube_r;
end

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
saveas(gcf,fullfile(figures_path,'frobenius_spectrum_non_rephasing.png'));
dlmwrite(fullfile(figures_path,'frobenius_spectrum_non_rephasing.txt'),[rZfaxis' frobeniusSpectrum],'\t');


%% Construct the fourier amplitude maps at peaks selected from Frobenius Spectrum
if ~exist(fullfile(figures_path,'Beat_Analysis_Fourier'),'dir');
    system(['mkdir' ' ' fullfile(figures_path,'Beat_Analysis_Fourier_NonRephasing')]);
end
make_exc_mat = @(v) abs(triu(repmat(v',[1 length(v)]) - repmat(v',[1 length(v)])' + diag(v)));
kristin_excitons = make_exc_mat([14661 14719 14778 14806 14832 14854 14874 14905 15085]);
novo_excitons = make_exc_mat([14663 14719 14778 14846 14888 14908 14925 14993 15110]);
renger_excitons = make_exc_mat([14599 14749 14815 14859 14970 15198]);
% peaks = [127.1];
peaks = [338.9];
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
%     if peaks(k)<max(max(kristin_excitons - diag(diag(kristin_excitons))));
%         line(w2wn(x_axis_w_roi)*1e-2,ones(size(x_axis_w_roi))*D_excitons(r)*1e-2,'linestyle','-.','linewidth',0.75,'Color','r');
%         line(w2wn(x_axis_w_roi)*1e-2,ones(size(x_axis_w_roi))*D_excitons(c)*1e-2,'linestyle','-.','linewidth',0.75,'Color','r');
%         line(ones(size(y_axis_w_roi))*D_excitons(r)*1e-2,w2wn(y_axis_w_roi)*1e-2,'linestyle','-.','linewidth',0.75,'Color','r');
%         line(ones(size(y_axis_w_roi))*D_excitons(c)*1e-2,w2wn(y_axis_w_roi)*1e-2,'linestyle','-.','linewidth',0.75,'Color','r');
%     end
    if peaks(k)<max(max(renger_excitons - diag(diag(renger_excitons))));
        line(w2wn(x_axis_w_roi)*1e-2,ones(size(x_axis_w_roi))*D_excitons_r(r1)*1e-2,'linestyle','-','linewidth',0.75,'Color','r');
        line(w2wn(x_axis_w_roi)*1e-2,ones(size(x_axis_w_roi))*D_excitons_r(c1)*1e-2,'linestyle','-','linewidth',0.75,'Color','r');
        line(ones(size(y_axis_w_roi))*D_excitons_r(r1)*1e-2,w2wn(y_axis_w_roi)*1e-2,'linestyle','-','linewidth',0.75,'Color','r');
        line(ones(size(y_axis_w_roi))*D_excitons_r(c1)*1e-2,w2wn(y_axis_w_roi)*1e-2,'linestyle','-','linewidth',0.75,'Color','r');
    end
    if peaks(k)<max(max(novo_excitons - diag(diag(novo_excitons))));
        line(w2wn(x_axis_w_roi)*1e-2,ones(size(x_axis_w_roi))*D_excitons_n(r2)*1e-2,'linestyle','--','linewidth',0.75,'Color','r');
        line(w2wn(x_axis_w_roi)*1e-2,ones(size(x_axis_w_roi))*D_excitons_n(c2)*1e-2,'linestyle','--','linewidth',0.75,'Color','r');
        line(ones(size(y_axis_w_roi))*D_excitons_n(r2)*1e-2,w2wn(y_axis_w_roi)*1e-2,'linestyle','--','linewidth',0.75,'Color','r');
        line(ones(size(y_axis_w_roi))*D_excitons_n(c2)*1e-2,w2wn(y_axis_w_roi)*1e-2,'linestyle','--','linewidth',0.75,'Color','r');
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
    title(['Population Frequency: ' num2str(peaks(k),str_precision) ' cm^{-1}'],'FontSize',14);
%     saveas(gcf,fullfile(figures_path,'Beat_Analysis_Fourier',['beats_at_' num2str(peaks(k)) '_wn.eps']));
    plot2svg(fullfile(figures_path,'Beat_Analysis_Fourier_NonRephasing',['beats_at_' num2str(peaks(k)) '_wn.svg']),fig);
    setenv('PATH','C:\Program Files (x86)\Inkscape');
    export_string = ['inkscape ' fullfile(figures_path,'Beat_Analysis_Fourier_NonRephasing',['beats_at_' num2str(peaks(k)) '_wn.svg'])  ...
        ' --export-png=' fullfile(figures_path,'Beat_Analysis_Fourier_NonRephasing',['beats_at_' num2str(peaks(k)) '_wn.png']) ...
        ' --export-dpi=600'];
    system(export_string);
    close all;
end

 
 %% Make a pretty 2D plot
 
%  clear data residual fit;
%  if (~exist('all_data_mat_r','var') || ~exist('T2vec','var'));
%     load(fullfile(data_path,'data_stack.mat'));
%  end
%  if (~exist('mycmap','var'))
%      load(fullfile(base,'2D','Frank','Matlab','2D Analysis','ColorMap42D.mat'));
%  end
% Trange = 1:185;
% T = -T2vec(Trange);
% 
% 
% ind = dsearchn(T,180);
% img_data = squeeze(all_data_mat_r(:,:,ind));
% 
% min_img = min(min(img_data));
% max_img = max(max(img_data));
% tol = mean(diff(linspace(min_img,max_img,100)));
% 
% img_data_pos = img_data;
% img_data_pos(img_data<=0) = NaN;
% img_data_neg = img_data;
% img_data_neg(img_data>0) = NaN;
% 
% v = linspace(min_img,max_img,200);
% v(dsearchn(v',0)) = [];
% fig = figure;
% [~,h1] = contour(w2wn(x_axis_w_roi)*1e-2,w2wn(y_axis_w_roi)*1e-2,img_data_pos,v,'linestyle','-');
% set(fig,'Colormap',mycmap)
% colorbar;
% hold all;
% [~,h2] = contour(w2wn(x_axis_w_roi)*1e-2,w2wn(y_axis_w_roi)*1e-2,img_data_neg,v,'linestyle','-');
% % [~,h3] = contour(w2wn(x_axis_w_roi_dec)*1e-2,w2wn(y_axis_w_roi_dec)*1e-2,img_data_neg_dec,v,'linestyle','-');
% line(w2wn(x_axis_w_roi)*1e-2,w2wn(x_axis_w_roi)*1e-2,'linewidth',0.75,'Color','k');
% xlabel('Excitation Frequency in cm^{-1}x10^{2}','FontSize',14,'FontName','Arial');
% ylabel('Detection Frequency in cm^{-1}x10^{2}','FontSize',14,'FontName','Arial');
% set(gca,'FontSize',14,'FontName','Arial')
% % title('2DES of PSII RC at t_{2} = 170 fs','FontSize',14,'FontName','Arial');
% set(h1,'LineWidth',.3);
% set(h2,'LineWidth',.3);
% 
% set(gcf,'Units','normalized');
% % Create textarrow
% [fx3, fy3] = axescoord2figurecoord([155,150.90],[144,148.40]);
% hAnnotation3 = handle(annotation(fig,'textarrow',fx3,...
%     fy3,'TextEdgeColor','none','FontSize',20,...
%     'FontName','Arial',...
%     'String',{'3'}));
% 
% % Create textarrow
% [fx2, fy2] = axescoord2figurecoord([148,149.10],[143, 147.70]);
% hAnnotation2 = handle(annotation(fig,'textarrow',fx2,...
%     fy2,'TextEdgeColor','none','FontSize',20,...
%     'FontName','Arial',...
%     'String',{'2'}));
% 
% % Create textarrow
% [fx4, fy4] = axescoord2figurecoord([145,148.59],[154, 148.59]);
% hAnnotation4 = handle(annotation(fig,'textarrow',fx4,...
%     fy4,'TextEdgeColor','none','FontSize',20,...
%     'FontName','Arial',...
%     'String',{'4'}));
% 
% % Create textarrow
% [fx1, fy1] = axescoord2figurecoord([144,148.06],[144, 147.18]);
% hAnnotation1 = handle(annotation(fig,'textarrow',fx1,...
%     fy1,'TextEdgeColor','none','FontSize',20,...
%     'FontName','Arial',...
%     'String',{'1'}));
% 
% % hFig = ancestor(Axes,'Figure');
% % pos = hgconvertunits(hFig, get(Axes, 'position'), get(Axes, 'units'), 'Normalized', hFig);
% % pos = pos(1:2)+pos(3:4)/2;%Setting initial annotation position within the figure bounds
% % har = annotation('doublearrow',[pos(1),pos(1)],[pos(2),pos(2)],'Color',[.6 .6 .6]);
% % localPinObject(har);%matlabroot\toolbox\matlab\scribe\@scribe\@scribeobject1D\createPinContextMenu.m\function localPinObject
% % hThis = handle(har);
% % hThis.Pin(1).DataPosition = [X1 Y1 Z1]; 
% % hThis.Pin(1).updateTarget;
% % hThis.Pin(2).DataPosition = [X2 Y2 Z2];
% % hThis.Pin(2).updateTarget;
% 
% % plot2svg(fullfile(figures_path,'Pretty_2D',['2D_at_' num2str(T(ind),3) '_fs_newer.svg']));
% 
% % set(h3,'LineWidth',0.5);
%  