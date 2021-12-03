close all; clear all; clc
if ispc;
    base = 'Z:';
else
    base = '/mnt/data';
end
c = 299.792458; %nm/fs
w2wn = @(w) 1E7*w/(2*pi*c);
wn2f = @(wn) (c*wn)./(1E7);
wn2w = @(wn) (2*pi*c*wn)./(1E7);
wn2T = @(wn) 1./((c*wn)./(1E7)); %wave number to period
T2wn = @(tau) 1./((c*tau)./(1E7));
w2l = @(w) 2*pi*c./w;
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','NLFittingNew'));
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','DOPBoxV1_7','DOPBoxV1-7','DOPbox'));
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis'));
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','plot2svg'));
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','NLFitting','cm_and_cb_utilities','cm_and_cb_utilities'));
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','freezeColors_v23_cbfreeze','freezeColors'));
%% Path inputs
data_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','14-03-11','PostProcessingData');
results_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','14-03-11','PostProcessingData','NewFits');
figures_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','14-03-11','FiguresNew');
data_file_name = 'data_stack_r.mat';
optModelName = 'Optimal_Model_Rephasing_and_NonRephasing';
%% create the basis functions for re-expanding the data
if (~exist('all_data_mat_r','var') || ~exist('T2vec','var'));
    data = load(fullfile(data_path,'data_stack_r.mat'));
end
s = fieldnames(data);
Nbf = 60;
dY = size(data.(s{1}),1);
dX = size(data.(s{1}),2);
Bx = dop( dX, Nbf );
By = dop( dY, Nbf );
mean_data = squeeze(mean(data.(s{1}),3));
mean_data = mean_data/max(max(mean_data));
%% Load the model
load(fullfile(results_path,optModelName));

%reconstruct the complex amplitude
cAmps = zeros(size(best_model.Fit.alpha,2),best_model.Nosc);
for k=1:best_model.Nosc;
    cAmps(:,k) = best_model.Fit.alpha(2*(k-1)+1,:) + 1i*best_model.Fit.alpha(2*(k-1)+2,:);
end
peaks = T2wn(1./best_model.Fit.beta(2:2:(2*best_model.Nosc)));
line_width_vec = 1./best_model.Fit.beta(1:2:(2*best_model.Nosc));
rmaps = zeros(dY,dX,best_model.Nosc);
nrmaps = zeros(dY,dX,best_model.Nosc);
for k=1:best_model.Nosc;
    rmaps(:,:,k) = By*reshape(cAmps(1:(Nbf*Nbf),k),[Nbf Nbf])*Bx';
    nrmaps(:,:,k) = By*reshape(cAmps((Nbf*Nbf+1):end,k),[Nbf Nbf])*Bx';
end

% %% plot the rephasing maps
% Ncontours = 50;
% ylim_low = dsearchn(data.(s{3})',2.7);
% xlim_low = dsearchn(data.(s{2})',2.7);
% ylim_high = dsearchn(data.(s{3})',3.05);
% xlim_high = dsearchn(data.(s{2})',3.05);
% max_bkg = max(max(mean_data));
% min_bkg = min(min(mean_data));
% line_diag_eq = @(x,b) x+b;
% for j=1:best_model.Nosc;
%     fig = figure;
%     timg = abs(squeeze(rmaps(:,:,j)));
%     timg = timg/max(max(timg(ylim_low:ylim_high,xlim_high:xlim_low)));
%     max_timg = max(max(abs(timg(ylim_low:ylim_high,xlim_high:xlim_low))));
%     min_timg = min(min(abs(timg(ylim_low:ylim_high,xlim_high:xlim_low))));
%     v = linspace(min_timg,max_timg,Ncontours);
%     contourf(w2wn(data.(s{2}))*1e-3,w2wn(data.(s{3}))*1e-3,timg,v,'linewidth',0.01);
%     a = makeColorMapWStops(linspace(1,Ncontours,6),[[1,1,1]; [1,1,1]; [0,0,1]; [0,1,0]; [1,1,0]; [1,0,0]]);
%     colormap(a)
%     colorbar;
%     hold all;
%     contour(w2wn(data.(s{2}))*1e-3,w2wn(data.(s{3}))*1e-3,abs(mean_data),linspace(min_bkg,0.5*max_bkg,8),'linewidth',0.5,'linecolor','k');
%     line(w2wn(data.(s{2}))*1e-3,w2wn(line_diag_eq(data.(s{2}),0))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
%     line(w2wn(data.(s{2}))*1e-3,w2wn(line_diag_eq(data.(s{2}),-wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
%     line(w2wn(data.(s{2}))*1e-3,w2wn(line_diag_eq(data.(s{2}),+wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
%     line(w2wn(data.(s{2}))*1e-3,w2wn(line_diag_eq(data.(s{2}),-2*wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
%     xlim([14.0 17.0]);
%     ylim([14.705 17.857]);
%     xlabel('Excitation frequency in cm^{-1}x10^{3}','FontSize',14);
%     ylabel('Detection frequency in cm^{-1}x10^{3}','FontSize',14);
%     set(gca,'XTick',[14.0 14.5 15.5 16.0])
%     set(gca,'YTick',[14.8   15.4   16.0   16.6])
%     set(gca,'FontSize',14)
%     axis('square');
%     title(['w2 = ' num2str(peaks(j)) ' with inv linewidth ' num2str(line_width_vec(j)) ' fs']);
%     svgname = fullfile(figures_path,['r_beats_at_' num2str(peaks(j),3) '_wn.svg']);
%     pdfname = fullfile(figures_path,['r_beats_at_' num2str(peaks(j),3) '_wn.pdf']);
%     plot2svg(svgname,fig);
%     export_string = ['inkscape ' svgname  ...
%         ' --export-pdf=' pdfname ...
%         ' --export-dpi=600'];
%     system(export_string);
%     close all;
% end
% %% Plot the non-rephasing maps
% if (~exist('all_data_mat_r','var') || ~exist('T2vec','var'));
%     data = load(fullfile(data_path,'data_stack_r.mat'));
% end
% s = fieldnames(data);
% mean_data = squeeze(mean(data.(s{1}),3));
% mean_data = mean_data/max(max(mean_data));
% 
% Ncontours = 50;
% ylim_low = dsearchn(data.(s{3})',2.7);
% xlim_low = dsearchn(data.(s{2})',2.7);
% ylim_high = dsearchn(data.(s{3})',3.05);
% xlim_high = dsearchn(data.(s{2})',3.05);
% max_bkg = max(max(mean_data));
% min_bkg = min(min(mean_data));
% line_diag_eq = @(x,b) x+b;
% for j=1:best_model.Nosc;
%     fig = figure;
%     timg = abs(squeeze(nrmaps(:,:,j)));
%     timg = timg/max(max(timg(ylim_low:ylim_high,xlim_high:xlim_low)));
%     max_timg = max(max(abs(timg(ylim_low:ylim_high,xlim_high:xlim_low))));
%     min_timg = min(min(abs(timg(ylim_low:ylim_high,xlim_high:xlim_low))));
%     v = linspace(min_timg,max_timg,Ncontours);
%     contourf(w2wn(data.(s{2}))*1e-3,w2wn(data.(s{3}))*1e-3,timg,v,'linewidth',0.01);
%     a = makeColorMapWStops(linspace(1,Ncontours,6),[[1,1,1]; [1,1,1]; [0,0,1]; [0,1,0]; [1,1,0]; [1,0,0]]);
%     colormap(a)
%     colorbar;
%     hold all;
%     contour(w2wn(data.(s{2}))*1e-3,w2wn(data.(s{3}))*1e-3,abs(mean_data),linspace(min_bkg,0.5*max_bkg,8),'linewidth',0.5,'linecolor','k');
%     line(w2wn(data.(s{2}))*1e-3,w2wn(line_diag_eq(data.(s{2}),0))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
%     line(w2wn(data.(s{2}))*1e-3,w2wn(line_diag_eq(data.(s{2}),-wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
%     line(w2wn(data.(s{2}))*1e-3,w2wn(line_diag_eq(data.(s{2}),+wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
%     line(w2wn(data.(s{2}))*1e-3,w2wn(line_diag_eq(data.(s{2}),-2*wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
%     xlim([14.0 17.0]);
%     ylim([14.705 17.857]);
%     xlabel('Excitation frequency in cm^{-1}x10^{3}','FontSize',14);
%     ylabel('Detection frequency in cm^{-1}x10^{3}','FontSize',14);
%     set(gca,'XTick',[14.0 14.5 15.5 16.0])
%     set(gca,'YTick',[14.8   15.4   16.0   16.6])
%     set(gca,'FontSize',14)
%     axis('square');
%     title(['w2 = ' num2str(peaks(j)) ' with inv linewidth ' num2str(line_width_vec(j)) ' fs']);
%     svgname = fullfile(figures_path,['nr_beats_at_' num2str(peaks(j),3) '_wn.svg']);
%     pdfname = fullfile(figures_path,['nr_beats_at_' num2str(peaks(j),3) '_wn.pdf']);
%     plot2svg(svgname,fig);
%     export_string = ['inkscape ' svgname  ...
%         ' --export-pdf=' pdfname ...
%         ' --export-dpi=600'];
%     system(export_string);
%     close all;
% end
% fprintf(1,'All Done with that crap\n');