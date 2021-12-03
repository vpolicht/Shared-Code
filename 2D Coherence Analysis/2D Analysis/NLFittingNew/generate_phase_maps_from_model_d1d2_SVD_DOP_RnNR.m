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
addpath(fullfile(base,'2D','Frank','Matlab','2D Analysis','2DPhaseUnwrapping','2D_SRNCP_phase_unwrapper_without_a_mask'));
%% Path inputs
data_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','PostProcessingData');
results_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','PostProcessingData','NewFits');
figures_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','HighFreqZoomedRectangular');
data_file_name_r = 'data_stack_r.mat';
data_file_name_nr = 'data_stack_nr.mat';
optModelName = 'uber_fit';
%% Load the rephasing data, compress it.

load(fullfile(results_path,optModelName));
best_model = uber_model;
peaks = T2wn(1./best_model.Fit.beta(2:2:(2*best_model.Nosc)));
line_width_vec = 1./best_model.Fit.beta(1:2:(2*best_model.Nosc));

load(fullfile(results_path,'ComplexAmplitudes.mat'));
PhaseM = angle(Cm);
% PhaseM = zeros(size(Cm));
% for k=1:best_model.Nosc; 
%     PhaseM(:,:,k) = Miguel_2D_unwrapper(single(angle(squeeze(Cm(:,:,k)))));
% end;

Ncontours = 50;
ylim_low = dsearchn(yaxis',wn2w(14000));
xlim_low = dsearchn(xaxis',wn2w(14000));
ylim_high = dsearchn(yaxis',wn2w(16000));
xlim_high = dsearchn(xaxis',wn2w(17000));
line_diag_eq = @(x,b) x+b;
for j=1:best_model.Nosc;
    fig = figure;
    timg = (squeeze(PhaseM(:,:,j)));
    min_timg = min(min((timg(ylim_low:ylim_high,xlim_high:xlim_low))));
    timg = timg + min_timg;
    timgabs = squeeze(abs(Cm(:,:,j)));
    timgabs = timgabs/max(max(timgabs));
    cutoff = 0.2;
    timg(timgabs<cutoff) = -100;
    timg = timg - min_timg;
    max_timg = max(max(timg(timgabs>cutoff)));
    min_timg = min(min(timg(timgabs>cutoff)));
    

    
    
    v = linspace(-pi,pi,Ncontours);
    contourf(w2wn(xaxis)*1e-3,w2wn(yaxis)*1e-3,timg,v,'linewidth',0.01);
    colormap('HSV')
%     colorbar;
    hold all;
    line(w2wn(xaxis)*1e-3,w2wn(line_diag_eq(xaxis,0))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
    line(w2wn(xaxis)*1e-3,w2wn(line_diag_eq(xaxis,-wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
    line(w2wn(xaxis)*1e-3,w2wn(line_diag_eq(xaxis,+wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
    line(w2wn(xaxis)*1e-3,w2wn(line_diag_eq(xaxis,-2*wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
    xlim([14.0 17.0]);
    ylim([14.0 16.0]);
    xlabel('Excitation frequency in cm^{-1}x10^{3}','FontSize',14);
    ylabel('Detection frequency in cm^{-1}x10^{3}','FontSize',14);
    set(gca,'XTick',[14.0 14.5 15.0 15.5 16.0 16.5 17.0])
    set(gca,'YTick',[14.0 14.5 15.0 15.5 16.0])
    set(gca,'FontSize',14)
%     grid
    title(['w2 = ' num2str(peaks(j)) ' with inv linewidth ' num2str(line_width_vec(j)) ' fs']);
    svgname = fullfile(figures_path,['r_phase_at_' num2str(peaks(j),3) '_wn.svg']);
    pdfname = fullfile(figures_path,['r_phase_at_' num2str(peaks(j),3) '_wn.pdf']);
    plot2svg(svgname,fig);
    export_string = ['inkscape ' svgname  ...
        ' --export-pdf=' pdfname ...
        ' --export-dpi=300'];
    system(export_string);
    close all;
end