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
figures_path = fullfile(base,'2D','Frank','LyX_Lab_Notebook','LaserLab','13-09-05','FiguresRephasingOnly');
data_file_name_r = 'data_stack_r.mat';
data_file_name_nr = 'data_stack_nr.mat';
optModelName = 'Optimal_Model_Rephasing_2';
%% Load the rephasing data, compress it.
if (~exist('all_data_mat_r','var') || ~exist('T2vec','var'));
    rdata = load(fullfile(data_path,data_file_name_r));
    rs = fieldnames(rdata);
end
Trange = 1:185;
T = -rdata.(rs{4})(Trange);
rdata.(rs{1}) = rdata.(rs{1})(:,:,Trange);
mean_data_r = squeeze(mean(rdata.(rs{1}),3));
mean_data_r = mean_data_r/max(max(real(mean_data_r)));
f = @(L,flen) [exp(-linspace(0,flen-1,flen)'.^4/(flen/2)^4); zeros(L-flen,1)];
L1 = 150;
L2 = 150;
f1 = f(size(rdata.(rs{1}),1),L1);
f2 = f(size(rdata.(rs{1}),2),L2);
Fmat = f1*f2';
Rcompressed = zeros(150,150,size(rdata.(rs{1}),3));
for k=1:size(rdata.(rs{1}),3);
Rcompressed(:,:,k)  = idct(idct(dct(dct(squeeze(rdata.(rs{1})(:,:,k)))')'.*Fmat,150)',150)';
end
R = zeros(L1*L2,size(rdata.(rs{1}),3));
for k=1:L2
    R(((k-1)*L1+1):(k*L1),:) = squeeze(Rcompressed(:,k,:));
end
[Ur, Sr, Vr] = svd(R,0);
Yr = (Sr*Vr')';
%% Load the non-rephasing data, compress it.
% if (~exist('all_data_mat_nr','var') || ~exist('T2vec','var'));
%     nrdata = load(fullfile(data_path,data_file_name_nr));
%     nrs = fieldnames(nrdata);
% end
% nrdata.(nrs{1}) = nrdata.(nrs{1})(:,:,Trange);
% mean_data_nr = squeeze(mean(nrdata.(nrs{1}),3));
% mean_data_nr = mean_data_nr/max(max(real(mean_data_nr)));
% L1 = 150;
% L2 = 150;
% f1 = f(size(rdata.(rs{1}),1),L1);
% f2 = f(size(rdata.(rs{1}),2),L2);
% Fmat = f1*f2';
% NRcompressed = zeros(150,150,size(nrdata.(nrs{1}),3));
% for k=1:size(nrdata.(nrs{1}),3);
% NRcompressed(:,:,k)  = idct(idct(dct(dct(squeeze(nrdata.(nrs{1})(:,:,k)))')'.*Fmat,150)',150)';
% end
% NR = zeros(L1*L2,size(nrdata.(nrs{1}),3));
% for k=1:L2
%     NR(((k-1)*L1+1):(k*L1),:) = squeeze(NRcompressed(:,k,:));
% end
% [Unr, Snr, Vnr] = svd(NR,0);
% Ynr = (Snr*Vnr')';
%% Load the model
load(fullfile(results_path,optModelName));


%% reconstruct the complex amplitude
cAmps = zeros(size(best_model.Fit.alpha,2),best_model.Nosc);
for k=1:best_model.Nosc;
    cAmps(:,k) = best_model.Fit.alpha(2*(k-1)+1,:) + 1i*best_model.Fit.alpha(2*(k-1)+2,:);
end
peaks = T2wn(1./best_model.Fit.beta(2:2:(2*best_model.Nosc)));
line_width_vec = 1./best_model.Fit.beta(1:2:(2*best_model.Nosc));
rmapsu = zeros(L1*L2,best_model.Nosc);
% nrmapsu = zeros(L1*L2,best_model.Nosc);
for k=1:best_model.Nosc;
    rmapsu(:,k) = Ur*cAmps(1:size(Ur,2),k);
%     nrmapsu(:,k) = Unr*cAmps(size(Ur,2)+1:end,k);
end
rmaps = zeros(L1,best_model.Nosc,L2);
% nrmaps = zeros(L1,best_model.Nosc,L2);
for k=1:L2
    rmaps(:,:,k) = rmapsu(((k-1)*L1+1):(k*L1),:);
%     nrmaps(:,:,k) = nrmapsu(((k-1)*L1+1):(k*L1),:);
end
rmaps = permute(rmaps,[1 3 2]);
% nrmaps = permute(nrmaps,[1 3 2]);

%% Construct the phase maps in the original pixel space

P = best_model.modelfun(best_model.Fit.beta);
vectorize = @(m) reshape(m,[size(m,1)*size(m,2) size(m,3)]);
matricize = @(m,X,Y,Z) reshape(m,[X Y Z]);
X = vectorize(rdata.(rs{1}));
a = P\X';
C = zeros(size(a,2),best_model.Nosc); for k=1:best_model.Nosc; C(:,k) = a(2*(k-1)+1,:) + 1i*a(2*(k-1)+2,:); end;
Cm = zeros(size(rdata.(rs{1}),1),size(rdata.(rs{1}),2),best_model.Nosc);
PhaseM = abs(angle(Cm));
for k=1:best_model.Nosc; Cm(:,:,k) = matricize(C(:,k),size(rdata.(rs{1}),1),size(rdata.(rs{1}),2),1); end;
xaxis = rdata.(rs{2});
yaxis = rdata.(rs{3});
T = rdata.(rs{4});
%save(fullfile(results_path,'ComplexAmplitudes.mat'),'Cm','xaxis','yaxis','T');

%% Re-expand the maps to full size
rmapsfull = zeros(size(rdata.(rs{1})));
rpmapsfull = zeros(size(rdata.(rs{1})));
% nrmapsfull = zeros(size(rdata.(rs{1})));
for k=1:best_model.Nosc;
    tempr = dct(dct(squeeze(rmaps(:,:,k)))')';
    ptempr = dct(dct(squeeze(rpmaps(:,:,k)))')';
%     tempnr = dct(dct(squeeze(nrmaps(:,:,k)))')';
    rmapsfull(1:L1,1:L2,k) = tempr;
    rpmapsfull(1:L1,1:L2,k) = ptempr;
%     nrmapsfull(1:L1,1:L2,k) = tempnr;
    rmapsfull(:,:,k) = idct(idct(rmapsfull(:,:,k))')';
    rpmapsfull(:,:,k) = idct(idct(rpmapsfull(:,:,k))')';
%     nrmapsfull(:,:,k) = idct(idct(nrmapsfull(:,:,k))')';
end

%% plot the rephasing amp maps
Ncontours = 50;
ylim_low = dsearchn(rdata.(rs{3})',wn2w(14000));
xlim_low = dsearchn(rdata.(rs{2})',wn2w(14000));
ylim_high = dsearchn(rdata.(rs{3})',wn2w(17000));
xlim_high = dsearchn(rdata.(rs{2})',wn2w(17000));
max_bkg = max(max(mean_data_r));
min_bkg = min(min(mean_data_r));
line_diag_eq = @(x,b) x+b;
% for j=1:best_model.Nosc;
%     fig = figure;
%     timg = abs(squeeze(rmapsfull(:,:,j)));
%     timg = timg/max(max(timg(ylim_low:ylim_high,xlim_high:xlim_low)));
%     max_timg = max(max(abs(timg(ylim_low:ylim_high,xlim_high:xlim_low))));
%     min_timg = min(min(abs(timg(ylim_low:ylim_high,xlim_high:xlim_low))));
%     v = linspace(min_timg,max_timg,Ncontours);
%     contourf(w2wn(rdata.(rs{2}))*1e-3,w2wn(rdata.(rs{3}))*1e-3,timg,v,'linewidth',0.01);
%     a = makeColorMapWStops(linspace(1,Ncontours,6),[[1,1,1]; [1,1,1]; [0,0,1]; [0,1,0]; [1,1,0]; [1,0,0]]);
%     colormap(a)
% %     colorbar;
%     hold all;
%     contour(w2wn(rdata.(rs{2}))*1e-3,w2wn(rdata.(rs{3}))*1e-3,abs(mean_data_r),linspace(min_bkg,0.5*max_bkg,8),'linewidth',0.5,'linecolor','k');
%     line(w2wn(rdata.(rs{2}))*1e-3,w2wn(line_diag_eq(rdata.(rs{2}),0))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
%     line(w2wn(rdata.(rs{2}))*1e-3,w2wn(line_diag_eq(rdata.(rs{2}),-wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
%     line(w2wn(rdata.(rs{2}))*1e-3,w2wn(line_diag_eq(rdata.(rs{2}),+wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
%     line(w2wn(rdata.(rs{2}))*1e-3,w2wn(line_diag_eq(rdata.(rs{2}),-2*wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
%     xlim([14.0 17.0]);
%     ylim([14.0 17.0]);
%     xlabel('Excitation frequency in cm^{-1}x10^{3}','FontSize',14);
%     ylabel('Detection frequency in cm^{-1}x10^{3}','FontSize',14);
%     set(gca,'XTick',[14.0 14.5 15.0 15.5 16.0 16.5 17.0])
%     set(gca,'YTick',[14.0 14.5 15.0 15.5 16.0 16.5 17.0])
%     set(gca,'FontSize',14)
%     title(['w2 = ' num2str(peaks(j)) ' with inv linewidth ' num2str(line_width_vec(j)) ' fs']);
%     svgname = fullfile(figures_path,['r_beats_at_' num2str(peaks(j),3) '_wn.svg']);
%     pdfname = fullfile(figures_path,['r_beats_at_' num2str(peaks(j),3) '_wn.pdf']);
%     plot2svg(svgname,fig);
%     export_string = ['inkscape ' svgname  ...
%         ' --export-pdf=' pdfname ...
%         ' --export-dpi=300'];
%     system(export_string);
%     close all;
% end

%% Plot the rephasing phase maps

Ncontours = 50;
ylim_low = dsearchn(yaxis',wn2w(14000));
xlim_low = dsearchn(xaxis',wn2w(14000));
ylim_high = dsearchn(yaxis',wn2w(17000));
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
    

    
    
    v = linspace(min_timg,max_timg,Ncontours);
    contourf(w2wn(xaxis)*1e-3,w2wn(yaxis)*1e-3,timg,v,'linewidth',0.01);
    a = makeColorMapWStops(linspace(1,Ncontours,5),[[1,1,1]; [0,0,1]; [0,1,0]; [1,1,0]; [1,0,0]]);
    colormap(a)
    colorbar;
    hold all;
    line(w2wn(xaxis)*1e-3,w2wn(line_diag_eq(xaxis,0))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
    line(w2wn(xaxis)*1e-3,w2wn(line_diag_eq(xaxis,-wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
    line(w2wn(xaxis)*1e-3,w2wn(line_diag_eq(xaxis,+wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
    line(w2wn(xaxis)*1e-3,w2wn(line_diag_eq(xaxis,-2*wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
    xlim([14.0 17.0]);
    ylim([14.0 17.0]);
    xlabel('Excitation frequency in cm^{-1}x10^{3}','FontSize',14);
    ylabel('Detection frequency in cm^{-1}x10^{3}','FontSize',14);
    set(gca,'XTick',[14.0 14.5 15.0 15.5 16.0 16.5 17.0])
    set(gca,'YTick',[14.0 14.5 15.0 15.5 16.0 16.5 17.0])
    set(gca,'FontSize',14)
    title(['w2 = ' num2str(peaks(j)) ' with inv linewidth ' num2str(line_width_vec(j)) ' fs']);
    return
    svgname = fullfile(figures_path,['r_phase_at_' num2str(peaks(j),3) '_wn.svg']);
    pdfname = fullfile(figures_path,['r_phase_at_' num2str(peaks(j),3) '_wn.pdf']);
    plot2svg(svgname,fig);
    export_string = ['inkscape ' svgname  ...
        ' --export-pdf=' pdfname ...
        ' --export-dpi=300'];
    system(export_string);
    close all;
end


% %% Plot the non-rephasing maps
% 
% Ncontours = 50;
% ylim_low = dsearchn(rdata.(rs{3})',wn2w(14000));
% xlim_low = dsearchn(rdata.(rs{2})',wn2w(14000));
% ylim_high = dsearchn(rdata.(rs{3})',wn2w(17000));
% xlim_high = dsearchn(rdata.(rs{2})',wn2w(17000));
% max_bkg = max(max(mean_data_nr));
% min_bkg = min(min(mean_data_nr));
% line_diag_eq = @(x,b) x+b;
% for j=1:best_model.Nosc;
%     fig = figure;
%     timg = abs(squeeze(nrmapsfull(:,:,j)));
%     timg = timg/max(max(timg(ylim_low:ylim_high,xlim_high:xlim_low)));
%     max_timg = max(max(abs(timg(ylim_low:ylim_high,xlim_high:xlim_low))));
%     min_timg = min(min(abs(timg(ylim_low:ylim_high,xlim_high:xlim_low))));
%     v = linspace(min_timg,max_timg,Ncontours);
%     contourf(w2wn(rdata.(rs{2}))*1e-3,w2wn(rdata.(rs{3}))*1e-3,timg,v,'linewidth',0.01);
%     a = makeColorMapWStops(linspace(1,Ncontours,6),[[1,1,1]; [1,1,1]; [0,0,1]; [0,1,0]; [1,1,0]; [1,0,0]]);
%     colormap(a)
% %     colorbar;
%     hold all;
%     contour(w2wn(rdata.(rs{2}))*1e-3,w2wn(rdata.(rs{3}))*1e-3,abs(mean_data_nr),linspace(min_bkg,0.5*max_bkg,8),'linewidth',0.5,'linecolor','k');
%     line(w2wn(rdata.(rs{2}))*1e-3,w2wn(line_diag_eq(rdata.(rs{2}),0))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
%     line(w2wn(rdata.(rs{2}))*1e-3,w2wn(line_diag_eq(rdata.(rs{2}),-wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
%     line(w2wn(rdata.(rs{2}))*1e-3,w2wn(line_diag_eq(rdata.(rs{2}),+wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
%     line(w2wn(rdata.(rs{2}))*1e-3,w2wn(line_diag_eq(rdata.(rs{2}),-2*wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
%     xlim([14.0 17.0]);
%     ylim([14.0 17.0]);
%     xlabel('Excitation frequency in cm^{-1}x10^{3}','FontSize',14);
%     ylabel('Detection frequency in cm^{-1}x10^{3}','FontSize',14);
%     set(gca,'XTick',[14.0 14.5 15.0 15.5 16.0 16.5 17.0])
%     set(gca,'YTick',[14.0 14.5 15.0 15.5 16.0 16.5 17.0])
%     set(gca,'FontSize',14)
%     title(['w2 = ' num2str(peaks(j)) ' with inv linewidth ' num2str(line_width_vec(j)) ' fs']);
%     svgname = fullfile(figures_path,['nr_beats_at_' num2str(peaks(j),3) '_wn.svg']);
%     pdfname = fullfile(figures_path,['nr_beats_at_' num2str(peaks(j),3) '_wn.pdf']);
%     plot2svg(svgname,fig);
%     export_string = ['inkscape ' svgname  ...
%         ' --export-pdf=' pdfname ...
%         ' --export-dpi=300'];
%     system(export_string);
%     close all;
% end
fprintf(1,'All Done with that crap\n');