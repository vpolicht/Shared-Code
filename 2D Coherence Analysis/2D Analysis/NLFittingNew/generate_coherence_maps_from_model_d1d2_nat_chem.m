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
optModelName = 'OldModel_for_NatChem';
%% Load the rephasing data, compress it.
if (~exist('all_data_mat_r','var') || ~exist('T2vec','var'));
    rdata = load(fullfile(data_path,data_file_name_r));
    rs = fieldnames(rdata);
end
Trange = 1:185;
T = -rdata.(rs{4})(Trange);
rdata.(rs{1}) = rdata.(rs{1})(:,:,Trange);
% mean_data_r = squeeze(mean(rdata.(rs{1}),3));
% mean_data_r = mean_data_r/max(max(real(mean_data_r)));
% Nbf = 60;
% Bx = dop( size(rdata.(rs{1}),2), Nbf );
% By = dop( size(rdata.(rs{1}),1), Nbf );
[SRX, SRY, SRZ] = size(rdata.(rs{1}));
% Rcompressed = zeros(Nbf,Nbf,SRZ);
% for k=1:SRZ; Rcompressed(:,:,k) = By' * squeeze(rdata.(rs{1})(:,:,k)) * Bx; end;
% R = zeros(Nbf*Nbf,SRZ);
% for k=1:Nbf
%     R(((k-1)*Nbf+1):(k*Nbf),:) = squeeze(Rcompressed(:,k,:));
% end
% [Ur, Sr, Vr] = svd(R,0);
% Yr = (Sr*Vr')';
%% Load the non-rephasing data, compress it.
% if (~exist('all_data_mat_nr','var') || ~exist('T2vec','var'));
%     nrdata = load(fullfile(data_path,data_file_name_nr));
%     nrs = fieldnames(nrdata);
% end
% nrdata.(nrs{1}) = nrdata.(nrs{1})(:,:,Trange);
% mean_data_nr = squeeze(mean(nrdata.(nrs{1}),3));
% mean_data_nr = mean_data_nr/max(max(real(mean_data_nr)));
% Nbf = 60;
% NRcompressed = zeros(Nbf,Nbf,SRZ);
% for k=1:SRZ; NRcompressed(:,:,k) = By' * squeeze(nrdata.(nrs{1})(:,:,k)) * Bx; end;
% NR = zeros(Nbf*Nbf,SRZ);
% for k=1:Nbf
%     NR(((k-1)*Nbf+1):(k*Nbf),:) = squeeze(NRcompressed(:,k,:));
% end
% [Unr, Snr, Vnr] = svd(NR,0);
% Ynr = (Snr*Vnr')';
%% Load the model
load(fullfile(results_path,optModelName));

%% Construct the phase maps in the original pixel space

% P = best_model.modelfun(best_model.Fit.beta);
% vectorize = @(m) reshape(m,[size(m,1)*size(m,2) size(m,3)]);
% matricize = @(m,X,Y,Z) reshape(m,[X Y Z]);
% X = vectorize(rdata.(rs{1}));
% [sX,sY,sZ] = size(rdata.(rs{1}));
% a = P\X';
% Yhat = matricize((P*a)',sX,sY,sZ);
% Z = rdata.(rs{1}) - Yhat;
% fprintf(1,'Saving linear coeffs...\n');
% save(fullfile(data_path,'linear_coeffs_old_model.mat'),'a');
% fprintf(1,'Saving residual...\n');
% save(fullfile(data_path,'residual_oldmodel.mat'),'Z');
load(fullfile(data_path,'residual_oldmodel.mat'),'Z');

filter_f = @(L,a) [0.5; subsref(tukeywin(2*L,a),struct('type','()','subs',{{(L+2):(2*L),1}}))];
filter_m = @(x,y,L,a) permute(repmat(filter_f(L,a),[1 x y]),[2 3 1]);
fZ = Z.*filter_m(sX,sY,sZ,0.3);
F = fft(fZ,2^10,3);
Faxis = w2wn(MakeFourierOmegaAxis(T,2^10));
frobeniusSpectrum = zeros(size(F,3),1);

for k=1:size(F,3);
    frobeniusSpectrum(k) = norm(F(:,:,k),'fro');
end
figure; plot(Faxis,frobeniusSpectrum);
xlim([0 1000]);
saveas(gcf,fullfile(figures_path,'frobenius_spectrum_nat_chem.png'));

%% plot the rephasing amp maps
% Ncontours = 50;
% ylim_low = dsearchn(rdata.(rs{3})',wn2w(14000));
% xlim_low = dsearchn(rdata.(rs{2})',wn2w(14000));
% ylim_high = dsearchn(rdata.(rs{3})',wn2w(16000));
% xlim_high = dsearchn(rdata.(rs{2})',wn2w(17000));
% max_bkg = max(max(mean_data_r));
% min_bkg = min(min(mean_data_r));
% line_diag_eq = @(x,b) x+b;
% peaks = [91.24 127.1 250.9 338.9 729.9 853.7 974.3];
% for j=1:best_model.Nosc;
%     tind = dsearchn(Faxis',peaks(k));
%     fig = figure;
%     timg = abs(squeeze(F(:,:,tind)));
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
%     xlim([14.00 17.00]);
%     ylim([14.00 16.00]);
%     xlabel('Excitation frequency in cm^{-1}x10^{3}','FontSize',14);
%     ylabel('Detection frequency in cm^{-1}x10^{3}','FontSize',14);
%     set(gca,'XTick',[14.0 14.5 15.0 15.5 16.0 16.5 17.0])
%     set(gca,'YTick',[14.0 14.5 15.0 15.5 16.0])
% %     set(gca,'XTick',[14.5 15.0 15.5])
% %     set(gca,'YTick',[14.5 15.0 15.5])
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

% % %% Plot the rephasing phase maps
% % 
% % Ncontours = 50;
% % ylim_low = dsearchn(yaxis',wn2w(14000));
% % xlim_low = dsearchn(xaxis',wn2w(14000));
% % ylim_high = dsearchn(yaxis',wn2w(16000));
% % xlim_high = dsearchn(xaxis',wn2w(16000));
% % line_diag_eq = @(x,b) x+b;
% % for j=1:best_model.Nosc;
% %     fig = figure;
% %     timg = (squeeze(PhaseM(:,:,j)));
% %     min_timg = min(min((timg(ylim_low:ylim_high,xlim_high:xlim_low))));
% %     timg = timg + min_timg;
% %     timgabs = squeeze(abs(Cm(:,:,j)));
% %     timgabs = timgabs/max(max(timgabs));
% %     cutoff = 0.2;
% %     timg(timgabs<cutoff) = -100;
% %     timg = timg - min_timg;
% %     max_timg = max(max(timg(timgabs>cutoff)));
% %     min_timg = min(min(timg(timgabs>cutoff)));
% %     
% % 
% %     
% %     
% %     v = linspace(min_timg,max_timg,Ncontours);
% %     contourf(w2wn(xaxis)*1e-3,w2wn(yaxis)*1e-3,timg,v,'linewidth',0.01);
% %     a = makeColorMapWStops(linspace(1,Ncontours,5),[[1,1,1]; [0,0,1]; [0,1,0]; [1,1,0]; [1,0,0]]);
% %     colormap(a)
% %     colorbar;
% %     hold all;
% %     line(w2wn(xaxis)*1e-3,w2wn(line_diag_eq(xaxis,0))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
% %     line(w2wn(xaxis)*1e-3,w2wn(line_diag_eq(xaxis,-wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
% %     line(w2wn(xaxis)*1e-3,w2wn(line_diag_eq(xaxis,+wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
% %     line(w2wn(xaxis)*1e-3,w2wn(line_diag_eq(xaxis,-2*wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
% %     xlim([14.0 17.0]);
% %     ylim([14.0 17.0]);
% %     xlabel('Excitation frequency in cm^{-1}x10^{3}','FontSize',14);
% %     ylabel('Detection frequency in cm^{-1}x10^{3}','FontSize',14);
% %     set(gca,'XTick',[14.0 14.5 15.0 15.5 16.0 16.5 17.0])
% %     set(gca,'YTick',[14.0 14.5 15.0 15.5 16.0 16.5 17.0])
% %     set(gca,'FontSize',14)
% %     title(['w2 = ' num2str(peaks(j)) ' with inv linewidth ' num2str(line_width_vec(j)) ' fs']);
% %     svgname = fullfile(figures_path,['r_phase_at_' num2str(peaks(j),3) '_wn.svg']);
% %     pdfname = fullfile(figures_path,['r_phase_at_' num2str(peaks(j),3) '_wn.pdf']);
% %     plot2svg(svgname,fig);
% %     export_string = ['inkscape ' svgname  ...
% %         ' --export-pdf=' pdfname ...
% %         ' --export-dpi=300'];
% %     system(export_string);
% %     close all;
% % end
% 
% 
% %% Plot the non-rephasing maps
% % 
% % Ncontours = 50;
% % ylim_low = dsearchn(rdata.(rs{3})',wn2w(14000));
% % xlim_low = dsearchn(rdata.(rs{2})',wn2w(14000));
% % ylim_high = dsearchn(rdata.(rs{3})',wn2w(17000));
% % xlim_high = dsearchn(rdata.(rs{2})',wn2w(17000));
% % max_bkg = max(max(mean_data_nr));
% % min_bkg = min(min(mean_data_nr));
% % line_diag_eq = @(x,b) x+b;
% % for j=1:best_model.Nosc;
% %     fig = figure;
% %     timg = abs(squeeze(nrmapsfull(:,:,j)));
% %     timg = timg/max(max(timg(ylim_low:ylim_high,xlim_high:xlim_low)));
% %     max_timg = max(max(abs(timg(ylim_low:ylim_high,xlim_high:xlim_low))));
% %     min_timg = min(min(abs(timg(ylim_low:ylim_high,xlim_high:xlim_low))));
% %     v = linspace(min_timg,max_timg,Ncontours);
% %     contourf(w2wn(rdata.(rs{2}))*1e-3,w2wn(rdata.(rs{3}))*1e-3,timg,v,'linewidth',0.01);
% %     a = makeColorMapWStops(linspace(1,Ncontours,6),[[1,1,1]; [1,1,1]; [0,0,1]; [0,1,0]; [1,1,0]; [1,0,0]]);
% %     colormap(a)
% % %     colorbar;
% %     hold all;
% %     contour(w2wn(rdata.(rs{2}))*1e-3,w2wn(rdata.(rs{3}))*1e-3,abs(mean_data_nr),linspace(min_bkg,0.5*max_bkg,8),'linewidth',0.5,'linecolor','k');
% %     line(w2wn(rdata.(rs{2}))*1e-3,w2wn(line_diag_eq(rdata.(rs{2}),0))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
% %     line(w2wn(rdata.(rs{2}))*1e-3,w2wn(line_diag_eq(rdata.(rs{2}),-wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
% %     line(w2wn(rdata.(rs{2}))*1e-3,w2wn(line_diag_eq(rdata.(rs{2}),+wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
% %     line(w2wn(rdata.(rs{2}))*1e-3,w2wn(line_diag_eq(rdata.(rs{2}),-2*wn2w(peaks(j))))*1e-3,'linestyle','--','linewidth',0.5,'Color','k');
% %     xlim([14.0 17.0]);
% %     ylim([14.0 17.0]);
% %     xlabel('Excitation frequency in cm^{-1}x10^{3}','FontSize',14);
% %     ylabel('Detection frequency in cm^{-1}x10^{3}','FontSize',14);
% %     set(gca,'XTick',[14.0 14.5 15.0 15.5 16.0 16.5 17.0])
% %     set(gca,'YTick',[14.0 14.5 15.0 15.5 16.0 16.5 17.0])
% %     set(gca,'FontSize',14)
% %     title(['w2 = ' num2str(peaks(j)) ' with inv linewidth ' num2str(line_width_vec(j)) ' fs']);
% %     svgname = fullfile(figures_path,['nr_beats_at_' num2str(peaks(j),3) '_wn.svg']);
% %     pdfname = fullfile(figures_path,['nr_beats_at_' num2str(peaks(j),3) '_wn.pdf']);
% %     plot2svg(svgname,fig);
% %     export_string = ['inkscape ' svgname  ...
% %         ' --export-pdf=' pdfname ...
% %         ' --export-dpi=300'];
% %     system(export_string);
% %     close all;
% % end
fprintf(1,'All Done with that crap\n');