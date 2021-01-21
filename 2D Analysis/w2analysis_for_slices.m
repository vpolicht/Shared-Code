clear all; close all; clc
load(['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-09-05/PostProcessingData/' 'data_stack']);
c = 299.792458;
w2wn = @(w) 1E7*w/(2*pi*c);
wn2f = @(wn) (c*wn)./(1E7);
wn2w = @(wn) (2*pi*c*wn)./(1E7);
data.data = all_data_mat;
clear all_data_mat;
data.data2 = squeeze(sum(reshape(data.data,[4 size(data.data,1)/4 size(data.data,2) size(data.data,3)]),1)); %bin 4 consecutive pixels
data.data = data.data(:,:,1:185);
data.x_axis_w_roi = x_axis_w_roi;
data.y_axis_w_roi = y_axis_w_roi;
% data.y_axis_w_roi = squeeze(mean(reshape(data.y_axis_w_roi,[4 length(data.y_axis_w_roi)/4]),1)); %bin 4 consecutive pixels
data.T = -T2vec;
data.T = data.T(1:185);
data.W2 = MakeFourierOmegaAxis(data.T,length(data.T));
data.W2 = data.W2(1:floor(length(data.T)/2));
data.WN = w2wn(data.W2);
data.waxis = wn2w(linspace(100,1000,200))';
data.wnaxis = w2wn(MakeFourierOmegaAxis(data.T,2^12));
data.basis_A = GenCompressedBasis(data.waxis,data.T);
segment_starts = [80];
segment_ends = [2050];
ss_this_data = segment_starts(segment_starts<=max(round(data.T)));
se_this_data = segment_ends(segment_ends<=max(round(data.T)));
data_starts = 0*ss_this_data; data_ends = 0*se_this_data;
for k=1:length(data_starts); data_starts(k) = dsearchn(data.T,ss_this_data(k)); end;
for k=1:length(data_ends); data_ends(k) = dsearchn(data.T,se_this_data(k)); end;
exciton_lines = sort(2*pi*c./([662.1 669.3 671.0 672.2 673.5 675 677.2 680.8 693.1]));

%%
load(['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-09-05/PostProcessingData/' 'bkg_cube']);
% bkg_cube = PointByPointExponentialFit(data.data2,data.T);
data.z = data.data(:,:,1:end) - bkg_cube(:,:,1:end);
% data.z = data.data2;
[data.Ts, data.Tsind] = sort(data.T);
data.z = data.z(:,:,data.Tsind);
data.zf = abs(fft(data.z,2^12,3));
data.zfaxis = w2wn(MakeFourierOmegaAxis(data.Ts,2^12));

%%

xmin = dsearchn(data.x_axis_w_roi',2*pi*c/700); xmax = dsearchn(data.x_axis_w_roi',(2*pi*c/630));
ymin = dsearchn(data.y_axis_w_roi',(2*pi*c/700)); ymax = dsearchn(data.y_axis_w_roi',(2*pi*c/630));
wnsearch = linspace(50,1000,100);
x_ticks = linspace(data.x_axis_w_roi(xmin),data.x_axis_w_roi(xmax),9);%exciton_lines;
y_ticks = linspace(data.y_axis_w_roi(ymin),data.y_axis_w_roi(ymax),9);%exciton_lines;
fig1 = contourf(data.x_axis_w_roi(xmin:xmax),data.y_axis_w_roi(ymin:ymax),abs(data.zf(ymin:ymax,xmin:xmax,dsearchn(data.zfaxis',wnsearch(1)))),10);
set(gca,'XTick',x_ticks,'XTickLabel',cellstr(num2str(((9:-1:1)'),3)));
set(gca,'YTick',x_ticks,'YTickLabel',cellstr(num2str((2*pi*c./y_ticks'),3)));
line_diag_eq = @(x,b) x+b;
x_sub = data.x_axis_w_roi(xmin:xmax);
y_sub = data.y_axis_w_roi(ymin:ymax);
max_z = max(max(max(abs(data.zf))));

contour_v = [0 linspace(0.006,0.2,20)];

for k=1:length(wnsearch);
    fig1 = contourf(x_sub,y_sub,abs(data.zf(ymin:ymax,xmin:xmax,dsearchn(data.zfaxis',wnsearch(k))))/max_z,contour_v);
    set(gca,'XTick',x_ticks,'XTickLabel',cellstr(num2str(((9:-1:1)'),3)));
    set(gca,'YTick',x_ticks,'YTickLabel',cellstr(num2str((2*pi*c./y_ticks'),3)));
    line(x_sub,line_diag_eq(x_sub,0),'linewidth',2,'Color','w');
    line(x_sub,line_diag_eq(x_sub,-wn2w(wnsearch(k))),'linewidth',2,'Color','w');
    line(x_sub,line_diag_eq(x_sub,+wn2w(wnsearch(k))),'linewidth',2,'Color','w');
    AddExcitonLinesToImage(gcf,x_sub,y_sub,exciton_lines);
    colorbar;
    drawnow;
    saveas(gcf,['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-09-05/Figures/Beat Analysis 1/' 'beats_at_' num2str(wnsearch(k)) '_wavenumbers' '.png']);
end

% %% This section for generating Chla pixel traces
% x = PickPointWithRadius(data,2*pi*c/650,2*pi*c/678);
% x = PickPointWithRadius(data,2*pi*c/656,2*pi*c/656);
% x = PickPointWithRadius(data,2*pi*c/670,2*pi*c/653);
% 
% x = PickPointWithRadius(data,2*pi*c/660,2*pi*c/663);
% x = PickPointWithRadius(data,2*pi*c/672,2*pi*c/654);
% 
% x = PickPointWithRadius(data,2*pi*c/653,2*pi*c/676);
% x = PickPointWithRadius(data,2*pi*c/670,2*pi*c/670);
% x = PickPointWithRadius(data,2*pi*c/670,2*pi*c/654);
% 
% 
% x = PickPointWithRadius(data,2*pi*c/648,2*pi*c/672);
% x = PickPointWithRadius(data,2*pi*c/672,2*pi*c/654);
% x = PickPointWithRadius(data,2*pi*c/643,2*pi*c/656);
% 
% x = PickPointWithRadius(data,2*pi*c/643,2*pi*c/674);
% x = PickPointWithRadius(data,2*pi*c/660,2*pi*c/660);
% x = PickPointWithRadius(data,2*pi*c/672,2*pi*c/654);
% x = PickPointWithRadius(data,2*pi*c/639,2*pi*c/656);
% 
% %% This section for generating d1d2 pixel traces
% 
% x = PickPointWithRadiusD1D2(data,2*pi*c/668,2*pi*c/676);
% x = PickPointWithRadiusD1D2(data,2*pi*c/680,2*pi*c/670);
% 
% x = PickPointWithRadiusD1D2(data,2*pi*c/666,2*pi*c/676);
% x = PickPointWithRadiusD1D2(data,2*pi*c/666,2*pi*c/685);
% 
% x = PickPointWithRadiusD1D2(data,2*pi*c/660,2*pi*c/676);
% x = PickPointWithRadiusD1D2(data,2*pi*c/681,2*pi*c/681);
% x = PickPointWithRadiusD1D2(data,2*pi*c/660,2*pi*c/690);
% 
% 
% x = PickPointWithRadiusD1D2(data,2*pi*c/676,2*pi*c/659);
% x = PickPointWithRadiusD1D2(data,2*pi*c/665,2*pi*c/688);
% x = PickPointWithRadiusD1D2(data,2*pi*c/657,2*pi*c/670);
% 
% x = PickPointWithRadiusD1D2(data,2*pi*c/657,2*pi*c/672);
% x = PickPointWithRadiusD1D2(data,2*pi*c/666,2*pi*c/666);
% x = PickPointWithRadiusD1D2(data,2*pi*c/678,2*pi*c/659);
% x = PickPointWithRadiusD1D2(data,2*pi*c/663,2*pi*c/690);
% 
% 
% x = PickPointWithRadiusD1D2(data,2*pi*c/640,2*pi*c/669);
% x = PickPointWithRadiusD1D2(data,2*pi*c/669,2*pi*c/669);
% x = PickPointWithRadiusD1D2(data,2*pi*c/672,2*pi*c/642);
% 
% x = PickPointWithRadiusD1D2(data,2*pi*c/636,2*pi*c/674);
% x = PickPointWithRadiusD1D2(data,2*pi*c/674,2*pi*c/636);
% x = PickPointWithRadiusD1D2(data,2*pi*c/662,2*pi*c/677);
% x = PickPointWithRadiusD1D2(data,2*pi*c/675,2*pi*c/675);
% x = PickPointWithRadiusD1D2(data,2*pi*c/665,2*pi*c/665);
% x = PickPointWithRadiusD1D2(data,2*pi*c/677,2*pi*c/662);
% 
% %% Blah, this is getting so fugly...
% clear all; close all; clc
% compress_data = @(D,L) squeeze(sum(reshape(D,[L size(D,1)/L size(D,2) size(D,3)]),1));
% compress_axis = @(a,L) squeeze(mean(reshape(a,[L length(a)/L]),1));
% c = 299.792458;
% w2wn = @(w) 1E7*w/(2*pi*c);
% wn2f = @(wn) (c*wn)./(1E7);
% wn2w = @(wn) (2*pi*c*wn)./(1E7);
% 
% load(['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-07-06/PostProcessingData/' 'data_stack']);
% data.data = compress_data(all_data_mat,5); clear all_data_mat;
% data.data = data.data(:,:,8:end);
% load(['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-07-06/PostProcessingData/' 'data_stack_r']);
% data.data_r = compress_data(all_data_mat_r,5); clear all_data_mat_r;
% data.data_r = data.data_r(:,:,8:end);
% load(['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-07-06/PostProcessingData/' 'data_stack_nr']);
% data.data_nr = compress_data(all_data_mat_nr,5);  clear all_data_mat_nr;
% data.data_nr = data.data_nr(:,:,8:end);
% data.x_axis_w_roi = x_axis_w_roi;
% data.y_axis_w_roi = compress_axis(y_axis_w_roi,5);
% data.T = -T2vec(8:66);
% 
% load(['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-07-06/PostProcessingData/' 'bkg_cube_phased']);
% data.bkg = bkg_cube; clear bkg_cube;
% load(['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-07-06/PostProcessingData/' 'bkg_cube_rephasing']);
% data.bkg_r = bkg_cube; clear bkg_cube;
% load(['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-07-06/PostProcessingData/' 'bkg_cube_nonrephasing']);
% data.bkg_nr = bkg_cube; clear bkg_cube;
% 
% data.z = data.data - data.bkg(:,:,8:66);
% data.zf = abs(fft(data.z,2^12,3));
% data.zfaxis = w2wn(MakeFourierOmegaAxis(data.T(8:66),2^12));
% 
% data.zr = data.data_r - data.bkg_r(:,:,8:66);
% data.zfr = abs(fft(data.z,2^12,3));
% 
% data.znr = data.data_nr - data.bkg_nr(:,:,8:66);
% data.zfnr = abs(fft(data.z,2^12,3));
% %%
% 
% % PickPointWithRadiusOverlayData(data,2*pi*c/650,2*pi*c/678);
% % PickPointWithRadiusOverlayData(data,2*pi*c/670,2*pi*c/653);
% % 
% % PickPointWithRadiusOverlayData(data,2*pi*c/660,2*pi*c/663);
% % PickPointWithRadiusOverlayData(data,2*pi*c/672,2*pi*c/654);
% % PickPointWithRadiusOverlayData(data,2*pi*c/653,2*pi*c/676);
% % PickPointWithRadiusOverlayData(data,2*pi*c/670,2*pi*c/654);
% % 
% % 
% % PickPointWithRadiusOverlayData(data,2*pi*c/648,2*pi*c/672);
% % PickPointWithRadiusOverlayData(data,2*pi*c/672,2*pi*c/654);
% % PickPointWithRadiusOverlayData(data,2*pi*c/643,2*pi*c/656);
% % PickPointWithRadiusOverlayData(data,2*pi*c/643,2*pi*c/674);
% % PickPointWithRadiusOverlayData(data,2*pi*c/672,2*pi*c/654);
% % PickPointWithRadiusOverlayData(data,2*pi*c/639,2*pi*c/656);
% 
% PickPointWithRadiusOverlayData(data,2*pi*c/660,2*pi*c/670);
% PickPointWithRadiusOverlayData(data,2*pi*c/655,2*pi*c/675);
% 
% %%
% xmin = dsearchn(data.x_axis_w_roi',2*pi*c/685); xmax = dsearchn(data.x_axis_w_roi',(2*pi*c/630));
% ymin = dsearchn(data.y_axis_w_roi',(2*pi*c/685)); ymax = dsearchn(data.y_axis_w_roi',(2*pi*c/630));
% wnsearch = linspace(50,1000,100);
% x_ticks = linspace(data.x_axis_w_roi(xmin),data.x_axis_w_roi(xmax),9);%exciton_lines;
% y_ticks = linspace(data.y_axis_w_roi(ymin),data.y_axis_w_roi(ymax),9);%exciton_lines;
% fig1 = contourf(data.x_axis_w_roi(xmin:xmax),data.y_axis_w_roi(ymin:ymax),abs(data.zf(ymin:ymax,xmin:xmax,dsearchn(data.zfaxis',wnsearch(1)))),10);
% set(gca,'XTick',x_ticks,'XTickLabel',cellstr(num2str((2*pi*c./x_ticks'),3)));
% set(gca,'YTick',x_ticks,'YTickLabel',cellstr(num2str((2*pi*c./y_ticks'),3)));
% line_diag_eq = @(x,b) x+b;
% x_sub = data.x_axis_w_roi(xmin:xmax);
% y_sub = data.y_axis_w_roi(ymin:ymax);
% 
% 
% for k=1:length(wnsearch);
%     fig1 = contourf(x_sub,y_sub,abs(data.zf(ymin:ymax,xmin:xmax,dsearchn(data.zfaxis',wnsearch(k)))),10);
%     set(gca,'XTick',x_ticks,'XTickLabel',cellstr(num2str((2*pi*c./x_ticks'),3)));
%     set(gca,'YTick',x_ticks,'YTickLabel',cellstr(num2str((2*pi*c./y_ticks'),3)));
%     line(x_sub,line_diag_eq(x_sub,0),'linewidth',2,'Color','w');
%     line(x_sub,line_diag_eq(x_sub,-wn2w(wnsearch(k))),'linewidth',2,'Color','w');
%     line(x_sub,line_diag_eq(x_sub,+wn2w(wnsearch(k))),'linewidth',2,'Color','w');
%     drawnow;
%     saveas(gcf,['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-07-06/Figures/Beats_Analysis/' 'phased_beats_at_' num2str(wnsearch(k)) '_wavenumbers' '.png']);end
% 
% for k=1:length(wnsearch);
%     fig1 = contourf(x_sub,y_sub,abs(data.zfr(ymin:ymax,xmin:xmax,dsearchn(data.zfaxis',wnsearch(k)))),10);
%     set(gca,'XTick',x_ticks,'XTickLabel',cellstr(num2str((2*pi*c./x_ticks'),3)));
%     set(gca,'YTick',x_ticks,'YTickLabel',cellstr(num2str((2*pi*c./y_ticks'),3)));
%     line(x_sub,line_diag_eq(x_sub,0),'linewidth',2,'Color','w');
%     line(x_sub,line_diag_eq(x_sub,-wn2w(wnsearch(k))),'linewidth',2,'Color','w');
%     line(x_sub,line_diag_eq(x_sub,+wn2w(wnsearch(k))),'linewidth',2,'Color','w');
%     drawnow;
%     saveas(gcf,['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-07-06/Figures/Beats_Analysis/' 'rephasing_beats_at_' num2str(wnsearch(k)) '_wavenumbers' '.png']);
% end
% 
% for k=1:length(wnsearch);
%     fig1 = contourf(x_sub,y_sub,abs(data.zfnr(ymin:ymax,xmin:xmax,dsearchn(data.zfaxis',wnsearch(k)))),10);
%     set(gca,'XTick',x_ticks,'XTickLabel',cellstr(num2str((2*pi*c./x_ticks'),3)));
%     set(gca,'YTick',x_ticks,'YTickLabel',cellstr(num2str((2*pi*c./y_ticks'),3)));
%     line(x_sub,line_diag_eq(x_sub,0),'linewidth',2,'Color','w');
%     line(x_sub,line_diag_eq(x_sub,-wn2w(wnsearch(k))),'linewidth',2,'Color','w');
%     line(x_sub,line_diag_eq(x_sub,+wn2w(wnsearch(k))),'linewidth',2,'Color','w');
%     drawnow;
%     saveas(gcf,['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-07-06/Figures/Beats_Analysis/' 'non_rephasing_beats_at_' num2str(wnsearch(k)) '_wavenumbers' '.png']);
% end
% 
% 
% 
% 
% 
% 
% 
% 
% % [slice_a, slice_d, slice_a_axis, slice_d_axis] = CutSlice(data,2*pi*c/(650),2*pi*c/(672));
% % slice_a_bkg = 0*slice_a;
% % slice_a_nobkg = 0*slice_a;
% % slice_a_freq = zeros(size(slice_a,1),length(data.basis_A));
% % slice_a_end_states = zeros(size(slice_a,1),1);
% % for k=1:size(slice_a,1);
% %     [~,slice_a_bkg(k,:),~] = single_exp_fit(data.T,slice_a(k,:)');
% %     slice_a_nobkg(k,:) = slice_a(k,:) - slice_a_bkg(k,:);
% % end
% % slice_d_bkg = 0*slice_d;
% % slice_d_nobkg = 0*slice_d;
% % slice_d_freq = zeros(size(slice_d,1),length(data.basis_A));
% % slice_d_end_states = zeros(size(slice_d,1),1);
% % for k=1:size(slice_d,1);
% %     [~,slice_d_bkg(k,:),~] = single_exp_fit(data.T,slice_d(k,:)');
% %     slice_d_nobkg(k,:) = slice_d(k,:) - slice_d_bkg(k,:);
% % end
% % % 
% % % for k=1:size(slice_a,1);
% % %     [slice_a_freq(k,:), slice_a_end_states(k)] = BPDN(data.basis_A,slice_a_nobkg(k,:)');
% % % end
% % 
% % figure; imagesc(w2wn(MakeFourierOmegaAxis(data.T,2^12)),slice_a_axis,abs(fft(slice_a_nobkg,2^12,2)))
% % figure; imagesc(w2wn(MakeFourierOmegaAxis(data.T,2^12)),slice_d_axis,abs(fft(slice_d_nobkg,2^12,2)))