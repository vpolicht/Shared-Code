function [point_trace_t,point_trace_f] = PickPointWithRadiusD1D2(data,exc_center,det_center)
c = 299.792458;
exc_center_ind = dsearchn(data.x_axis_w_roi',exc_center);
det_center_ind = dsearchn(data.y_axis_w_roi',det_center);
point_trace_t = squeeze(sum(sum(data.z((det_center_ind-1):(det_center_ind+1),(exc_center_ind-1):(exc_center_ind+1),:),2),1));
point_trace_f = abs(fft(point_trace_t,2^12));
fig1 = figure; plot(data.T(1:end),point_trace_t,'LineWidth',3);
title_string_t = sprintf('D1D2 Exc at %g, Probed at %g vs T',2*pi*c/exc_center,2*pi*c/det_center);
xlabel('femtoseconds');
title(title_string_t,'FontSize',20);
saveas(fig1,['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-07-24/Figures/Pixel_Traces/' 'time_trace_at' num2str(2*pi*c/exc_center,3) '-' num2str(2*pi*c/det_center,3) '.png']);
fig2 = figure; plot(data.wnaxis,abs(point_trace_f),'LineWidth',3); xlim([0 950]);
title_string_f = sprintf('D1D2 a Exc at %g, Probed at %g vs W2',2*pi*c/exc_center,2*pi*c/det_center);
title(title_string_f,'FontSize',20);
xlabel('cm^{-1}');
saveas(fig2,['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-07-24/Figures/Pixel_Traces/' 'freq_trace_at' num2str(2*pi*c/exc_center,3) '-' num2str(2*pi*c/det_center,3) '.png']);
end