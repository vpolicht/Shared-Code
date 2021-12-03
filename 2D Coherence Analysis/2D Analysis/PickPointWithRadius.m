function [point_trace_t,point_trace_f] = PickPointWithRadius(data,exc_center,det_center)
c = 299.792458;
exc_center_ind = dsearchn(data.x_axis_w_roi',exc_center);
det_center_ind = dsearchn(data.y_axis_w_roi',det_center);
point_trace_t = squeeze(sum(sum(data.z((det_center_ind-1):(det_center_ind+1),(exc_center_ind-1):(exc_center_ind+1),:),2),1));
point_trace_f = abs(fft(point_trace_t,2^12));
fig1 = figure; plot(data.T(1:end),point_trace_t,'g');
title_string_t = sprintf('d1d2 a Excited at %g, Probed at %g vs Population Time',2*pi*c/exc_center,2*pi*c/det_center);
xlabel('femtoseconds');
title(title_string_t);
saveas(fig1,['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-09-04/Figures/Pixel_Traces/' 'nonrephasing_time_trace_at' num2str(2*pi*c/exc_center,3) '-' num2str(2*pi*c/det_center,3) '.png']);
fig2 = figure; plot(data.wnaxis,abs(point_trace_f),'g'); xlim([0 950]);
title_string_f = sprintf('d1d2 a Excited at %g, Probed at %g vs Population Frequency',2*pi*c/exc_center,2*pi*c/det_center);
title(title_string_f);
xlabel('cm^{-1}');
saveas(fig2,['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-09-04/Figures/Pixel_Traces/' 'nonrephasing_freq_trace_at' num2str(2*pi*c/exc_center,3) '-' num2str(2*pi*c/det_center,3) '.png']);
end