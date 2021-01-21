function [] = PickPointWithRadiusOverlayData(data,freq_1,freq_2)
assert(freq_1 ~= freq_2);
if freq_1 < freq_2; 
    freq_low = freq_1; freq_high = freq_2;
else
    freq_low = freq_2; freq_high = freq_1;
end
c = 299.792458;
freq_low_ind = dsearchn(data.x_axis_w_roi',freq_low);
freq_high_ind = dsearchn(data.y_axis_w_roi',freq_high);
point_trace_t_r_DP = squeeze(sum(sum(data.zr((freq_high_ind-1):(freq_high_ind+1),(freq_low_ind-1):(freq_low_ind+1),:),2),1));
point_trace_t_nr_DP = squeeze(sum(sum(data.znr((freq_high_ind-1):(freq_high_ind+1),(freq_low_ind-1):(freq_low_ind+1),:),2),1));
point_trace_t_r_DM = squeeze(sum(sum(data.zr((freq_low_ind-1):(freq_low_ind+1),(freq_high_ind-1):(freq_high_ind+1),:),2),1));
point_trace_t_nr_DM = squeeze(sum(sum(data.znr((freq_low_ind-1):(freq_low_ind+1),(freq_high_ind-1):(freq_high_ind+1),:),2),1));
point_trace_t_r_DH = squeeze(sum(sum(data.zr((freq_high_ind-1):(freq_high_ind+1),(freq_high_ind-1):(freq_high_ind+1),:),2),1));
point_trace_t_nr_DH = squeeze(sum(sum(data.znr((freq_high_ind-1):(freq_high_ind+1),(freq_high_ind-1):(freq_high_ind+1),:),2),1));
point_trace_t_r_DL = squeeze(sum(sum(data.zr((freq_low_ind-1):(freq_low_ind+1),(freq_low_ind-1):(freq_low_ind+1),:),2),1));
point_trace_t_nr_DL = squeeze(sum(sum(data.znr((freq_low_ind-1):(freq_low_ind+1),(freq_low_ind-1):(freq_low_ind+1),:),2),1));
fig1 = figure; 
plot(data.T,point_trace_t_r_DP,'g'); hold all; plot(data.T,point_trace_t_nr_DP,'r'); %%%
offset1 = 4*max(std(point_trace_t_r_DP,[],1),std(point_trace_t_nr_DP,[],1));
plot(data.T,point_trace_t_r_DH + offset1,'k-'); plot(data.T,point_trace_t_nr_DH + offset1,'k-.'); %%%
offset2 = 4*max(std(point_trace_t_r_DH,[],1),std(point_trace_t_nr_DH,[],1)) + offset1;
plot(data.T,point_trace_t_r_DL + offset2,'b-'); plot(data.T,point_trace_t_nr_DL + offset2,'b-.'); %%%
offset3 = 4*max(std(point_trace_t_r_DL,[],1),std(point_trace_t_nr_DL,[],1)) + offset2;
plot(data.T,point_trace_t_r_DM + offset3,'g'); plot(data.T,point_trace_t_nr_DM + offset3,'r'); %%%
hold off;
title_string_t = sprintf('Chlorophyll a at %g x %g vs Population Time',2*pi*c/freq_low,2*pi*c/freq_high);
xlabel('femtoseconds');
title(title_string_t);
set(gca,'XGrid','on')
legend({'R upper CP','NR upper CP','R diag upper','NR diag upper','R diag lower','NR diag lower','R lower CP','NR lower CP'},'Location','NorthEastOutside');
saveas(fig1,['/mnt/data/2D/Frank/LyX_Lab_Notebook/LaserLab/13-07-06/Figures/Pixel_Traces/' 'composite_time_trace_at' num2str(2*pi*c/freq_low,3) '-' num2str(2*pi*c/freq_high,3) '.png']);

end
