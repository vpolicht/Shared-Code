function [a_slice, d_slice, slice_a_axis, slice_d_axis] = CutSlice(data,exc_freq,det_freq)
assert(exc_freq<max(data.x_axis_w_roi) && exc_freq>min(data.x_axis_w_roi));
assert(det_freq<max(data.y_axis_w_roi) && det_freq>min(data.y_axis_w_roi));


%given a center point detemined by exc_freq, det_freq, find the
%diagonal and anti-diagonal lines which pass through this point.
line_eq = @(x,x1,y1,m) m*(x-x1)+y1;
det_line_diag = line_eq(data.x_axis_w_roi,exc_freq,det_freq,1);
det_line_adiag = line_eq(data.x_axis_w_roi,exc_freq,det_freq,-1);
exc_line_diag_inds = 1:length(data.x_axis_w_roi);
exc_line_adiag_inds = 1:length(data.x_axis_w_roi);
%search for the closest detected pixel to the desired detection point
det_line_diag_inds = 0*data.x_axis_w_roi; 
det_line_adiag_inds = 0*data.x_axis_w_roi;
for k=1:length(data.x_axis_w_roi);
    det_line_diag_inds(k) = dsearchn(data.y_axis_w_roi',det_line_diag(k));
    det_line_adiag_inds(k) = dsearchn(data.y_axis_w_roi',det_line_adiag(k));
end  
%figure out which requested points are out of bounds on the detection axis
det_line_diag_oob = ~(det_line_diag<max(data.y_axis_w_roi) & det_line_diag>min(data.y_axis_w_roi));
det_line_adiag_oob = ~(det_line_adiag<max(data.y_axis_w_roi) & det_line_adiag>min(data.y_axis_w_roi));
det_line_diag_inds(det_line_diag_oob) = [];
exc_line_diag_inds(det_line_diag_oob) = [];
det_line_adiag_inds(det_line_adiag_oob) = [];
exc_line_adiag_inds(det_line_adiag_oob) = [];
d_slice = zeros(length(det_line_diag_inds),size(data.data,3));
a_slice = zeros(length(det_line_adiag_inds),size(data.data,3));
for k=1:length(det_line_diag_inds); d_slice(k,:) = squeeze(data.data(det_line_diag_inds(k),exc_line_diag_inds(k),:)); end;
for k=1:length(det_line_adiag_inds); a_slice(k,:) = squeeze(data.data(det_line_adiag_inds(k),exc_line_adiag_inds(k),:)); end;
% figure; plot(data.x_axis_w_roi(exc_line_diag_inds),data.y_axis_w_roi(det_line_diag_inds));
% figure; plot(data.x_axis_w_roi(exc_line_adiag_inds),data.y_axis_w_roi(det_line_adiag_inds));
slice_a_axis = data.x_axis_w_roi(exc_line_adiag_inds);
slice_d_axis = data.x_axis_w_roi(exc_line_diag_inds);
end