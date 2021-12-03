function [ trace ] = RetrieveAveragedPixelTrace( D, x_axis, y_axis, point_vec, box_size )
exc_ind = dsearchn(x_axis',point_vec(1));
det_ind = dsearchn(y_axis',point_vec(2));
box_size = box_size - mod(box_size,2);
trace = squeeze(mean(mean(D((det_ind-box_size/2):(det_ind+box_size/2),(exc_ind-box_size/2):(exc_ind+box_size/2),:),1),2));
end

