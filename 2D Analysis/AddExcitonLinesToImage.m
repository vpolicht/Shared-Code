function [] = AddExcitonLinesToImage(fig_handle,x_axis,y_axis,line_locations)
x_line_locs = 0*line_locations;
y_line_locs = 0*line_locations;
for k=1:length(line_locations); 
    x_line_locs(k) = dsearchn(x_axis',line_locations(k));
    y_line_locs(k) = dsearchn(y_axis',line_locations(k));
end
figure(fig_handle);
dummy_x = 0*x_axis+1;
dummy_y = 0*y_axis+1;
for k=1:length(line_locations);
    line(x_axis(x_line_locs(k))*dummy_y,y_axis,'linewidth',2,'Color','r');
    line(x_axis,y_axis(y_line_locs(k))*dummy_x,'linewidth',2,'Color','r');
end