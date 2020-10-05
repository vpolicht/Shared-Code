function [Data_average_Free_windowed,Map2D,Real_cursor,Imag_cursor]=FT_C(handles,pump_index,FFTsize,Data_in,Filter,Phase,Abs_disp_mag,spike_thrsh,activate_pr)
Zero_Padding=str2double(get(handles.editZeroPadding,'String')); 
Data_in=Data_in-repmat(mean(Data_in,1),size(Data_in,1),1);

for i=1:size(Data_in,2)
SpikeInd=find(abs(Data_in(:,i))>spike_thrsh);
Data_in(SpikeInd,i)=0;
Data_inFilter(:,i)=Data_in(:,i).*Filter; %%multiply each column for the filter
end
ZeroPadding=str2double(get(handles.editZeroPadding,'String')); %%Input User
padsize=[ZeroPadding*2^(nextpow2(size(Data_inFilter,1)))-size(Data_inFilter,1)];
Data_inFilterPad= padarray(Data_inFilter,padsize,'post');

Map2D=fft(Data_inFilterPad,FFTsize);%.*Phase;
Data_average_Free_windowed=Data_in;
RealFFTData_inFilterPad=real(Map2D);
Real_cursor=real(Map2D(pump_index,:));
Imag_cursor=imag(Map2D(pump_index,:));

% 