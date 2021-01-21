function [wprCutaver, DataCutaver]=probe_filter(handles,hObject,wpr, DataOut)
ProbeAverageNsquare=str2double(get(handles.editNumAverProbe,'String')); %from user
MinProbe=str2double(get(handles.editMinProbe,'String')); 
MaxProbe=str2double(get(handles.editMaxProbe,'String'));  
%%Check MaxProbe>MinProbe already in handles dynamic text%%%
A=find(wpr>=MinProbe,1); 
B=find(wpr>=MaxProbe,1);
DataCut=DataOut(:,A:B); %Data cut according to probe wavelength
%% Call function average_N.vi
wprCut=wpr(A:B);
NumCol=size(DataCut,2);
wprCutaver=zeros(length(wprCut)/2^ProbeAverageNsquare,1);
DataCutaver=zeros(size(DataCut,1),length(wprCut)/2^ProbeAverageNsquare);
for i=1:round(length(wprCut)/2^ProbeAverageNsquare)
    wprCutaver(i)=mean(wprCut(i*2^ProbeAverageNsquare-2^ProbeAverageNsquare+1:i*2^ProbeAverageNsquare));
    DataCutaver(:,i)=mean(DataCut(:,i*2^ProbeAverageNsquare-2^ProbeAverageNsquare+1:i*2^ProbeAverageNsquare),2);
end
