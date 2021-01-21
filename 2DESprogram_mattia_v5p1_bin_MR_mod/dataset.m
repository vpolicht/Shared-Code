function [ Map2DSave,T2In,wlpump,wprCutaver ] = dataset( hObject, eventdata, handles )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
T2In=getappdata(0, 'T2');
DataIn=getappdata(0, 'Data');
NumScans=getappdata(0, 'NumScans');
wlpump=getappdata(0, 'wlpump');
mStep=getappdata(0, 'mStep');
IndexPump=getappdata(0, 'IndexPump');
Bin0=str2double(get(handles.editBin0,'String'));
PhIntIn=getappdata(0, 'PhInt');
if (get(handles.checkboxNM,'Value') == get(handles.checkboxNM,'Max'))
wpr=getappdata(0, 'wpr');
else
    wpr=converter./getappdata(0, 'wpr');
end
ZeroPadding=str2double(get(handles.editZeroPadding,'String'));
    InterceptCal=str2double(get(handles.textCalibInter,'String'));
    SlopeCal=str2double(get(handles.CalibS,'String'));
    UnwrapPhase=1; %%user button from handles
    InvertPhase=0; %%%User button from handles
    padsize=[2^(nextpow2(length(mStep))+ZeroPadding)-length(mStep)];
    FTsize=length(mStep)+padsize;
    XX=linspace(1,FTsize,FTsize);
    XX=XX/FTsize;
    XX=fliplr(XX);
    if (get(handles.checkboxNM,'Value') == get(handles.checkboxNM,'Max'))
wlpump=(300000./(SlopeCal.*XX+InterceptCal));
else
   wlpump=converter./(300000./(SlopeCal.*XX+InterceptCal));
end
    
ProbeAverageNsquare=str2double(get(handles.editNumAverProbe,'String')); %from user
MinProbe=str2double(get(handles.editMinProbe,'String')); 
MaxProbe=str2double(get(handles.editMaxProbe,'String'));

if MinProbe<MaxProbe
    IndexProbe=find(wpr>=MinProbe & wpr<=MaxProbe);
else
    IndexProbe=find(wpr<=MinProbe & wpr>=MaxProbe);
end
 
wprCut=wpr(IndexProbe);
wprCutaver=zeros(floor(length(wprCut)/2^ProbeAverageNsquare),1);
DevStdInput=str2double(get(handles.editFWHMfilter,'String')); 

MaxPump=str2double(get(handles.editMaxPump,'String'));
MinPump=str2double(get(handles.ediMinPump,'String'));

if MaxPump>MinPump
IndexPump=find(wlpump>MinPump & wlpump<MaxPump);
else
    IndexPump=find(wlpump<MinPump & wlpump>MaxPump);
end

wlpump=wlpump(IndexPump);

Map2DSave=cell(length(T2In),1);
PhIntTrunc=cell(length(T2In),1);

    
for k=1:length(T2In)
%DataOut=DataIn{k}(Bin0+3:end,:); %%%selecting part of the matrix corrisponding to t1>0
DataOut=DataIn{k}(Bin0:end,:);
PhInt=PhIntIn{k}; %%%selected interferogram
spike_thrsh=str2double(get(handles.editSpikeThreshold,'String'));

%% FFT Unwrap.vi%%
    
    %PhIntPad= padarray(PhInt,padsize,'post');
    PhIntPad= padpost2(PhInt,padsize);
    PhIntPadRot=circshift(PhIntPad,-Bin0); 
    FFTPhIntPadRot=fft(PhIntPadRot, FTsize);
     if UnwrapPhase==1
        phase=unwrap(angle(FFTPhIntPadRot));
    end
    phase=angle(FFTPhIntPadRot);
    %phase=fliplr(phase);
    phase=flipud(phase);
    amplFFT=abs(FFTPhIntPadRot);
    if InvertPhase==1
        phase=(-1)*phase;
    end
    expIphase=exp(1i*phase);

    PhIntPadTrunc=PhIntPad(Bin0:end);
    FFTPhIntPadTrunc=fft(PhIntPadTrunc,FTsize).*expIphase;
    
    %% Probe Filter %%
Data=DataOut;
DataCut=Data(:,IndexProbe); %Data cut according to probe wavelength

%% Call function average_N.vi
DataCutaver=zeros(size(DataCut,1),floor(length(wprCut)/2^ProbeAverageNsquare));
for i=1:floor(length(wprCut)/2^ProbeAverageNsquare)
    wprCutaver(i)=mean(wprCut(i*2^ProbeAverageNsquare-2^ProbeAverageNsquare+1:i*2^ProbeAverageNsquare));
    DataCutaver(:,i)=mean(DataCut(:,i*2^ProbeAverageNsquare-2^ProbeAverageNsquare+1:i*2^ProbeAverageNsquare),2);
end
prwlFilt=wprCutaver;

%% Gaussian Filter

nOrder=str2double(get(handles.editOrderG,'String'));
DevStdInput_rel=DevStdInput*size(DataCutaver,1); %the program works with indexes so std deviation must be used realtive to the employed vector
    x=linspace(1,size(DataCutaver,1)*2,size(DataCutaver,1)*2);
    gaussianF=(exp(-(x.^2./(2*DevStdInput_rel^2)).^nOrder))';
    UniGaussWindow=ones(size(DataCutaver,1)*2,1).*gaussianF;
    UniGaussWindow=UniGaussWindow(length(UniGaussWindow)/2:length(UniGaussWindow)-1);
    Phase_Spectra=1; %from user
    if Phase_Spectra==1
        phase=exp(1i*phase);
    end
   
NoPhase=1; %from user
if NoPhase==1;
    phase=ones(length(phase),1);
end;

%% Ft loop

Data_in=DataCutaver;
Data_in=Data_in-repmat(mean(Data_in,1),size(Data_in,1),1);
Filter=UniGaussWindow;
for i=1:size(Data_in,2)
SpikeInd=find(abs(Data_in(:,i))>spike_thrsh);
Data_in(SpikeInd,i)=0;
Data_inFilter(:,i)=Data_in(:,i).*Filter; %%multiply each column for the filter

end


%Data_inFilterPad= padarray(Data_inFilter,padsize,'post');
Data_inFilterPad= padpost2(Data_inFilter,padsize);
%%
lp=size(Data_inFilterPad,2);
support=expIphase;
phase_padding=expIphase;
for w=1:lp-1
    phase_padding=cat(2,phase_padding,support);
end
%phase_padding=flipud(phase_padding);
phase_padding=-1.*phase_padding;
%%
FFTData_inFilterPad=fft(Data_inFilterPad,FTsize);%.*phase;
Data_average_Free_windowed=Data_in;
Map2D=FFTData_inFilterPad;
%% Pump Filter
DataCut2D=Map2D(IndexPump,:); %Data cut according to probe wavelength
DataCut2D=DataCut2D.*phase_padding(IndexPump,:);

%% Call function average_N.vi  %%%Nel programma non media sul pump

Pump_spectrum=amplFFT(IndexPump);
square_root=0;
NormalizedBySpectrum=0;
if square_root==1
    DataCut2D=sqrt(DataCut2D);
end
if (get(handles.checkboxNBP,'Value') == get(handles.checkboxNBP,'Max'))
    for i=1:size(DataCut2D,2)
     DataCut2D(:,i)=DataCut2D(:,i)./(Pump_spectrum);
    end
end
%%%Application of Filtering %%%
% DataCut2D_filt=sgolayfilt(real(DataCut2D)',3,7);
% DataCut2Dfilt=sgolayfilt(DataCut2D_filt',3,7);
DataCut2Dfilt=real(DataCut2D);


Map2DSave{k}=DataCut2Dfilt;
PhIntTrunc{k}=PhIntPadTrunc;


end

end

