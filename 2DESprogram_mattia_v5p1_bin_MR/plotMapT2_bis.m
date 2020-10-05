function plotMapT2(hObject,handles,converter)
T2value=str2double(get(handles.editT2,'String')); 
T2In=getappdata(0, 'T2');
DataIn=getappdata(0, 'Data');
%%
%DataIn=DataIn';
%%
indexT2=find(T2In>=T2value,1);
textLabel = sprintf('%.1f', T2In(indexT2));
set(handles.textT22D, 'String', textLabel);
NScans=getappdata(0, 'NumScans');
textLabel = sprintf('%.1f', NScans(indexT2));
set(handles.textNumScan, 'String', textLabel);
Bin0=str2double(get(handles.editBin0,'String'));
%DataOut=DataIn{indexT2}(Bin0+3:end,:); %%%selecting part of the matrix corrisponding to t1>0
DataOut=DataIn{indexT2}(Bin0:end,:);
%DataOut=-1.*DataOut;
PhIntIn=getappdata(0, 'PhInt');
PhInt=PhIntIn{indexT2}; %%%selected interferogram

if (get(handles.checkboxNM,'Value') == get(handles.checkboxNM,'Max'))
wpr=getappdata(0, 'wpr');
else
    wpr=converter./getappdata(0, 'wpr');
end
Bin0=str2double(get(handles.editBin0,'String'));
%% FFT Unwrap.vi%%
ZeroPadding=str2double(get(handles.editZeroPadding,'String'));
    InterceptCal=str2double(get(handles.textCalibInter,'String'));
    SlopeCal=str2double(get(handles.CalibS,'String'));
    UnwrapPhase=1; %%user button from handles
    InvertPhase=0; %%%User button from handles
    
    padsize=[2^(nextpow2(length(PhInt))+ZeroPadding)-length(PhInt)];
    PhIntPad= padpost2(PhInt,padsize);
    %PhIntPad= padarray(PhInt,padsize,'post');
    
    FTsize=length(PhIntPad);
    A=linspace(1,FTsize,FTsize);
    A=A/FTsize;
    A=fliplr(A);
    if (get(handles.checkboxNM,'Value') == get(handles.checkboxNM,'Max'))
wlpump=(300000./(SlopeCal.*A+InterceptCal));
else
   wlpump=converter./(300000./(SlopeCal.*A+InterceptCal));
    end
   

    
%PhIntPadRot=circshift(PhIntPad,-Bin0-3);
PhIntPadRot=circshift(PhIntPad,-Bin0);

    
    FFTPhIntPadRot=fft(PhIntPadRot, FTsize);
    %FFTPhIntPadRot=fliplr(FFTPhIntPadRot);
    %%%%%modifica
    %
%     FFTPhIntPadRot=fft(PhIntPadRot,2*FTsize);
%     FFTPhIntPadRot=FFTPhIntPadRot(1:FTsize+1);
%     FFTPhIntPadRot=FourierDir([1:FTsize],PhIntPadRot',[1:FTsize]);
%     FFTPhIntPadRot=FFTPhIntPadRot';
    %phase=angle(FFTPhIntPadRot);
    phase=angle(FFTPhIntPadRot);
     if UnwrapPhase==1
        phase=unwrap(angle(FFTPhIntPadRot));
        
     end
    %phase=angle(FFTPhIntPadRot);
    amplFFT=abs(FFTPhIntPadRot);
    phase=flipud(phase);
    if InvertPhase==1
        phase=(-1)*phase;
    end
     %phase=fliplr(phase);
    expIphase=exp(1i*phase);

    PhIntPadTrunc=PhIntPad(Bin0:end);
    FFTPhIntPadTrunc=fft(PhIntPadTrunc,FTsize).*expIphase;
    %FFTPhIntPadTrunc=fliplr( FFTPhIntPadTrunc);
%     FFTPhIntPadTrunc=fft(PhIntPadTrunc,2*FTsize);
%     FFTPhIntPadTrunc=FFTPhIntPadTrunc(1:FTsize+1).*expIphase;
    %FFTPhIntPadTrunc=FourierDir([1:FTsize],PhIntPadTrunc,[1:FTsize]).*expIphase;
    %% Probe Filter %%
ProbeAverageNsquare=str2double(get(handles.editNumAverProbe,'String')); %from user
MinProbe=str2double(get(handles.editMinProbe,'String')); 
MaxProbe=str2double(get(handles.editMaxProbe,'String'));



Data=DataOut;
if MinProbe<MaxProbe
    IndexProbe=find(wpr>=MinProbe & wpr<=MaxProbe);
else
    IndexProbe=find(wpr<=MinProbe & wpr>=MaxProbe);
end
DataCut=Data(:,IndexProbe); 
wprCut=wpr(IndexProbe);
wprCutaver=zeros(floor(length(wprCut)/2^ProbeAverageNsquare),1);
DataCutaver=zeros(size(DataCut,1),floor(length(wprCut)/2^ProbeAverageNsquare));
for i=1:floor(length(wprCut)/2^ProbeAverageNsquare)
    wprCutaver(i)=mean(wprCut(i*2^ProbeAverageNsquare-2^ProbeAverageNsquare+1:i*2^ProbeAverageNsquare));
    DataCutaver(:,i)=mean(DataCut(:,i*2^ProbeAverageNsquare-2^ProbeAverageNsquare+1:i*2^ProbeAverageNsquare),2);
end
setappdata(0, 'wprCutaver' ,wprCutaver);
prwlFilt=wprCutaver;
Probe_resolution=prwlFilt(2)-prwlFilt(1); %%%print on GUI
%% Gaussian Filter
if get(handles.gau,'Value') == get(handles.gau,'Max')
DevStdInput=str2double(get(handles.editFWHMfilter,'String')); 
nOrder=str2double(get(handles.editOrderG,'String'));    
DevStdInput_rel=DevStdInput*size(DataCutaver,1); %the program works with indexes so std deviation must be used realtive to the employed vector
    x=linspace(1,size(DataCutaver,1)*2,size(DataCutaver,1)*2);
   gaussianF=(exp(-(x.^2./(2*DevStdInput_rel^2)).^nOrder))';
    UniGaussWindow=ones(size(DataCutaver,1)*2,1).*gaussianF;
    UniGaussWindow=UniGaussWindow(length(UniGaussWindow)/2:length(UniGaussWindow)-1);
end 
   if get(handles.bh,'Value') == get(handles.bh,'Max')
        x=linspace(1,size(DataCutaver,1)*2,size(DataCutaver,1)*2);
        UniGaussWindow=blackmanharris(length(x));
        UniGaussWindow=UniGaussWindow(length(UniGaussWindow)/2:length(UniGaussWindow)-1);
   end

   if get(handles.han,'Value') == get(handles.han,'Max')
       DevStdInput=str2double(get(handles.editFWHMfilter,'String'));
        x=linspace(1,size(DataCutaver,1)*2,size(DataCutaver,1)*2);
        UniGaussWindow=tukeywin(length(x));
        tfun=bartlett(round(length(x)*DevStdInput));
        if DevStdInput<=1
        tfun=padarray(tfun,abs(length(x)-length(tfun)),'post');
        else
                tfun=tfun(1:length(x));
        end
         UniGaussWindow=UniGaussWindow.*tfun;
        UniGaussWindow=UniGaussWindow(length(UniGaussWindow)/2:length(UniGaussWindow)-1);
   end
   
    Phase_Spectra=0; %from user
    if Phase_Spectra==1
        phase=exp(1i*phase);
    end
   
NoPhase=0; %from user
if NoPhase==1;
    phase=ones(length(phase),1);
end;

%% Ft loop

Data_in=DataCutaver;
Data_in=Data_in-repmat(mean(Data_in,1),size(Data_in,1),1);
setappdata(0, 'DataAverageFree' ,Data_in);
Xcursor=str2double(get(handles.editXcursor,'String'));
IndexXcursor=find(wprCutaver>=Xcursor,1);
axes(handles.axesXcursor)
if (get(handles.checkboxFS,'Value') == get(handles.checkboxFS,'Max'))
    x=linspace(0,length(Data_in(:,IndexXcursor)),length(Data_in(:,IndexXcursor)))*250/length(Data_in(:,IndexXcursor));
    plot(x,Data_in(:,IndexXcursor)*100,'LineWidth',1.5)  
hold on
else
plot(Data_in(:,IndexXcursor)*100,'LineWidth',1.5)
hold on
end
spike_thrsh=str2double(get(handles.editSpikeThreshold,'String'));
Filter=UniGaussWindow./max(UniGaussWindow);
index=50;
for i=1:size(Data_in,2)
SpikeInd=find(abs(Data_in(:,i))>spike_thrsh);
Data_in(SpikeInd,i)=0;
Data_inFilter(:,i)=Data_in(:,i).*Filter; %%multiply each column for the filter

end

%% Smooth over t1%
nSmooth=str2double(get(handles.editNSmooth,'String'));
if (get(handles.checkboxSmoothing,'Value') == get(handles.checkboxSmoothing,'Max'))
%Data_inFilter=medfilt1(Data_inFilter,nSmooth);
Data_inFilter=medfiltmr(Data_inFilter,nSmooth);
end
%%
setappdata(0, 'DataAverageFreeWindowed' ,Data_inFilter);

axes(handles.axesXcursor)
if (get(handles.checkboxFS,'Value') == get(handles.checkboxFS,'Max'))
    plot(x,Data_inFilter(:,IndexXcursor)*100,'r','LineWidth',1.5)
    xlabel('t_1 [fs]');
hold on
else
plot(Data_inFilter(:,IndexXcursor)*100,'r','LineWidth',1.5)
xlabel('t_1 [Motor Steps]');
end
ylabel('\Delta T/T [%]');
legend('Signal average free','Signal average free and windowed')
set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
hold off

padsize=[2^(nextpow2(size(Data_inFilter,1))+ZeroPadding)-size(Data_inFilter,1)];
%Data_inFilterPad= padarray(Data_inFilter,padsize,'post');
Data_inFilterPad= padpost2(Data_inFilter,padsize);
%%%%
lp=size(Data_inFilterPad,2);
support=expIphase;
phase_padding=expIphase;
for w=1:lp-1
    phase_padding=cat(2,phase_padding,support);
end
%phase_padding=flipud(phase_padding);
phase_padding=-1.*phase_padding;
%%%%



%FFTData_inFilterPad=fft(Data_inFilterPad,FTsize);%.*phase;
FFTData_inFilterPad=fft(Data_inFilterPad,FTsize);
% FFTData_inFilterPad=fft(Data_inFilterPad,2*FTsize);
% FFTData_inFilterPad=FFTData_inFilterPad(1:FTsize+1);
%FFTData_inFilterPad=FourierDir([1:FTsize],Data_inFilterPad,[1:FTsize]);
Data_average_Free_windowed=Data_in;
RealFFTData_inFilterPad=real(FFTData_inFilterPad);
Real_cursor=real(FFTData_inFilterPad(index,:));
Imag_cursor=imag(FFTData_inFilterPad(index,:));
Map2D=FFTData_inFilterPad;
%% Pump Filter

MaxPump=str2double(get(handles.editMaxPump,'String'));
MinPump=str2double(get(handles.ediMinPump,'String'));
if MaxPump>MinPump
IndexPump=find(wlpump>MinPump & wlpump<MaxPump);
else
    IndexPump=find(wlpump<MinPump & wlpump>MaxPump);
end
DataCut2D=Map2D(IndexPump,:); %Data cut according to probe wavelength
%%%%%%%%% moltiplico per la fase qui
DataCut2D=DataCut2D.*phase_padding(IndexPump,:);
%setappdata(0, 'IndexPump' ,IndexPump);
%% Call function average_N.vi  %%%Nel programma non media sul pump
wlpump=wlpump(IndexPump);
%setappdata(0, 'wlpump' ,wlpump);

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
if get(handles.buttondis,'Value') == get(handles.buttondis,'Max')
%DataCut2D_filt=sgolayfilt(imag(DataCut2D)',3,7);
DataCut2D_filt=imag(DataCut2D);
end
if get(handles.areal,'Value') == get(handles.areal,'Max')
%DataCut2D_filt=sgolayfilt(real(DataCut2D)',3,7);
DataCut2D_filt=real(DataCut2D);
end
if get(handles.buttonmag,'Value') == get(handles.buttonmag,'Max')
%DataCut2D_filt=sgolayfilt(abs(DataCut2D)',3,7);
DataCut2D_filt=abs(DataCut2D);
end
%DataCut2Dfilt=sgolayfilt(DataCut2D_filt',3,7);
DataCut2Dfilt=DataCut2D_filt;
%%%%%%%%%%%%%%%%%%%%% here the shift
speed=300000000;
sh=str2double(get(handles.edit36,'String'));
sh=sh*1000000000000;
f1=(speed./wlpump).*1000000000;
f2=f1+sh;
newpump=(speed./f2)*1000000000;
wlpump=newpump;
%%%%%%%%%%%%%%%%%%%%%%
if (get(handles.checkboxContour,'Value') == get(handles.checkboxContour,'Max'))
    axes(handles.axes2Dmap)
contourf(wprCutaver,wlpump,DataCut2Dfilt,20)
colormap bone
%shading flat
hold on
v = caxis;
textLabel = sprintf('%.3f', v(1));
set(handles.editCMmin, 'String', textLabel);
textLabel = sprintf('%.3f', v(2));
set(handles.editCMmax, 'String', textLabel);
%[C,hfigc] = contour(wprCutaver, wlpump, DataCut2Dfilt);

% if (get(handles.checkboxNM,'Value') == get(handles.checkboxNM,'Max'))
% g=line([300 1000],[300 1000],'Color','k','LineStyle','--','LineWidth',3);
% else
%   g=line([converter/300 converter/1000],[converter/300 converter/1000],'Color','k','LineStyle','--','LineWidth',3);
% end

%set(hfigc, ...
%    'LineWidth',1.0, ...
 %   'Color', [1 1 1]);


xlabel('\lambda_{probe} [nm]');
ylabel('\lambda_{pump} [nm]'); 
set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
colorbar('southoutside')
hold off
else
axes(handles.axes2Dmap)
pcolor(wprCutaver,wlpump,DataCut2Dfilt)
colormap jet
shading interp
hold on
v = caxis;
textLabel = sprintf('%.3f', v(1));
set(handles.editCMmin, 'String', textLabel);
textLabel = sprintf('%.3f', v(2));
set(handles.editCMmax, 'String', textLabel);
[C,hfigc] = contour(wprCutaver, wlpump, DataCut2Dfilt);

if (get(handles.checkboxNM,'Value') == get(handles.checkboxNM,'Max'))
g=line([300 1000],[300 1000],'Color','k','LineStyle','--','LineWidth',3);
else
  g=line([converter/300 converter/1000],[converter/300 converter/1000],'Color','k','LineStyle','--','LineWidth',3);
end

set(hfigc, ...
    'LineWidth',1.0, ...
    'Color', [1 1 1]);


xlabel('\lambda_{probe} [nm]');
ylabel('\lambda_{pump} [nm]'); 
set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
colorbar('southoutside')
hold off
end

PP2D=mean(DataCut2Dfilt,1);
axes(handles.axesPP)
hold on
plot(wprCutaver,PP2D./max(abs(PP2D)),'m','LineWidth',1.5);
set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
box on



