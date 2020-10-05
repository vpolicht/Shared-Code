function calibration(hObject,handles,pat)

a=regexp(pat,'\');
ind=pat(1:a(end));
file =dir(pat);
pati=[ind,file.name];
dass=dlmread(pati);
DataCal=dass(4:end,3:end);
PhIntCal=dass(3:end,2);
wlCal=dass(1,3:end);
PumpSpec=dass(2,3:end); %#ok<*NASGU>
mStepCal=dass(3:end,1);
set(hObject,'BackgroundColor','white')
%% Smooth of each column%%%%
ntimes=3;  %#number of points for the smoothing
for i=1:size(DataCal,2)
   %DataCal(:,i)=smooth(DataCal(:,i),ntimes);
   DataCal(:,i)=smoothmr(DataCal(:,i),ntimes);
end
%% Nsquare_Average.vi%%
CaliAverageNsquare=str2double(get(handles.editNAverCal,'String')); %from user
wlCalaver=zeros(floor(length(wlCal)/2^CaliAverageNsquare),1);
PumpSpecaver=zeros(floor(length(wlCal)/2^CaliAverageNsquare),1);
DataCalaver=zeros(size(DataCal,1),floor(length(wlCal)/2^CaliAverageNsquare));
for i=1:floor(length(wlCal)/2^CaliAverageNsquare)
    wlCalaver(i)=mean(wlCal(i*2^CaliAverageNsquare-2^CaliAverageNsquare+1:i*2^CaliAverageNsquare));
    DataCalaver(:,i)=mean(DataCal(:,i*2^CaliAverageNsquare-2^CaliAverageNsquare+1:i*2^CaliAverageNsquare),2);
    PumpSpecaver(i)=mean(PumpSpec(i*2^CaliAverageNsquare-2^CaliAverageNsquare+1:i*2^CaliAverageNsquare));
end
%% Resize of Data in order to reduce noise%%%%
%%%Resize Motor step%%%%
A=length(mStepCal)*0.1; 
B= length(mStepCal)*0.8; 
mStepCal=mStepCal(round(A)+1:round(A)+round(B)-1);
%%%Resize Wavelength
Threshold=str2double(get(handles.editThreshold,'String')); %%express in percentage
KK=(max(wlCal)-min(wlCal))*Threshold+min(wlCal);
index_min=find(PumpSpec>KK,1);   %%%%%Inserire il controllo
index_max=find(fliplr(PumpSpec)>KK,1);
 axes(handles.axesPumpSpectrometer)
plot(wlCal,PumpSpec,'r','LineWidth',1.5)
if (get(handles.uv,'Value') == get(handles.uv,'Max'))
xlim([300 420])
else
xlim([480 750])
end
xlabel('\lambda_{spectrometer} [nm]');
          ylabel('Counts');
         set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
         set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
          set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
PumpSpec=PumpSpec(index_min:length(PumpSpec)-index_max);          
wlCal=wlCal(index_min:length(PumpSpec)-index_max);
%%%Resize Data
DataCal=DataCal(round(A)+1:round(A)+round(B)-1,index_min:length(PumpSpec)-index_max);
%% FFT Data%%%%%
DataCal(isnan(DataCal)) = 0;
DataCal=DataCal-ones(size(DataCal,1),1)*mean(DataCal,1);  %%%Possibilità moltiplicare con finestra rettagolare
%sizeFFT=2^nextpow2(size(DataCal,1));
sizeFFT=10000;
FFTDataCal=abs(fft(DataCal,sizeFFT,1));
VarRem=0.02;
FFTDataCal=FFTDataCal(round(VarRem*sizeFFT)+1:round(sizeFFT/2)+round(VarRem*sizeFFT),:);
[A,B]=max(FFTDataCal);
MaxCorr=B+VarRem*sizeFFT;
OptFreq=300000./wlCal;
LL=abs(mStepCal(1)-mStepCal(end));
PP=1/(LL/length(mStepCal)*sizeFFT); %frequency resolution
PseudoFreq=MaxCorr*PP; % X data [step^-1]
%OptFreq=300000./wlCal; %%Y data [ThZ]
%%%%%%%%%%%%%%%%%modifica
if (get(handles.uv,'Value') == get(handles.uv,'Max'))
[OptFreq,PseudoFreq ]=rescaleCali( wlCal,MaxCorr,PP);
else
    %[OptFreq,PseudoFreq ]=rescaleCali2( wlCal,MaxCorr,PP);
    OptFreq=300000./wlCal; %%Y data [ThZ]
end
%%%%%%%%%%%%%%%%%
axes(handles.axesCalibration)
plot(PseudoFreq,OptFreq,'o')
P = polyfit(PseudoFreq,OptFreq,1);
hold on
plot(PseudoFreq, P(1)*PseudoFreq+P(2),'r','LineWidth',1.5) 
         xlabel('Pseudo Frequencies [steps^{-1}]');
          ylabel('Optical Frequencies [THz]');
         set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
         set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
          set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
          %%% --> Inserire un filtro Nutterworth su PhilCal con 0.008 e 0.2 cut off Frequencies
 
          textLabel = sprintf('%.1f', P(1));
set(handles.CalibS, 'String', textLabel);
textLabel = sprintf('%.1f', P(2));
set(handles.textCalibInter, 'String', textLabel);
axes(handles.axesCorrMap)
pcolor(DataCal')
shading flat
colormap winter
xlabel('t_1[Wedges steps]');
          ylabel('Counts');
         set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
         set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
          set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
axes(handles.axesPhotodiode)
    plot(PhIntCal)
    xlabel('t_1 [Wedges Steps]');
          ylabel('Amplitude');
          ylim([-2000 2000])
         set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
         set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
          set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
%% Phase Button%%
 PhaseButton=1; %%%Input user
 DevStdInput=str2double(get(handles.editGaussCal,'String'));  %%%%Input user
 DevStdInput_rel=DevStdInput*length(PhIntCal);
    %% Gaussian Filter%%%%
        x=linspace(-length(PhIntCal),length(PhIntCal),length(PhIntCal)*2);
    %UniGaussWindow=ones(length(PhIntCal)*2,1).*gaussmf(x,[DevStdInput_rel 0])';
    gauss=exp(-(x-0).^2/(2*DevStdInput_rel^2));
    UniGaussWindow=ones(length(PhIntCal)*2,1).*gauss';
    UniGaussWindow=UniGaussWindow(length(UniGaussWindow)/2:length(UniGaussWindow)-1);    
    UniGaussWindow=UniGaussWindow./max(UniGaussWindow);
    yfilt= UniGaussWindow.*PhIntCal;
    PhIntCalFilt=yfilt;
     %% Butterworth filter design order=2, bandpass (0.008-0.2)   
%  [b,a] = butter(2,[0.008 0.2],'bandpass');
b=[0.0630 0 -0.1259 0 0.0630];
a=[1 -3.1497 3.7363 -2.0139 0.4274];
PhIntCalFilt = filter(b,a,PhIntCalFilt);
    %% FFT phasing unwrap.vi%%%%
    ZeroPadding=str2double(get(handles.edit0Pad,'String')); %%Input User
    padsize=[2^(nextpow2(length(PhIntCalFilt))+ZeroPadding)-length(PhIntCalFilt)];
    %PhIntCalFiltPad= padarray(PhIntCalFilt,padsize,'post');
    PhIntCalFiltPad= padpost(PhIntCalFilt,padsize);
    
    %% Now Call pump_axis.vi%%%
     InterceptCal=P(2);
     SlopeCal=P(1);
    FTsize=length(PhIntCalFiltPad);
    A=linspace(1,FTsize,FTsize);
    A=A/FTsize;
    A=fliplr(A);
    wlpump=300000./(SlopeCal.*A+InterceptCal); %%%PumpAxis [nm]
    save('pippo.mat','wlpump')

    if PhaseButton==1
        %% FT_mag_phase_unwrap.vi%%%%
    UnwrapPhase=1; %%user button
    Bin0=str2double(get(handles.editBin0,'String')); %%user choice
    InvertPhase=0; %%%User button
    
    PhIntCalFiltPadRot=circshift(PhIntCalFiltPad,-Bin0);
    axes(handles.axesPhotodiode)
    plot(PhIntCalFiltPadRot)
    xlabel('t_1 [Wedges Steps]');
          ylabel('Amplitude');
          ylim([-2000 2000])
         set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
         set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
          set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
    FFTPhIntCalFiltPadRot=fft(PhIntCalFiltPadRot, FTsize);
     if UnwrapPhase==1
        phase=unwrap(angle(FFTPhIntCalFiltPadRot));
    end
    phase=angle(FFTPhIntCalFiltPadRot);
    phase=fliplr(phase);
    amplFFT=abs(FFTPhIntCalFiltPadRot);
    if InvertPhase==1
        phase=(-1)*phase;
    end
    phase=phase;
    expIphase=exp(1i*phase);
    axes(handles.axesAmplPh)
    %plotyy(wlpump,amplFFT,wlpump,phase)
    plot(wlpump,amplFFT,'LineWidth',1.5)
    if (get(handles.uv,'Value') == get(handles.uv,'Max'))
    xlim([300 420])
    else
    xlim([480 750])
    end
    xlabel('\lambda [nm]');
          ylabel('Amplitude');
         set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
         set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
          set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
    axes(handles.axesPhasePh)
   plot(wlpump,phase,'LineWidth',1.5)
    if (get(handles.uv,'Value') == get(handles.uv,'Max'))
    xlim([300 420])
    else
    xlim([480 750])
    end
    xlabel('\lambda [nm]');
          ylabel('Phase');
          
         set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
         set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
          set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
    %% FT_ABS_DISP%%%%
    PhIntCalFiltPadTrunc=PhIntCalFiltPad(Bin0:end);
    KK=fft(PhIntCalFiltPadTrunc,FTsize).*expIphase;

%% Calculation of Resolution%%% ---> Fare un check
    [AverageFreq,GroupDelay, Sampling]=Calc_Reso(OptFreq,InterceptCal, SlopeCal,length(PhIntCalFiltPadRot));
 %%nm Output  Fare un check su sampling --> fare un check di analisi dimensinale dovrebbe venire Sampling=3.3 circa AverageWl=650 Resolution=11 T1max=121                                                                                       
Average_Wavelength=300000/AverageFreq; %%nm Output
NumpointScan=length(PhIntCalFilt)-Bin0;
[Resolution,t1max, Reso]=Calc_Reso(OptFreq,InterceptCal, SlopeCal,NumpointScan);

%% Plot in main GUI%%%
HighWlFit=700; %%input
LowWlFit=560; %%input
[wlpumpCut, interceptPhase, SlopePhase, phasecut]=linearphaseFit(wlpump,amplFFT, phase,HighWlFit,LowWlFit);
    end
    %% Create global%%
    setappdata(0, 'PhIntCalFiltPad' ,PhIntCalFiltPad);
    setappdata(0, 'FTsize' ,FTsize);%storing matrix in GUI variables
    setappdata(0, 'wlpump' ,wlpump);%storing matrix in GUI variables