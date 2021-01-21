function calibration(hObject,handles,pat)
% Modified by VP
tic
if ispc
    a=regexp(pat,'\');
else
    a=regexp(pat,'/');
end

ind = pat(1:a(end));
file = dir(pat);
pati = [ind,file.name];
dass = dlmread(pati);
DataCal = dass(4:end,3:end);
PhIntCal = dass(4:end,2); %changed from 3:end
wlCal = dass(1,3:end); 
PumpSpec = dass(2,3:end); %#ok<*NASGU>
mStepCal = dass(4:end,1); %changed from 3:end
set(hObject,'BackgroundColor','w')

%% Smooth of each column%%%%
ntimes = 3;  %#number of points for the smoothing

for m = 1:size(DataCal,2)
   DataCal(:,m) = smooth(DataCal(:,m),ntimes);
end

%% removed averaging over w1

wlCalaver = wlCal;
PumpSpecaver = PumpSpec;
DataCalaver = DataCal;
%% Resize of Data in order to reduce noise

% Remove noisy edges of calibration scan
A = round(length(mStepCal)*0.035); % Changed from 0.1
B = round(length(mStepCal)*0.95); % Changed from 0.8, seems pretty arbitrary
mStepCal = mStepCal(A+1:A+B-1);

% Resize Wavelength
Threshold = str2double(get(handles.editThreshold,'String')); %%express in percentage
KK = max(PumpSpec)*Threshold;
index_min = find(PumpSpec > KK,1);   %%%%%Inserire il controllo
index_max = find(fliplr(PumpSpec) > KK,1);

axes(handles.axesPumpSpectrometer);
plot(wlCal,PumpSpec,'r','LineWidth',1.5)

if (get(handles.uv,'Value') == get(handles.uv,'Max'))
    xlim([300 420])
else
   if (get(handles.but_OPA3,'Value') == get(handles.but_OPA3,'Max'))
        xlim([650 1000])
    else 
        xlim([450 800])
   end
end

xlabel('\lambda_{spectrometer} [nm]');
ylabel('Counts');
set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)

PumpSpec = PumpSpec(index_min:length(PumpSpec)-index_max);          
wlCal = wlCal(index_min:length(PumpSpec)-index_max);

% Resize Data
DataCal = DataCal(A+1:A+B-1,index_min:length(PumpSpec)-index_max);
DataCal(isnan(DataCal)) = 0;
DataCal = DataCal - mean(DataCal,1);% Subtract pump spectrum

%% FFT Data

sizeFFT = 2^14; % Changed from 10000;
FFTDataCala = (fft(DataCal,sizeFFT,1));

FFTDataCal = abs(FFTDataCala(1:sizeFFT/2+1,:));
[A,B] = max(FFTDataCal); % Find peak along pseudo frequency as function of optical frequency

MaxCorr = B;% + VarRem*sizeFFT;
LL = abs(mStepCal(1)-mStepCal(end));
PP = 1/(LL/length(mStepCal)*sizeFFT); % pseudo frequency resolution
PseudoFreq = MaxCorr*PP; % X data [step^-1]

difdif = find(diff(MaxCorr) > mean(abs(diff(MaxCorr))));
if ~isempty(difdif)
    fitind = difdif(end):length(MaxCorr);
else
    fitind = 1:length(MaxCorr);
end

dx = mean(abs(diff(mStepCal))); % Motor Step resolution
f_fft = ((1/dx)/2).*linspace(0,1,sizeFFT/2+1); % Pseudo frequency FFT axis

OptFreq = 299792./wlCal; % in THz
if (get(handles.uv,'Value') == get(handles.uv,'Max'))
    [OptFreq,PseudoFreq ] = rescaleCali( wlCal,MaxCorr,PP);
else
    OptFreq = 299792./wlCal; %% Y data [ThZ]
end

% Fit Pseudo frequency against optical frequency
axes(handles.axesCalibration)
plot(PseudoFreq(fitind),OptFreq(fitind),'o')

P = polyfit(PseudoFreq(fitind),OptFreq(fitind),1);

hold on
plot(PseudoFreq, P(1)*PseudoFreq+P(2),'r','LineWidth',1.5) 
xlabel('Pseudo Frequencies [steps^{-1}]');
ylabel('Optical Frequencies [THz]');
set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)

vm = P(1).*PseudoFreq + P(2); % in THz

%% Try to calculate the t1 axis from corrected motor frequency axis

P2 = polyfit(f_fft(B),OptFreq,1);
vmtot = (P2(1).*f_fft + P2(2)).*1e12; % in Hz

df = mean(abs(diff(vmtot))); % Frequency resolution in Hz
t_ifft = ((1e15./(2*df)).*linspace(0,1,sizeFFT/2+1)); % t1 ifft Time axis in s

ifftdata = ifft(FFTDataCala,sizeFFT,1);

t_ifft = t_ifft(1:length(mStepCal));
ifftdata = ifftdata(1:length(mStepCal),:);

toc
%% --> Inserire un filtro Nutterworth su PhilCal con 0.008 e 0.2 cut off Frequencies
 
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
         
%% integrazione della mappa raccolta dallo spettrometro

PhIntSpec = sum(DataCal')'./max(sum(DataCal')');

%% plot del fotodiodo

axes(handles.axesPhotodiode)
plot(PhIntCal)

xlabel('t_1 [Wedges Steps]');
ylabel('Amplitude');
ylim([-1 1])
set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)

%% Phase Button%%

PhaseButton = 1; 
DevStdInput = str2double(get(handles.editGaussCal,'String')); 
DevStdInput_rel = DevStdInput*length(PhIntCal);

%% Gaussian Filter%%%%

x = linspace(-length(PhIntCal),length(PhIntCal),length(PhIntCal)*2);

UniGaussWindow = ones(length(PhIntCal)*2,1).*gaussmf(x,[DevStdInput_rel 0])';
UniGaussWindow = UniGaussWindow(length(UniGaussWindow)/2:length(UniGaussWindow)-1);    
UniGaussWindow = UniGaussWindow./max(UniGaussWindow);
yfilt = UniGaussWindow.*PhIntCal;
PhIntCalFilt = yfilt;
    
%% gaussian filter for spectra
     
PhaseButton = 1; %%%Input user
DevStdInput = str2double(get(handles.editGaussCal,'String'));  %%%%Input user
DevStdInput_rel_Spec = DevStdInput*length(PhIntSpec);

%% Gaussian Filter%%%%

x = linspace(-length(PhIntSpec),length(PhIntSpec),length(PhIntSpec)*2);

gauss = exp(-(x-0).^2/(2*DevStdInput_rel_Spec^2));
UniGaussWindow = ones(length(PhIntSpec)*2,1).*gauss';
UniGaussWindow = UniGaussWindow(length(UniGaussWindow)/2:length(UniGaussWindow)-1);    
UniGaussWindow = UniGaussWindow./max(UniGaussWindow);

ySpecfilt = UniGaussWindow.*PhIntSpec;

PhIntSpecFilt = ySpecfilt;

%% Butterworth filter design order=2, bandpass (0.008-0.2)   

% [b,a] = butter(2,[0.008 0.2],'bandpass');
% PhIntCalFilt = filter(b,a,PhIntCalFilt);
% PhIntSpecFilt = filter(b,a,PhIntSpecFilt);

% code to change bin0 reference
PhIntCalFilt = PhIntSpecFilt;

% FFT phasing unwrap.vi%%%%
ZeroPadding = str2double(get(handles.edit0Pad,'String')); %%Input User
padsize = [2^(nextpow2(length(PhIntCalFilt))+ZeroPadding)-length(PhIntCalFilt)];
PhIntCalFiltPad = padarray(PhIntCalFilt,padsize,'post');

%% Now Call pump_axis.vi%%%

InterceptCal = P(2);
SlopeCal = P(1);
FTsize = length(PhIntCalFiltPad);
A = linspace(1,FTsize,FTsize);
A = A/FTsize;
A = fliplr(A);
wlpump = 299792./(SlopeCal.*A+InterceptCal); %%%PumpAxis [nm]
savefile = pat(1:regexp(pat,'2D')-1);

if PhaseButton == 1
    %% FT_mag_phase_unwrap.vi%%%%
    
    UnwrapPhase = 0; %%user button
    Bin0 = str2double(get(handles.editBin0,'String')); %%user choice
    InvertPhase = 0; %%%User button

    PhIntCalFiltPadRot = circshift(PhIntCalFiltPad,-Bin0);
    
    axes(handles.axesPhotodiode)
    plot(PhIntCalFiltPadRot)
    xlabel('t_1 [Wedges Steps]');
    ylabel('Amplitude');
    ylim([-1 1])
    set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
    set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
    set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
    
    FFTPhIntCalFiltPadRot = fft(PhIntCalFiltPadRot, FTsize);
    
    if UnwrapPhase == 1
        phase = unwrap(angle(FFTPhIntCalFiltPadRot));
    else    
        phase = angle(FFTPhIntCalFiltPadRot);
    end
    
    phase = fliplr(phase);
    amplFFT = abs(FFTPhIntCalFiltPadRot);
    
    if InvertPhase == 1
        phase = (-1)*phase;
    end
    
    expIphase = exp(1i*phase);
    
    axes(handles.axesAmplPh)
    plot(wlpump,amplFFT,'LineWidth',1.5)
    
    if (get(handles.uv,'Value') == get(handles.uv,'Max'))
        xlim([300 420])
    else
        if (get(handles.but_OPA3,'Value') == get(handles.but_OPA3,'Max'))
            xlim([650 1000])
        else 
            xlim([450 800])
        end
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
        if (get(handles.but_OPA3,'Value') == get(handles.but_OPA3,'Max'))
            xlim([650 1000])
        else 
            xlim([450 800])
        end
    end
    
    xlabel('\lambda [nm]');
    ylabel('Phase');

    set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
    set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
    set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
% FT_ABS_DISP%%%%

    PhIntCalFiltPadTrunc = PhIntCalFiltPad(Bin0:end);
    KK = fft(PhIntCalFiltPadTrunc,FTsize).*expIphase;

% Calculation of Resolution%%% ---> Fare un check

    [AverageFreq,GroupDelay, Sampling] = Calc_Reso(OptFreq,InterceptCal,...
        SlopeCal,length(PhIntCalFiltPadRot));
    %%nm Output  Fare un check su sampling --> fare un check di analisi dimensinale dovrebbe venire Sampling=3.3 circa AverageWl=650 Resolution=11 T1max=121                                                                                       
    Average_Wavelength = 299792/AverageFreq; %%nm Output
    NumpointScan = length(PhIntCalFilt)-Bin0;
    [Resolution,t1max, Reso] = Calc_Reso(OptFreq,InterceptCal, SlopeCal,NumpointScan);

% Plot in main GUI%%%
 
    HighWlFit = 700; %%input
    LowWlFit = 560; %%input
    [wlpumpCut, interceptPhase, SlopePhase, phasecut] = linearphaseFit(wlpump,amplFFT,...
        phase,HighWlFit,LowWlFit);
end

%% Create global%%
setappdata(0, 'PhIntCalFiltPad' ,PhIntCalFiltPad);
setappdata(0, 'FTsize' ,FTsize);%storing matrix in GUI variables
setappdata(0, 'wlpump' ,wlpump);%storing matrix in GUI variables
setappdata(0, 'PumpSpecaver',PumpSpecaver);
setappdata(0, 'wlCalaver',wlCalaver);   