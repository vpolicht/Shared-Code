function [wlpump, expIphase, phase, amplFFT]=FFT_phasing_unwrap(handles,hObject,PhInt)
%% FFT phasing unwrap.vi%%%%
    %%%Input from GUI%%%%
    ZeroPadding = str2double(get(handles.editZeroPadding,'String'));%%Input User from handles
    Bin0 = str2double(get(handles.editBin0,'String')); 
    InterceptCal = str2double(get(handles.textCalibInter,'String'));
    SlopeCal = str2double(get(handles.CalibS,'String'));
    UnwrapPhase = 1; %%user button from handles
    InvertPhase = 0; %%%User button from handles
    
    padsize = [ZeroPadding*2^(nextpow2(length(PhInt)))-length(PhInt)];
    PhIntPad = padarray(PhInt,padsize,'post');
    
    %% Now Call pump_axis.vi%%%

    FTsize = length(PhIntPad);
    A = linspace(1,FTsize,FTsize);
    A = A/FTsize;
    A = fliplr(A);
    wlpump = 300000./(SlopeCal.*A+InterceptCal); %%%PumpAxis [nm]

        %% FT_mag_phase_unwrap.vi%%%%
    PhIntPadRot = circshift(PhIntPad,-Bin0); 
    FFTPhIntPadRot = fft(PhIntPadRot,FTsize);
     if UnwrapPhase == 1
        phase = unwrap(angle(FFTPhIntPadRot));
     end
    
    phase = angle(FFTPhIntPadRot);
    phase = fliplr(phase);
    amplFFT = abs(FFTPhIntPadRot);
    
    if InvertPhase == 1
        phase = (-1)*phase;
    end
    expIphase = exp(1i*phase);
    %% FT_ABS_DISP%%%%
    PhIntPadTrunc=PhIntPad(Bin0:end);
    FFTPhIntPadTrunc=fft(PhIntPadTrunc,FTsize).*expIphase;