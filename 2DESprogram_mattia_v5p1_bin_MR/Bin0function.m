function Bin0function(handles)
FTsize=getappdata(0, 'FTsize');
wlpump=getappdata(0, 'wlpump');
PhIntCalFiltPad=getappdata(0, 'PhIntCalFiltPad');
    Bin0=str2double(get(handles.editBin0,'String')); %%user choice
    InvertPhase=0; %%%User button
    UnwrapPhase=1;
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