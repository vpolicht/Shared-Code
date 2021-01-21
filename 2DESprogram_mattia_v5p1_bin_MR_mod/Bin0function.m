function Bin0function(handles)

FTsize = getappdata(0, 'FTsize');
wlpump = getappdata(0, 'wlpump');
PhIntCalFiltPad = getappdata(0, 'PhIntCalFiltPad');
Bin0 = str2double(get(handles.editBin0,'String')); %%user choice

InvertPhase = 0; %%%User button
UnwrapPhase = 0;

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

if InvertPhase==1
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
        xlim([500 800])
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
        xlim([500 800])
    end
end

xlabel('\lambda [nm]');
ylabel('Phase');

set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)