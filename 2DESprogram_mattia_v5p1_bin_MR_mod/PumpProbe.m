function PumpProbe(hObject,handles,pat,converter)
a=regexp(pat,'\');
ind=pat(1:a(end));
file =dir(pat);
pati=[ind,file.name];
dass=dlmread(pati);
%%%%Attenzione la matrice è trasposta rispetto a quanto scritto nel file di
%%%%descrizione
t2PP=dass(1,2:end);
if (get(handles.checkboxNM,'Value') == get(handles.checkboxNM,'Max'))
wprPP=dass(2:end,1);
else
    wprPP=converter./dass(2:end,1);
end
DataPP=dass(2:end,2:end);
set(hObject,'BackgroundColor','white')
T2=str2double(get(handles.editT2,'String')); %%%Selected by the user
Shift=str2double(get(handles.editShiftTime,'String')); %%%Selected by the user
t2PP=t2PP+Shift;
ReverseSign=0; %%%Selected by the user
if ReverseSign==1
    DataPP=-DataPP;
end
IndT2=find(t2PP>=T2,1);
textLabel = sprintf('%.0f', t2PP(IndT2));
set(handles.textT2, 'String', textLabel);

MinProbe=str2double(get(handles.editMinProbe,'String')); 
MaxProbe=str2double(get(handles.editMaxProbe,'String'));  
PumpAverage=str2double(get(handles.editNumAverProbe,'String')); 
if MinProbe<MaxProbe
IndexProbeSel=find(wprPP>MinProbe & wprPP<MaxProbe);
else
    IndexProbeSel=find(wprPP<MinProbe & wprPP>MaxProbe);
end
wprPPcut=wprPP(IndexProbeSel);
DataPPcut=DataPP(IndexProbeSel,IndT2);
%wprPPcutaver=zeros(length(wprPPcut)/2^PumpAverage,1);
%DataPPcutaver=zeros(length(wprPPcut)/2^PumpAverage,1);
wprPPcutaver=zeros(round(length(wprPPcut)/2^PumpAverage),1);
DataPPcutaver=zeros(round(length(wprPPcut)/2^PumpAverage),1);
for i=1:length(wprPPcut)/2^PumpAverage
    wprPPcutaver(i)=mean(wprPPcut(i*2^PumpAverage-2^PumpAverage+1:i*2^PumpAverage));
    DataPPcutaver(i)=mean(DataPPcut(i*2^PumpAverage-2^PumpAverage+1:i*2^PumpAverage));
end
DataPPcutaver=DataPPcutaver./max(abs(DataPPcutaver));
DataPPcut=DataPPcut./max(abs(DataPPcut));
max(abs(DataPPcut))
min(abs(DataPPcut))
axes(handles.axesPP)
plot(wprPPcut,DataPPcut,'LineWidth',1.5)
hold on 
plot(wprPPcutaver,DataPPcutaver,'LineWidth',1.5,'Color','r')
hold off
xlabel('\lambda_{probe} [nm]');
ylabel('\DeltaT/T'); 
xlim([MinProbe MaxProbe])
legend('Pump-Probe','Averaged Pump-Probe','Location','best')
set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)

    