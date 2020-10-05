function [ ] = film( hObject, eventdata, handles)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%load('2D_CdSe_NRs_10nJ_2DMap.mat')
Map2DSave=handles.map;
 T2In=handles.t;
 wlpump=handles.pump;
 wprCutaver=handles.probe;
set(handles.text58,'String','movie');
N=length(T2In);
[mx,mn]=findmxmn( N,Map2DSave );
for i=1:N
    ss=get(handles.text58,'String');

switch ss
    case'stop'
        i=N;
    case'movie'
        
    

    A=cell2mat(Map2DSave(i));
    c=size(A);
    textLabel=num2str(T2In(i));
    textLabel=strcat(textLabel,'fs');
    set(handles.text57,'String',textLabel);
axes(handles.axes12)
pcolor(wprCutaver,wlpump,A)
colormap jet
shading flat
hold on
[C,hfigc] = contour(wprCutaver, wlpump, A);

set(hfigc, ...
    'LineWidth',1.0, ...
    'Color', [1 1 1]);

 xlabel('\lambda_{probe} [nm]');
 ylabel('\lambda_{pump} [nm]');  
%  set(findall(gcf,'-property','TickLength'),'LineWidth',1.5)
%  set(findall(gcf,'-property','TickLength'),'TickLength',[0.02,0.035])
%  set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
  caxis([mn,mx])
  %colorbar('southoutside')
hold off
pause(0.2)

 

end

end

